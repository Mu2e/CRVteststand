from __future__ import print_function
import sys, os
import gc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from multiprocessing import Pool

import ROOT
from ROOT import gROOT, gStyle, gDirectory, gPad

from array import array

import constants
import crv_event
import crv_spill
import utils
import geometry
import geometry_constants
import filepath

import pandas as pd
import pickle
from scipy import stats

stdFEBTemp = 45 # degC
stdCMBTemp = 20 # degC
SPEperFEBTempVover = 0.51 # SPE/degC, Vover contribution
SPEperFEBTempAFEGain = -1.46 # SPE/degC, AFE gain contribution (at 0x384)
SPEperFEBTemp = SPEperFEBTempVover+SPEperFEBTempAFEGain # SPE/degC, both AFE and Vover contribution
SPEperFEBTempVover = -0.95 # SPE/degC, both AFE and Vover contribution
SPEperCMBTemp = -6.95 # SPE/degC

gROOT.Reset()
gROOT.SetBatch(1)
gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")
gStyle.SetOptStat(111110)
gStyle.SetOptFit(0)
gStyle.SetLineScalePS(0.3)

# -35 dB / 1 V; 4096 ticks spanning 1.54V
def dBratio (setting, ref = 0x384):
    settingV = float(setting)/4096.*1.54
    refV = float(ref)/4096.*1.54
    diffdB = (settingV - refV) * -35.
    diffAtten = pow(10. , (diffdB/20.))
    return diffAtten

def eventCalibWrap(tEvent, pedestal, i):
    try:
        pulseSizeList = tEvent.FillCalibrationHistogramsStock(pedestal, None, None, True)
        return pulseSizeList
    except Exception as e:
        print(f"Exception in eventCalibWrap, %i: {e}"%i)
        return []

firstPass = False

def main():
    df = pd.DataFrame()

    fileList = filepath.getfilelist(["gainAFEScan"],"recoROOT",4)[0]
    topdir = os.path.dirname(os.path.abspath("__file__"))
    analysis_dir = os.path.join(topdir,"analysis_root_files")
    foutname = os.path.join(analysis_dir,"gainAFEScanBuffer.root")
    pdfname = os.path.join(analysis_dir,"gainAFEScanRefit.pdf")
    pklname = os.path.join(analysis_dir,"gainAFEScanDF.pkl")

    if firstPass:
        fout = ROOT.TFile(foutname, "UPDATE")
        # fout = ROOT.TFile(foutname, "RECREATE")
    else:
        fout = ROOT.TFile(foutname, "READ")
    fC = ROOT.TCanvas("c0", "c0", 1200*2, 900*2)
    fC.Divide(2,2,0.0001,0.001)
    fC.Print(pdfname+"[", "pdf")

    for k, filename in enumerate(fileList):
        print(filename)
        runNum = filepath.filenameparser(filename, 'run')
        subrunNum = filepath.filenameparser(filename, 'subrun')
        ii = filepath.datatag['gainAFEScan']['run#'].index(runNum)
        AFEsetting = filepath.datatag['gainAFEScan']['gainAFE'][ii]

        fFile = ROOT.TFile(filename, "READ")
        runSummaryTree = fFile.Get("runSummary")
        spillTree = fFile.Get("spills")
        spillTree.GetEntry(0)
        nFEB = spillTree.spill_number_of_febs
        runSummaryTree.GetEntry(0)

        pedestal = np.reshape(np.array(runSummaryTree.pedestals, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        # calibRaw = np.reshape(np.array(runSummaryTree.calibConstants, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        # calibAdj = np.reshape(np.array(runSummaryTree.calibConstantsTemperatureCorrected, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        FEBtemp = np.array(runSummaryTree.febTemperaturesAvg, dtype=np.float32)
        CMBtemp = np.reshape(np.array(runSummaryTree.meanTemperatures, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        biasV = np.reshape(np.array(runSummaryTree.biasVoltagesAvg, dtype=np.float32), (nFEB, int(geometry_constants.nChannelPerFEB/8)))
        PEs = np.reshape(np.array(runSummaryTree.PEs, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        PEsCorrected = np.reshape(np.array(runSummaryTree.PEsTemperatureCorrected, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))

        calib_hist = [[None]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)]
        calibcorr1_hist = [[None]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)]
        calibcorr2_hist = [[None]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)]
        calibcorr3_hist = [[None]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)]
        fitfunc = [[[None]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)] for i in range(4)]

        if firstPass:
            for iFEB in range(nFEB):
                for iCh in range(geometry_constants.nChannelPerFEB):
                    calibName = "hist_calib_run_%04i_%03i_FEB%i_Ch%i"%(runNum, subrunNum, iFEB, iCh)
                    calib_hist[iFEB][iCh] = ROOT.TH1F(calibName, "Calibration (Raw) run_%04i_%03i_FEB%i_Ch%i; ADC*ns; counts"%(runNum, subrunNum, iFEB, iCh),
                                                      300, 0., 3000.) 
                    calibcorr1_hist[iFEB][iCh] = ROOT.TH1F(calibName+'_corrected1', "Calib. (stock temp. correction); ADC*ns; counts",
                                                          300, 0., 3000.) 
                    calibcorr2_hist[iFEB][iCh] = ROOT.TH1F(calibName+'_corrected2', "Calib. (CMB temp. corrected only); ADC*ns; counts",
                                                          300, 0., 3000.)
                    calibcorr3_hist[iFEB][iCh] = ROOT.TH1F(calibName+'_corrected3', "Calib. (w/ correction gain scaled); ADC*ns; counts",
                                                          300, 0., 3000.) 
            runTree = fFile.Get("run")
            iEvent = 0
            nEvent = runTree.GetEntries()
            nSpill = spillTree.GetEntries()
            nSpillUsed = 0

            for iSpill in range(nSpill):
                tSpill = crv_spill.crv_spill(spillTree, iSpill, False,  False) # not raw, not cosmic run
                nEventInSpill = tSpill.nEventsActual

                if nEventInSpill > 0:
                    tSpill.getTempCMB(runTree, nFEB, iEvent, False)

                # if iSpill != 0 and iSpill%5 == 0:
                    # print("Run%04i_%03i processing %i/%i spills..."%(runNum,subrunNum,iSpill,nSpill))

                tSpill.checkDQM(False) # dqmVerbose = False
                if tSpill.tsEpoch == None:
                    if iSpill == 0: # first spills of subruns >= 001 have timestamp processing issues... refer to next spill for timestamp
                        try:
                            nextSpillTsEpoch = crv_spill.crv_spill(spillTree, 1, False, False).tsEpoch
                            tSpill.tsEpoch = nextSpillTsEpoch - 120
                            tnbit = tSpill.dqmRedFlag.bit_length()
                            tmask = (((1<<tnbit) - 1) ^ 0x10) # this gives 0b1...101111
                            tSpill.dqmRedFlag &= tmask
                            print ("*** spill %04i timestamp amended: %i-%i"%(tSpill.spillNumber, nextSpillTsEpoch, 120))
                        except:
                            sys.exit("Error: %s reading first 2 spill timestamps failed"%(filename))
                    else:
                        print ("!!! spill %04i skipped; the spill has no valid timestamp; 0x%x"%(tSpill.spillNumber, tSpill.dqmRedFlag))
                        continue
                DQMgood = (tSpill.dqmRedFlag==0x0)
                if DQMgood:
                    # populate a list of events
                    argumentList = []
                    for index in range(nEventInSpill):
                        argumentList.append((crv_event.crv_event(runTree, iEvent+index, 0b1100, nFEB), pedestal,index))
                    with Pool(processes = 8) as pool:
                        results = pool.starmap(eventCalibWrap, argumentList)
                    
                    # print(results)
                    
                    for index in range(nEventInSpill):
                        pulseSizeList = results[index]
                        
                        if not pulseSizeList:
                            continue

                        for iFEB in range(nFEB):
                            for iCh in range(geometry_constants.nChannelPerFEB):
                                tcorrection1 = (stdFEBTemp - tSpill.temperatureFEB[iFEB]) * SPEperFEBTemp + (stdCMBTemp - tSpill.temperatureCMB[iFEB, iCh]) * SPEperCMBTemp
                                tcorrection2 = (stdCMBTemp - tSpill.temperatureCMB[iFEB, iCh]) * SPEperCMBTemp
                                tcorrection3 = tcorrection1 * dBratio(AFEsetting)
                                for val in pulseSizeList[iFEB][iCh]:
                                    calib_hist[iFEB][iCh].Fill(val)
                                    calibcorr1_hist[iFEB][iCh].Fill(val+tcorrection1)
                                    calibcorr2_hist[iFEB][iCh].Fill(val+tcorrection2)
                                    calibcorr3_hist[iFEB][iCh].Fill(val+tcorrection3)
                        
                        # del pulseSizeList
                        # tEvent.adc = None
                        # tEvent = None
                        # del tEvent
                        
                        if index % 500 == 499:
                            gc.collect()
                    
                iEvent += nEventInSpill
                print("Run%04i_%03i processed %i/%i events..."%(runNum,subrunNum,iEvent,nEvent))
                if iEvent > nEvent:
                    sys.exit("Error: %s nEvent = %i, Spill # %i is reading iEvent = %i"%(filename, nEvent, iSpill, iEvent))

            fout.cd()
            for iFEB in range(nFEB):
                for iCh in range(geometry_constants.nChannelPerFEB):
                    calib_hist[iFEB][iCh].Write()
                    calibcorr1_hist[iFEB][iCh].Write()
                    calibcorr2_hist[iFEB][iCh].Write()
                    calibcorr3_hist[iFEB][iCh].Write()
            print("Done filling histograms.")
        else:        
            for iFEB in range(nFEB):
                for iCh in range(geometry_constants.nChannelPerFEB):
                    calibName = "hist_calib_run_%04i_%03i_FEB%i_Ch%i"%(runNum, subrunNum, iFEB, iCh)
                    calib_hist[iFEB][iCh] = fout.Get(calibName).Clone()
                    calibcorr1_hist[iFEB][iCh] = fout.Get(calibName+'_corrected1').Clone()
                    calibcorr2_hist[iFEB][iCh] = fout.Get(calibName+'_corrected2').Clone()
                    calibcorr3_hist[iFEB][iCh] = fout.Get(calibName+'_corrected3').Clone()

        # fit calib and plot
        for iFEB in range(nFEB):
            for iCh in range(geometry_constants.nChannelPerFEB):
                fC.Clear("D")
                peak1PE_list = [0.]*4
                peak1PEerr_list = [0.]*4
                for j in range(4):
                    tPad = fC.cd(j+1)
                    tPad.SetLogy()
                    tHist = None
                    if j == 0:
                        tHist = calib_hist[iFEB][iCh]
                    elif j == 1:
                        tHist = calibcorr1_hist[iFEB][iCh]
                    elif j == 2:
                        tHist = calibcorr2_hist[iFEB][iCh]
                    elif j == 3:
                        tHist = calibcorr3_hist[iFEB][iCh]
                    tHist.Draw("")

                    rawFitLow = 400. * dBratio(AFEsetting, 0x300)
                    rawBinLow = tHist.FindBin(rawFitLow)
                    rawFitHigh = 600. * dBratio(AFEsetting, 0x300)
                    rawBinHigh = tHist.FindBin(rawFitHigh)
                    maxbin1 = 0
                    maxbin2 = 0
                    maxbin3 = 0
                    maxbinContent = 0
                    for ibin in range(rawBinLow, rawBinHigh+1):
                        binContent = tHist.GetBinContent(ibin)
                        if binContent > maxbinContent:
                            maxbin1=ibin
                            maxbin2=maxbin1
                            maxbin3=maxbin2
                            maxbinContent=binContent
                    peak1PE = tHist.GetBinCenter(tHist.FindBin((tHist.GetBinCenter(maxbin1)+tHist.GetBinCenter(maxbin2)+tHist.GetBinCenter(maxbin3))/3.))
                    fitfuncmin = peak1PE*0.85
                    fitfuncmax = peak1PE*1.15
                    if fitfuncmax-fitfuncmin<70:
                        fitfuncmin = peak1PE-34.9
                        fitfuncmax = peak1PE+34.9
                    fitfunc[j][iFEB][iCh] = ROOT.TF1("f%04i_%03i_FEB%i_Ch%i"%(runNum, subrunNum, iFEB, iCh), "gaus", fitfuncmin, fitfuncmax)
                    fitfunc[j][iFEB][iCh].SetParameter(1, peak1PE)
                    fitfunc[j][iFEB][iCh].SetLineWidth(1)
                    fitfunc[j][iFEB][iCh].SetLineColor(2) # kRed
                    tHist.Fit(fitfunc[j][iFEB][iCh], "QR")
                    peak1PE_list[j] = fitfunc[j][iFEB][iCh].GetParameter(1)
                    peak1PEerr_list[j] = fitfunc[j][iFEB][iCh].GetParError(1)
                    fitfunc[j][iFEB][iCh].DrawClone("SAME")
                    
                    t2 = ROOT.TPaveText(.4, .65, .75, .89, "NDC")
                    t2.SetFillColor(0)
                    t2.SetFillColorAlpha(0,0)
                    t2.AddText("1 PE = %.1f #pm %.1f ADC*ns"%(fitfunc[j][iFEB][iCh].GetParameter(1),fitfunc[j][iFEB][iCh].GetParError(1)))
                    t2.AddText("#sigma = %.1f #pm %.1f ADC*ns"%(fitfunc[j][iFEB][iCh].GetParameter(2),fitfunc[j][iFEB][iCh].GetParError(2)))
                    t2.AddText("#chi^{2}/NDF = %.2f / %i"%(fitfunc[j][iFEB][iCh].GetChisquare(),fitfunc[j][iFEB][iCh].GetNDF()))
                    if j >= 1:
                        t2.AddText("ref CMB T %.0f deg C"%stdCMBTemp)
                        if j!=2:
                            t2.AddText("ref FEB T %.0f deg C"%stdFEBTemp)
                    t2.DrawClone("SAME")

                fC.Print(pdfname, "Title: Run%04i_%03i FEB%i Ch%i"%(runNum, subrunNum, iFEB, iCh))         

                df_ = pd.DataFrame()
                df_['run'] = [runNum]
                df_['subrun'] = [subrunNum]
                df_['gainAFE'] = [AFEsetting]
                df_['FEB'] = [iFEB]
                df_['ch'] = [iCh]
                df_['pedestal'] = [pedestal[iFEB][iCh]]
                # df_['calibRaw'] = [calibRaw[iFEB][iCh]]
                df_['calibRaw'] = [peak1PE_list[0]]
                df_['calibRawErr'] = [peak1PEerr_list[0]]
                # df_['calibAdj'] = [calibAdj[iFEB][iCh]]
                df_['calibAdj1'] = [peak1PE_list[1]]
                df_['calibAdj1Err'] = [peak1PEerr_list[1]]
                df_['calibAdj2'] = [peak1PE_list[2]]
                df_['calibAdj2Err'] = [peak1PEerr_list[2]]
                df_['calibAdj3'] = [peak1PE_list[3]]
                df_['calibAdj3Err'] = [peak1PEerr_list[3]]            
                df_['FEBtemp'] = [FEBtemp[iFEB]]
                df_['CMBtemp'] = [CMBtemp[iFEB][iCh]]
                df_['biasV'] = [biasV[iFEB][int(iCh/8)]]
                # df_['PEs'] = [PEs[iFEB][iCh]]
                # df_['PEsCorrected'] = [PEsCorrected[iFEB][iCh]]

                # display(df_)

                df = pd.concat([df,df_], ignore_index=True)

        fFile.Close()

    fC.Print(pdfname+"]", "pdf")
    fout.Close()
    # display(df)

    with open(pklname, 'wb') as f:
        pickle.dump(df, f)
    
    return
    
if __name__ == "__main__":
    main()
