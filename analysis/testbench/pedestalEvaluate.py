from __future__ import print_function
import sys, json, os
import subprocess
import multiprocessing as mp
import numpy as np
from array import array

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import ROOT
# %jsroot on
from ROOT import TCanvas, TH1F, TH2F, TF1, TMath, TGraph, TFile, TSpectrum, TPaveText, TMultiGraph, TGraphErrors, TLine
from ROOT import gStyle, gROOT, gDirectory, gPad

import constants
import crv_event
import crv_spill
import utils
import geometry
import geometry_constants
import filepath

lock = mp.Manager().Lock()

gROOT.Reset()
gROOT.SetBatch(1)
gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")
gStyle.SetOptStat(111110)
gStyle.SetOptFit(0)
gStyle.SetLineScalePS(0.3)

topdir = os.path.dirname(os.path.abspath("__file__"))
analysis_dir = os.path.join(topdir,"analysis_root_files")

tgeometry = geometry.testBenchGeometry(geometry_constants.setup_dict['crvled-001'][0])
filelist = filepath.getfilelist(["LED_low_PE"])[0]
nFEB = geometry_constants.setup_dict['crvled-001'][0]['nFEB']
nCalibPeak = 2

summary_pass1_branch_list = [("FEB", "i"), ("ch", "i"),
                             ("pedestal", "vdouble"), ("calibration","vstring")]
summary_pass2_branch_list = [("spectrum","vstring"), ("spe","vdouble")]
event_pass1_branch_list = [("fitStatus", "vint"), 
                           ("area", "vdouble"), ("pulseHeight", "vdouble"), ("beta", "vdouble"), 
                           ("time", "vdouble"), ("LEtime", "vdouble"), 
                           ("recoStartBin", "vint"), ("recoEndBin", "vint")]
# event_pass2_branch_list = [("PEs", "vdouble")] # use spectrum directly. no need to save again

def getPedestals(filename):
    # get pedestals from the json file
    pedFile = open(os.path.join(analysis_dir,"LED_low_PE_pedestal_run_%04i_%03i.txt"%(runNum, subrunNum)),"r")
    with open(pedFile, 'r') as fin:
        line = fin.readline()
        ped_dict = json.loads(line)    
    pedestalTags = ["stock", "stockRebin", "allPtsExceptP5ADC", "stockRebinM0P15"]
    pedestalList = []
    pedestalList.append(ped_dict["stock"])
    pedestalList.append(ped_dict["stockRebin0"])
    pedestalList.append(ped_dict["allPtsExceptP5ADC"])
    temp_ped = [[None]*geometry_constants.nChannelPerFEB for i in range(nFEB)]
    for iFEB in range(nFEB):
        for iCh in range(geometry_constants.nChannelPerFEB):
            if ped_dict["stockRebin0"][iFEB][iCh]:
                temp_ped[iFEB][iCh] = ped_dict["stockRebin0"][iFEB][iCh] - 0.15
    pedestalList.append(ped_dict[temp_ped])    
    nPed = len(pedestalList)
    return nPed, pedestalTags, pedestalList

def traceProcess(filename, nPed, pedestalTags, pedestalList):
    lock.acquire()
    print("Processing %s..."%filename)
    sys.stdout.flush()
    lock.release()

    runNum = filepath.filenameparser(filename, 'run')
    subrunNum = filepath.filenameparser(filename, 'subrun')
    
    foutname = os.path.join(analysis_dir,"LED_low_PE_pulseReco_run_%04i_%03i.root"%(runNum, subrunNum))
    fout = ROOT.TFile(foutname, "RECREATE")
    fout.mkdir("calibration")
    fout.mkdir("spectrum")
    fout.cd()
    
    summary_tree = ROOT.TTree('summary','summary')
    summary_dict = root_utils.treeInitialization(summary_tree, summary_pass1_branch_list)
    event_tree = ROOT.TTree('event','event')
    event_dict = root_utils.treeInitialization(event_tree, event_pass1_branch_list)
    
    calibPlot = [[[None]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)] for iPed in range(nPed)]
    
    for iFEB in range(nFEB):
        for iCh in range(geometry_constants.nChannelPerFEB):
            root_utils.treeEntryDictPurge(summary_pass1_branch_list, summary_dict)
            summary_dict["FEB"][0] = iFEB
            summary_dict["ch"][0] = iCh
            for iPed in range(nPed):
                summary_dict["pedestal"].push_back(pedestalList[iPed][iFEB][iCh])
                calibName = "calibration/hist_calib_run_%04i_%03i_FEB%i_Ch%i_ped%i"%(runNum, subrunNum, iFEB, iCh, iPed)
                calibPlot[iPed][iFEB][iCh] = ROOT.TH1F(calibName, "Calibration (%s pedestal); ADC*ns; counts"%(pedestalTags[iPed]),
                                                       300, 0., 3000.) 
                
                summary_dict["calibration"].push_back(calibName)
            summary_tree.Fill()
    fout.cd()
    summary_tree.Write()
    
    fFile = ROOT.TFile(filename, "READ")
    fTree = fFile.Get("run")
    nEntries = fTree.GetEntries()
    
    # get entries and fill in...
    for iEntry in range(nEntries):
        
        # for test
        if iEntry == 4000:
            break
            
        if iEntry % 20000 == 0:
            lock.acquire()
            print("Run%04i_%03i fitted %i/%i entries..."%(runNum, subrunNum, iEntry, nEntries))
            sys.stdout.flush()
            lock.release()
        tEvent = crv_event.crv_event(fTree, iEntry, 0b1100, nFEB)
        
        root_utils.treeEntryDictPurge(event_pass1_branch_list, event_dict)
        for iPed in range(nPed):
            # Fill calib plot
            tEvent.FillCalibrationHistogramsStock(pedestalList[iPed], calibPlot[iPed])
            # find peak
            tempFactor = [[1.]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)] # fill tEvent.PEs with area
            tEvent.peakFitter(pedestalList[iPed], tempFactor, tempFactor, False)
            for iFEB in range(nFEB):
                for iCh in range(geometry_constants.nChannelPerFEB):
                    event_dict["fitStatus"].push_back(tEvent.fitStatus[iFEB][iCh])
                    event_dict["area"].push_back(tEvent.PEs[iFEB][iCh])
                    event_dict["pulseHeight"].push_back(tEvent.pulseHeight[iFEB][iCh])
                    event_dict["beta"].push_back(tEvent.beta[iFEB][iCh])
                    event_dict["time"].push_back(tEvent.time[iFEB][iCh])
                    event_dict["LEtime"].push_back(tEvent.LEtime[iFEB][iCh])
                    event_dict["recoStartBin"].push_back(tEvent.recoStartBin[iFEB][iCh])
                    event_dict["recoEndBin"].push_back(tEvent.recoEndBin[iFEB][iCh])
                    
        event_tree.Fill()
    
    fout.cd()
    event_tree.Write("",ROOT.TObject.kWriteDelete)
    
    fout.cd("calibration")
    for iFEB in range(nFEB):
        for iCh in range(geometry_constants.nChannelPerFEB):
            for iPed in range(nPed):
                calibPlot[iPed][iFEB][iCh].Write()
    
    fout.Close()     
    
    lock.acquire()
    print("Run%04i_%03i trace fitting completed."%(runNum, subrunNum))
    sys.stdout.flush()
    # make a backup copy
    os.system("cp %s %s.bak"%(foutname,foutname))
    lock.release()
    
    return
        
def spectrumGen(filename, nPed, pedestalTags, pedestalList):
    runNum = filepath.filenameparser(filename, 'run')
    subrunNum = filepath.filenameparser(filename, 'subrun')
    rootname = "LED_low_PE_pulseReco_run_%04i_%03i.root"%(runNum, subrunNum)
    pdfname = "LED_low_PE_calibration_run_%04i_%03i.pdf"%(runNum, subrunNum)
    pdfname = os.path.join(analysis_dir, pdfname)
    
    fC = TCanvas("c0", "c0", 1200*2, 900)
    fC.Divide(2,1,0.0001,0.001)
    fC.Print(pdfname+"[", "pdf")
    
    lock.acquire()
    print("Generating spectrum for Run%04i_%03i..."%(runNum, subrunNum))
    sys.stdout.flush()
    lock.release()
    
    foutname = os.path.join(analysis_dir, rootname)
    fout = ROOT.TFile(foutname, "UPDATE")
    summary_tree = fout.Get('summary')
    event_tree = fout.Get('event')
    nEntries = event_tree.GetEntries()
    
    summary_pass2_dict, summary_pass2_branches = root_utils.treeInitialization(summary_tree, summary_pass2_branch_list, True)
    
    spectrum = [[[[None]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)] for iCalibPeak in range(nCalibPeak)] for iPed in range(nPed)]
    calibConstant = [[[[np.nan]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)] for iCalibPeak in range(nCalibPeak)] for iPed in range(nPed)]
    
    # fit calib plots and extract calib constants. Compare. Save
    for iFEB in range(nFEB):
        for iCh in range(geometry_constants.nChannelPerFEB):
            index = iFEB * geometry_constants.nChannelPerFEB + iCh
            summary_tree.GetEntry(index)
            root_utils.treeEntryDictPurge(summary_pass2_branch_list, summary_pass2_dict)
            for iPed in range(nPed):
                
                calibName = summary_tree.calibration[iPed]
                tCalibHist = fout.Get(calibName).Clone()
                fC.Clear("D")
                nn, c2, peakLocations = crv_event.CalculateCalibrationConstantSingleChannel(tCalibHist, 2, iFEB, iCh, fC, 1, 2, 0, pedestalTags[iPed], 50)
                c1 = peakLocations[0][0] # for single peak, spe is the location of the first peak 
                calibConstant[iPed][1][iFEB][iCh] = c2
                calibConstant[iPed][0][iFEB][iCh] = c1 
                
                fC.cd(2)
                fl = ROOT.TLine(0., 0., 1., c1)
                fl.SetLineColor(9)
                fl.SetLineWidth(1)
                fl.SetLineStyle(2)
                fl.Draw("SAME")
                t4 = ROOT.TPaveText(.7, .65, .9, .88, "NDC")
                t4.SetLineColorAlpha(0,0)
                t4.SetFillStyle(0)
                t4.AddText("1 PE fit, SPE = %.3f"%c1)
                t4.AddText("2 PE fit, SPE = %.3f"%c2)
                t4.AddText("Difference %.3f%%"%((c2-c1)/c1*100.))
                t4.SetTextAlign(12)
                t4.DrawClone("SAME") 
                
                fC.Print(pdfname, "Title: Calibration FEB %i Ch %i (%s pedestal)"%(iFEB, iCh, pedestalTags[iPed]))
                
                for iCalibPeak in range(nCalibPeak):
                    spectName = "spectrum/hist_spect_run_%04i_%03i_FEB%i_Ch%i_ped%i_%ipeak"%(runNum, subrunNum, iFEB, iCh, iPed, iCalibPeak)
                    specTitle = "PEs (%s pedestal, %i peak, FEB %i Ch %i); nPE (reco); counts"%(pedestalTags[iPed], iCalibPeak+1, iFEB, iCh)
                    spectrum[iPed][iCalibPeak][iFEB][iCh] = ROOT.TH1F(spectName, specTitle, 600, 0., 30.)
                    summary_pass2_dict["spectrum"].push_back(spectName)
                    summary_pass2_dict["spe"].push_back(calibConstant[iPed][iCalibPeak][iFEB][iCh])
                    
            for tbranch in summary_pass2_branches:
                tbranch.Fill()
                
    fout.cd()
    summary_tree.Write("",ROOT.TObject.kWriteDelete)
    
    fC.Clear("D")
    fC.cd(1)
    fMultiGraph = ROOT.TMultiGraph("tCalibByChRatio_run_%04i_%03i"%(runNum, subrunNum), 
                                   "Calibration / Stock Calibration; Channel; Calibration Ratio")
    TG = [None]*8
    TH = [None]*8
    legText = ["stock,1PE", "stock,2PE", "stockRebin,1PE", "stockRebin,2PE", "allPtsExcept5+ADC,1PE", "allPtsExcept5+ADC,2PE", "stockRebin-0.15,1PE", "stockRebin-0.15,2PE"]
    channels = np.array([ii for ii in range(nFEB*geometry_constants.nChannelPerFEB)])
    fLegend = ROOT.TLegend(0.70, 0.5, 0.95, 0.9)
    for iPed in range(nPed):
        for iCalibPeak in range(nCalibPeak):
            i = iPed * nCalibPeak + iCalibPeak
            tCalibRatio = np.array([[calibConstant[iPed][iCalibPeak][iFEB][iCh]/calibConstant[0][0][iFEB][iCh] for iCh in range(nFEB*geometry_constants.nChannelPerFEB)] for iFEB in range(nFEB)]).flatten()
            tfilter = np.isnan(tCalibRatio)
            TG[i] = ROOT.TGraph(tCalibRatio[tfilter].size, array('d',channels[tfilter].tolist()), array('d',tCalibRatio[tfilter]).tolist())
            TG[i].SetName("tCalibByChRatio%i_run_%04i_%03i"%(i, runNum, subrunNum))
            TG[i].SetMarkerStyle(constants.rootmarkers[i])
            TG[i].SetMarkerColor(constants.rootcolors[i])
            TG[i].SetLineColor(constants.rootcolors[i])
            TG[i].SetMarkerSize(1)
            fMultiGraph.Add(TG[i])
            fLegend.AddEntry(TG[i], legText[i])
            TH[i] = ROOT.TH1F("tHStackCalibByChRatio_%irun_%04i_%03i"%(i, runNum, subrunNum),
                              "Calibration constants compared to the stock method; pedestal/pedestal_{0}; count",
                              50, 0.98, 1.08)
            for v in tCalibRatio[tfilter]:
                TH[i].Fill(v)
    fMultiGraph.GetXaxis().SetLimits(-0.5,150.)
    fMultiGraph.Draw("ALP")
    fLegend.Draw("SAME")
    
    fC.cd(2)
    fStack = ROOT.THStack("tHStackCalibByChRatio_run_%04i_%03i"%(runNum, subrunNum),
                          "Calibration constants compared to the stock method; pedestal/pedestal_{0}; count")
    fLegend = ROOT.TLegend(0.70, 0.3, 0.95, 0.9)
    for iPed in range(nPed):
        for iCalibPeak in range(nCalibPeak):
            i = iPed * nCalibPeak + iCalibPeak
            TH[i].SetLineColor(constants.rootcolors[i])
            TH[i].SetLineWidth(2)
            if i>0:
                fStack.Add(TH[i])
                fLegend.AddEntry(TH[i], "#split{%s}{#mu = %.3lf, #sigma = %.3lf}"%(legText[i], TH[i].GetMean(), TH[i].GetRMS()))
    fStack.Draw("NOSTACK")
    fLegend.Draw("SAME")
    
    fC.Print(pdfname+')', "Title: Calibration Constants Compared to Stock Value")
    lock.acquire()
    print("Run%04i_%03i calibration constants calculated."%(runNum, subrunNum))
    sys.stdout.flush()
    lock.release()
                    
    # calculate PEs, save.
    for iEntry in range(nEntries):
        
        # for test
        if iEntry == 4000:
            break
            
        if iEntry % 20000 == 0:
            lock.acquire()
            print("Run%04i_%03i processed %i/%i entries for spectrum generation..."%(runNum, subrunNum, iEntry, nEntries))
            sys.stdout.flush()
            lock.release()
        
        event_tree.GetEntry(iEntry)
        for iPed in range(nPed):
            tFitStatus = event_tree.fitStatus[iPed*geometry_constants.nChannelPerFEB*nFEB : (iPed+1)*geometry_constants.nChannelPerFEB*nFEB]
            tArea = event_tree.area[iPed*geometry_constants.nChannelPerFEB*nFEB : (iPed+1)*geometry_constants.nChannelPerFEB*nFEB]
            for iFEB in range(nFEB):
                for iCh in range(geometry_constants.nChannelPerFEB):
                    index = iFEB*geometry_constants.nChannelPerFEB+iCh
                    if tFitStatus[index] == 1:
                        if np.isnan(calibConstant[iPed][0][iFEB][iCh])==False:
                            tPE1 = tArea[index]/calibConstant[iPed][0][iFEB][iCh]
                            spectrum[iPed][0][iFEB][iCh].Fill(tPE1)
                        if np.isnan(calibConstant[iPed][1][iFEB][iCh])==False:
                            tPE2 = tArea[index]/calibConstant[iPed][1][iFEB][iCh]
                            spectrum[iPed][1][iFEB][iCh].Fill(tPE2)
                        
    fout.cd("spectrum")
    for iFEB in range(nFEB):
        for iCh in range(geometry_constants.nChannelPerFEB):
            for iPed in range(nPed):
                for iCalibPeak in range(nCalibPeak):
                    spectrum[iPed][iCalibPeak][iFEB][iCh].Write()
    
    fout.Close()     
    
    lock.acquire()
    print("Run%04i_%03i spectrum generation completed."%(runNum, subrunNum))
    sys.stdout.flush()
    lock.release()
    
    return     

def compareSpectrum(filename, nPed, pedestalTags, pedestalList)        
    runNum = filepath.filenameparser(filename, 'run')
    subrunNum = filepath.filenameparser(filename, 'subrun')
    rootname = "LED_low_PE_pulseReco_run_%04i_%03i.root"%(runNum, subrunNum)
    pdfname = "LED_low_PE_spectrum_run_%04i_%03i.pdf"%(runNum, subrunNum)
    pdfname = os.path.join(analysis_dir, pdfname)
    
    fC1 = TCanvas("c1", "c1", 600*nPed, 450*nCalibPeak)
    fC1.Divide(nPed,nCalibPeak,0.0001,0.001)
    fC1.Print(pdfname+"[", "pdf")
    fC2 = TCanvas("c2", "c1", 600*nCalibPeak, 450)
    fC2.Divide(nCalibPeak,1,0.0001,0.001)
    
    lock.acquire()
    print("Plotting spectrum for Run%04i_%03i..."%(runNum, subrunNum))
    sys.stdout.flush()
    lock.release()
    
    finname = os.path.join(analysis_dir, rootname)
    fin = ROOT.TFile(finname, "READ")
    summary_tree = fout.Get('summary')
    
    spectrum = [[[[None]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)] for iCalibPeak in range(nCalibPeak)] for iPed in range(nPed)]
    peakGraph = [[[[None]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)] for iCalibPeak in range(nCalibPeak)] for iPed in range(nPed)]
    peakMG = [[[None]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)] for iCalibPeak in range(nCalibPeak)]
    fLegend = [[[None]*geometry_constants.nChannelPerFEB for iFEB in range(nFEB)] for iCalibPeak in range(nCalibPeak)]
    legText = ["stock,1PE", "stockRebin,1PE", "allPtsExcept5+ADC,1PE", "stockRebin-0.15,1PE", 
               "stock,2PE", "stockRebin,2PE", "allPtsExcept5+ADC,2PE", "stockRebin-0.15,2PE"]
    
    for iFEB in range(nFEB):
        for iCh in range(geometry_constants.nChannelPerFEB):
            index = iFEB * geometry_constants.nChannelPerFEB + iCh
            fC1.Clear("D")
            fC1.SetLogy()
            fC2.Clear("D")
            fC2.SetLogy(0)
            summary_tree.GetEntry(index)
            for iPed in range(nPed):
                for iCalibPeak in range(nCalibPeak):                   
                    tSpecName = summary_tree.spectrum[iPed*nCalibPeak+iCalibPeak]
                    spectrum[iPed][iCalibPeak][iFEB][iCh] = fout.Get(tSpecName).Clone()
            
            # draw spectrum
            for iCalibPeak in range(nCalibPeak):
                peakMG[iCalibPeak][iFEB][iCh] = ROOT.TMultiGraph("peakMG_FEB%i_Ch%i_%i"%(iFEB, iCh, iCalibPeak), "Reconstruction Quality FEB%i Ch%i (calibration with %iPE); PE; nPE(true)/nPE(reco)"%(iFEB, iCh, 1+iCalibPeak))
                fLegend[iCalibPeak][iFEB][iCh] = ROOT.TLegend(0.75, 0.15, 0.95, 0.6)
                for iPed in range(nPed):
                    tpad = fC1.cd(1+iCalibPeak*nPed+iPed)
                    spectrum[iPed][iCalibPeak][iFEB][iCh].SetLineColor(4)
                    spectrum[iPed][iCalibPeak][iFEB][iCh].Draw()
                    tpad.Update()
                    spectrum[iPed][iCalibPeak][iFEB][iCh].GetXaxis().SetNdivisions(30)
                                        
                    tSpec = ROOT.TSpectrum(30)
                    nfoundraw = tSpec.Search(spectrum[iPed][iCalibPeak][iFEB][iCh])
                    peposraw = tSpec.GetPositionX()
                    peposraw.sort()
                    
                    pereco = []
                    petrue = []
                    peratio = []
                    
                    for iPE in range(nfoundraw):
                        if len(petrue)==0:
                            if peposraw[iPE]<0.8:
                                continue
                            else:
                                petrue.append(peposraw[iPE])
                                pereco.append(round(peposraw[iPE], 0))
                        else:
                            if peposraw[iPE]-petrue[-1]<0.7:
                                continue
                            elif peposraw[iPE]-petrue[-1]<2.5:
                                petrue.append(peposraw[iPE])
                                pereco.append(pereco[-1]+round(peposraw[iPE]-petrue[-1], 0))
                            else:
                                break
                    for ii in range(len(pereco)):
                        
                        func = ROOT.TF1("f%i%i%i%i%i"%(iFEB,iCh,iPed,iCalibPeak,ii), "gaus", petrue[ii]-0.4, petrue[ii]+0.4)
                        func.SetParameter(1, petrue[ii])
                        func.SetParLimits(0, 0., 10.*spectrum[iPed][iCalibPeak][iFEB][iCh].GetEntries())
                        func.SetLineWidth(1)
                        func.SetLineColor(2) # kRed
                        spectrum[iPed][iCalibPeak][iFEB][iCh].Fit(func, "QR0")
                        func.DrawClone("SAME")
                        petrue[ii] = func.GetParameter(1)  
                        peratio.append(petrue[ii]/pereco[ii])
                    
                    peakGraph[iPed][iCalibPeak][iFEB][iCh] = ROOT.TGraph(len(peratio), array('d', pereco), array('d', peratio))
                    peakGraph[iPed][iCalibPeak][iFEB][iCh].SetName("peakGraph_run_%04i_%03i_FEB%i_Ch%i_%i_%i"%(runNum, subrunNum, iFEB, iCh, iPed, iCalibPeak))
                    peakGraph[iPed][iCalibPeak][iFEB][iCh].SetTitle("Reconstruction Quality FEB%i Ch%i (calibration with %iPE); PE; nPE(true)/nPE(reco)"%(iFEB, iCh, 1+iCalibPeak))
                    peakGraph[iPed][iCalibPeak][iFEB][iCh].SetMarkerStyle(constants.rootmarkers[i])
                    peakGraph[iPed][iCalibPeak][iFEB][iCh].SetMarkerColor(constants.rootcolors[i])
                    peakGraph[iPed][iCalibPeak][iFEB][iCh].SetLineColor(constants.rootcolors[i])
                    peakGraph[iPed][iCalibPeak][iFEB][iCh].SetLineWidth(1)
                    peakGraph[iPed][iCalibPeak][iFEB][iCh].SetMarkerSize(1)
                    peakMG[iCalibPeak][iFEB][iCh].Add(peakGraph[iPed][iCalibPeak][iFEB][iCh])
                    fLegend[iCalibPeak][iFEB][iCh].AddEntry(peakGraph[iPed][iCalibPeak][iFEB][iCh], legText[iCalibPeak*nPed+iPed])      
                    
            f1.Print(pdfname, "Title: Spectrum FEB %i Ch %i"%(iFEB, iCh))
            
            # fit peak positions, compare
            for iCalibPeak in range(nCalibPeak):
                fC2.cd(1+iCalibPeak)
                peakMG[iCalibPeak][iFEB][iCh].GetXaxis().SetLimits(0.,30.)
                peakMG[iCalibPeak][iFEB][iCh].Draw("ALP")
                fl = ROOT.TLine(0., 1., 30., 1)
                fl.SetLineColor(1)
                fl.SetLineWidth(1)
                fl.SetLineStyle(2)
                fl.Draw("SAME")
                fLegend[iCalibPeak][iFEB][iCh].Draw("SAME")
            f2.Print(pdfname+(')' if (iFEB==nFEB-1 and iCh==geometry_constants.nChannelPerFEB-1 ) else ''), "Title: Spectrum peak positions FEB %i Ch %i"%(iFEB, iCh))

    fin.Close()
    
    lock.acquire()
    print("Run%04i_%03i done."%(runNum, subrunNum))
    sys.stdout.flush()
    lock.release()
    
    return

def routine(filename):
    nPed, pedestalTags, pedestalList = getPedestals(filename)
    traceProcess(filename, nPed, pedestalTags, pedestalList)
    spectrumGen(filename, nPed, pedestalTags, pedestalList)
    compareSpectrum(filename, nPed, pedestalTags, pedestalList)
    return
    
def main():
    nWorker = 7
    p = mp.Pool(processes=nWorker)
    # for filename in filelist[:2]:
        # p.apply_async(routine, (filename,))
    # p.close()
    # p.join()
    routine(filelist[0])
    
    return 0
    
if __name__ == "__main__":
    main()