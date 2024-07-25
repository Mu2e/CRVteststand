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

# for each subrun file, draw pedstal history, 3 pedestal comparisons
def checkPedestal(filename):
    lock.acquire()
    print("Processing %s..."%filename)
    sys.stdout.flush()
    lock.release()
    runNum = filepath.filenameparser(filename, 'run')
    subrunNum = filepath.filenameparser(filename, 'subrun')
    fout = open(os.path.join(analysis_dir,"LED_low_PE_pedestal_run_%04i_%03i.txt"%(runNum, subrunNum)),"w")
    pdfname = os.path.join(analysis_dir,"LED_low_PE_pedestal_run_%04i_%03i.pdf"%(runNum, subrunNum))
    pdfpages = PdfPages(os.path.join(analysis_dir,"LED_low_PE_pedestal_outlier_run_%04i_%03i.pdf"%(runNum, subrunNum)))   
    pdfStart = True
    
    # n_outlier_plotted = 0
    last_local_ped = None
    
    fFile = ROOT.TFile(filename, "READ")
    fTree = fFile.Get("run")
    nEntries = fTree.GetEntries()
    
    graph_pedHistory = [[None]*geometry_constants.nChannelPerFEB for i in range(2)]
    hist_pedStock = [[None]*geometry_constants.nChannelPerFEB for i in range(2)]
    hist_pedStockRebin = [[None]*geometry_constants.nChannelPerFEB for i in range(2)]
    hist_pedAllPts = [[None]*geometry_constants.nChannelPerFEB for i in range(2)]
    hist_pedAllPtsExceptP5ADC = [[None]*geometry_constants.nChannelPerFEB for i in range(2)]
    for iFEB in range(nFEB):
        for iCh in range(geometry_constants.nChannelPerFEB):
            graph_pedHistory[iFEB][iCh] = ROOT.TGraph()
            graph_pedHistory[iFEB][iCh].SetName("graph_pedHistory_run_%04i_%03i_FEB%i_Ch%i"%(runNum, subrunNum, iFEB, iCh))
            graph_pedHistory[iFEB][iCh].SetTitle("Local pedestal estimation; Time [ns]; Pedestal[ADC]")
        
            hist_pedStock[iFEB][iCh] = ROOT.TH1F("hist_pedStock_run_%04i_%03i_FEB%i_Ch%i"%(runNum, subrunNum, iFEB, iCh),
                                                 "Pedestal (stock method) Run_%04i_%03i_FEB%i_Ch%i; ADC; counts"%(runNum, subrunNum, iFEB, iCh),
                                                 1000, -50., 50.) 
            hist_pedStockRebin[iFEB][iCh] = ROOT.TH1F("hist_pedStockRebin_run_%04i_%03i_FEB%i_Ch%i"%(runNum, subrunNum, iFEB, iCh),
                                                      "Pedestal (stock method rebinned) Run_%04i_%03i_FEB%i_Ch%i; ADC; counts"%(runNum, subrunNum, iFEB, iCh),
                                                      2001, -50.025, 50.025) 
            hist_pedAllPts[iFEB][iCh] = ROOT.TH1F("hist_pedAllPts_run_%04i_%03i_FEB%i_Ch%i"%(runNum, subrunNum, iFEB, iCh),
                                                  "Pedestal (stock method rebinned) Run_%04i_%03i_FEB%i_Ch%i; ADC; counts"%(runNum, subrunNum, iFEB, iCh),
                                                  101, -50.5, 50.5) 
            hist_pedAllPtsExceptP5ADC[iFEB][iCh] = ROOT.TH1F("hist_pedAllPtsExceptP5ADC_run_%04i_%03i_FEB%i_Ch%i"%(runNum, subrunNum, iFEB, iCh),
                                                             "Pedestal (stock method rebinned) Run_%04i_%03i_FEB%i_Ch%i; ADC; counts"%(runNum, subrunNum, iFEB, iCh),
                                                             101, -50.5, 50.5) 
    
    # get entries and fill in...
    for iEntry in range(nEntries):
        
        # for test
        # if iEntry == 4000:
            # break
        
        if iEntry % 20000 == 0:
            lock.acquire()
            print("Run%04i_%03i processed %i/%i entries..."%(runNum, subrunNum, iEntry, nEntries))
            sys.stdout.flush()
            lock.release()
        tEvent = crv_event.crv_event(fTree, iEntry, 0b1100, nFEB)
        
        tEvent.FillPedestalHistogramsStock(hist_pedStock)
        tEvent.FillPedestalHistogramsStock(hist_pedStockRebin)
        tEvent.FillPedestalHistogramsAllPts(hist_pedAllPts)
        local_ped = tEvent.FillPedestalHistogramsAllPtsExceptP5ADC(hist_pedAllPtsExceptP5ADC)
        for iFEB in range(nFEB):
            for iCh in range(geometry_constants.nChannelPerFEB):
                t = tEvent.timeSinceSpill[iFEB][iCh]
                v = local_ped[iFEB][iCh]
                if v:
                    graph_pedHistory[iFEB][iCh].SetPoint(graph_pedHistory[iFEB][iCh].GetN(), t, v)
                
                if iEntry >0:
                    if v:
                        if last_local_ped[iFEB][iCh]:
                            if (v-last_local_ped[iFEB][iCh]<-3 or v-last_local_ped[iFEB][iCh]>3):
                                # if n_outlier_plotted<5:
                                if True:
                                    fig, ax = plt.subplots()
                                    tEvent.plotTrace(ax, iFEB, iCh, "Event %i, local ped. = %.2f"%(tEvent.eventNumber,v))
                                    pdfpages.savefig(fig)
                                    plt.close(fig)
                                    # n_outlier_plotted += 1
                                local_ped[iFEB][iCh] = None                               
                
        last_local_ped = local_ped
        
    # draw pedestal histories
    fC = TCanvas("c", "c", 1200, 900)
    
    for iFEB in range(nFEB):
        for iPlot in range(int(geometry_constants.nChannelPerFEB/4)):
            fC.Clear()
            fMultiGraph = ROOT.TMultiGraph("tMulti_run_%04i_%03i_FEB%i_Plot%i"%(runNum, subrunNum, iFEB, iPlot), 
                                           "Local pedestal estimation; Time [ns]; Pedestal[ADC]")
            fLegend = ROOT.TLegend(0.80, 0.7, 0.95, 0.9)
            for i in range(4):
                iCh = iPlot*4+i
                graph_pedHistory[iFEB][iCh].SetMarkerStyle(constants.rootmarkers[i])
                graph_pedHistory[iFEB][iCh].SetMarkerColor(constants.rootcolors[i])
                graph_pedHistory[iFEB][iCh].SetLineColor(constants.rootcolors[i])
                graph_pedHistory[iFEB][iCh].SetMarkerSize(1)
                fMultiGraph.Add(graph_pedHistory[iFEB][iCh])
                fLegend.AddEntry(graph_pedHistory[iFEB][iCh], "FEB%i Ch%02i"%(iFEB,iCh))
    
            # fMultiGraph.GetHistogram().SetMaximum(50.)
            # fMultiGraph.GetHistogram().SetMinimum(-50.)  
            fMultiGraph.Draw("ALP")
            fLegend.Draw("SAME")
            
            if ( (runNum==1361 and subrunNum==0)  or  ((runNum!=1361 or subrunNum!=0) and iPlot == 0) ):
                fC.Print(pdfname+("(" if pdfStart else ""), "Title: Pedestal History FEB%i Ch%02i-%02i"%(iFEB,iPlot*4,iPlot*4+3))
                pdfStart = False
    
    # fit pedestals
    pedStock, _, funcPedStock = crv_event.CalculatePedestalStock(nFEB, hist_pedStock)
    pedStockRebin0, _, funcPedStockRebin0 = crv_event.CalculatePedestalStock(nFEB, hist_pedStockRebin)
    pedStockRebin1, _, funcPedStockRebin1 = crv_event.CalculatePedestalStock(nFEB, hist_pedStockRebin, 1)
    pedAllPts, _, funcPedAllPts = crv_event.CalculatePedestalStock(nFEB, hist_pedAllPts, 3)
    pedAllPtsExceptP5ADC, _, funcPedAllPtsExceptP5ADC = crv_event.CalculatePedestalStock(nFEB, hist_pedAllPtsExceptP5ADC, 3)
    
    ped_dict = {"stock": pedStock,
                "stockRebin0": pedStockRebin0,
                "stockRebin1": pedStockRebin1,
                "allPts": pedAllPts,
                "allPtsExceptP5ADC": pedAllPtsExceptP5ADC}
    json.dump(ped_dict,fout)
      
    # draw pedestal comparisons
    hist_diff = [None]*4
    for i in range(4):
        hist_diff[i] = ROOT.TH1F("hist_diff_%i_run_%04i_%03i"%(i, runNum, subrunNum), 
                                 "Pedestal difference compared to the stock method; #Delta pedestal; count", 
                                 150, -0.1, 0.05)
    for iFEB in range(nFEB):
        for iCh in range(geometry_constants.nChannelPerFEB):
            fC.Clear()
            fC.SetLogy()
            
            fStack = ROOT.THStack("tHStack_run_%04i_%03i_FEB%i_Ch%i"%(runNum, subrunNum, iFEB, iCh),
                                  "Pedestal comparison Run_%04i_%03i_FEB%i_Ch%i; ADC; counts"%(runNum, subrunNum, iFEB, iCh))
            fLegend = ROOT.TLegend(0.6, 0.5, 0.95, 0.9)
            
            if pedStock[iFEB][iCh]:
                hist_pedStock[iFEB][iCh].SetLineColor(constants.rootcolors[0])
                hist_pedStock[iFEB][iCh].SetLineWidth(1)
                funcPedStock[iFEB][iCh].SetLineColor(constants.rootcolors[0])
                fStack.Add(hist_pedStock[iFEB][iCh])
                fLegend.AddEntry(hist_pedStock[iFEB][iCh], "#splitline{stock method}{#mu = %.3lf, #sigma = %.3lf}"%(funcPedStock[iFEB][iCh].GetParameter(1), funcPedStock[iFEB][iCh].GetParameter(2)))
            if pedStockRebin0[iFEB][iCh]:
                hist_pedStockRebin[iFEB][iCh].SetLineColor(constants.rootcolors[1])
                hist_pedStockRebin[iFEB][iCh].SetLineWidth(1)
                funcPedStockRebin0[iFEB][iCh].SetLineColor(constants.rootcolors[1])
                funcPedStockRebin1[iFEB][iCh].SetLineColor(constants.rootcolors[4])
                fStack.Add(hist_pedStockRebin[iFEB][iCh])
                fLegend.AddEntry(hist_pedStockRebin[iFEB][iCh], 
                                 "#splitline{stock method rebinned}{#splitline{#mu = %.3lf, #sigma = %.3lf (range #pm 4)}{#mu = %.3lf, #sigma = %.3lf (range #pm 1)}}"%
                                 (funcPedStockRebin0[iFEB][iCh].GetParameter(1), funcPedStockRebin0[iFEB][iCh].GetParameter(2), funcPedStockRebin1[iFEB][iCh].GetParameter(1), funcPedStockRebin1[iFEB][iCh].GetParameter(2)))
            if pedAllPts[iFEB][iCh]:
                hist_pedAllPts[iFEB][iCh].SetLineColor(constants.rootcolors[2])
                hist_pedAllPts[iFEB][iCh].SetLineWidth(1)
                funcPedAllPts[iFEB][iCh].SetLineColor(constants.rootcolors[2])
                fStack.Add(hist_pedAllPts[iFEB][iCh])
                fLegend.AddEntry(hist_pedAllPts[iFEB][iCh], "#splitline{all points}{#mu = %.3lf, #sigma = %.3lf}"%(funcPedAllPts[iFEB][iCh].GetParameter(1), funcPedAllPts[iFEB][iCh].GetParameter(2)))
            if pedAllPtsExceptP5ADC[iFEB][iCh]:
                hist_pedAllPtsExceptP5ADC[iFEB][iCh].SetLineColor(constants.rootcolors[3])
                hist_pedAllPtsExceptP5ADC[iFEB][iCh].SetLineWidth(1)
                funcPedAllPtsExceptP5ADC[iFEB][iCh].SetLineColor(constants.rootcolors[3])
                fStack.Add(hist_pedAllPtsExceptP5ADC[iFEB][iCh])
                fLegend.AddEntry(hist_pedAllPtsExceptP5ADC[iFEB][iCh], "#splitline{all pts. (peak-removed)}{#mu = %.3lf, #sigma = %.3lf}"%(funcPedAllPtsExceptP5ADC[iFEB][iCh].GetParameter(1), funcPedAllPtsExceptP5ADC[iFEB][iCh].GetParameter(2)))

            gStyle.SetOptStat(0)
            gROOT.ForceStyle()
            fStack.Draw("NOSTACK")
            if pedStock[iFEB][iCh]:
                fStack.GetXaxis().SetLimits(pedStock[iFEB][iCh]-10., pedStock[iFEB][iCh]+20.)
            
            if pedStock[iFEB][iCh]:
                funcPedStock[iFEB][iCh].Draw("SAME")
            if pedStockRebin0[iFEB][iCh]:
                funcPedStockRebin0[iFEB][iCh].Draw("SAME")
            if pedStockRebin1[iFEB][iCh]:
                funcPedStockRebin1[iFEB][iCh].Draw("SAME")
            if pedAllPts[iFEB][iCh]:
                funcPedAllPts[iFEB][iCh].Draw("SAME")
            if pedAllPtsExceptP5ADC[iFEB][iCh]:
                funcPedAllPtsExceptP5ADC[iFEB][iCh].Draw("SAME")

            fLegend.Draw("SAME")
            fC.Print(pdfname, "Title: Pedestal Comparison FEB%i Ch%02i"%(iFEB,iCh))

            if pedStock[iFEB][iCh]:
                if pedStockRebin0[iFEB][iCh]:
                    hist_diff[0].Fill(pedStockRebin0[iFEB][iCh]-pedStock[iFEB][iCh])
                if pedStockRebin1[iFEB][iCh]:
                    hist_diff[1].Fill(pedStockRebin1[iFEB][iCh]-pedStock[iFEB][iCh])
                if pedAllPts[iFEB][iCh]:
                    hist_diff[2].Fill(pedAllPts[iFEB][iCh]-pedStock[iFEB][iCh])
                if pedAllPtsExceptP5ADC[iFEB][iCh]:
                    hist_diff[3].Fill(pedAllPtsExceptP5ADC[iFEB][iCh]-pedStock[iFEB][iCh])
            
    # difference summary
    fC.Clear()
    fC.SetLogy(0)

    pedestals = np.array([pedStock, pedStockRebin0, pedAllPts, pedAllPtsExceptP5ADC, pedStockRebin1])
    channels = np.array([ii for ii in range(nFEB*geometry_constants.nChannelPerFEB)])
    tags = ["stock method", "stock rebinned #pm 4", "all points", "all pts. (peak-removed)", "stock rebinned #pm 1"]
    fMultiGraph = ROOT.TMultiGraph("tPedByCh_run_%04i_%03i_FEB%i_Plot%i"%(runNum, subrunNum, iFEB, iPlot), 
                                   "Pedestals; Channel; Pedestal[ADC]")
    TG = [None]*5
    fLegend = ROOT.TLegend(0.80, 0.7, 0.95, 0.9)
    for i in range(5):
        tped = pedestals[i].flatten()
        tfilter = tped!=None
        TG[i] = ROOT.TGraph(tped[tfilter].size, array('d',channels[tfilter].tolist()), array('d',tped[tfilter].tolist()))
        TG[i].SetName("tPedByCh%i_run_%04i_%03i_FEB%i_Plot%i"%(i, runNum, subrunNum, iFEB, iPlot))
        TG[i].SetMarkerStyle(constants.rootmarkers[i])
        TG[i].SetMarkerColor(constants.rootcolors[i])
        TG[i].SetLineColor(constants.rootcolors[i])
        TG[i].SetMarkerSize(1)
        fMultiGraph.Add(TG[i])
        fLegend.AddEntry(TG[i], tags[i])
    fMultiGraph.GetXaxis().SetLimits(-0.5,150.)
    fMultiGraph.Draw("ALP")
    fLegend.Draw("SAME")
    fC.Print(pdfname, "Title: Pedestals by Channel")
    
    fC.Clear()
    fC.SetLogy(0)
    fMultiGraph = ROOT.TMultiGraph("tPedByChDiff_run_%04i_%03i_FEB%i_Plot%i"%(runNum, subrunNum, iFEB, iPlot), 
                                   "Pedestals-Stock Value; Channel; #Delta Pedestal[ADC]")
    TG = [None]*4
    fLegend = ROOT.TLegend(0.80, 0.7, 0.95, 0.9)
    for i in range(4):
        tped0 = pedestals[0].flatten()
        tped = pedestals[i+1].flatten()
        tfilter = (tped0!=None) & (tped!=None)
        TG[i] = ROOT.TGraph(tped[tfilter].size, array('d',channels[tfilter].tolist()), array('d',(tped[tfilter]-tped0[tfilter]).tolist()))
        TG[i].SetName("tPedByChDiff%i_run_%04i_%03i_FEB%i_Plot%i"%(i+1, runNum, subrunNum, iFEB, iPlot))
        TG[i].SetMarkerStyle(constants.rootmarkers[i+1])
        TG[i].SetMarkerColor(constants.rootcolors[i+1])
        TG[i].SetLineColor(constants.rootcolors[i+1])
        TG[i].SetMarkerSize(1)
        fMultiGraph.Add(TG[i])
        fLegend.AddEntry(TG[i], tags[i+1])
    fMultiGraph.GetXaxis().SetLimits(-0.5,150.)
    fMultiGraph.Draw("ALP")
    fLegend.Draw("SAME")
    fl = ROOT.TLine(-0.5, 0., 128., 0.)
    fl.SetLineColor(1)
    fl.SetLineWidth(1)
    fl.SetLineStyle(2)
    fl.Draw("SAME")
    fC.Print(pdfname, "Title: Pedestals Compared to Stock Value")
    
    fC.Clear()
    fStack = ROOT.THStack("tHStack_diff_run_%04i_%03i"%(runNum, subrunNum),
                          "Pedestal difference compared to the stock method; #Delta pedestal; count")
    fLegend = ROOT.TLegend(0.7, 0.55, 0.95, 0.9)    
    hist_diff[0].SetLineColor(constants.rootcolors[1])
    hist_diff[0].SetLineWidth(2)
    fStack.Add(hist_diff[0])
    fLegend.AddEntry(hist_diff[0], "#splitline{#splitline{stock method rebinned}{fit peak #pm 4}}{#mu = %.3lf, #sigma = %.3lf}"%(hist_diff[0].GetMean(), hist_diff[0].GetRMS()))
    hist_diff[1].SetLineColor(constants.rootcolors[4])
    hist_diff[1].SetLineWidth(2)
    fStack.Add(hist_diff[1])
    fLegend.AddEntry(hist_diff[1], "#splitline{#splitline{stock method rebinned}{fit peak #pm 1}}{#mu = %.3lf, #sigma = %.3lf}"%(hist_diff[1].GetMean(), hist_diff[1].GetRMS()))
    hist_diff[2].SetLineColor(constants.rootcolors[2])
    hist_diff[2].SetLineWidth(2)
    fStack.Add(hist_diff[2])
    fLegend.AddEntry(hist_diff[2], "#splitline{all points}{#mu = %.3lf, #sigma = %.3lf}"%(hist_diff[2].GetMean(), hist_diff[2].GetRMS()))
    hist_diff[3].SetLineColor(constants.rootcolors[3])
    hist_diff[3].SetLineWidth(2)
    fStack.Add(hist_diff[3])
    fLegend.AddEntry(hist_diff[3], "#splitline{all pts. (peak-removed)}{#mu = %.3lf, #sigma = %.3lf}"%(hist_diff[3].GetMean(), hist_diff[3].GetRMS()))
    fStack.Draw("NOSTACK")
    fLegend.Draw("SAME")
    fC.Print(pdfname+")", "Title: Difference between Methods")
                                                                              
    fout.close()
    fFile.Close()
    pdfpages.close()
    
    lock.acquire()
    print("Run%04i_%03i done."%(runNum, subrunNum))
    sys.stdout.flush()
    lock.release()
    
    return

def main():
    nWorker = 5
    p = mp.Pool(processes=nWorker)
    for filename in filelist[10:]:
        p.apply_async(checkPedestal, (filename,))
    p.close()
    p.join()
    
    
if __name__ == "__main__":
    main()
