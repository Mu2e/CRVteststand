import ROOT
from ROOT import TCanvas, TH1F, TFile, TPaveText
from ROOT import gStyle, gROOT, gDirectory, gPad
import numpy as np
import argparse

nFEB = 4 # Just to get a handle for adjusting FEB numbers
nPE_per_layer_cutoff = 20
filename = "rec.mu2e.CRV_wideband_cosmics.crvaging-004.001082_000.root"

parser = argparse.ArgumentParser(description='Optional command arguments',formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-n','--nFEB', type=int, nargs='?', default=nFEB, 
                    help='number of FEBs present in the data')
parser.add_argument('-c','--cutoff', type=int, nargs='?', default=nPE_per_layer_cutoff, 
                    help='single layer trigger PE cutoff')
parser.add_argument('-f','--filename', type=str, nargs='?', default=filename, 
                    help='file name to process (after reco)')
args = parser.parse_args()
nFEB = args.nFEB
nPE_per_layer_cutoff = args.cutoff
filename = args.filename
pdfname = filename.split("/")[-1].split(".root")[0]+".efficiency_demo.pdf"
                    
gROOT.Reset()
gROOT.SetBatch(1)
gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")
gStyle.SetOptStat(111110)
gStyle.SetOptFit(0)
gStyle.SetLineScalePS(0.3)

fFile = TFile(filename, "READ")
fTree = fFile.Get("run") # Tree name is "run"

nEvent = fTree.GetEntries() # Get the number of events
print nEvent, "events present in the file"

nMuon = 0 # Tally the number of events with layer 1 and layer 2 (FEB0) triggered
nTriggered = [0]*5 # Tally the number of layers triggered 

hist_nPE_per_channel = TH1F("hist_nPE_per_channel", "nPE per channel after calibration; n_{PE}; count", 205, -5., 200.)
hist_nPE_per_layer = TH1F("hist_nPE_per_layer", "nPE per layer after calibration; n_{PE}; count", 150, 0., 300.)
hist_n_triggered_layer= TH1F("hist_n_triggered_layer", "n triggered layer; number of triggered layers; count", 7, -0.5, 6.5)

for iEvent in range(nEvent):
    fTree.GetEntry(iEvent)
    
    if iEvent > 0 and (iEvent % 5000) == 0:
        print "  ...", iEvent, "processed"
    
    # Data to analyze are: fTree.PEsTemperatureCorrected[nFEB][64]
    # type float[nFEB][64]
    #
    # Layers:0.2 PEsTemperatureCorrected[0][0:31]
    #        0.3 PEsTemperatureCorrected[0][32:63]
    #
    #        1.0 PEsTemperatureCorrected[1][0:31]
    #        1.1 PEsTemperatureCorrected[1][32:63]
    #        1.2 PEsTemperatureCorrected[2][0:31]
    #        1.3 PEsTemperatureCorrected[2][32:63]
    #
    #        2.0 PEsTemperatureCorrected[3][0:31]
    #        2.1 PEsTemperatureCorrected[3][32:63]
    #
    # We want to see when there are muons going through the modules 
    # (i.e., there is a track going through layers 0.2, 0.3, 2.0, and 2.1),
    # if the central module (with layers 1.x) triggers at >99.99% of the time
    # (satisfying 3/4 layers criteria)
    #
    # For a start, ignoring the hit timing and the track geometry (these are
    # things for you to play with), let's say when there are more than nPE_per_layer_cutoff
    # PEs in a layer, the layer is triggered. Then we can loop through the 
    # layers and determine if each layer is triggered 
    
    PE_array = np.reshape(fTree.PEsTemperatureCorrected, (nFEB,64)) 
    # just some pyroot stuff, go with float[nFEB][64] array if you code in C++
        
    # first select the events with top (and bottom) modules fired
    exist_muon_track = 1
    # list of starting channel of a layer: (FEB#, ch#)
    list_starting_channel = [(0, 0), (0, 32)]
    if nFEB == 4:
        list_starting_channel += [(3, 0), (3, 32)] # if 4 FEB, also add bottom two layers
    
    for item in list_starting_channel:
        nPE_in_layer = 0
        for i in range(32):
            # count # of PEs in layers
            nPE_in_layer += max(0, PE_array[item[0]][item[1]+i])
            # nPE_in_layer += PE_array[item[0]][item[1]+i]
        if nPE_in_layer <= nPE_per_layer_cutoff:
            exist_muon_track *= 0 
            
    if exist_muon_track == 0:
        continue # if some layers not hit in top/bottom layers, maybe not a muon track
    else:
        nMuon += 1
        n_triggered_layer = 0
        for iFEB in range(1,3): # iFEB values are 1 and 2
            for iLayer in range(2):
                nPE_this_layer = 0
                for i in range(32): # loop around channels in this layer
                    nPE_this_channel = PE_array[iFEB][iLayer*32+i]
                    hist_nPE_per_channel.Fill(nPE_this_channel)
                    
                    if nPE_this_channel < 0:
                        nPE_this_channel = 0
                    # I got a sense that the negative PEs are just place holders. 
                    # Definitely include these if that is confirmed 
                        
                    nPE_this_layer += nPE_this_channel
                
                hist_nPE_per_layer.Fill(nPE_this_layer)
                if nPE_this_layer > nPE_per_layer_cutoff:
                    n_triggered_layer += 1
        
        hist_n_triggered_layer.Fill(n_triggered_layer)
        nTriggered[n_triggered_layer] +=1

print nMuon, "events had cosmic muon tracks"
print nTriggered[3], "events triggered 3 layers"
print nTriggered[4], "events triggered 4 layers"
effic = float(nTriggered[3]+nTriggered[4])/float(nMuon)

fC0 = TCanvas("c0", "c0", 1200, 900)
fC0.SetLogy()
hist_nPE_per_channel.Draw()
fC0.Print(pdfname+"(", "Title:nPE_per_channel")
hist_nPE_per_layer.Draw()
fC0.Print(pdfname, "Title:nPE_per_layer")
fC0.SetLogy(0)
hist_n_triggered_layer.Draw()
pav=TPaveText(.78, .7, .98, .25,"NDC")
pav.AddText("Events present: "+str(nEvent))
pav.AddText("Cosmic muon tracks: "+str(nMuon))
for i in range(5):
    pav.AddText(str(i)+" layers triggered: "+str(nTriggered[i]))
pav.AddText("Efficiency = %.5f"%effic)
pav.Draw("SAME")
fC0.Print(pdfname+")", "Title:n_triggered_layer")

fFile.Close()
    
