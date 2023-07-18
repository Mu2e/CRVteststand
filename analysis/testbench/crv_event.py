from __future__ import print_function

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import collections as mc
import ROOT
from ROOT import TCanvas, TH1F, TH2F, TF1, TMath, TGraph, TFile, TSpectrum, TPaveText, TMultiGraph, TGraphErrors, TLine
from ROOT import gStyle, gROOT, gDirectory, gPad
import numpy as np
from array import array
import os, sys, json, glob, datetime 
import math
import argparse, textwrap
import constants, geometry_constants
import root_utils
from root_utils import treeInitialization
import pprint

nPreSignalRegions = [[0]*geometry_constants.nChannelPerFEB for i in range(10)] # container
nNoiseHits = [[0]*geometry_constants.nChannelPerFEB for i in range(10)] # container

class crv_event:
    def __init__(self, runtree, iEntry, detailLevel = 0b1010, nFEB = 4, isRaw = False):

        # initializer to load CRV_event from ROOT tree.
        # input: runtree: ROOT TTree object that corresponds to "run" tree in the crvparsed or crvreco ROOT files.
        #        iEntry: index of the entry to load.
        #        detailLevel (default 0b1010): controls how much information to load: 
        #                     0b1000: basic event info 
        #                     0b0100: expanded info, with adc traces
        #                     0b0010: main peak fitting info
        #                     0b0001: reflected peak fitting info
        #        nFEB (default 4): number of FEBs present
        #        isRaw: indicates if the input file is among the crvparsed without fitting.

        # in the case of reading raw file, force detailLevel to "basic"+trace
        if isRaw:
            detailLevel = 0b1100
        # always load at least the basic event info
        if (detailLevel & 0b1111)==0b0: 
            detailLevel = 0b1000
            print("*** WARNING: CRV_event: Attribute 'detailLevel' unspecified. Changed to '0b1000' to load only basic event information.")

        try:
            runtree.GetEntry(iEntry)
        except Exception as e:
            print("Error: CRV_event: reading entry #%i unsuccessful!"%(iEntry))
            print("                ", type(e), e)
            return

        # prepare the containers.
        attributeList = [# tree stuff
                         "spillNumber", "eventNumber", "temperature", "pedestal", "tdcSinceSpill", "timeSinceSpill", 
                         "fitStatus", "fitStatusReflectedPulse", "adc", 
                         "PEs", "PEsTemperatureCorrected", "pulseHeight", 
                         "beta", "time", "LEtime", 
                         "recoStartBin", "recoEndBin",
                         "PEsReflectedPulse", "PEsTemperatureCorrectedReflectedPulse", "pulseHeightReflectedPulse", 
                         "betaReflectedPulse", "timeReflectedPulse", "LEtimeReflectedPulse", 
                         "recoStartBinReflectedPulse", "recoEndBinReflectedPulse",
                         # operation containers
                         #"crvSpill", # container to hold a linked CRV_spill object
                         "counterTriggered", # matrix containing counter triggering info
                         "layerTriggered", # boolean list indicating layer triggering
                         "trackExist" # boolean object t
                         ]
        for item in attributeList:
            exec("self."+item+"=None")

        self.isRaw = isRaw 
        self.detailLevel = detailLevel
        self.nFEB = nFEB
        self.temp_dict = {}
        self.iEntry = iEntry
        self.nSample = int(np.size(runtree.runtree_adc if isRaw else runtree.adc)/nFEB/(geometry_constants.nChannelPerFEB))
            
        # start filling in the contents using the trees 
        # FIXME: v3 also has spillIndex/I (runtree_spill_index/I), 
        #                    spillTimestamp/G (runtree_spillTimestamp/G), 
        #                    boardStatus[6][22]/I (runtree_boardStatus[6][22]/I) 
        if (detailLevel & 0b1000):
            self.spillNumber = runtree.runtree_spill_num if isRaw else runtree.spillNumber #runtree_spill_num/I or spillNumber/I
            self.eventNumber = runtree.runtree_event_num if isRaw else runtree.eventNumber #runtree_event_num/I or eventNumber/I
            self.temperature = np.reshape(np.array(runtree.runtree_temperature if isRaw else runtree.temperature, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #runtree_temperature[4][64]/L or temperature[4][64]/F
            self.pedestal = None if isRaw else np.reshape(np.array(runtree.pedestal, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #pedestal[4][64]/F
            self.tdcSinceSpill = np.reshape(np.array(runtree.runtree_tdc_since_spill if isRaw else runtree.tdcSinceSpill, dtype=np.int64), (nFEB, geometry_constants.nChannelPerFEB)) #runtree_tdc_since_spill[4][64]/L or tdcSinceSpill[4][64]/L
            self.timeSinceSpill = np.reshape(np.array(runtree.runtree_time_since_spill if isRaw else runtree.timeSinceSpill, dtype=np.float64), (nFEB, geometry_constants.nChannelPerFEB)) #runtree_time_since_spill[4][64]/D or timeSinceSpill[4][64]/D
            self.fitStatus = None if isRaw else np.reshape(np.array(runtree.fitStatus, dtype=np.int32), (nFEB, geometry_constants.nChannelPerFEB)) #fitStatus[4][64]/I
            self.fitStatusReflectedPulse = None if isRaw else np.reshape(np.array(runtree.fitStatusReflectedPulse, dtype=np.int32), (nFEB, geometry_constants.nChannelPerFEB)) #fitStatusReflectedPulse[4][64]/I
            # FIXME: global timestamp (from Spill time + time since spill)??? 
        if (detailLevel & 0b0100):
            self.adc = np.reshape(np.array(runtree.runtree_adc if isRaw else runtree.adc, dtype=np.int16), (nFEB, geometry_constants.nChannelPerFEB, -1)) #adc[4][64][127]/S
        if (detailLevel & 0b0010):
            self.PEs = np.reshape(np.array(runtree.PEs, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #PEs[4][64]/F
            self.PEsTemperatureCorrected = np.reshape(np.array(runtree.PEsTemperatureCorrected, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #PEsTemperatureCorrected[4][64]/F
            self.pulseHeight = np.reshape(np.array(runtree.pulseHeight, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #pulseHeight[4][64]/F
            self.beta = np.reshape(np.array(runtree.beta, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #beta[4][64]/F
            self.time = np.reshape(np.array(runtree.time, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #time[4][64]/F
            self.LEtime = np.reshape(np.array(runtree.LEtime, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #LEtime[4][64]/F
            self.recoStartBin = np.reshape(np.array(runtree.recoStartBin, dtype=np.int32), (nFEB, geometry_constants.nChannelPerFEB)) #recoStartBin[4][64]/I
            self.recoEndBin = np.reshape(np.array(runtree.recoEndBin, dtype=np.int32), (nFEB, geometry_constants.nChannelPerFEB)) #recoEndBin[4][64]/I
        if (detailLevel & 0b0001):
            self.PEsReflectedPulse = np.reshape(np.array(runtree.PEsReflectedPulse, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #PEsReflectedPulse[4][64]/F
            self.PEsTemperatureCorrectedReflectedPulse = np.reshape(np.array(runtree.PEsTemperatureCorrectedReflectedPulse, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #PEsTemperatureCorrectedReflectedPulse[4][64]/F
            self.pulseHeightReflectedPulse = np.reshape(np.array(runtree.pulseHeightReflectedPulse, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #pulseHeightReflectedPulse[4][64]/F
            self.betaReflectedPulse = np.reshape(np.array(runtree.betaReflectedPulse, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #betaReflectedPulse[4][64]/F
            self.timeReflectedPulse = np.reshape(np.array(runtree.timeReflectedPulse, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #timeReflectedPulse[4][64]/F
            self.LEtimeReflectedPulse = np.reshape(np.array(runtree.LEtimeReflectedPulse, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB)) #LEtimeReflectedPulse[4][64]/F
            self.recoStartBinReflectedPulse = np.reshape(np.array(runtree.recoStartBinReflectedPulse, dtype=np.int32), (nFEB, geometry_constants.nChannelPerFEB)) #recoStartBinReflectedPulse[4][64]/I
            self.recoEndBinReflectedPulse = np.reshape(np.array(runtree.recoEndBinReflectedPulse, dtype=np.int32), (nFEB, geometry_constants.nChannelPerFEB)) #recoEndBinReflectedPulse[4][64]/I
    
    #FIXME: 
    # def linkCRVSpill(spilltree, ):
        # self.crvSpill = CRV_Spill()

    #FIXME:
    # def isTrack(optional draw)

    def geomAttributeFill(self, benchGeometryList, attribute):
        tAttribute = None
        tAttribute = eval("self.%s"%attribute)
        # print (tAttribute)
        if tAttribute is None:
            sys.exit("ERROR: CRV_event: %s is not loaded"%(attribute))
        geomAttribute = []
        for iGeometry in benchGeometryList:
            tGeomAttribute = np.zeros_like(iGeometry.mappingFEBCh, dtype=float)
            with np.nditer(iGeometry.mappingFEBCh, flags=['multi_index']) as it:
                for x in it:
                    tFEB = int(x.real)
                    tCh = int(x.imag)
                    if tFEB < 0:
                        pass
                    else:
                        tGeomAttribute[it.multi_index] = tAttribute[tFEB][tCh]
            geomAttribute.append(tGeomAttribute)
        return geomAttribute

    def triggerCounter(self, benchGeometryList, counterThreshold = 8):
        if "geomPECorrected" not in self.temp_dict.keys():
            geomPECorrected = self.geomAttributeFill(benchGeometryList, "PEsTemperatureCorrected")
            self.temp_dict.update({"geomPECorrected":geomPECorrected})
        triggeredCounter = []
        for iGeomPEC in self.temp_dict["geomPECorrected"]:
            triggeredCounter.append((iGeomPEC[:, :, 1::2]+iGeomPEC[:, :, ::2])>counterThreshold)
        self.temp_dict.update({"triggeredCounter":triggeredCounter})
        return # list of triggered counter mask, boolean type 

    def triggerLayer(self, benchGeometryList, layerThreshold = 60):
        if "geomPECorrected" not in self.temp_dict.keys():
            geomPECorrected = self.geomAttributeFill(benchGeometryList, "PEsTemperatureCorrected")
            self.temp_dict.update({"geomPECorrected":geomPECorrected})
        triggeredLayer = []
        for iGeomPEC in self.temp_dict["geomPECorrected"]:
            triggeredLayer.append(iGeomPEC.sum(axis=2)>layerThreshold)
        self.temp_dict.update({"triggeredLayer":triggeredLayer})
        return # list of triggered layer mask, boolean type 

    def plotTrace(self, ax, FEB, Channel, AdditionalTitle): # Channel can be a range or a number 
        if self.adc is None:
            sys.exit("ERROR: CRV_event: adc is not loaded")
        strChannel = Channel
        if isinstance(strChannel, int):
            strChannel = "%i:%i"%(strChannel,strChannel+1)

        ax.set_title("FEB %i Ch %s "%(FEB, Channel)+AdditionalTitle)
        ax.set_xlabel('ns')
        ax.set_ylabel('ADC counts')
        x=[constants.sampleRate*i for i in range(np.size(self.adc, 2))]
        label = ""        
        label_list = ["'ch %i', "%(i) for i in range(geometry_constants.nChannelPerFEB)]
        exec("label_list = label_list["+strChannel+"]")
        for iLabel in label_list:
            label += iLabel
        label=label[:-2]
        exec("ax.legend(ax.plot(x, self.adc[%i"%FEB+","+strChannel+"].T), ["+label+"])")
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        return
        
    def plotEvent(self, pdfpages, benchGeometryList, ifPlotTriggeredCounter = True, ifPlotTrace = False):
        if "geomPECorrected" not in self.temp_dict.keys():
            geomPECorrected = self.geomAttributeFill(benchGeometryList, "PEsTemperatureCorrected")
            self.temp_dict.update({"geomPECorrected":geomPECorrected})
        if "geomLEtime" not in self.temp_dict.keys():
            geomLEtime = self.geomAttributeFill(benchGeometryList, "LEtime")
            self.temp_dict.update({"geomLEtime":geomLEtime})
        if ifPlotTriggeredCounter:
            if "triggeredCounter" not in self.temp_dict.keys():
                self.triggerCounter(benchGeometryList) # default counterThreshold = 8
        
        canvasSize = benchGeometryList[0].canvasSize
        nX, nY = 0, 0
        if canvasSize[0]>canvasSize[1]:
            nX = len(benchGeometryList)
            nY = 1
        else:
            nX = 1
            nY = len(benchGeometryList)
        fig, axes = plt.subplots(nX,nY,figsize=(nY*canvasSize[0], nX*canvasSize[1]))
        if nX*nY == 1:
            axes = [axes]
        
        for index, iGeometry in enumerate(benchGeometryList):
            iGeometry.plotGeometry(axes[index], 'Event #%i'%self.iEntry, False, True)
            plotArea = (np.ma.masked_where((iGeometry.idleMask|iGeometry.badChMask), self.temp_dict["geomPECorrected"][index])).flatten()
            plotColor = self.temp_dict["geomLEtime"][index].flatten()
            scatTriggeredCh = axes[index].scatter(iGeometry.coordinates["channelX"].flatten(), 
                                                  iGeometry.coordinates["channelY"].flatten(),
                                                  s=plotArea*0.5, marker='o', c=plotColor, 
                                                  #cmap='viridis', 
                                                  cmap='inferno', 
                                                  vmin = 0., vmax = constants.sampleRate*self.nSample,
                                                  zorder = 10)
            axes[index].legend([scatTriggeredCh],['PECorrected'])
            cbar = plt.colorbar(scatTriggeredCh)
            #cbar.ax.get_yaxis().set_ticks([])
            cbar.ax.get_yaxis().labelpad = 15
            cbar.ax.set_ylabel('time [ns]', rotation=270)
            if ifPlotTriggeredCounter:
                linesCoord_triggered = iGeometry.plotCntrCoordGen(self.temp_dict["triggeredCounter"][index])
                if linesCoord_triggered:
                    lc_triggered = mc.LineCollection(linesCoord_triggered, linewidths=1.5, linestyles='solid', 
                                                     colors='#113285',alpha=1.,zorder=6)
                    axes[index].add_collection(lc_triggered)
        pdfpages.savefig(fig)
        plt.close(fig)
        
        #FIXME: add straight line fitting for tracks

        if ifPlotTrace:
            for iFEB in range(self.nFEB):
                fig, axes = plt.subplots(8, 8, figsize=(8*6, 8*4))
                for iCh in range(geometry_constants.nChannelPerFEB):
                    self.plotTrace(axes[int(iCh/8), iCh%8], iFEB, iCh, '')
                plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
                pdfpages.savefig(fig)
                plt.close(fig)
                    
    #FIXME:
    # get nFEB from geometry_constants setup_dict / shape of TTree
    
    def FillPedestalHistogramsStock(self, pedestalHist):
    # adaption of calibCrv::CrvEvent::FillCalibrationHistograms for testing purposes
        if self.adc is None:
            sys.exit("ERROR: CRV_event: adc is not loaded")
        for iFEB in range(self.nFEB):
            for iCh in range(geometry_constants.nChannelPerFEB):
                if np.isnan(self.timeSinceSpill[iFEB][iCh]):
                # missing FEB/channel in raw data
                    continue
                data = self.adc[iFEB][iCh]
                if data[0]==0 and data[1]==0 and data[3]==0:
                # FIXME temporary check for bad events
                # where other channels work, so that timSinceSpill wasn't marked as NAN
                    continue
                # divide prespill region into three parts
                # using only 1 region (i.e. don't devide the prespill region) returns a pedestal which is too high
                numberOfRegions=3
                for i in range(numberOfRegions):
                    average = 0
                    minADC = float('nan')
                    maxADC = float('nan')
                    for j in range(int(constants.numberOfPreSignalSamples/numberOfRegions)):
                        tadc = data[i*int(constants.numberOfPreSignalSamples/numberOfRegions)+j]
                        average += tadc
                        if tadc < minADC or math.isnan(minADC):
                            minADC = tadc
                        if tadc > maxADC or math.isnan(maxADC):
                            maxADC = tadc
                    if maxADC-minADC > 5:
                        continue
                    average/=float(int(constants.numberOfPreSignalSamples/numberOfRegions))
                    pedestalHist[iFEB][iCh].Fill(average)
        return
    
    def FillPedestalHistogramsAllPts(self, pedestalHist):
        if self.adc is None:
            sys.exit("ERROR: CRV_event: adc is not loaded")
        for iFEB in range(self.nFEB):
            for iCh in range(geometry_constants.nChannelPerFEB):
                if np.isnan(self.timeSinceSpill[iFEB][iCh]):
                # missing FEB/channel in raw data
                    continue
                data = self.adc[iFEB][iCh]
                if data[0]==0 and data[1]==0 and data[3]==0:
                # FIXME temporary check for bad events
                # where other channels work, so that timSinceSpill wasn't marked as NAN
                    continue
                for i in range(constants.numberOfPreSignalSamples):
                    pedestalHist[iFEB][iCh].Fill(data[i])
        return
    
    def FillPedestalHistogramsAllPtsExceptP5ADC(self, pedestalHist):
        local_ped = [[None]*geometry_constants.nChannelPerFEB for i in range(self.nFEB)]
        
        if self.adc is None:
            sys.exit("ERROR: CRV_event: adc is not loaded")
        for iFEB in range(self.nFEB):
            for iCh in range(geometry_constants.nChannelPerFEB):
                if np.isnan(self.timeSinceSpill[iFEB][iCh]):
                # missing FEB/channel in raw data
                    continue
                data = self.adc[iFEB][iCh]
                if data[0]==0 and data[1]==0 and data[3]==0:
                # FIXME temporary check for bad events
                # where other channels work, so that timSinceSpill wasn't marked as NAN
                    continue
                # divide prespill region into three parts
                # using only 1 region (i.e. don't devide the prespill region) returns a pedestal which is too high
                
                # first do a "local estimation"
                numberOfRegions=3
                crude_pedestal = []
                for i in range(numberOfRegions):
                    average = 0
                    minADC = float('nan')
                    maxADC = float('nan')
                    for j in range(int(constants.numberOfPreSignalSamples/numberOfRegions)):
                        tadc = data[i*int(constants.numberOfPreSignalSamples/numberOfRegions)+j]
                        average += tadc
                        if tadc < minADC or math.isnan(minADC):
                            minADC = tadc
                        if tadc > maxADC or math.isnan(maxADC):
                            maxADC = tadc
                    if maxADC-minADC > 5:
                        continue
                    average/=float(int(constants.numberOfPreSignalSamples/numberOfRegions))
                    crude_pedestal.append(average)
                
                if len(crude_pedestal) == 0:
                    continue # too noisy local ped. ignore this sector
                crude_pedestal = float(sum(crude_pedestal))/float(len(crude_pedestal))
                local_ped[iFEB][iCh] = crude_pedestal
                
                pedfilter = [True]*constants.numberOfPreSignalSamples
                for i in range(constants.numberOfPreSignalSamples):
                    if data[i] > crude_pedestal+5:
                        for j in range(max(0,i-4), min(constants.numberOfPreSignalSamples,i+11)):
                            pedfilter[j] = False
                    
                for i in range(constants.numberOfPreSignalSamples):
                    if pedfilter[i]:
                        pedestalHist[iFEB][iCh].Fill(data[i])
        return local_ped

    def FillCalibrationHistogramsStock(self, all_pedestals, rawCalibHist, temperatureCorrectedCalibHist = None):
    # adaption of calibCrv::CrvEvent::FillCalibrationHistograms for testing purposes       
        global nPreSignalRegions, nNoiseHits
        if self.adc is None:
            sys.exit("ERROR: CRV_event: adc is not loaded")
            
        for iFEB in range(self.nFEB):
            for iCh in range(geometry_constants.nChannelPerFEB):
                if np.isnan(self.timeSinceSpill[iFEB][iCh]):
                # missing FEB/channel in raw data
                    continue
                data = self.adc[iFEB][iCh]
                if data[0]==0 and data[1]==0 and data[3]==0:
                # FIXME temporary check for bad events
                # where other channels work, so that timSinceSpill wasn't marked as NAN
                    continue
                pedestal = all_pedestals[iFEB][iCh]
                
                # in the stock codes, example pulses were drawn, omitted here.
                
                nBins = constants.numberOfPreSignalSamples
                # remove the pedestal and find the local maxima
                waveform = []
                peaks = []
                # sum=0
                for ibin in range(nBins):
                    waveform.append(data[ibin]-pedestal)
                    if ibin>1 and ibin<nBins-3: #don't search for peaks too close to the sample start or end
                        if (data[ibin-1]<data[ibin] and data[ibin]>data[ibin+1] and data[ibin]-pedestal>constants.noiseThreshold):
                            peaks.append((ibin,False))
                        if (data[ibin-1]<data[ibin] and data[ibin]==data[ibin+1] and data[ibin+1]>data[ibin+2] and data[ibin]-pedestal>constants.noiseThreshold):
                            peaks.append((ibin,True))
                    # sum += math.fabs(data[ibin]-pedestal)
                # don't use noisy events  # doesn't seem to do much
                # if (float(sum)/float(nBins)>1.0):
                    # return  //FIXME
                    
                nPreSignalRegions[iFEB][iCh] += 1
                nNoiseHits[iFEB][iCh] += len(peaks)
                
                # fit all peaks
                for i in range(len(peaks)):
                    # select a range of up to 4 points before and after the maximum point
                    # -find up to 5 points before and after the maximum point for which the waveform is stricly decreasing
                    # -remove 1 point on each side. this removes potentially "bad points" belonging to a second pulse (i.e. in double pulses) 
                    maxBin = peaks[i][0]
                    startBin = maxBin
                    endBin = maxBin
                    for ibin in range(maxBin-1, max(-1, maxBin-5-1), -1):
                        if waveform[ibin]<=wavform[ibin+1]:
                            startBin = ibin
                        else:
                            break
                    for ibin in range(maxBin+1, min(nBins, maxBin+5+1)):
                        if waveform[ibin]<=wavform[ibin-1]:
                            endBin = ibin
                        else:
                            break
                    if maxBin - startBin > 1:
                        startBin += 1
                    if endBin - maxBin > 1:
                        endBin -= 1
                        
                    # fill the graph
                    binWidth = 1.0/constants.RATE
                    t = [ibin*binWidth for ibin in range(startBin, endBin+1)]
                    v = waveform[startBin:(endBin+1)]
                    g = ROOT.TGraph(len(t), array('d', t), array('d', v))
                    
                    # set the fit function
                    f = ROOT.TF1("peakfinder","[0]*(TMath::Exp(-(x-[1])/[2]-TMath::Exp(-(x-[1])/[2])))")
                    f.SetParameter(0, waveform[maxBin]*ROOT.TMath.E())
                    f.SetParameter(1, (maxBin+0.5)*binWidth if peaks[i][1] else maxBin*binWidth)
                    f.SetParameter(2, 12.6)
                    
                    # do the fit
                    fr = g.Fit(f,"NQS")

                    if not fr.IsValid():
                        continue
                    if fr.Parameter(0)<=0 or fr.Parameter(2)<4:
                        continue # probably misreconstructed
                    if fr.Parameter(2)>25:
                        continue # probably not a noise hit
                    if math.fabs(fr.Parameter(1)-maxBin*binWidth)>30:
                        continue
                    if fr.Parameter(0)/(waveform[maxBin]*ROOT.TMath.E())>1.5:
                        continue
                    
                    pulseArea = fr.Parameter(0)*fr.Parameter(2)
                    rawCalibHist[iFEB][iCh].Fill(pulseArea)
                    
                    # temperature correction of noise pulse area
                    if temperatureCorrectedCalibHist is not None:
                        if self.temperature[iFEB][iCh] != 0:
                            pulseArea *= (constants.calibTemperatureIntercept-constants.referenceTemperature)/(constants.calibTemperatureIntercept-self.temperature[iFEB][iCh])
                            temperatureCorrectedCalibHist[iFEB][iCh].Fill(pulseArea)

                    # plotting of example 1PE / 2PE plots, temperature plot omitted
        return
    
    def peakFitter(self, all_pedestals, calibrationFactor, calibrationFactorTemperatureCorrected, fitReflected = True):
    # adaption of recoCrv::CrvRecoEvent::PeakFitter
        
        if self.adc is None:
            sys.exit("ERROR: CRV_event: adc is not loaded")
            
        # prepare containers. All initialized to nan
        self.pedestal = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
        
        self.fitStatus = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
        self.recoStartBin = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
        self.recoEndBin = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
        self.PEs = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
        self.pulseHeight = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
        self.time = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
        self.beta = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
        self.LEtime = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
        self.PEsTemperatureCorrected = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
        
        if fitReflected:
            self.fitStatusReflectedPulse = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
            self.recoStartBinReflectedPulse = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
            self.recoEndBinReflectedPulse = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
            self.PEsReflectedPulse = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
            self.pulseHeightReflectedPulse = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
            self.timeReflectedPulse = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
            self.betaReflectedPulse = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
            self.LEtimeReflectedPulse = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
            self.PEsTemperatureCorrectedReflectedPulse = numpy.full((self.nFEB, geometry_constants.nChannelPerFEB), np.nan, dtype=np.float32)
        
        for iFEB in range(self.nFEB):
            for iCh in range(geometry_constants.nChannelPerFEB):
                if np.isnan(self.timeSinceSpill[iFEB][iCh]):
                # missing FEB/channel in raw data
                    continue
                data = self.adc[iFEB][iCh]
                if data[0]==0 and data[1]==0 and data[3]==0:
                # FIXME temporary check for bad events
                # where other channels work, so that timSinceSpill wasn't marked as NAN
                    continue
                pedestal = all_pedestals[iFEB][iCh]
                self.pedestal[iFEB][iCh] = pedestal
                
                # remove the pedestal and find the maxima in the signal region
                waveform = []
                peaks = []
                peakBinsStart=0
                peakBinsEnd=0
                
                for ibin in range(constants.signalRegionStart, min(len(data),constants.signalRegionEnd+1)):
                    waveform.append(data[ibin]-pedestal)
                    if ibin <= constants.signalRegionStart:
                        continue
                    
                    if data[ibin-1]<data[ibin]: # rising edge
                        peakBinsStart=ibin
                        peakBinsEnd=ibin
                        
                    if data[ibin-1]==data[ibin]: #potentially at a peak with consecutive ADC values which are equal
                        peakBinsEnd=ibin
                        
                    if data[ibin-1]>data[ibin]:
                        if peakBinsStart>0: # found a peak
                            if data[peakBinsStart]-pedestal>6: # ignores fluctuations of the baseline
                                peaks.append((data[peakBinsStart]-pedestal, peakBinsStart, peakBinsEnd))
                        peakBinsStart = 0
                        
                # ignore events without pulses
                if len(peaks)==0:
                    self.fitStatus[iFEB][iCh] = 0
                    self.fitStatusReflectedPulse[iFEB][iCh] = 0
                    continue
                
                # only look at the largest peak (and if available the peak after that for a potential reflected pulse)
                iPeak = 0
                maxPeak = 0
                for i in range(len(peaks)):
                    if peaks[i][0]>maxPeak:
                        maxPeak = peaks[i][0]
                        iPeak = i
                
                binWidth=1.0/constants.RATE
                    
                for i in range(2):
                    
                    if fitReflected == False and i == 1:
                        continue
                    
                    if iPeak == len(peaks)-1:
                        self.fitStatusReflectedPulse[iFEB][iCh] = 0
                        continue
                        
                    peakBinsStart = peaks[iPeak+i][1]-constants.signalRegionStart
                    peakBinsEnd   = peaks[iPeak+i][2]-constants.signalRegionStart
                    averagePeakBin = 0.5*(peakBinsStart+peakBinsEnd)
                    
                    #select a range of up to 4 points before and after the peak
                    #-find up to 5 points before and after the peak for which the waveform is stricly decreasing
                    #-remove 1 point on each side. this removes potentially "bad points" belonging to a second pulse (i.e. in double pulses)
                    nBins = len(waveform)
                    tmp_recoStartBin=peakBinsStart-1
                    tmp_recoEndBin=peakBinsEnd+1
                    for ibin in range(peakBinsStart-1, max(-1, peakBinsStart-5-1), -1):
                        if waveform[ibin]<=wavform[ibin+1]:
                            tmp_recoStartBin = ibin
                        else:
                            break
                    for ibin in range(peakBinsEnd+1, min(nBins, peakBinsEnd+5+1)):
                        if waveform[ibin]<=wavform[ibin-1]:
                            tmp_recoEndBin = ibin
                        else:
                            break
                    if peakBinsStart - tmp_recoStartBin > 1:
                        tmp_recoStartBin += 1
                    if tmp_recoEndBin - peakBinsEnd > 1:
                        tmp_recoEndBin -= 1
                     
                    if i == 0:
                        self.recoStartBin[iFEB][iCh] = tmp_recoStartBin
                        self.recoEndBin[iFEB][iCh] = tmp_recoEndBin
                    else:
                        self.recoStartBinReflectedPulse[iFEB][iCh] = tmp_recoStartBin
                        self.recoEndBinReflectedPulse[iFEB][iCh] = tmp_recoEndBin
                        
                    # fill the graph
                    g = ROOT.TGraph()
                    for ibin in range(tmp_recoStartBin, tmp_recoEndBin+1):
                        t = ibin*binWidth
                        v = waveform[ibin]
                        g.SetPoint(g.GetN(), t, v)
                        
                    # set the fit function
                    f = ROOT.TF1("peakfitter","[0]*(TMath::Exp(-(x-[1])/[2]-TMath::Exp(-(x-[1])/[2])))")
                    f.SetParameter(0, waveform[peakBinsStart]*ROOT.TMath.E())
                    f.SetParameter(1, averagePeakBin*binWidth)
                    f.SetParameter(2, constants.DEFAULT_BETA)
                    bound_low = [waveform[peakBinsStart]*ROOT.TMath.E()*0.7, averagePeakBin*binWidth-15.0, 5.0]
                    bound_high = [waveform[peakBinsStart]*ROOT.TMath.E()*1.5, averagePeakBin*binWidth+15.0, 40.0]
                    f.SetParLimits(0, bound_low[0], bound_high[0])
                    f.SetParLimits(1, bound_low[1], bound_high[1])
                    f.SetParLimits(2, bound_low[2], bound_high[2])
                    
                    # do the fit, ignore draws
                    fr = g.Fit(f,"NQS")
                    invalidFit=False
                    if int(fr)!=0: 
                        invalidFit = True
                    if not fr.IsValid():
                        invalidFit = True
                    tolerance=0.01
                    for ipar in range(3):
                        v = fr.Parameter(ipar)
                        if (v-bound_low[ipar])/(bound_high[ipar]-bound_low[ipar])<tolerance:
                            invalidFit = True
                        if (bound_high[ipar]-v)/(bound_high[ipar]-bound_low[ipar])<tolerance:
                            invalidFit = True
                        
                    if invalidFit:
                        if i == 0:
                            self.PEs[iFEB][iCh] = waveform[peakBinsStart]*ROOT.TMath.E()*constants.DEFAULT_BETA/calibrationFactor[iFEB][iCh] #using maximum ADC value of this pulse and a typical value of beta
                            self.pulseHeight[iFEB][iCh] = waveform[peakBinsStart] 
                            self.time[iFEB][iCh] = averagePeakBin*binWidth
                            self.beta[iFEB][iCh] = constants.DEFAULT_BETA
                            self.LEtime[iFEB][iCh] = self.time[iFEB][iCh]-0.985*constants.DEFAULT_BETA # time-0.985*beta for 50% pulse height
                            self.fitStatus[iFEB][iCh] = 2                        
                        else:
                            self.PEsReflectedPulse[iFEB][iCh] = waveform[peakBinsStart]*ROOT.TMath.E()*constants.DEFAULT_BETA/calibrationFactor[iFEB][iCh] #using maximum ADC value of this pulse and a typical value of beta
                            self.pulseHeightReflectedPulse[iFEB][iCh] = waveform[peakBinsStart]
                            self.timeReflectedPulse[iFEB][iCh] = averagePeakBin*binWidth
                            self.betaReflectedPulse[iFEB][iCh] = constants.DEFAULT_BETA
                            self.LEtimeReflectedPulse[iFEB][iCh] = self.time[iFEB][iCh]-0.985*constants.DEFAULT_BETA # time-0.985*beta for 50% pulse height
                            self.fitStatusReflectedPulse[iFEB][iCh] = 2 
                    else:
                        if i == 0:
                            self.PEs[iFEB][iCh] = fr.Parameter(0)*fr.Parameter(2) / calibrationFactor[iFEB][iCh]
                            self.pulseHeight[iFEB][iCh] = fr.Parameter(0)/ROOT.TMath.E()
                            self.time[iFEB][iCh] = fr.Parameter(1)
                            self.beta[iFEB][iCh] = fr.Parameter(2)
                            self.LEtime[iFEB][iCh] = fr.Parameter(1)-0.985*fr.Parameter(2) # at 50% of pulse height
                            self.fitStatus[iFEB][iCh] = 1
                        else:
                            self.PEsReflectedPulse[iFEB][iCh] = fr.Parameter(0)*fr.Parameter(2) / calibrationFactor[iFEB][iCh]
                            self.pulseHeightReflectedPulse[iFEB][iCh] = fr.Parameter(0)/ROOT.TMath.E()
                            self.timeReflectedPulse[iFEB][iCh] = fr.Parameter(1)
                            self.betaReflectedPulse[iFEB][iCh] = fr.Parameter(2)
                            self.LEtimeReflectedPulse[iFEB][iCh] = fr.Parameter(1)-0.985*fr.Parameter(2) # at 50% of pulse height
                            self.fitStatusReflectedPulse[iFEB][iCh] = 1
                    
                self.time[iFEB][iCh] += constants.signalRegionStart*binWidth
                self.LEtime[iFEB][iCh] += constants.signalRegionStart*binWidth
                self.timeReflectedPulse[iFEB][iCh] += constants.signalRegionStart*binWidth
                self.LEtimeReflectedPulse[iFEB][iCh] += constants.signalRegionStart*binWidth

                self.PEsTemperatureCorrected[iFEB][iCh] = -1
                self.PEsTemperatureCorrectedReflectedPulse[iFEB][iCh] = -1
                if self.temperature[iFEB][iCh] != 0:
                    temperatureCorrection=(constants.PETemperatureIntercept-constants.referenceTemperature)/(constants.PETemperatureIntercept-self.temperature[iFEB][iCh])
                    if calibrationFactorTemperatureCorrected[iFEB][iCh]!=0:
                        calibTemperatureCorrection=(constants.calibTemperatureIntercept-constants.referenceTemperature)/(constants.calibTemperatureIntercept-self.temperature[iFEB][iCh])
                        temperatureCorrection*=calibTemperatureCorrection*calibrationFactor[iFEB][iCh]/calibrationFactorTemperatureCorrected[iFEB][iCh]
                    self.PEsTemperatureCorrected[iFEB][iCh] = self.PEs[iFEB][iCh]*temperatureCorrection
                    self.PEsTemperatureCorrectedReflectedPulse[iFEB][iCh] = self.PEsReflectedPulse[iFEB][iCh]*temperatureCorrection
        return   
    
def CalculatePedestalStock(nFEB, pedestalHist, fitRange = 4, canvas = None):
# adaption of calibCrv::CrvEvent::FillCalibrationHistograms for testing purposes
    pedestals = [[None]*geometry_constants.nChannelPerFEB for i in range(nFEB)]
    noiseStdDev = [[None]*geometry_constants.nChannelPerFEB for i in range(nFEB)]
    funcPedestal = [[None]*geometry_constants.nChannelPerFEB for i in range(nFEB)]
    for iFEB in range(nFEB):
        for iCh in range(geometry_constants.nChannelPerFEB):
            if canvas:
                canvas.cd(1+iFEB*geometry_constants.nChannelPerFEB+iCh)
                gPad.cd(1)
                pedestalHist[iFEB][iCh].SetLineColor(1) # kblack
                pedestalHist[iFEB][iCh].DrawClone()
                
            pedestals[iFEB][iCh] = 0.
            if pedestalHist[iFEB][iCh].GetEntries()<200:
                continue
            n = pedestalHist[iFEB][iCh].GetNbinsX()
            overflow = pedestalHist[iFEB][iCh].GetBinContent(0)+pedestalHist[iFEB][iCh].GetBinContent(n+1)
            if overflow/float(pedestalHist[iFEB][iCh].GetEntries())>0.1:
                continue
                
            maxbinPedestal = pedestalHist[iFEB][iCh].GetMaximumBin()   
            peakPedestal = pedestalHist[iFEB][iCh].GetBinCenter(maxbinPedestal)
            # funcPedestal[iFEB][iCh] = ROOT.TF1("f_stock_%i_%i"%(iFEB,iCh), "gaus", peakPedestal-fitRange, peakPedestal+fitRange)
            funcPedestal[iFEB][iCh] = ROOT.TF1("f_stock_%i_%i"%(iFEB,iCh), "gaus", pedestalHist[iFEB][iCh].GetMean()-fitRange, pedestalHist[iFEB][iCh].GetMean()+fitRange)
            funcPedestal[iFEB][iCh].SetLineWidth(2)
            funcPedestal[iFEB][iCh].SetLineColor(2) # kRed
            pedestalHist[iFEB][iCh].Fit(funcPedestal[iFEB][iCh], "QR")
            if canvas:
                funcPedestal[iFEB][iCh].DrawClone("SAME")
            pedestals[iFEB][iCh] = funcPedestal[iFEB][iCh].GetParameter(1)
            noiseStdDev[iFEB][iCh] = pedestalHist[iFEB][iCh].GetStdDev()
            
            if canvas:
                t1 = ROOT.TPaveText(.15, .7, .50, .8, "NDC")
                t1.SetFillColor(0);
                t1.AddText("Pedestal = %.2f", pedestals[iFEB][iCh])
                t1.AddText("Noise StdDev = %.2f", noiseStdDev[iFEB][iCh])
                t1.SetTextAlign(12)
                t1.DrawClone("SAME")

    return pedestals, noiseStdDev, funcPedestal

def CalculateCalibrationConstantSingleChannel(calibHist, nPEpeaksToFit, iFEB, iCh, canvas = None, iPad1 = None, iPad2 = None, k = 0, tag = '', peakEventLimit = 50): 
    # k is whether temp corrected
    global nPreSignalRegions
    calibHist.SetLineColor(4) # kBlue
    if canvas and iPad1:
        tPad = canvas.cd(iPad1)
        tPad.SetLogy()
        calibHist.DrawClone("")

    maxbin = 0
    maxbinContent = 0
    for ibin in range(1, calibHist.GetNbinsX()):
        if calibHist.GetBinCenter(ibin)<250:
            continue # find 1PE maximum only between 350 and 1100
        if calibHist.GetBinCenter(ibin)>1500:
            break # 1100
        binContent = calibHist.GetBinContent(ibin)
        if binContent > maxbinContent:
            maxbin=ibin
            maxbinContent=binContent

    peak1PE = calibHist.GetBinCenter(maxbin)
    func1 = ROOT.TF1("f1", "gaus", peak1PE*0.8, peak1PE*1.2)
    func1.SetParameter(1, peak1PE)
    func1.SetLineWidth(1)
    func1.SetLineColor(2) # kRed
    calibHist.Fit(func1, "QR")
    peak1PE = func1.GetParameter(1)
    peak1PEerr = func1.GetParError(1)

    peaks = []
    peakserr = []
    if peak1PE>200.0 and peak1PE<1800.0:
        peaks.append(peak1PE)
        peakserr.append(peak1PEerr)
        if canvas and iPad1:
            func1.DrawClone("SAME")
        for iPeak in range(2, nPEpeaksToFit+1):
            peakTmp = iPeak*peak1PE
            rangeStart = peakTmp-0.35*peak1PE
            rangeEnd = peakTmp+0.35*peak1PE;
            if rangeEnd>calibHist.GetXaxis().GetXmax():
                break
            if calibHist.Integral(calibHist.FindFixBin(rangeStart),calibHist.FindFixBin(rangeEnd))<peakEventLimit:
                break
            func2 = ROOT.TF1("f%i"%iPeak, "gaus", rangeStart, rangeEnd)
            func2.SetParameter(1,peakTmp)
            func2.SetLineWidth(1)
            func2.SetLineColor(2) # kRed
            if (int(calibHist.Fit(func2, "QRL"))!=0): # use least log-likelihood fit
                break
            peakTmp = func2.GetParameter(1)
            peakTmperr = func2.GetParError(1)
            if (peakTmp/peak1PE>iPeak*0.95 and peakTmp/peak1PE<iPeak*1.05):
                peaks.append(peakTmp)
                peakserr.append(peakTmperr)
                if canvas and iPad1:
                    func2.DrawClone("SAME")
            else:
                break
    nPEpeaks=len(peaks)

    if canvas:
        t2 = ROOT.TPaveText(.50, .65, .89, .89, "NDC")
        t2.SetFillColor(0)
        for iPeak in range(len(peaks)):
            t2.AddText("%iPE = %4.0f #pm %4.0f ADC*ns"%(iPeak+1,peaks[iPeak],peakserr[iPeak]))
        if k == 1:
            t2.AddText("Temperature corrected values")
            t2.AddText("refT %f deg C, calibT0 %f deg C"%(constants.referenceTemperature,constants.calibTemperatureIntercept))
    if len(peaks) == 0:
        if canvas and iPad1:
            t2.DrawClone("SAME")
        return nPEpeaks, np.nan, None

    # cross-talk and noise
    events1PEandMore=0
    events2PEandMore=0
    for ibin in range(1, calibHist.GetNbinsX()):
        area = calibHist.GetBinCenter(ibin)
        binContent = calibHist.GetBinContent(ibin)
        if(area>0.5*peak1PE): 
            events1PEandMore += binContent
        if(area>1.5*peak1PE): 
            events2PEandMore += binContent
    
    noiseRate = events1PEandMore/(nPreSignalRegions[iFEB][iCh]*constants.numberOfPreSignalSamples/(constants.RATE*1e9))
    xtalkProbability = events2PEandMore/events1PEandMore
    if canvas and iPad1:
        t2.AddText("Noise rate = %4.2f MHz"%(noiseRate/1.0e6))
        t2.AddText("Xtalk prob. = %4.2f"%(xtalkProbability))
        t2.DrawClone("SAME")

    # Fit plot to determine calibration constants
    if canvas and iPad2:
        tPad = canvas.cd(iPad2)
        tPad.SetLogy(0)

    graph = ROOT.TGraphErrors()
    graph.SetPoint(0,0,0)
    graph.SetPointError(0,0,0)
    for iPeak in range(len(peaks)):
        graph.SetPoint(iPeak+1,iPeak+1,peaks[iPeak])
        graph.SetPointError(iPeak+1,0,peakserr[iPeak])
    if k == 0:
        graph.SetTitle("FEB%i Ch%i Calib. Fit; n PE peak; Area [ADC*ns]"%(iFEB,iCh))
    else:
        graph.SetTitle("FEB%i Ch%i Temp. corrected calib fit; n PE peak; Area [ADC*ns]"%(iFEB,iCh))
    graph.GetXaxis().SetRangeUser(0, len(peaks)+0.5)
    graph.GetYaxis().SetRangeUser(0, peaks[0]*(len(peaks)+0.5))
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(1) # kBlack
    if canvas and iPad2:
        graph.DrawClone("AP")

    funcFit = ROOT.TF1("calibration%i%i%i%s"%(iFEB,iCh,k,tag),"[0]*x", -0.5, len(peaks)+0.5)
    funcFit.SetLineColor(2) # kRed
    graph.Fit(funcFit, "QR")
    funcFit.DrawClone("SAME")
    if k == 0:
        rawCalibrationConstant = funcFit.GetParameter(0)
    else:
        temperatureCorrectedCalibrationConstant = funcFit.GetParameter(0)

    if canvas and iPad2:
        t3 = ROOT.TPaveText(.15, .65, .7, .88, "NDC")
        t3.SetLineColorAlpha(0,0)
        t3.SetFillStyle(0)
        t3.AddText("Calibration constant")
        t3.AddText("%.0f ADC*ns/PE"%(funcFit.GetParameter(0)))
        t3.SetTextAlign(12)
        t3.DrawClone("SAME")
    return nPEpeaks, rawCalibrationConstant if k == 0 else temperatureCorrectedCalibrationConstant, (peaks, peakserr)

def CalculateCalibrationConstantsStock(nFEB, pedestals, rawCalibHist, temperatureCorrectedCalibHist, nPEpeaksToFit, canvas = None):
# adaption of calibCrv::CrvEvent::CalculateCalibrationConstants for testing purposes
    deadChannels = [[False]*geometry_constants.nChannelPerFEB for i in range(nFEB)]
    noPedestal = [[False]*geometry_constants.nChannelPerFEB for i in range(nFEB)]
    noCalibration = [[False]*geometry_constants.nChannelPerFEB for i in range(nFEB)]
    nPEpeaks = [[0]*geometry_constants.nChannelPerFEB for i in range(nFEB)]
    rawCalibrationConstant = [[np.nan]*geometry_constants.nChannelPerFEB for i in range(nFEB)]
    temperatureCorrectedCalibrationConstant = [[np.nan]*geometry_constants.nChannelPerFEB for i in range(nFEB)]
    
    calibHist = [rawCalibHist, temperatureCorrectedCalibHist]
            
    for iFEB in range(nFEB):
        for iCh in range(geometry_constants.nChannelPerFEB):
            
            # temperatureHist omitted

            if pedestals[iFEB][iCh] is None:
                noPedestal[iFEB][iCh] = True
                continue
            
            if rawCalibHist[iFEB][iCh].GetEntries()<200:
                deadChannels[iFEB][iCh] = True
                continue

            for k in range(2):
                # omitted noiseRate, xtalkProbability
                if calibHist[k] is None:
                    continue
                
                if k == 1 and temperatureCorrectedCalibHist[iFEB][iCh].GetEntries()<200:
                    noCalibration[iFEB][iCh] = True
                    continue
                
                # calibration plot
                nn, cc, _ = CalculateCalibrationConstantSingleChannel(calibHist[k][iFEB][iCh], nPEpeaksToFit, iFEB, iCh, canvas, 1+4*(iFEB*geometry_constants.nChannelPerFEB+iCh)+2*k, 2+4*(iFEB*geometry_constants.nChannelPerFEB+iCh)+2*k, k)
                nPEpeaks[iFEB][iCh] = nn
                if k == 0:
                    rawCalibrationConstant[iFEB][iCh] = cc
                else:
                    temperatureCorrectedCalibrationConstant[iFEB][iCh] = cc
    
    return rawCalibrationConstant, temperatureCorrectedCalibrationConstant
                    
