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
            print("WARNING: CRV_event: Attribute 'detailLevel' unspecified. Changed to '0b1000' to load only basic event information.")

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
            self.temperature = np.reshape(runtree.runtree_temperature if isRaw else runtree.temperature, (nFEB, geometry_constants.nChannelPerFEB)) #runtree_temperature[4][64]/L or temperature[4][64]/F
            self.pedestal = None if isRaw else np.reshape(runtree.pedestal, (nFEB, geometry_constants.nChannelPerFEB)) #pedestal[4][64]/F
            self.tdcSinceSpill = np.reshape(runtree.runtree_tdc_since_spill if isRaw else runtree.tdcSinceSpill, (nFEB, geometry_constants.nChannelPerFEB)) #runtree_tdc_since_spill[4][64]/L or tdcSinceSpill[4][64]/L
            self.timeSinceSpill = np.reshape(runtree.runtree_time_since_spill if isRaw else runtree.timeSinceSpill, (nFEB, geometry_constants.nChannelPerFEB)) #runtree_time_since_spill[4][64]/D or timeSinceSpill[4][64]/D
            self.fitStatus = None if isRaw else np.reshape(runtree.fitStatus, (nFEB, geometry_constants.nChannelPerFEB)) #fitStatus[4][64]/I
            self.fitStatusReflectedPulse = None if isRaw else np.reshape(runtree.fitStatusReflectedPulse, (nFEB, geometry_constants.nChannelPerFEB)) #fitStatusReflectedPulse[4][64]/I
            # FIXME: global timestamp (from Spill time + time since spill)??? 
        if (detailLevel & 0b0100):
            self.adc = np.reshape(runtree.runtree_adc if isRaw else runtree.adc, (nFEB, geometry_constants.nChannelPerFEB, -1)) #adc[4][64][127]/I
        if (detailLevel & 0b0010):
            self.PEs = np.reshape(runtree.PEs, (nFEB, geometry_constants.nChannelPerFEB)) #PEs[4][64]/F
            self.PEsTemperatureCorrected = np.reshape(runtree.PEsTemperatureCorrected, (nFEB, geometry_constants.nChannelPerFEB)) #PEsTemperatureCorrected[4][64]/F
            self.pulseHeight = np.reshape(runtree.pulseHeight, (nFEB, geometry_constants.nChannelPerFEB)) #pulseHeight[4][64]/F
            self.beta = np.reshape(runtree.beta, (nFEB, geometry_constants.nChannelPerFEB)) #beta[4][64]/F
            self.time = np.reshape(runtree.time, (nFEB, geometry_constants.nChannelPerFEB)) #time[4][64]/F
            self.LEtime = np.reshape(runtree.LEtime, (nFEB, geometry_constants.nChannelPerFEB)) #LEtime[4][64]/F
            self.recoStartBin = np.reshape(runtree.recoStartBin, (nFEB, geometry_constants.nChannelPerFEB)) #recoStartBin[4][64]/I
            self.recoEndBin = np.reshape(runtree.recoEndBin, (nFEB, geometry_constants.nChannelPerFEB)) #recoEndBin[4][64]/I
        if (detailLevel & 0b0001):
            self.PEsReflectedPulse = np.reshape(runtree.PEsReflectedPulse, (nFEB, geometry_constants.nChannelPerFEB)) #PEsReflectedPulse[4][64]/F
            self.PEsTemperatureCorrectedReflectedPulse = np.reshape(runtree.PEsTemperatureCorrectedReflectedPulse, (nFEB, geometry_constants.nChannelPerFEB)) #PEsTemperatureCorrectedReflectedPulse[4][64]/F
            self.pulseHeightReflectedPulse = np.reshape(runtree.pulseHeightReflectedPulse, (nFEB, geometry_constants.nChannelPerFEB)) #pulseHeightReflectedPulse[4][64]/F
            self.betaReflectedPulse = np.reshape(runtree.betaReflectedPulse, (nFEB, geometry_constants.nChannelPerFEB)) #betaReflectedPulse[4][64]/F
            self.timeReflectedPulse = np.reshape(runtree.timeReflectedPulse, (nFEB, geometry_constants.nChannelPerFEB)) #timeReflectedPulse[4][64]/F
            self.LEtimeReflectedPulse = np.reshape(runtree.LEtimeReflectedPulse, (nFEB, geometry_constants.nChannelPerFEB)) #LEtimeReflectedPulse[4][64]/F
            self.recoStartBinReflectedPulse = np.reshape(runtree.recoStartBinReflectedPulse, (nFEB, geometry_constants.nChannelPerFEB)) #recoStartBinReflectedPulse[4][64]/I
            self.recoEndBinReflectedPulse = np.reshape(runtree.recoEndBinReflectedPulse, (nFEB, geometry_constants.nChannelPerFEB)) #recoEndBinReflectedPulse[4][64]/I
    
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

