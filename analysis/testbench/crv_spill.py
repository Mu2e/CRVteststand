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
import crv_event
import pprint
from datetime import datetime
import time 

class crv_spill:
    def __init__(self, spilltree, iSpill, isRaw = False):

        # initializer to load CRV_spill from ROOT tree.
        # input: spilltree: ROOT TTree object that corresponds to "spill" tree in the crvparsed or crvreco ROOT files.
        #        iSpill: index of the spill to load.
        
        try:
            spilltree.GetEntry(iSpill)
        except Exception as e:
            print("Error: CRV_spill: reading spill #%i unsuccessful!"%(iSpill))
            print("                ", type(e), e)
            return
        
        # prepare the containers.
        attributeList = [# tree stuff
                         "spillIndex", "spillNumber", "nEvents", "nEventsActual", "spillStored",
                         "nFEB", "nChannel", "nSample", 
                         # parsed from boardStatus
                         "boardID", "recentSpillNumber", "temperatureFEB",
                         "supply15V", "supply10V", "supply5V", "supplyN5V", 
                         "supply3V3", "supply2V5", "supply1V8", "supply1V2", 
                         "busSiPMBias", "TrigCtrlReg", "settingPipelineLen", "settingSampleLen",
                         "tsEpoch",
                         # operation containers
                         "temperatureCMB", # getTempCMB provides array of [nFEB][nChannelPerFEB]
                         "dqmRedFlag"
                         ]
        for item in attributeList:
            exec("self."+item+"=None")

        self.temp_dict = {}
        self.iSpill = iSpill
        self.isRaw = isRaw

        self.spillIndex = spilltree.spill_index # spill_index/I
        self.spillNumber = spilltree.spill_num # spill_num/I
        self.nEvents = spilltree.spill_nevents # spill_nevents/I
        self.nEventsActual = spilltree.spill_neventsActual # spill_neventsActual/I
        self.spillStored = spilltree.spill_stored # spill_stored/O, boolean 
        self.nFEB = spilltree.spill_number_of_febs # spill_number_of_febs/I
        self.nChannel = spilltree.spill_channels_per_feb # spill_channels_per_feb/I
        self.nSample = spilltree.spill_number_of_samples # spill_number_of_samples/I

        # board status from docdb-37635
        boardStatus = np.reshape(spilltree.spill_boardStatus, (self.nFEB, 22)) # spill_boardStatus[4][22]/I
        self.boardID = np.array([boardStatus[i][0] for i in range(self.nFEB)])
        self.recentSpillNumber = np.array([boardStatus[i][1] for i in range(self.nFEB)])
        self.temperatureFEB = np.array([boardStatus[i][2]*0.01 for i in range(self.nFEB)]) # typo in docdb-37635        
        self.supply1V2 = np.array([boardStatus[i][3]*0.001 for i in range(self.nFEB)])
        self.supply1V8 = np.array([boardStatus[i][4]*0.001 for i in range(self.nFEB)])
        self.supply5V = np.array([boardStatus[i][5]*0.002 for i in range(self.nFEB)])
        self.supply10V = np.array([boardStatus[i][6]*0.004 for i in range(self.nFEB)])
        self.supply2V5 = np.array([boardStatus[i][7]*0.001 for i in range(self.nFEB)])
        self.supplyN5V = np.array([boardStatus[i][8]*(-0.002) for i in range(self.nFEB)])
        self.supply15V = np.array([boardStatus[i][9]*0.006 for i in range(self.nFEB)])
        self.supply3V3 = np.array([boardStatus[i][10]*0.001 for i in range(self.nFEB)])        
        self.busSiPMBias = np.array([boardStatus[i][11:19] for i in range(self.nFEB)])/50.
        self.TrigCtrlReg = np.array([boardStatus[i][19] for i in range(self.nFEB)])
        self.settingPipelineLen = np.array([boardStatus[i][20] for i in range(self.nFEB)])
        self.settingSampleLen = np.array([boardStatus[i][21] for i in range(self.nFEB)])        
        # Trigger Control Register read
        # Bit 0: Always reads 0
        # Bit 1: Selects the trigger input type as a pulse or an FM data stream
        # The assumption is that a trigger pulse comes from the LEMO connector, the FM encoded trigger
        # message comes from the RJ-45 connector. The microprocessor controls the multiplexer that routes
        # either the LEMO or the RJ-45 signal to the trigger input on the FPGAs.
        # Bit 2: Trigger Inhibit. If trigger inhibit is enabled, this bit goes to one in response to a trigger.
        # Bit 3: Trigger Inhibit enable.
        # Bit 4: Spill Inhibit. If spill inhibit is enabled, this bit goes to one in response to end of spill.
        # Bit 5: Spill Inhibit clear request. Returns a ‘1’ from the time clear inhibit was written until spill inhibit clears.
        # Bit 6: Spill Inhibit enable.
        # Bit 7: Not used.
        # Bit 8: On card test pulser enabled
        # Bit 9: Test pulser enabled for one spill only

        if hasattr(spilltree, "spill_timestamp"):
            self.tsEpoch = spilltree.spill_timestamp
            if self.tsEpoch == 0:
                self.tsEpoch = None

        if self.tsEpoch == None:
            if spilltree.spill_timestamp_year != 0:
                raw_dt = time.struct_time([1900 + spilltree.spill_timestamp_year, # tm_year
                                           1 + spilltree.spill_timestamp_mon, # tm_mon
                                           spilltree.spill_timestamp_mday, # tm_mday
                                           spilltree.spill_timestamp_hour, # tm_hour
                                           spilltree.spill_timestamp_min, # tm_min
                                           spilltree.spill_timestamp_sec, # tm_sec
                                           (spilltree.spill_timestamp_wday-1)%7, # tm_wday
                                           1 + spilltree.spill_timestamp_yday, # tm_yday
                                           spilltree.spill_timestamp_isdst])  # tm_isdst

                self.tsEpoch = time.mktime(raw_dt)
                # use "datetime.fromtimestamp(self.tsEpoch)" to convert to datetime objects
                
                root_utils.lastTsEpoch = (self.spillNumber, self.tsEpoch)
            else: 
                if root_utils.lastTsEpoch and self.spillNumber>root_utils.lastTsEpoch[0]:
                    self.tsEpoch = root_utils.lastTsEpoch[1]+float(self.spillNumber-root_utils.lastTsEpoch[0])*(constants.tDataTaking+constants.tDataTransfer)
                else:
                    sys.exit("Error: CRV_spill: spill # %i reading spill timestamp unsuccessful! \n"%self.spillNumber +
                             "                  last good timestamp at spill # %i: %i"%(root_utils.lastTsEpoch if root_utils.lastTsEpoch else (-1,-1)))

    def getTempCMB(self, runtree, nFEBFile, iEntryStart, doAverage = True):
        averageTemp = None
        nEvent = self.nEventsActual            
        if nEvent == 0:
            print("WARNING: CRV_spill: getTempCMB: Spill # %i has no actual event"%self.spillNumber)
            averageTemp = np.empty((nFEBFile,geometry_constants.nChannelPerFEB,))
            averageTemp[:] = np.nan
        else:
            if doAverage:                
                for index in range(self.nEventsActual):
                    tEvent = crv_event.crv_event(runtree, iEntryStart+index, 0b1000, nFEBFile, self.isRaw)
                    if tEvent.spillNumber != self.spillNumber:
                        nEvent = index
                        print("WARNING: CRV_spill: getTempCMB: Spill # %i reporting %i actual events, containing %i events"%(self.spillNumber, self.nEventsActual, nEvent))
                        break
                    if index == 0:                        
                        averageTemp = tEvent.temperature
                    else:                        
                        averageTemp = (averageTemp*index+tEvent.temperature)/(index+1)
            else:
                tEvent = crv_event.crv_event(runtree, iEntryStart, 0b1000, nFEBFile, self.isRaw)
                averageTemp = tEvent.temperature
        self.temperatureCMB = averageTemp
        return nEvent

    #??? active FEB mask???

    def checkDQM(self, verbose = True):
        # filling in dqmRedFlag.
        # a spill that has everything normal should see dqmRedFlag == 0x0
        # bit masks for different criteria. 
        # 0x01: spillIndex != spillNumber
        # 0x02: spillStored == 0
        # 0x04: nEventsActual == 0
        # 0x08: boardStatus values (other than bulk biases) out of reasonable range
        # 0xf0: reserved
        # all the higher bits: bulk bias voltages < 50 V for each AFE (nFEB*8)
        self.dqmRedFlag = 0x0
        if verbose:
            print ("***", self.spillIndex, self.spillNumber)
        if (self.spillIndex != self.spillNumber):
            self.dqmRedFlag |= 0x01
        if verbose:
            print ("***", self.spillStored)
        if (self.spillStored == 0):
            self.dqmRedFlag |= 0x02
        if verbose:
            print ("***", self.nEventsActual)
        if (self.nEventsActual == 0):
            self.dqmRedFlag |= 0x04
        if verbose:
            print ("***", self.supply15V, self.supply10V, self.supply5V, self.supplyN5V, self.supply3V3, self.supply2V5, self.supply1V8, self.supply1V2)
        if ((self.supply15V < 14.5).any() or (self.supply15V > 15.5).any() or
            # (self.supply10V < 9.5).any() or (self.supply10V > 10.5).any() or # One of the FEB seems to be consistently around 10.8V
            (self.supply10V < 9.0).any() or (self.supply10V > 11.0).any() or
            (self.supply5V < 4.5).any() or (self.supply5V > 5.5).any() or
            (self.supplyN5V < -5.5).any() or (self.supplyN5V > -4.5).any() or
            (self.supply3V3 < 2.8).any() or (self.supply3V3 > 3.8).any() or
            (self.supply2V5 < 2.).any() or (self.supply2V5 > 3.).any() or
            (self.supply1V8 < 1.3).any() or (self.supply1V8 > 2.3).any() or
            (self.supply1V2 < 0.7).any() or (self.supply1V2 > 1.7).any() ):
            self.dqmRedFlag |= 0x08
        for i in range(self.nFEB):
            if verbose:
                print ("***", self.busSiPMBias[i])
            for j in range(8):
                if self.busSiPMBias[i][j] < 50:
                    self.dqmRedFlag |= (0b1 << (8+i*8+j))
        if verbose and self.dqmRedFlag != 0:
            print ("Spill index %i failed DQM check: %x"%(self.spillIndex, self.dqmRedFlag))
        return

        # FIXME: temperature DQM???
            