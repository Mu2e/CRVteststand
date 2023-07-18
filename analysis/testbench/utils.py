from __future__ import print_function

import sys
import re
import os
import json
import datetime
from datetime import datetime
import threading
import readline
import gc
from copy import deepcopy
import numpy as np
import ROOT
from ROOT import TFile, TCanvas, TH1F, TH2F, TF1, TMath, TGraph
from ROOT import gStyle, gROOT, gDirectory
from array import array
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import math
import argparse, textwrap
import root_utils
import constants, geometry_constants
import crv_event, crv_spill
import filepath
import manualDQC

################################## UTILS FOR TIME PLOT ######################################

def ts2datetime(ts): # inheritance from g-2 field analysis code
    try:
        if len(ts) == 0:
             return ts
    except:
        pass
    return np.vectorize(datetime.fromtimestamp)(ts)

def plot_ts(*args, **kwargs): # inheritance from g-2 field analysis code
    from matplotlib.dates import DateFormatter
    if 'dformat' in kwargs:
        dformat = kwargs['dformat']
        del kwargs['dformat']
    else:
        if np.amax(args[0])-np.amin(args[0])<3600*24*10: # shorter than 10 days by default include hr/min
            dformat = '%m/%d\n%H:%M'
        else:                                            # otherwise only use yr/mon/day
            dformat = '%Y\n%m/%d'
    formatter = DateFormatter(dformat)
    ax = plt.plot_date(ts2datetime(args[0]), *args[1:], **kwargs)
    plt.gca().xaxis.set_major_formatter(formatter)
    return ax

def smoothing(x, nSmooth, axis): # running average using nSmooth points
    if axis == 'x':
        x_smoothed = np.array(x[(nSmooth-1):])
    elif axis == 'y':
        x_smoothed = np.array([np.mean(x[i:(i+nSmooth)]) for i in range(len(x[(nSmooth-1):]))])
    else:
        sys.exit("Error: utils.smoothing: bad axis argument")
    return x_smoothed

def plot_dqm(filename_list, plot_dict, dqmFilter = '', dqmVerbose = False, nSmooth = 1, show = False, title = ';;', isRaw = False): 
    
    # plot_dict is a dictionary of {"keyword":[["attribute[slicing]"],[next plot]]}; no omission is allowed if slicing is used
    # keyword is the selection criteria in the file name; if keyword is not an attribute name, 
    # it would be searched for in crv_spill.temp_dict.
    #
    # if an attribute has more than 1 element and the slicing is not specified, all elements 
    # are plotted together
    
    plotted_attribute = []
    full_plot_dict = {}
    plotData = {}
    fig_list = []

    for k in plot_dict.keys():
        full_plot_dict.update({k:[]})
        plotData.update({k:{}})

    for k,plotlist in plot_dict.items():
        for i, vl in enumerate(plotlist):
            full_plot_dict[k].append([])
            for v in vl:
                tAttribute = v.split('[')[0] if '[' in v else v
                if tAttribute not in plotted_attribute:
                    plotted_attribute.append(tAttribute)
                nColon = v.count(':')
                if nColon == 0:
                    full_plot_dict[k][i].append(v)
                else:
                    tNameStem = '['.join(v.split('[')[:-1])
                    tSlicingIndex = v.split('[')[-1].split(']')[0].split(':')
                    tStart = int(tSlicingIndex[0])
                    tEnd = int(tSlicingIndex[1])
                    tSkip = 1
                    if nColon > 1:
                        tSkip = int(tSlicingIndex[2])
                    for j in range(tStart, tEnd, tSkip):
                        full_plot_dict[k][i].append(tNameStem+'[%i]'%(j))
    
    # print (plotted_attribute)
    # print (full_plot_dict)
    
    for k,plotlist in full_plot_dict.items():
        # treat each group of files matching keyword 'k' together
        plotData[k].update({"ts":[], "dqm":[]})
        for vl in plotlist:
            for v in vl:
                plotData[k].update({v:[]})
        for filename in filename_list:
            if k in filename.split('/')[-1] or k == '*':
                fFile = TFile(filename, "READ")
                if 'crvaging' in filename:
                    isCosmic = True
                else:
                    isCosmic = False
                print ("Reading file: "+filename.split('/')[-1])
                runtree = fFile.Get("run")
                spilltree = fFile.Get("spills")
                iEvent = 0
                nEvent = runtree.GetEntries()
                nSpill = spilltree.GetEntries()
                for iSpill in range(nSpill):
                    tSpill = crv_spill.crv_spill(spilltree, iSpill, isRaw, isCosmic)
                    if "temperatureCMB" in plotted_attribute:
                        tSpill.getTempCMB(runtree, tSpill.nFEB, iEvent, False) # populate tSpill.temperatureCMB[nFEB][nChannelPerFEB]; no need to do avg, copies of same values per spill
                        iEvent += tSpill.nEventsActual
                        if iEvent > nEvent:
                            sys.exit("Error: utils.plot_dqm: %s nEvent = %i, Spill # %i is reading iEvent = %i"%(filename, nEvent, iSpill, iEvent))
                    # print (tSpill.temperatureCMB)
                    if dqmFilter:
                        tSpill.checkDQM(dqmVerbose) # populate tSpill.dqmRedFlag
                    if tSpill.tsEpoch == None:
                        if iSpill == 0: # first spills of subruns >= 001 have timestamp processing issues... refer to next spill for timestamp
                            try:
                                nextSpillTsEpoch = crv_spill.crv_spill(spilltree, 1, isRaw, isCosmic).tsEpoch
                                tSpill.tsEpoch = nextSpillTsEpoch - (constants.timeDataSpill_cosmic if isCosmic else constants.timeDataSpill_led)
                                tnbit = tSpill.dqmRedFlag.bit_length()
                                tmask = (((1<<tnbit) - 1) ^ 0x10) # this gives 0b1...101111
                                tSpill.dqmRedFlag &= tmask
                                print ("*** spill %04i timestamp amended: %i-%i"%(tSpill.spillNumber, nextSpillTsEpoch, (constants.timeDataSpill_cosmic if isCosmic else constants.timeDataSpill_led)))
                            except:
                                sys.exit("Error: utils.plot_dqm: %s reading first 2 spill timestamps failed"%(filename))
                        else:
                            print ("!!! spill %04i skipped; the spill has no valid timestamp; 0x%x"%(tSpill.spillNumber, tSpill.dqmRedFlag))
                            continue
                    if dqmFilter:
                        if not eval("tSpill.dqmRedFlag"+dqmFilter):
                            if (tSpill.dqmRedFlag & 0x6) != 0x6: # not empty spill
                                print ("!!! spill %i, %04i DQM = 0x%x"%(iSpill, tSpill.spillNumber, tSpill.dqmRedFlag))
                    plotData[k]["ts"].append(tSpill.tsEpoch)
                    if dqmFilter:
                        plotData[k]['dqm'].append(eval("tSpill.dqmRedFlag"+dqmFilter))
                    else:
                        plotData[k]['dqm'].append(True)
                    for vl in plotlist:
                        for v in vl:
                            if len(plotData[k][v]) < len(plotData[k]["ts"]): # prevent filling in twice in case of repeated attribute in different figures
                                try:
                                    plotData[k][v].append(eval("tSpill.%s"%v))
                                except:
                                    try:
                                        tval = eval("tSpill.temp_dict['"+v.split('[')[0]+"']"+("["+'['.join(v.split('[')[1:]) if len(v.split('['))>1 else ''))
                                        plotData[k][v].append(tval)
                                    except:
                                        sys.exit("Error: utils.plot_dqm: %s does not exist in the spill attribute"%(v))
                fFile.Close()
        
        if len(plotData[k]["ts"]) == 0:
            print ("!!! No data found matching keyword %s"%(k))
            continue
        
        plotData[k]["ts"] = smoothing(plotData[k]["ts"], nSmooth, 'x')
        plotData[k]["dqm"] = smoothing(plotData[k]["dqm"], nSmooth, 'x')
        for vl in plotlist:
            for v in vl:
                plotData[k][v] = smoothing(plotData[k][v], nSmooth, 'y')
        
        # plot the list in corresponding figures
        for vl in plotlist:
            fig = plt.figure()
            for i, v in enumerate(vl):
                if dqmFilter:
                    # print (k, v)
                    # print (plotData[k]['ts'].shape, plotData[k]['ts'])
                    # print (plotData[k][v].shape, plotData[k][v])
                    plot_ts(plotData[k]['ts'], plotData[k][v], '.', markersize=2,  label="", rasterized=True, color='gray')
                tdqc = plotData[k]['dqm']
                plot_ts(plotData[k]['ts'][tdqc], plotData[k][v][tdqc], '.', markersize=1.2,  label="%s"%(v), rasterized=True, color=constants.colors[i%10])
            ax = plt.gca()
            ax.set_title(title.split(';')[0]+('' if k == '*' else k))
            ax.set_xlabel(title.split(';')[1])
            ax.set_ylabel(title.split(';')[2])
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            lgd_h = 0.045*np.ceil(len(vl)/2.+1)
            fig.subplots_adjust(bottom=lgd_h+0.1)
            plt.legend(markerscale=4, ncol=2, bbox_to_anchor=(0.5,0), loc="lower center", bbox_transform=fig.transFigure, labelspacing=0.2, handletextpad=0.1, columnspacing=0.8)
            plt.tight_layout() 

            if show:
                plt.show()
            fig_list.append(fig)

    return fig_list

std_plotAttribute = [["temperatureFEB[0:2]", "temperatureCMB[0][0]"],
                     ["temperatureFEB[2:4]", "temperatureCMB[2][0]"],
                     ["temperatureFEB[4:6]", "temperatureCMB[4][0]"],
                     ["supply15V[0:2]", "supply10V[0:2]"],
                     ["supply5V[0:2]", "supplyN5V[0:2]"],
                     ["supply3V3[0:2]", "supply2V5[0:2]"],
                     ["supply1V8[0:2]", "supply1V2[0:2]"],
                     ["supply15V[2:4]", "supply10V[2:4]"],
                     ["supply5V[2:4]", "supplyN5V[2:4]"],
                     ["supply3V3[2:4]", "supply2V5[2:4]"],
                     ["supply1V8[2:4]", "supply1V2[2:4]"],
                     ["supply15V[4:6]", "supply10V[4:6]"],
                     ["supply5V[4:6]", "supplyN5V[4:6]"],
                     ["supply3V3[4:6]", "supply2V5[4:6]"],
                     ["supply1V8[4:6]", "supply1V2[4:6]"],
                     ["busSiPMBias[0][0:4]"],
                     ["busSiPMBias[0][4:8]"],
                     ["busSiPMBias[1][0:4]"],
                     ["busSiPMBias[1][4:8]"],
                     ["busSiPMBias[2][0:4]"],
                     ["busSiPMBias[2][4:8]"],
                     ["busSiPMBias[3][0:4]"],
                     ["busSiPMBias[3][4:8]"],
                     ["busSiPMBias[4][0:4]"],
                     ["busSiPMBias[4][4:8]"],
                     ["busSiPMBias[5][0:4]"],
                     ["busSiPMBias[5][4:8]"],
                     ["settingPipelineLen[0:2]", "settingSampleLen[0:2]"],
                     ["settingPipelineLen[2:4]", "settingSampleLen[2:4]"],
                     ["settingPipelineLen[4:6]", "settingSampleLen[4:6]"],
                     ["temperatureCMB[0][0:16:4]"], # channels on the same CMB read a single CMB temperature
                     ["temperatureCMB[0][16:32:4]"],
                     ["temperatureCMB[0][32:48:4]"],
                     ["temperatureCMB[0][48:64:4]"],
                     ["temperatureCMB[1][0:16:4]"],
                     ["temperatureCMB[1][16:32:4]"],
                     ["temperatureCMB[1][32:48:4]"],
                     ["temperatureCMB[1][48:64:4]"],
                     ["temperatureCMB[2][0:16:4]"],
                     ["temperatureCMB[2][16:32:4]"],
                     ["temperatureCMB[2][32:48:4]"],
                     ["temperatureCMB[2][48:64:4]"],
                     ["temperatureCMB[3][0:16:4]"],
                     ["temperatureCMB[3][16:32:4]"],
                     ["temperatureCMB[3][32:48:4]"],
                     ["temperatureCMB[3][48:64:4]"],
                     ["temperatureCMB[4][0:16:4]"],
                     ["temperatureCMB[4][16:32:4]"],
                     ["temperatureCMB[4][32:48:4]"],
                     ["temperatureCMB[4][48:64:4]"],
                     ["temperatureCMB[5][0:16:4]"],
                     ["temperatureCMB[5][16:32:4]"],
                     ["temperatureCMB[5][32:48:4]"],
                     ["temperatureCMB[5][48:64:4]"]]

##################################   OTHER UTILS   ######################################

def calibExtract_txt(filename):
    # OBSOLETE METHOD
    # localfilename = filename.split('/')[-1]
    # if localfilename.split('.')[0] != 'cal':
        # tpath = '/'.join(filename.split('/')[:-2])+'/'+filepath.file_calib_subdir
        # localfilename = 'cal.'+'.'.join(localfilename.split('.')[1:-1])+'.txt'
        # filename = tpath + localfilename
    # print ("Loading calibration information from", localfilename)
    
    if localfilename.split('.')[-1] != 'txt' and localfilename.split('.')[-1] != 'TXT':
        localfilename = filepath.findlinked(localfilename, "caliTXT")
    print ("Loading calibration information from", localfilename)
    
    try:
        with open(filename, 'r') as fin:
            lines = fin.readlines()
            lines = lines[2:] # rid of comment line and header, data in the order of  
                              #'FEB' 'Channel' 'Pedestal' 'Calib' 'CalibT'
            max_FEB = 0
            max_Ch = 0
            
            for i in range(len(lines)):
                fields = lines[i].split()
                tFEB = int(fields[0])
                tCh = int(fields[1])
                tPed = float(fields[2])
                tcalib = float(fields[3])
                tCalibT = float(fields[4])
                max_FEB = max(max_FEB, tFEB)
                max_Ch = max(max_Ch, tCh)
                lines[i] = (tFEB, tCh, tPed, tcalib, tCalibT)

            nFEB = max_FEB + 1
            nCh = max_Ch + 1

            pedestal = np.empty((nFEB, nCh), dtype=np.float32)
            calibRaw = np.empty((nFEB, nCh), dtype=np.float32)
            calibAdj = np.empty((nFEB, nCh), dtype=np.float32)
            
            for i in range(len(lines)):
                ttpl = lines[i]
                pedestal[ttpl[0], ttpl[1]] = ttpl[2]
                calibRaw[ttpl[0], ttpl[1]] = ttpl[3]
                calibAdj[ttpl[0], ttpl[1]] = ttpl[4]

        return pedestal, calibRaw, calibAdj

    except:
        print ("*** SEVERE WARNING: utils.calibExtract: %s does not exit."%(filename))
        return None, None, None
    
def calibExtract_root(filename):
    localfilename = filename.split('/')[-1]
    print ("Loading calibration information from", localfilename)
    
    try:
        fFile = TFile(filename, "READ")
        runSummarytree = fFile.Get("runSummary")
        spilltree = fFile.Get("spills")
        spilltree.GetEntry(0)
        nFEB = spilltree.spill_number_of_febs
        runSummarytree.GetEntry(0)
        
        pedestal = np.reshape(np.array(runSummarytree.pedestals, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        calibRaw = np.reshape(np.array(runSummarytree.calibConstants, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        calibAdj = np.reshape(np.array(runSummarytree.calibConstantsTemperatureCorrected, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        
        return pedestal, calibRaw, calibAdj

    except:
        print ("*** SEVERE WARNING: utils.calibExtract: error during extracting calibrations from %s."%(filename))
        return None, None, None
    
spill_summary_branch_list = [("filename", "vstring"), ("runNum", "i"), ("subrunNum", "i"), ("runType", "vstring"), 
                             ("nFEB", "i"), ("nChannel", "i"), ("iSpill", "i"), ("spillIndex", "i"), 
                             ("eventStart", "i"), ("eventNum", "i"), ("dqcRedFlag", "vstring"), # dqc flag is saved as string for possible expansion
                             ("tsEpoch", "L"), ("busSiPMBias", "vdouble"), ("temperatureFEB", "vdouble"), ("temperatureCMB", "vdouble")]
subrun_summary_branch_list = [("filename", "vstring"), ("runNum", "i"), ("subrunNum", "i"), ("runType", "vstring"), 
                              ("nFEB", "i"), ("nChannel", "i"), ("spillStart", "i"), ("spillNum", "i"), 
                              ("pedestal", "vdouble"), ("calibRaw", "vdouble"), ("calibAdj", "vdouble"),
                              ("tsEpochStart", "L"), ("tsEpochEnd", "L"),
                              ("busSiPMBias", "vdouble"), ("temperatureFEB", "vdouble"), ("temperatureCMB", "vdouble"), # averaged
                              ("spectrum", "vstring"),("spectrumCorrected", "vstring"),("spectrumCustomized", "vstring")]
spec_nbin = 200
spec_xmin = 0.
spec_xmax = 200.

def processSubRun(filename, fROOT_output, spill_summary_tree, subrun_summary_tree, spill_summary_entry_dict, subrun_summary_entry_dict, 
                  applyManualDQC = True, dqmFilter = '== 0x0', dqmVerbose = False, spectrumTempCorrection = "110", weightedAvg = True):
    
    # spectrumTempCorrection controls the temperature correction used for the spectrum
    # 1st digit: 0/1 whether generate uncorrected spectrum
    # 2nd digit: 0/1 whether generate using temp. corrected values in reco files
    # 3nd digit: 0/1 whether generate using customized correction
    
    localfilename = filename.split('/')[-1]
    print ("Reading file:", localfilename)
    if 'crvaging' in localfilename:
        isCosmic = True
    else:
        isCosmic = False
    if localfilename[:3] == 'rec':
        isRaw = False
    else:
        isRaw = True
    runNum = filepath.filenameparser(filename, 'run')
    subrunNum = filepath.filenameparser(filename, 'subrun')
    runType = "cosmics" if isCosmic else "LED"
    spectrumNameStem = "run_%06i_%03i_"%(runNum,subrunNum)+runType+"_"
    
    root_utils.treeEntryDictPurge(spill_summary_branch_list, spill_summary_entry_dict)
    root_utils.treeEntryDictPurge(subrun_summary_branch_list, subrun_summary_entry_dict)
    spill_summary_entry_dict["filename"].push_back(filename)
    subrun_summary_entry_dict["filename"].push_back(filename)
    spill_summary_entry_dict["runNum"][0] = runNum
    subrun_summary_entry_dict["runNum"][0] = runNum
    spill_summary_entry_dict["subrunNum"][0] = subrunNum
    subrun_summary_entry_dict["subrunNum"][0] = subrunNum
    spill_summary_entry_dict["runType"].push_back(runType)
    subrun_summary_entry_dict["runType"].push_back(runType)    
    
    if applyManualDQC:
        if (runNum, subrunNum) in manualDQC.exclude_subrun:
            return
    
    fFile = TFile(filename, "READ")
    runtree = fFile.Get("run")
    spilltree = fFile.Get("spills")
    iEvent = 0
    nEvent = runtree.GetEntries()
    nSpill = spilltree.GetEntries()
    nSpillUsed = 0
    nFEB = None
    
    temperatureFEB_avg = None
    temperatureCMB_avg = None
    busSiPMBias_avg = None
    divider = 0
    
    for iSpill in range(nSpill):
        root_utils.treeEntryDictPurge(spill_summary_branch_list, spill_summary_entry_dict, ["filename", "runType"])
        tSpill = crv_spill.crv_spill(spilltree, iSpill, isRaw,  isCosmic)
        nFEB = tSpill.nFEB
        nEventInSpill = tSpill.nEventsActual
        if nEventInSpill > 0:
            tSpill.getTempCMB(runtree, nFEB, iEvent, False)
        
        if iSpill != 0 and iSpill%200 == 0:
            print("Run%04i_%03i processed %i/%i spills..."%(runNum,subrunNum,iSpill,nSpill))
        
        if iSpill == 0:
            temperatureFEB_avg = np.zeros(nFEB, dtype=np.float64)
            temperatureCMB_avg = np.zeros((nFEB, geometry_constants.nChannelPerFEB), dtype=np.float64)
            busSiPMBias_avg = np.zeros((nFEB, geometry_constants.nAFEPerFEB), dtype=np.float64)
            hspectrum = [None]*(nFEB*geometry_constants.nChannelPerFEB)
            hspectrum_corrected = [None]*(nFEB*geometry_constants.nChannelPerFEB)
            hspectrum_customized = [None]*(nFEB*geometry_constants.nChannelPerFEB)
            spill_summary_entry_dict["nFEB"][0] = nFEB
            subrun_summary_entry_dict["nFEB"][0] = nFEB
            spill_summary_entry_dict["nChannel"][0] = geometry_constants.nChannelPerFEB
            subrun_summary_entry_dict["nChannel"][0] = geometry_constants.nChannelPerFEB
            subrun_summary_entry_dict["spillStart"][0] = spill_summary_tree.GetEntries()
            subrun_summary_entry_dict["spillNum"][0] = nSpill
            
            fROOT_output.cd()
            for iFEB in range(nFEB):
                for iCh in range(geometry_constants.nChannelPerFEB):
                    tSpecName = spectrumNameStem+"FEB%02i_Ch%02i"%(iFEB, iCh)
                    if spectrumTempCorrection[0] == '1':
                        hspectrum[iFEB*geometry_constants.nChannelPerFEB+iCh] = TH1F(tSpecName+"_uncorrected", tSpecName+";PE (uncorrected);count", spec_nbin, spec_xmin, spec_xmax)
                        subrun_summary_entry_dict["spectrum"].push_back("spectrum/"+tSpecName+"_uncorrected")
                    else:
                        subrun_summary_entry_dict["spectrum"].push_back("")
                    if spectrumTempCorrection[1] == '1':
                        hspectrum_corrected[iFEB*geometry_constants.nChannelPerFEB+iCh] = TH1F(tSpecName+"_corrected", tSpecName+";PE (corrected);count", spec_nbin, spec_xmin, spec_xmax)
                        subrun_summary_entry_dict["spectrumCorrected"].push_back("spectrum/"+tSpecName+"_corrected")
                    else:
                        subrun_summary_entry_dict["spectrumCorrected"].push_back("")
                    if spectrumTempCorrection[2] == '1':
                        hspectrum_customized[iFEB*geometry_constants.nChannelPerFEB+iCh] = TH1F(tSpecName+"_customized", tSpecName+";PE (standalone correction);count", spec_nbin, spec_xmin, spec_xmax)
                        subrun_summary_entry_dict["spectrumCustomized"].push_back("spectrum/"+tSpecName+"_customized")
                    else:
                        subrun_summary_entry_dict["spectrumCustomized"].push_back("")
        
        if applyManualDQC:
            if (runNum, subrunNum) in manualDQC.manual_reset_spillNum_offset_dict:
                if iSpill in manualDQC.manual_reset_spillNum_offset_dict[(runNum, subrunNum)]:
                    root_utils.startSpill = tSpill.spillNumber-tSpill.spillIndex
                    print ("*** offset between spillNumber and spillIndex realigned.")
        tSpill.checkDQM(dqmVerbose) # populate tSpill.dqmRedFlag
        if tSpill.tsEpoch == None:
            if iSpill == 0: # first spills of subruns >= 001 have timestamp processing issues... refer to next spill for timestamp
                try:
                    nextSpillTsEpoch = crv_spill.crv_spill(spilltree, 1, isRaw, isCosmic).tsEpoch
                    tSpill.tsEpoch = nextSpillTsEpoch - (constants.timeDataSpill_cosmic if isCosmic else constants.timeDataSpill_led)
                    tnbit = tSpill.dqmRedFlag.bit_length()
                    tmask = (((1<<tnbit) - 1) ^ 0x10) # this gives 0b1...101111
                    tSpill.dqmRedFlag &= tmask
                    print ("*** spill %04i timestamp amended: %i-%i"%(tSpill.spillNumber, nextSpillTsEpoch, (constants.timeDataSpill_cosmic if isCosmic else constants.timeDataSpill_led)))
                except:
                    sys.exit("Error: utils.processSubRun: %s reading first 2 spill timestamps failed"%(filename))
            else:
                print ("!!! spill %04i skipped; the spill has no valid timestamp; 0x%x"%(tSpill.spillNumber, tSpill.dqmRedFlag))
                continue
        if applyManualDQC:
            if (runNum, subrunNum) in manualDQC.manaulDQC_dict:
                if iSpill in manualDQC.manaulDQC_dict[(runNum, subrunNum)]:
                    tSpill.dqmRedFlag = manualDQC.manaulDQC_dict[(runNum, subrunNum)][iSpill]
                    print ("*** spill %i, %04i dqmRedFlag manually set to 0x%x"%(iSpill, tSpill.spillNumber, tSpill.dqmRedFlag))
        DQMgood = True
        if dqmFilter:
            DQMgood = eval("tSpill.dqmRedFlag"+dqmFilter)
        
        if DQMgood:
            for index in range(nEventInSpill):
                tEvent = crv_event.crv_event(runtree, iEvent+index, 0b1010, nFEB, tSpill.isRaw)
                if tEvent.spillNumber != tSpill.spillNumber:
                    nEventInSpill = index
                    print("*** SEVERE WARNING: utils: processSubRun: Spill # %i reporting %i actual events, containing %i events"%(tSpill.spillNumber, tSpill.nEventsActual, nEventInSpill))
                    break
                for iFEB in range(nFEB):
                    for iCh in range(geometry_constants.nChannelPerFEB):
                        if spectrumTempCorrection[0] == '1':
                            hspectrum[iFEB*geometry_constants.nChannelPerFEB+iCh].Fill(tEvent.PEs[iFEB, iCh])
                        if spectrumTempCorrection[1] == '1':
                            hspectrum_corrected[iFEB*geometry_constants.nChannelPerFEB+iCh].Fill(tEvent.PEsTemperatureCorrected[iFEB, iCh])
                        if spectrumTempCorrection[2] == '1':
                            pass
                            # FIXME
                            # hspectrum_customized[iFEB*geometry_constants.nChannelPerFEB+iCh].Fill(tEvent.PEs[iFEB, iCh])
                        
            tmutiplicity = nEventInSpill if weightedAvg else 1
            temperatureFEB_avg += (tSpill.temperatureFEB * tmutiplicity)
            temperatureCMB_avg += (tSpill.temperatureCMB * tmutiplicity)
            busSiPMBias_avg += (tSpill.busSiPMBias * tmutiplicity)
            divider += tmutiplicity
        
        spill_summary_entry_dict["iSpill"][0] = iSpill
        spill_summary_entry_dict["spillIndex"][0] = tSpill.spillIndex
        spill_summary_entry_dict["eventStart"][0] = iEvent
        spill_summary_entry_dict["eventNum"][0] = nEventInSpill
        spill_summary_entry_dict["dqcRedFlag"].push_back("0x%x"%(tSpill.dqmRedFlag))
        spill_summary_entry_dict["tsEpoch"][0] = tSpill.tsEpoch
        if iSpill == 0:
            subrun_summary_entry_dict["tsEpochStart"][0] = tSpill.tsEpoch
        subrun_summary_entry_dict["tsEpochEnd"][0] = tSpill.tsEpoch
        
        for iFEB in range(nFEB):
            spill_summary_entry_dict["temperatureFEB"].push_back(tSpill.temperatureFEB[iFEB])
            for iCh in range(geometry_constants.nChannelPerFEB):
                spill_summary_entry_dict["temperatureCMB"].push_back(np.nan if (tSpill.temperatureCMB is None) else tSpill.temperatureCMB[iFEB, iCh])
            for iAFE in range(geometry_constants.nAFEPerFEB):
                spill_summary_entry_dict["busSiPMBias"].push_back(tSpill.busSiPMBias[iFEB, iAFE])
        spill_summary_tree.Fill()
        
        iEvent += nEventInSpill
        if iEvent > nEvent:
            sys.exit("Error: utils.processSubRun: %s nEvent = %i, Spill # %i is reading iEvent = %i"%(filename, nEvent, iSpill, iEvent))
                
    fFile.Close()
    
    if divider != 0:
        temperatureFEB_avg /= float(divider)
        temperatureCMB_avg /= float(divider)
        busSiPMBias_avg /= float(divider)
    
    pedestal, calibRaw, calibAdj = calibExtract_root(filename)
    # pedestal_txt, calibRaw_txt, calibAdj_txt = calibExtract_txt(filename)

    isEmpty = True if (pedestal is None) else False
    for iFEB in range(nFEB):
        subrun_summary_entry_dict["temperatureFEB"].push_back(np.nan if (divider == 0) else temperatureFEB_avg[iFEB])
        for iCh in range(geometry_constants.nChannelPerFEB):
            subrun_summary_entry_dict["pedestal"].push_back(np.nan if isEmpty else pedestal[iFEB, iCh])
            subrun_summary_entry_dict["calibRaw"].push_back(np.nan if isEmpty else calibRaw[iFEB, iCh])
            subrun_summary_entry_dict["calibAdj"].push_back(np.nan if isEmpty else calibAdj[iFEB, iCh])
            subrun_summary_entry_dict["temperatureCMB"].push_back(np.nan if (divider == 0) else temperatureCMB_avg[iFEB, iCh])
        for iAFE in range(geometry_constants.nAFEPerFEB):
            subrun_summary_entry_dict["busSiPMBias"].push_back(np.nan if (divider == 0) else busSiPMBias_avg[iFEB, iAFE])
    subrun_summary_tree.Fill()
    
    fROOT_output.cd()
    fROOT_output.cd("spectrum")
    for iFEB in range(nFEB):
        for iCh in range(geometry_constants.nChannelPerFEB):
            ii = iFEB*geometry_constants.nChannelPerFEB+iCh
            if spectrumTempCorrection[0] == '1':
                hspectrum[ii].Write()
            if spectrumTempCorrection[1] == '1':
                hspectrum_corrected[ii].Write()
            if spectrumTempCorrection[2] == '1':
                hspectrum_customized[ii].Write()
    return 
