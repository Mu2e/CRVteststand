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
import constants, geometry_constants
import crv_event, crv_spill
import filepath

##################################  GLOBAL CONTAINER   ######################################

lastTsEpoch = None

################################## UTILS FOR ROOT TREE ######################################

def type_convert_array_ROOT(arraytype):
    convert_dict = {"b":"B", "B":"b", "u":"L", "h":"S",
                    "H":"s", "i":"I", "I":"i", "l":"L",
                    "L":"l", "q":"G", "Q":"g", "f":"F",
                    "d":"D"}
    if arraytype in convert_dict:
        return convert_dict[arraytype]
    else:
        print ("WARNING: utils.type_convert_array_ROOT: array type ",arraytype," not allowed in array.array; changed to 64-bit floating point.")
        return "D"

# input branchlist is a list of 2-tuples. Each 2-tuple contains branch name and data type.
# data types are consistent with data types in array.array;
# type "v/Vxx" indicates ROOT.vector.
# e.g. a type vint8 means ROOT vector of int8_t
# function returns a container dictionary, keys are the branch names, values are arrays/vector
# to store entry data
def treeInitialization(tree, branchlist, exportBranches = False): 
    ROOT_entry_dict = {}
    branches = []
    for branch in branchlist: # branch is (name, type)
        branchname = branch[0]
        branchtype = branch[1]
        tcontainer = None
        if branchtype[0] == 'v' or branchtype[0] == 'V':
            tcontainer = ROOT.vector(branchtype[1:])()
        else: 
            tcontainer = array(branchtype,[0])
        ROOT_entry_dict.update({branchname:tcontainer})
        if branchtype[0] == 'v' or branchtype[0] == 'V':
            tbranch = tree.Branch(branchname, ROOT_entry_dict[branchname])
            branches.append(tbranch)
        else:
            roottype = type_convert_array_ROOT(branchtype)
            tbranch = tree.Branch(branchname,ROOT_entry_dict[branchname],branchname+"/"+roottype)
            branches.append(tbranch)
    if exportBranches:
        return ROOT_entry_dict, branches
    else:
        return ROOT_entry_dict

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
        dformat = '%m/%d\n%H:%M'
    formatter = DateFormatter(dformat)
    ax = plt.pyplot.plot_date(ts2datetime(args[0]), *args[1:], **kwargs)
    plt.pyplot.gca().xaxis.set_major_formatter(formatter)
    return ax

def smoothing(x, nSmooth, axis): # running average using nSmooth points
    if axis == 'x':
        x_smoothed = np.array(x[(nSmooth-1):])
    elif axis == 'y':
        x_smoothed = np.array([np.mean(x[i:(i+nSmooth)]) for i in range(len(x[(nSmooth-1):]))])
    else:
        sys.exit("Error: utils.smoothing: bad axis argument")
    return x_smoothed

def plot_dqm(filename_list, plot_dict, dqmFilter = '', nSmooth = 1, show = False, title = ';;', isRaw = False): 
    
    # plot_dict is a dictionary of {"keyword":["attribute[slicing]"]}; no omission is allowed if slicing is used
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

    for k,vl in plot_dict.iteritems():
        for v in vl:
            tAttribute = v.split('[')[0] if '[' in v else v
            if tAttribute not in plotted_attribute:
                plotted_attribute.append(tAttribute)
            nColon = v.count(':')
            if nColon == 0:
                full_plot_dict[k].append(v)
            else:
                tNameStem = '['.join(v.split('[')[:-1])
                tSlicingIndex = v.split('[')[-1].split(']')[0].split(':')
                tStart = int(tSlicingIndex[0])
                tEnd = int(tSlicingIndex[1])
                tSkip = 1
                if nColon > 1:
                    tSkip = int(tSlicingIndex[2])
                for i in range(tStart, tEnd, tSkip):
                    full_plot_dict[k].append(tNameStem+'[%i]'%(i))

    for k,vl in full_plot_dict.iteritems():
        # treat each group of files matching keyword 'k' together
        plotData[k].update({"ts":[], "dqm":[]})
        for v in vl:
            plotData[k].update({v:[]})
        for filename in filename_list:
            if k in filename.split('/')[-1] or k == '*':
                fFile = TFile(filename, "READ")
                print ("Reading file: "+filename.split('/')[-1])
                runtree = fFile.Get("run")
                spilltree = fFile.Get("spills")
                iEvent = 0
                nEvent = runtree.GetEntries()
                nSpill = spilltree.GetEntries()
                for iSpill in range(nSpill):
                    tSpill = crv_spill.crv_spill(spilltree, iSpill, isRaw)
                    if "temperatureCMB" in plotted_attribute:
                        tSpill.getTempCMB(runtree, tSpill.nFEB, iEvent, True) # populate tSpill.temperatureCMB[nFEB][nChannelPerFEB]
                        iEvent += tSpill.nEventsActual
                        if iEvent > nEvent:
                            sys.exit("Error: utils.plot_dqm: %s nEvent = %i, Spill # %i is reading iEvent = %i"%(filename, nEvent, iSpill, iEvent))
                    if dqmFilter:
                        tSpill.checkDQM() # populate tSpill.dqmRedFlag
                    plotData[k]["ts"].append(tSpill.tsEpoch)
                    if dqmFilter:
                        exec("plotData[k]['dqm'].append(tSpill.dqmRedFlag"+dqmFilter+")")
                    else:
                        plotData[k]['dqm'].append(True)
                    for v in vl:
                        try:
                            exec("plotData[k][v].append(tSpill."+v+")")
                        except:
                            try:
                                exec("plotData[k][v].append(tSpill.temp_dict['"+v.split('[')[0]+"']"+('['+'['.join(v.split('[')[1:]) if len(v.split('['))>1 else '')+")")
                            except:
                                sys.exit("Error: utils.plot_dqm: %s does not exist as a spill attribute"%(v))
                fFile.Close()
        
        plotData[k]["ts"] = smoothing(plotData[k]["ts"], nSmooth, 'x')
        plotData[k]["dqm"] = smoothing(plotData[k]["dqm"], nSmooth, 'x')
        for v in vl:
            plotData[k][v] = smoothing(plotData[k][v], nSmooth, 'y')
        
        # plot the list in the same figure
        fig = plt.figure()
        for i, v in enumerate(vl):
            if dqmFilter:
                plot_ts(plotData[k]['ts'], plotData[k][v], '.', markersize=2,  label="", rasterized=True, color='gray', dformat="%m/%d\n%H:%M")
            tdqc = plotData[k]['dqm']
            plot_ts(plotData[k]['ts'][tdqc], plotData[k][v][tdqc], '.', markersize=1.2,  label="%s"%(v), rasterized=True, color=constants.colors[i%10], dformat="%m/%d\n%H:%M")
        ax = plt.gca()
        ax.set_title(title.split(';')[0]+('' if k == '*' else k))
        ax.set_xlabel(title.split(';')[1])
        ax.set_ylabel(title.split(';')[2])
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        lgd_h = 0.045*int(len(vl)/4+1)
        fig.subplots_adjust(bottom=lgd_h+0.08)
        plt.legend(markerscale=4, ncol=4, bbox_to_anchor=(0.5,0), loc="lower center", bbox_transform=fig.transFigure, labelspacing=0.2, handletextpad=0.1, columnspacing=0.8)
        plt.tight_layout() 
    
        if show:
            plt.show()
        fig_list.append(fig)
        
    return fig_list

##################################   OTHER UTILS   ######################################

def SpillAttributeAvg(filename, attributeName, dqmFilter = '== 0x0'):
    localfilename = filename.split('/')[-1]
    print ("Loading FEB temperature from", localfilename)

    avgAttribute = None
    fFile = TFile(filename, "READ")
    runtree = fFile.Get("run")
    spilltree = fFile.Get("spills")
    iEvent = 0
    nEvent = runtree.GetEntries()
    nSpill = spilltree.GetEntries()
    nSpillUsed = 0
    for iSpill in range(nSpill):
        tSpill = crv_spill.crv_spill(spilltree, iSpill, isRaw)
        if iSpill == 0:
            nFEB = tSpill.nFEB
            if attributeName == 'temperatureFEB':
                avgAttribute = np.zeros(nFEB)
            elif attributeName == 'busSiPMBias':
                avgAttribute = np.zeros((nFEB, 8))
            elif attributeName == 'temperatureCMB':
                avgAttribute = np.zeros((nFEB, geometry_constants.nChannelPerFEB))
            else:
                avgAttribute = np.zeros(nFEB)
        if attributeName == 'temperatureCMB':
            tSpill.getTempCMB(runtree, tSpill.nFEB, iEvent, True) # populate tSpill.temperatureCMB[nFEB][nChannelPerFEB]
            iEvent += tSpill.nEventsActual
            if iEvent > nEvent:
                sys.exit("Error: utils.plot_dqm: %s nEvent = %i, Spill # %i is reading iEvent = %i"%(filename, nEvent, iSpill, iEvent))
        tSpill.checkDQM() # by default requires good DQM flag to be counted towards the average
        DQMgood = None
        exec("DQMgood = tSpill.dqmRedFlag"+dqmFilter)
        if DQMgood:
            nSpillUsed += 1
            exec("avgAttribute += tSpill."+attributeName)
    fFile.Close()
    avgAttribute = avgAttribute/float(nSpillUsed)
    return avgAttribute

def calibExtract(filename):
    localfilename = filename.split('/')[-1]
    if localfilename.split('.')[0] != 'cal':
        tpath = '/'.join(filename.split('/')[:-2])+'/'+filepath.file_calib_subdir
        localfilename = 'cal.'+'.'.join(localfilename.split('.')[1:-1])+'.txt'
        filename = tpath + localfilename
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

            pedestal = np.empty((nFEB, nCh))
            calibRaw = np.empty((nFEB, nCh))
            calibAdj = np.empty((nFEB, nCh))
            
            for i in range(len(lines)):
                ttpl = lines[i]
                pedestal[ttpl[0], ttpl[1]] = ttpl[2]
                calibRaw[ttpl[0], ttpl[1]] = ttpl[3]
                calibAdj[ttpl[0], ttpl[1]] = ttpl[4]

        return pedestal, calibRaw, calibAdj

    except:
        sys.exit("Error: utils.calibExtract: %s does not exit."%(filename))