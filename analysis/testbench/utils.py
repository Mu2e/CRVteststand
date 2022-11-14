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
import constants
import crv_event, crv_spill

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

def smoothing(x, y, nSmooth): # running average using nSmooth points
    x_smoothed = x[(nSmooth-1):]
    y_smoothed = np.array([np.mean(y[i:(i+nSmooth)]) for i in range(len(x_smoothed))])
    return x_smoothed, y_smoothed

def plot_dqm(plot_dict, nFEBFile, isRaw = False): 
    
    # plot_dict is a dictionary of {keyword:attribute[slicing]}
    # FIXME: Smoothing??? filter???

    tslist = []
    attribute_full_list = []

    attribute_list = [[] for i in range(len(attribute_list))]
    for filename in filename_list:
        fFile = TFile(filename, "READ")
        runtree = fFile.Get("run")
        spilltree = fFile.Get("spills")
        iEvent = 0
        nEvent = runtree.GetEntries()
        nSpill = spilltree.GetEntries()
        for iSpill in range(nSpill):
            tSpill = crv_spill.crv_spill(spilltree, iSpill, isRaw)
            tSpill.getTempCMB(runtree, nFEBFile, iEvent, True)
            iEvent += tSpill.nEventsActual
            if iEvent > nEvent:
                sys.exit("Error: utils.plot_dqm: %s nEvent = %i, Spill # %i is reading iEvent = %i"%(filename, nEvent, iSpill, iEvent))
        fFile.Close()

    
    #FIXME: Add display for removed spills