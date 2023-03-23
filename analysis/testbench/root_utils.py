from __future__ import print_function

import sys
import re
import os
import json
import datetime
from datetime import datetime
import threading
import numpy as np
import ROOT
from ROOT import TFile, TCanvas, TH1F, TH2F, TF1, TMath, TGraph
from ROOT import gStyle, gROOT, gDirectory
from array import array
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import math

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
