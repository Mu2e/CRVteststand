from __future__ import print_function
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import ROOT
from ROOT import gROOT, gStyle, gDirectory, gPad

from array import array

import constants
import crv_event
import crv_spill
import utils
import geometry
import geometry_constants
import filepath

import pandas as pd

def PortToDF(fileList, additional_dict = None, verbose = False):
    df = pd.DataFrame()

    for index, filename in enumerate(fileList):
        if verbose:
            print("Opening "+filename)
        runNum = filepath.filenameparser(filename, 'run') 
        subrunNum = filepath.filenameparser(filename, 'subrun') 

        fFile = ROOT.TFile(filename, "READ")
        runSummarytree = fFile.Get("runSummary")
        spilltree = fFile.Get("spills")
        spilltree.GetEntry(0)
        nFEB = spilltree.spill_number_of_febs
        runSummarytree.GetEntry(0)
        
        pedestal = np.reshape(np.array(runSummarytree.pedestals, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        calibRaw = np.reshape(np.array(runSummarytree.calibConstants, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        calibAdj = np.reshape(np.array(runSummarytree.calibConstantsTemperatureCorrected, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        FEBtemp = np.array(runSummarytree.febTemperaturesAvg, dtype=np.float32)
        CMBtemp = np.reshape(np.array(runSummarytree.meanTemperatures, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        biasV = np.reshape(np.array(runSummarytree.biasVoltagesAvg, dtype=np.float32), (nFEB, int(geometry_constants.nChannelPerFEB/8)))
        PEs = np.reshape(np.array(runSummarytree.PEs, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        PEsCorrected = np.reshape(np.array(runSummarytree.PEsTemperatureCorrected, dtype=np.float32), (nFEB, geometry_constants.nChannelPerFEB))
        
        for iFEB in range(nFEB):
            for iCh in range(geometry_constants.nChannelPerFEB):
                df_ = pd.DataFrame()
                df_['run'] = [runNum]
                df_['subrun'] = [subrunNum]
                df_['FEB'] = [iFEB]
                df_['ch'] = [iCh]
                df_['pedestal'] = [pedestal[iFEB][iCh]]
                df_['calibRaw'] = [calibRaw[iFEB][iCh]]
                df_['calibAdj'] = [calibAdj[iFEB][iCh]]
                df_['FEBtemp'] = [FEBtemp[iFEB] if nFEB>1 else FEBtemp]
                df_['CMBtemp'] = [CMBtemp[iFEB][iCh]]
                df_['biasV'] = [biasV[iFEB][int(iCh/8)]]
                df_['PEs'] = [PEs[iFEB][iCh]]
                df_['PEsCorrected'] = [PEsCorrected[iFEB][iCh]]

                if additional_dict:
                    for item in additional_dict.keys():
                        if isinstance(additional_dict[item][index], list): # redquire each run has a 2-D list representing the values of the item
                            df_[item] = [additional_dict[item][index][iFEB][iCh]]
                        else:
                            df_[item] = [additional_dict[item][index]] # each run corresponds to a single value
                
                df = pd.concat([df,df_], ignore_index=True)

        fFile.Close()

    return df