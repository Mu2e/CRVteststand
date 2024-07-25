from __future__ import print_function

import numpy as np
from array import array
import os, sys, json, glob, datetime 
import math
import copy

MV_PER_TRIM = -2 # mV/DAC
MV_PER_BULK = 20 # mV/DAC

CMBTempCoeff = -55.4 # mV/C, overvoltage vs CMB temp
FEBTempCoeff = 4.09 # mV/C, overvoltage vs FEB temp

class settingsFEB:
    def __init__(self, strBulk, strTrim, strGain, strCMBTemp, strFEBTemp, CMBTemp = None, FEBTemp = None):
        self.bulks = [float(int(k, 16)) for k in strBulk.split()]
        self.trims = [float(int(k, 16)) for k in strTrim.split()]
        self.gains = [float(int(k, 16)) for k in strGain.split()]
        self.CMBTemp = None
        if strCMBTemp:
            self.CMBTemp = [float(k) for k in strCMBTemp.split()]
        if CMBTemp:
            self.CMBTemp = [CMBTemp for i in range(16)] # each CMB
        if strFEBTemp:
            self.FEBTemp = float(strFEBTemp)
        if FEBTemp:
            self.FEBTemp = FEBTemp

        if len(self.bulks) != 8 or len(self.trims) != 64 or len(self.gains) != 8 or len(self.CMBTemp) != 16:
            print("Error dump:")
            print(strBulk, strTrim, strGain, strCMBTemp, strFEBTemp, sep = "\n")
            sys.exit("ERROR: settingsFEB: Unable to parse the strings.")

    def BalanceBiasTrimBulk(self, trimGoal = 0x800):
        newBulk = np.copy(self.bulks)
        newTrim = np.copy(self.trims)
        for iAFE in range(8):
            sumTrims = 0
            for i in range(8):
                chNum = iAFE * 8 + i
                sumTrims += self.trims[chNum]
            sumDiff = sumTrims - trimGoal * 8
            moveBulk = int(round(sumDiff / 8 / (MV_PER_BULK / MV_PER_TRIM)))
            moveTrim = moveBulk * int(MV_PER_BULK / MV_PER_TRIM)

            newBulk[iAFE] += moveBulk
            for i in range(8):
                chNum = iAFE * 8 + i
                newTrim[chNum] -= moveTrim
        self.bulks = newBulk
        self.trims = newTrim
        return

    def Print(self):
        print("    FEB TEMP = %.2f"% self.FEBTemp)
        print("BULK| TRIM                           |GAIN|CMB TEMP")
        print("----+--------------------------------+----+-------------")
        for i in range(8):
            output = ""
            output += "{:03X}".format(self.bulks[i])+" |"
            for j in range(8):
                output += ("{:03X}".format(self.trims[i*8+j])+" ")
            output += "|"
            output += "{:03X} |".format(self.gains[i])
            output += "%6.2f %6.2f"%(self.CMBTemp[i*2], self.CMBTemp[i*2+1])
            print(output)
        return
    
    def AdjustBulkByTemp(self, newCMBTempStr, newFEBTempStr):
        newSettingsFEB = copy.deepcopy(self)
        newCMBTemp = [float(k) for k in newCMBTempStr.split()]
        if len(newCMBTemp) != 16:
            print("Error dump:")
            print(newCMBTempStr)
            sys.exit("ERROR: settingsFEB: AdjustBulkByTemp: CMB temperature string has wrong length.")
        newFEBTemp = float(newFEBTempStr)
        for i in range(8):
            meanCMBTempdiff = (newCMBTemp[i*2]+newCMBTemp[i*2+1]-self.CMBTemp[i*2]-self.CMBTemp[i*2+1])/2.
            mVToShift = - CMBTempCoeff * meanCMBTempdiff - FEBTempCoeff * (newFEBTemp - self.FEBTemp)
            # bulkToShift = int(round(mVToShift / MV_PER_BULK))
            bulkToShift = mVToShift / float(MV_PER_BULK)
            newSettingsFEB.bulks[i] += bulkToShift
        # trim, gain untouched
        newSettingsFEB.FEBTemp = newFEBTemp
        newSettingsFEB.CMBTemp = newCMBTemp
        return newSettingsFEB

    def AdjustBulkByTemp_SiDetRun2(self, newCMBTempStr, newFEBTempStr):
        # bulk 0 only trace cmb0; no FEB tracing
        newSettingsFEB = copy.deepcopy(self)
        newCMBTemp = [float(k) for k in newCMBTempStr.split()]
        if len(newCMBTemp) != 16:
            print("Error dump:")
            print(newCMBTempStr)
            sys.exit("ERROR: settingsFEB: AdjustBulkByTemp: CMB temperature string has wrong length.")
        newFEBTemp = float(newFEBTempStr)
        for i in range(8):
            if i == 0:
                meanCMBTempdiff = newCMBTemp[0]-self.CMBTemp[0]
            else:
                meanCMBTempdiff = (newCMBTemp[i*2]+newCMBTemp[i*2+1]-self.CMBTemp[i*2]-self.CMBTemp[i*2+1])/2.
            mVToShift = - CMBTempCoeff * meanCMBTempdiff #- FEBTempCoeff * (newFEBTemp - self.FEBTemp)
            # bulkToShift = int(round(mVToShift / MV_PER_BULK))
            bulkToShift = mVToShift / float(MV_PER_BULK)
            newSettingsFEB.bulks[i] += bulkToShift
        # trim, gain untouched
        newSettingsFEB.FEBTemp = newFEBTemp
        newSettingsFEB.CMBTemp = newCMBTemp
        return newSettingsFEB
