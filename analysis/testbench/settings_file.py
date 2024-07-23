from __future__ import print_function

import numpy as np
from array import array
import os, sys, json, glob, datetime 
import math

MV_PER_TRIM = -2 # mV/DAC
MV_PER_BULK = 20 # mV/DAC

bulkRegs = [('{:x}'.format(4*i) if i > 0 else '')+'4'+('4' if j == 0 else '5') for i in range(4) for j in range(2)]
trimRegs = [('{:x}'.format(4*i) if i > 0 else '')+'3'+'{:x}'.format(j) for i in range(4) for j in range(16)]
gainRegs = [('{:x}'.format(4*i) if i > 0 else '')+'4'+('6' if j == 0 else '7') for i in range(4) for j in range(2)]

class settingsFile:
    def __init__(self, txtfile):
        self.lines = None
        self.settingStartIndex = None
        self.settingEndIndex = None
        self.bulks = np.empty(8) * np.nan
        self.trims = np.empty(64) * np.nan
        self.gains = np.empty(8) * np.nan
        try:
            with open(txtfile, 'r') as fin:
                self.lines = fin.readlines()
        except:
            sys.exit("ERROR: settings_file: Unable to open file "+txtfile)
                
        for index, line in enumerate(self.lines):
            fields = line.strip().split(" ")

            isSetting = False
            if fields[0]=="wr" and fields[1] in (bulkRegs+trimRegs+gainRegs):
                isSetting = True
            
            if self.settingStartIndex is None and isSetting:
                self.settingStartIndex = index
            if (self.settingStartIndex is not None) and (self.settingEndIndex is None) and (not isSetting):
                self.settingEndIndex = index

            if isSetting:
                if fields[1] in bulkRegs:
                    self.bulks[bulkRegs.index(fields[1])] = int(fields[2], 16)
                elif fields[1] in trimRegs:
                    self.trims[trimRegs.index(fields[1])] = int(fields[2], 16)
                elif fields[1] in gainRegs:
                    self.gains[gainRegs.index(fields[1])] = int(fields[2], 16)

        if np.any(np.isnan(self.bulks)) or np.any(np.isnan(self.trims)) or np.any(np.isnan(self.gains)):
            print(self.bulks)
            print(self.trims)
            print(self.gains)
            sys.exit("ERROR: settings_file: Incomplete settings in "+txtfile)
        else:
            self.bulks = self.bulks.astype(int)
            self.trims = self.trims.astype(int)
            self.gains = self.gains.astype(int)
        if self.settingEndIndex is None:
            self.settingEndIndex = len(self.lines)
        # print(self.lines)

    def BalanceBiasTrimBulk(self, trimGoal = 0x800):
        newBulk = np.copy(self.bulks)
        newTrim = np.copy(self.trims)
        for iAFE in range(8):
            sumTrims = 0
            division = 8
            for i in range(8):
                chNum = iAFE * 8 + i
                if self.trims[chNum] != 0xdead:
                    sumTrims += self.trims[chNum]
                else:
                    division -= 1
            sumDiff = sumTrims - trimGoal * division
            moveBulk = int(round(sumDiff / division / (MV_PER_BULK / MV_PER_TRIM)))
            moveTrim = moveBulk * int(MV_PER_BULK / MV_PER_TRIM)

            newBulk[iAFE] += moveBulk
            for i in range(8):
                chNum = iAFE * 8 + i
                if self.trims[chNum] != 0xdead:
                    newTrim[chNum] -= moveTrim
                else:
                    newTrim[chNum] = 0xdead
        self.bulks = newBulk
        self.trims = newTrim
        return

    def Print(self):
        print("BULK| TRIM                           |GAIN")
        print("----+--------------------------------+----")
        for i in range(8):
            output = ""
            output += "{:03X}".format(self.bulks[i])+" |"
            for j in range(8):
                if self.trims[i*8+j]!=0xdead:
                    output += ("{:03X}".format(self.trims[i*8+j])+" ")
                else:
                    output += "Nan "
            output += "|"
            output += "{:03X}".format(self.gains[i])
            print(output)
        return
    
    def DumpToTxt(self, txtfile, AFEmask = None):
        if txtfile != None:
            try:
                with open(txtfile, 'w') as fout:
                    for i in range(self.settingStartIndex):
                        fout.write(self.lines[i])
                    for i, bulk in enumerate(self.bulks):
                        fout.write("wr "+bulkRegs[i]+" "+"{:x}".format(bulk)+"\n")
                    for i, trim in enumerate(self.trims):
                        if trim == 0xdead:
                            fout.write("# wr "+trimRegs[i]+" dead\n")
                        else:
                            fout.write("wr "+trimRegs[i]+" "+"{:x}".format(trim)+"\n")
                    for i, gain in enumerate(self.gains):
                        fout.write("wr "+gainRegs[i]+" "+"{:x}".format(gain)+"\n")
                    for i in range(self.settingEndIndex, len(self.lines)):
                        fout.write(self.lines[i])
            except:
                sys.exit("ERROR: settings_file: Unable to open file "+txtfile)
        else:
            for i, bulk in enumerate(self.bulks):
                if (0b1 << i) & AFEmask:
                    print("wr "+bulkRegs[i]+" "+"{:x}".format(bulk))
            for i, trim in enumerate(self.trims):
                if (0b1 << int(i/8)) & AFEmask:
                    if trim == 0xdead:
                        print("# wr "+trimRegs[i]+" dead")
                    else:
                        print("wr "+trimRegs[i]+" "+"{:x}".format(trim))
            for i, gain in enumerate(self.gains):
                if (0b1 << i) & AFEmask:
                    print("wr "+gainRegs[i]+" "+"{:x}".format(gain))                        
        return