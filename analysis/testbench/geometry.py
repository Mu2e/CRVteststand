from __future__ import print_function
import sys
import numpy as np
import geometry_constants
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import collections as mc

# counter geometry: top left ([0][0]) is at (0,0) in the dimension 
# viewed facing the indicated end
# "trans" configure:
#  __________________
# |_____|_____|_____|___...   
#    |_____|_____|_____|___...   
#       |_____|_____|_____|___...   
#          |_____|_____|_____|   ...   
#
# "cis" configure:
#              __________________
#          ___|_____|_____|_____|...   
#      ___|_____|_____|_____|...   
#  ___|_____|_____|_____|...   
# |_____|_____|_____|... 
# FIXME: Further 'dimensions' can be added for vertically placed modules

def dimensionGen(type, ndc = geometry_constants.dicounterPerLayer, nl = geometry_constants.nLayer):
    if type != 'trans' and type != 'cis':
        sys.exit("Error: dimensionGen: 'type' cannot be " + type)
    offsetSign = 1 if type == 'trans' else -1
    _chanY = np.array([[-geometry_constants.scintillatorBarThickness*j-sum(geometry_constants.gapBetweenLayers[:j])for i in range(4*ndc)] for j in range(nl)])
    _chanX = np.array([[geometry_constants.xFiberPosition[i%4]+geometry_constants.dicounterTranslationalX*int(i/4)+geometry_constants.layerOffset*j*offsetSign for i in range(4*ndc)] for j in range(nl)])
    dimension = {"channelX": _chanX,
                 "channelY": _chanY,
                 "cntrX1":_chanX[:, ::2]-geometry_constants.channelToEdge,
                 "cntrX2":_chanX[:, 1::2]+geometry_constants.channelToEdge,
                 "cntrY1":_chanY[:, ::2]-geometry_constants.scintillatorBarThickness/2.,
                 "cntrY2":_chanY[:, 1::2]+geometry_constants.scintillatorBarThickness/2.}
    return dimension

class testBenchGeometry:
    # object containing coordinates of channels, extrusions, mapping 

    def __init__(self, geom_dict):

        # self.tag: a string used as a label
        # self.FEBList = [] list of relevant FEBs
        # self.coordinates = {"channelX":..., "channelY":..., "cntrX1":..., "cntrX2":..., "cntrY1":..., "cntrY2":...}
        # self.mappingFEBCh = ...
        # self.idleMask = ... ; boolean type indicating no FEB connection
        # self.badChMask = ... ; boolean type, bad channels with FEB connection
        # self.analysisMask = ... ; boolean type, channels with FEB connection, designated for analysis
        # self.triggerOnlyMask = ... ; boolean type, channels with FEB connection, designated for triggering only
        # self.canvasSize = (x,y), 2-tuple in inches
        # ...: 3-D array of dimension nModule*nLayer*nFiber(or mCounter)
        
        if not isinstance(geom_dict, dict):
            sys.exit("Error: testBenchGeometry: initialization requires a dict type argument")
        if not {'module','moduleOffsetX','moduleOffsetY','FEB'}.issubset(set(geom_dict.keys())):
            sys.exit("Error: testBenchGeometry: initialization requires a dict containing {'module','moduleOffsetX','moduleOffsetY','FEB'}")
        tLen = len(geom_dict['module'])
        for item in ['module','moduleOffsetX','moduleOffsetY','FEB']:
            if len(geom_dict[item])!=tLen:
                sys.exit("Error: testBenchGeometry: initialization requires a dict with values of same length")
        self.tag = geom_dict['tag']
        self.FEBList = []
        for iList in geom_dict['FEB']:
            for iFEB in iList:
                if iFEB not in self.FEBList and iFEB!=-1:
                    self.FEBList.append(iFEB)

        dimensionNameList = ["channelX", "channelY", "cntrX1", "cntrX2", "cntrY1", "cntrY2"]
        self.coordinates = {}
        for item in dimensionNameList:
            self.coordinates.update({item: []})
        self.mappingFEBCh = []
        for index, module in enumerate(geom_dict['module']):
            if module not in geometry_constants.module_geom_dict.keys():
                sys.exit("Error: testBenchGeometry: %s is not a defined module geometry"%(module))
            tDimension = dimensionGen(geometry_constants.module_geom_dict[module]['type'])
            for item in dimensionNameList:
                tOffset = geom_dict['moduleOffsetX'][index] if 'X' in item else geom_dict['moduleOffsetY'][index]
                self.coordinates[item].append(tOffset + tDimension[item])
            tMap = geometry_constants.module_geom_dict[module]['map'].copy()
            with np.nditer(tMap, op_flags=['readwrite']) as handle:
                for x in handle:
                    if int(x.real) not in range(len(geom_dict['FEB'][index])):
                        sys.exit("Error: testBenchGeometry: 'FEB'[%i] needs to have a length of at least %i; fill in the blank with -1"%(index, 1+int(x.real)))
                    x[...] = 1j*x.imag + geom_dict['FEB'][index][int(x.real)]
            self.mappingFEBCh.append(tMap)
        
        for item in dimensionNameList:
            self.coordinates[item] = np.array(self.coordinates[item])
        self.mappingFEBCh = np.array(self.mappingFEBCh)
        
        self.idleMask = (self.mappingFEBCh.real == -1.)

        self.badChMask = (self.mappingFEBCh == None)
        if 'badChannels' in geom_dict.keys():
            if geom_dict['badChannels']:
                for tTuple in geom_dict['badChannels']:
                    self.badChMask[tuple(np.argwhere(self.mappingFEBCh==tTuple[0]+1j*tTuple[1])[0])] = True
        
        self.analysisMask = (~self.idleMask)
        if 'triggerOnlyChs' in geom_dict.keys():
            if geom_dict['triggerOnlyChs']:
                for tTuple in geom_dict['triggerOnlyChs']:
                    self.analysisMask[tuple(np.argwhere(self.mappingFEBCh==tTuple[0]+1j*tTuple[1])[0])] = False
        self.triggerOnlyMask = ((~self.idleMask)&(~self.analysisMask))

        self.canvasSize = (8,6)
        if 'canvasSize' in geom_dict.keys():
            self.canvasSize = geom_dict['canvasSize']

    def chToFiberIndex(self, iFEB, iCh):
        # converts (iFEB, iCh) to a tuple that can access the coordinates of the corresponding fiber
        return tuple(np.argwhere(self.mappingFEBCh==iFEB+1j*iCh)[0])

    def chToCntrIndex(self, iFEB, iCh):
        # converts (iFEB, iCh) to a tuple that can access the coordinates of the corresponding counter
        tTuple = tuple(np.argwhere(self.mappingFEBCh==iFEB+1j*iCh)[0])
        tTuple = tuple(tTuple[0], tTuple[1], int(tTuple[2]/2))
        return tTuple

    def plotCntrCoordGen(self, mask):
        tchtrX1_masked = self.coordinates['cntrX1'][mask] # this returns flattened array
        tchtrX2_masked = self.coordinates['cntrX2'][mask]
        tchtrY1_masked = self.coordinates['cntrY1'][mask]
        tchtrY2_masked = self.coordinates['cntrY2'][mask]
        #tchtrX1_masked = np.ma.masked_where(~mask, self.coordinates['cntrX1']) # converting masked elements to NaN causes warnings
        #tchtrX2_masked = np.ma.masked_where(~mask, self.coordinates['cntrX2'])
        #tchtrY1_masked = np.ma.masked_where(~mask, self.coordinates['cntrY1'])
        #tchtrY2_masked = np.ma.masked_where(~mask, self.coordinates['cntrY2'])
        
        lineSegList = []
        if tchtrX1_masked.any():
            with np.nditer(tchtrX1_masked, flags=['multi_index']) as it:
                for x in it:
                    lineSegList.append([(tchtrX1_masked[it.multi_index], tchtrY1_masked[it.multi_index]),
                                        (tchtrX1_masked[it.multi_index], tchtrY2_masked[it.multi_index]),
                                        (tchtrX2_masked[it.multi_index], tchtrY2_masked[it.multi_index]),
                                        (tchtrX2_masked[it.multi_index], tchtrY1_masked[it.multi_index]),
                                        (tchtrX1_masked[it.multi_index], tchtrY1_masked[it.multi_index])])
        return lineSegList

    def plotGeometry(self, ax, titleExtra = '', dotChannel = False, verbose = True):
        idleCounterMask = (self.idleMask[:, :, 1::2])&(self.idleMask[:, :, ::2])
        linesCoord_solid = self.plotCntrCoordGen(~idleCounterMask)
        linesCoord_dotted = self.plotCntrCoordGen(idleCounterMask)
        if linesCoord_solid:
            lc_solid = mc.LineCollection(linesCoord_solid, linewidths=0.5, linestyles='solid', 
                                         colors='#080808',alpha=0.8)
            ax.add_collection(lc_solid)
        if linesCoord_dotted:
            lc_dotted = mc.LineCollection(linesCoord_dotted, linewidths=0.5, linestyles='dotted', 
                                         colors='#434343',alpha=0.8)
            ax.add_collection(lc_dotted)

        if dotChannel:
            activeChMask = (~self.idleMask)&(~self.badChMask)
            plt.scatter(self.coordinates['channelX'][activeChMask],
                        self.coordinates['channelY'][activeChMask],
                        c='#080808', alpha=0.8, marker='.')

        crossMask = self.badChMask|(self.idleMask!=np.repeat(idleCounterMask, 2, axis=2))
        plt.scatter(self.coordinates['channelX'][crossMask],
                    self.coordinates['channelY'][crossMask],
                    c='#434343', alpha=0.8, marker='x')
        
        ax.set_title(self.tag+' '+titleExtra)
        ax.set_xlabel('x [mm]')
        ax.set_ylabel('y [mm]')

        if verbose:
            with np.nditer(self.mappingFEBCh, flags=['multi_index']) as it:
                for x in it:
                    if (self.mappingFEBCh[it.multi_index].imag in [0,31,32, 63]) and self.mappingFEBCh[it.multi_index].real>=0:
                        ax.text(self.coordinates['channelX'][it.multi_index],self.coordinates['channelY'][it.multi_index]-2,
                                '%i_%i'%(self.mappingFEBCh[it.multi_index].real,self.mappingFEBCh[it.multi_index].imag),
                                color='black',verticalalignment='top', horizontalalignment='center',zorder = 5,fontsize = 3)

        ax.autoscale()
        plt.axis('equal')
        return



