from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import ROOT
from ROOT import TCanvas, TH1F, TH2F, TF1, TMath, TGraph, TFile, TSpectrum, TPaveText, TMultiGraph, TGraphErrors, TLine
from ROOT import gStyle, gROOT, gDirectory, gPad

import constants
import crv_event
import geometry
import geometry_constants

pdfpages = PdfPages("test.pdf")

geometry1 = geometry.testBenchGeometry(geometry_constants.geom_dict_crvaging001)
geometry4 = geometry.testBenchGeometry(geometry_constants.geom_dict_crvaging004)
geometry6 = geometry.testBenchGeometry(geometry_constants.geom_dict_crvaging006)
geometry7 = geometry.testBenchGeometry(geometry_constants.geom_dict_crvaging007)

fig,ax = plt.subplots(1,1,figsize=geometry1.canvasSize)
geometry1.plotGeometry(ax, '', True, True)
pdfpages.savefig(fig)
plt.close(fig)

fig,ax = plt.subplots(1,1,figsize=geometry4.canvasSize)
geometry4.plotGeometry(ax, '', True, True)
pdfpages.savefig(fig)
plt.close(fig)

fig,ax = plt.subplots(1,1,figsize=geometry6.canvasSize)
geometry6.plotGeometry(ax, '', True, True)
pdfpages.savefig(fig)
plt.close(fig)

fig,ax = plt.subplots(1,1,figsize=geometry7.canvasSize)
geometry7.plotGeometry(ax, '', True, True)
pdfpages.savefig(fig)
plt.close(fig)

print("Done plotting geometry")

fFile = TFile("/home/yongyiwu/Mu2e_CRV/Wideband/rec.mu2e.CRV_wideband_cosmics.crvaging-004.001134_000.root", "READ")
fTree = fFile.Get("run")

for iEntry in range(30):
    tEvent = crv_event.crv_event(fTree, iEntry, 0b1110)
    tEvent.plotEvent(pdfpages, [geometry4], True, True if iEntry<1 else False)

pdfpages.close()
