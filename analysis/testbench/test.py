from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import constants
import crv_event
import geometry
import geometry_constants

pdfpages = PdfPages("test.pdf")

geometry1 = geometry.testBenchGeometry(geometry_constants.geom_dict_crvaging001)
geometry4 = geometry.testBenchGeometry(geometry_constants.geom_dict_crvaging004)
geometry6 = geometry.testBenchGeometry(geometry_constants.geom_dict_crvaging001)

fig,ax = plt.subplots()
geometry1.plotGeometry(ax, '', True, True)
pdfpages.savefig(fig)
plt.close(fig)

fig,ax = plt.subplots()
geometry4.plotGeometry(ax, '', True, True)
pdfpages.savefig(fig)
plt.close(fig)

fig,ax = plt.subplots()
geometry6.plotGeometry(ax, '', True, True)
pdfpages.savefig(fig)
plt.close(fig)

fFile = TFile("/home/yongyiwu/Mu2e_CRV/Wideband/rec.mu2e.CRV_wideband_cosmics.crvaging-004.001134_000.root", "READ")
fTree = fFile.Get("run")

for iEntry in range(10):
    tEvent = crv_event.crv_event(fTree, iEntry, 0b1110)
    tEvent.plotEvent(pdfpages, [geometry4], True, True if iEntry<3 else False):

pdfpages.close()