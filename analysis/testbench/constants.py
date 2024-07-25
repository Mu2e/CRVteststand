sampleRate = 12.55 #ns

CRV_TDC_RATE = 159.324e6 # Hz
RATE = (CRV_TDC_RATE/2.0)/1.0e9 # GHZ

DEFAULT_BETA = 16.0

tDataTaking_cosmic = 120 #s
tDataTransfer_cosmic = 60 #s
timeDataSpill_cosmic = tDataTaking_cosmic + tDataTransfer_cosmic

tDataTaking_led = 20 #s
tDataTransfer_led = 220 #s
timeDataSpill_led = tDataTaking_led + tDataTransfer_led

numberOfPreSignalSamples = 60
signalRegionStart = 60
signalRegionEnd = 127
noiseThreshold = 5.0
calibTemperatureIntercept = 69.00
referenceTemperature = 25
PETemperatureIntercept = 85.12

import matplotlib as plt
cmap = plt.cm.get_cmap('tab10')
colors = [cmap(i) for i in range(cmap.N)]

rootcolors = [12, 46, 30, 38, 28, 40, 8, 9, 2, 4, 3, 6]
rootmarkers = [20, 21, 22, 23, 29, 33 ,34, 47, 43, 45, 41, 39]