sampleRate = 12.55 #ns

tDataTaking_cosmic = 120 #s
tDataTransfer_cosmic = 60 #s
timeDataSpill_cosmic = tDataTaking_cosmic + tDataTransfer_cosmic

tDataTaking_led = 20 #s
tDataTransfer_led = 220 #s
timeDataSpill_led = tDataTaking_led + tDataTransfer_led

import matplotlib as plt
cmap = plt.cm.get_cmap('tab10')
colors = [cmap(i) for i in range(cmap.N)]

rootcolors = [12, 46, 30, 38, 28, 40, 8, 9, 2, 4, 3, 6]