from __future__ import print_function
import numpy as np

# ===== SINGLE MODULE GEOMETRY =============================
# module dimensions, in [mm]
scintillatorBarThickness = 19.8 
scintillatorBarWidth     = 51.3
layerOffset              = 42      # offset for staggered layers
gapLarge                 = 0.5     # between di-counters
gapSmall                 = 0.0875  # between counters in a di-counter
gapBetweenModules        = 3.0    
gapBetweenLayers         = [9.525, 9.525, 9.525] 
aluminumSheetThickness   = 3.175
strongBackThickness      = 12.7

channelToEdge            = 12.65
channelToChannel         = 26.
dicounterTranslationalX  = 2*scintillatorBarWidth + gapSmall + gapLarge

moduleThickness          = strongBackThickness + scintillatorBarThickness*4 + sum(gapBetweenLayers) +aluminumSheetThickness

# layers of a module:
# 12.7mm aluminum strong back
# 19.8mm scintillator
# 9.525mm aluminum absorber
# 19.8mm scintillator
# 9.525mm aluminum absorber
# 19.8mm scintillator
# 9.525mm aluminum absorber
# 19.8mm scintillator
# 3.175mm aluminum cover sheet

dicounterPerLayer = 8
nLayer            = 4
nChannelPerFEB    = 64
xFiberPosition    = [0., channelToChannel, scintillatorBarWidth + gapSmall, channelToChannel + scintillatorBarWidth + gapSmall]

# ===== MODULE FEB/CH MAP =============================
# real part is the FEB index in this module, imaginary part is the channel number of that FEB
# allows high robustness of mapping (one FEB across multiple modules, multiple FEB in one module etc.)

# numbering scheme 'tr'
# FEB0 [31, ...,  0]
#      [63, ..., 32]
# FEB1 [31, ...,  0]
#      [63, ..., 32]
map_tr = np.array([[(0 if j<2 else 1)+1j*((j%2+1)*4*dicounterPerLayer-1-i)
                    for i in range(4*dicounterPerLayer)]
                    for j in range(nLayer)])  

# numbering scheme 'tl'
# FEB0 [ 0, ..., 31]
#      [32, ..., 63]
# FEB1 [ 0, ..., 31]
#      [32, ..., 63]
map_tl = np.array([[(0 if j<2 else 1)+1j*((j%2)*4*dicounterPerLayer+i)
                    for i in range(4*dicounterPerLayer)]
                    for j in range(nLayer)])  

# ===== MODULE GEOMETRY =============================
module_geom_dict = {'CRV_T_A_HORSTD': {'type':'cis', 'map':map_tr},
                    'CRV_T_B_HORSTD': {'type':'trans', 'map':map_tl},
                    'CRV_L_A_HORSTD': {'type':'trans', 'map':map_tl}, 
                    'CRV_L_B_HORSTD': {'type':'cis', 'map':map_tr}, 
                    'CRV_R_A_HORSTD': {'type':'cis', 'map':map_tr},
                    'CRV_R_B_HORSTD': {'type':'trans', 'map':map_tl}
                    }

# ===== TEST BENCH GEOMETRY =============================
# -1 in 'FEB' indicates unplugged FEB region
# 'badChannels' is a list of [(iFEB, iCh)] 
# If both ends of a module are tested, create a geom_dict for each side and record both in setup_dict
geom_dict_crvaging001 = {
    'tag'           : 'config_001_CRV_L_A',
    'module'        : ['CRV_L_A_HORSTD'],
    'moduleOffsetX' : [0.],
    'moduleOffsetY' : [0.],
    'FEB'           : [[1,0]],
    'badChannels'   : [],
    'triggerOnlyChs': [],
    'canvasSize'    : (8,2.4)
}

geom_dict_crvaging004 = {
    'tag'           : 'config_004_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(3)],
    'FEB'           : [[-1,0], [1,2], [3,-1]],
    'badChannels'   : [],
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)] + [(3, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (8,6)
}

geom_dict_crvaging006 = {
    'tag'           : 'config_006_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[-1,0], [1,2], [3,-1], [4,5]],
    'badChannels'   : [(4, 12),(4, 13),(4, 14),(4, 15)], # CMB reports no temperature; allegedly good data
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)] + [(3, i) for i in range(nChannelPerFEB)] + 
                      [(4, i) for i in range(nChannelPerFEB)] + [(5, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (8,8)
}

# setup_dict can be used to retrieve geom dict (and number of FEBs) from the file names
setup_dict = {'crvaging-001': [geom_dict_crvaging001],
              'crvaging-004': [geom_dict_crvaging004],
              'crvaging-006': [geom_dict_crvaging006]}
