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
moduleWidth              = dicounterTranslationalX*8+layerOffset*3

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
nCMBPerFEB        = 16
nAFEPerFEB        = 8
nChannelPerCMB    = 4
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
    'nFEB'          : 2,
    'badChannels'   : [],
    'triggerOnlyChs': [],
    'canvasSize'    : (8,2.4)
}

geom_dict_crvaging003 = {
    'tag'           : 'config_003_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(2)],
    'FEB'           : [[-1,0], [1,2]],
    'nFEB'          : 3,
    'badChannels'   : [],
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (8,4)
}

geom_dict_crvaging004 = {
    'tag'           : 'config_004_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(3)],
    'FEB'           : [[-1,0], [1,2], [3,-1]],
    'nFEB'          : 4,
    'badChannels'   : [],
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)] + [(3, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (8,6)
}

geom_dict_crvaging005 = {
    'tag'           : 'config_005_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[-1,0], [1,2], [3,-1], [-1,4]],
    'nFEB'          : 5,
    'badChannels'   : [],
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)] + [(3, i) for i in range(nChannelPerFEB)] + 
                      [(4, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (8,8)
}

geom_dict_crvaging006 = {
    'tag'           : 'config_006_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[-1,0], [1,2], [3,-1], [4,5]],
    'nFEB'          : 6,
    'badChannels'   : [(4, 12),(4, 13),(4, 14),(4, 15)], # CMB reports no temperature; allegedly good data
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)] + [(3, i) for i in range(nChannelPerFEB)] + 
                      [(4, i) for i in range(nChannelPerFEB)] + [(5, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (8,8)
}

geom_dict_crvaging007 = {
    'tag'           : 'config_007_CRV_L_A',
    'module'        : ['CRV_L_A_HORSTD', 'CRV_L_A_HORSTD', 'CRV_L_A_HORSTD', 'CRV_L_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[1,0], [-1,2], [3,-1], [4,5]],
    'nFEB'          : 6,
    'badChannels'   : [], # [(2, 4),(2, 5),(2, 6),(2, 7),(3, 0),(3, 1),(3, 2),(3, 3)], 
    'triggerOnlyChs': [(2, i) for i in range(nChannelPerFEB)] + [(3, i) for i in range(nChannelPerFEB)] + 
                      [(4, i) for i in range(nChannelPerFEB)] + [(5, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (8,8)
}

geom_dict_crvaging008 = {
    'tag'           : 'config_008_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[0, -1], [1,2], [3,4], [-1,5]],
    'nFEB'          : 6,
    'badChannels'   : [],
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)] + [(3, i) for i in range(nChannelPerFEB)] + 
                      [(4, i) for i in range(nChannelPerFEB)] + [(5, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (8,8)
}

geom_dict_crvaging009 = {
    'tag'           : 'config_009_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[0, -1], [1,2], [3,4], [-1,5]],
    'nFEB'          : 6,
    'badChannels'   : [],
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)] + [(3, i) for i in range(nChannelPerFEB)] + 
                      [(4, i) for i in range(nChannelPerFEB)] + [(5, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (8,8)
}

geom_dict_crvaging010 = {
    'tag'           : 'config_010_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[0, -1], [1,2], [3,4], [-1,5]],
    'nFEB'          : 6,
    'badChannels'   : [],
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)] + [(3, i) for i in range(nChannelPerFEB)] + 
                      [(4, i) for i in range(nChannelPerFEB)] + [(5, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (8,8)
}

geom_dict_crvaging011 = {
    'tag'           : 'config_011_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[0,6], [1,2], [3,4], [7,5]],
    'nFEB'          : 8,
    'badChannels'   : [],
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)] + [(3, i) for i in range(nChannelPerFEB)] + 
                      [(4, i) for i in range(nChannelPerFEB)] + [(5, i) for i in range(nChannelPerFEB)] + 
                      [(6, i) for i in range(nChannelPerFEB)] + [(7, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (8,8)
}

geom_dict_crvaging016 = {
    'tag'           : 'config_016_CRV_T_A/B',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_B_HORSTD', 'CRV_T_B_HORSTD'],
    'moduleOffsetX' : [0., 0., 0.+moduleWidth*1.2, 0.+moduleWidth*1.2],
    'moduleOffsetY' : [-(i%2)*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[0,1], [2,3], [4,5], [6,7]],
    'nFEB'          : 8,
    'badChannels'   : [],
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)] + [(1, i) for i in range(nChannelPerFEB)] + 
                      [(4, i) for i in range(nChannelPerFEB)] + [(5, i) for i in range(nChannelPerFEB)] + 
                      [(6, i) for i in range(nChannelPerFEB)] + [(7, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (16,5)
}

geom_dict_crvaging019 = { # not real config, placeholder
    'tag'           : 'config_019_CRV_T_A',
    'module'        : ['CRV_L_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_B_HORSTD', 'CRV_L_B_HORSTD'],
    'moduleOffsetX' : [0., 0., 0.+moduleWidth*1.2, 0.+moduleWidth*1.2],
    'moduleOffsetY' : [-(i%2)*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[0,1], [2,3], [4,5], [6,7]],
    'nFEB'          : 8,
    'badChannels'   : [],
    'triggerOnlyChs': [(0, i) for i in range(nChannelPerFEB)] + [(1, i) for i in range(nChannelPerFEB)] + 
                      [(4, i) for i in range(nChannelPerFEB)] + [(5, i) for i in range(nChannelPerFEB)] + 
                      [(6, i) for i in range(nChannelPerFEB)] + [(7, i) for i in range(nChannelPerFEB)],
                      # 2-tuple list if channels are designated in the triggering layers 
    'canvasSize'    : (16,5)
}

geom_dict_crvled001 = {
    'tag'           : 'config_LED001_CRV_L_A',
    'module'        : ['CRV_L_A_HORSTD'],
    'moduleOffsetX' : [0.],
    'moduleOffsetY' : [0.],
    'FEB'           : [[1,0]],
    'nFEB'          : 2,
    'badChannels'   : [],
    'triggerOnlyChs': [],
    'canvasSize'    : (8,2.4)
}


geom_dict_crvled002 = {
    'tag'           : 'config_LED002_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0.],
    'moduleOffsetY' : [0.],
    'FEB'           : [[0,1]],
    'nFEB'          : 2,
    'badChannels'   : [],
    'triggerOnlyChs': [],
    'canvasSize'    : (8,2.4)
}

geom_dict_crvled003 = {
    'tag'           : 'config_LED003_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[0, -1], [1,2], [3,4], [-1,5]],
    'nFEB'          : 6,
    'badChannels'   : [],
    'triggerOnlyChs': [],
    'canvasSize'    : (8,8)
}

geom_dict_crvled004 = {
    'tag'           : 'config_LED004',
    'module'        : [],
    'moduleOffsetX' : [],
    'moduleOffsetY' : [],
    'FEB'           : [],
    'nFEB'          : 1,
    'badChannels'   : [],
    'triggerOnlyChs': [],
    'canvasSize'    : (8,2)
}

geom_dict_crvled005 = {
    'tag'           : 'config_LED005_CRV_T_A',
    'module'        : ['CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD', 'CRV_T_A_HORSTD'],
    'moduleOffsetX' : [0., 0., 0., 0.],
    'moduleOffsetY' : [-i*(moduleThickness+4*25.4) for i in range(4)],
    'FEB'           : [[0,6], [1,2], [3,4], [7,5]],
    'nFEB'          : 8,
    'badChannels'   : [],
    'triggerOnlyChs': [],
    'canvasSize'    : (8,8)
}

# setup_dict can be used to retrieve geom dict (and number of FEBs) from the file names
setup_dict = {'crvaging-001': [geom_dict_crvaging001],
              'crvaging-003': [geom_dict_crvaging003],
              'crvaging-004': [geom_dict_crvaging004],
              'crvaging-005': [geom_dict_crvaging005],
              'crvaging-006': [geom_dict_crvaging006],
              'crvaging-007': [geom_dict_crvaging007],
              'crvaging-008': [geom_dict_crvaging008],
              'crvaging-009': [geom_dict_crvaging009],
              'crvaging-010': [geom_dict_crvaging010],
              'crvaging-011': [geom_dict_crvaging011],
              'crvaging-016': [geom_dict_crvaging016],
              'crvaging-019': [geom_dict_crvaging019],
              'crvled-001'  : [geom_dict_crvled001],
              'crvled-002'  : [geom_dict_crvled002],
              'crvled-003'  : [geom_dict_crvled003],
              'crvled-004'  : [geom_dict_crvled004],
              'crvled-005'  : [geom_dict_crvled005]}
