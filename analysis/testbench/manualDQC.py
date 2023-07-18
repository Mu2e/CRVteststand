# changes after inspecting the DQM outputs

# locations where offsets between spillNum and spillIndex need to be reset
# format : {(run, subrun):[iSpill]}
manual_reset_spillNum_offset_dict = {(1308,  0): [ 26], #seemingly spillNumber skipped one
                                     (1316,  1): [  1], #spill 0 has bad spillNumber
                                     (1321, 39): [  1], #seemingly spillNumber skipped one
                                     (  94,  8): [723], 
                                     ( 105,  1): [875], 
                                     ( 105,  6): [778], 
                                     (1034,  1): [ 13], 
                                     (1248,  1): [653], 
                                     (1271,  0): [ 24],
                                     (1161,  0): [434, 441, 466, 481], 
                                     (1208,  0): [34]
                                    }

# dqmRedFlag that needs to be manually changed
# format : {(run, subrun):{iSpill: new DQC flag}}
manaulDQC_dict = {(1316, 1):{
                             0: (0x6|0x21)} # spill number mismatch 
                 }

# temperature freezes. format (run, subrun):[FEB]
manual_CMB_temp_exclusion_dict = {(  94,  9):[1],
                                  (  94, 10):[1],
                                  (  94, 11):[1],
                                  ( 105,  0):[1],
                                  ( 105,  1):[1],
                                  ( 105,  2):[1],
                                  ( 105,  3):[1],
                                  ( 105,  4):[1],
                                  ( 105,  5):[1],
                                  ( 105,  6):[1],
                                  ( 105,  7):[1],
                                  ( 105,  8):[1],
                                  ( 105,  9):[1],
                                  ( 105, 10):[1], # short subrun, exclude anyways
                                  (1022,  2):[0, 1], # short subrun, exclude anyways
                                  (1247,  0):[0], 
                                  (1263,  0):[0], # short subrun, exclude anyways
                                  (1341,  5):[0],
                                  (1352,  0):[0],
                                  (1354,  0):[0],
                                  (1354,  1):[0]
                                 }

# list of subruns to exclude
exclude_subrun = [(94, 12), (105, 10), (1022, 2), (1034, 4), # short
                  (1341, 1), (1116, 0), (1137, 0), (1138, 0)] # cannot open