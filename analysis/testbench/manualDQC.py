# changes after inspecting the DQM outputs

# locations where offsets between spillNum and spillIndex need to be reset
# format : {(run, subrun):[iSpill]}
manual_reset_spillNum_offset_dict = {(1308,  0): [26], #seemingly spillNumber skipped one
                                     (1316,  1): [ 1], #spill 0 has bad spillNumber
                                     (1321, 39): [ 1]  #seemingly spillNumber skipped one
                                    }

# dqmRedFlag that needs to be manually changed
# format : {(run, subrun):{iSpill: new DQC flag}}
manaulDQC_dict = {(1316, 1):{
                             0: (0x6|0x21)} # spill number mismatch 
                 }