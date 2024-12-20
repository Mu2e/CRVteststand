# OBSOLETE METHODS TO ACCESS FILES SCATTERED AROUND THE SCRATCH AREA

import glob
import geometry_constants

file_parsed_subdir = "crvparsed/"
file_calib_subdir = "crvcalib/"
file_reco_subdir = "crvreco/"

data_location = {
                # LED
                "1303-1306"      : {# bad data due to long spills / external HD
                    "type":"LED", 
                    "config":"crvled-001",
                    "run#":[1303, 1304, 1305, 1306],
                    "location": "/pnfs/mu2e/scratch/outstage/ehrlich/wideband24LED/"}, 
                "1314"           : {# start of 20s/220s time window
                    "type":"LED", 
                    "config":"crvled-001",
                    "run#":[1314],
                    "location": "/pnfs/mu2e/scratch/outstage/ehrlich/wideband25LED/"}, 
                "1317"           : {# start temp drop
                    "type":"LED", 
                    "config":"crvled-001",
                    "run#":[1317],
                    "location": "/pnfs/mu2e/scratch/users/ehrlich/wideband26LED/"}, 
                "1316,1318-1321" : {
                    "type":"LED", 
                    "config":"crvled-001",
                    "run#":[1316, 1318, 1319, 1320, 1321],
                    "location": "/pnfs/mu2e/scratch/users/ehrlich/wideband27LED/"}, 
                "1328-1329"      : {
                    "type":"LED", 
                    "config":"crvled-001",
                    "run#":[1328, 1329],
                    "location": "/pnfs/mu2e/scratch/users/ehrlich/wideband28LED/"},
                # cosmic
                "1299-1312"      : {# cosmic data correstpond to temp scan
                    "type":"cosmics", 
                    "config":"crvaging-007",
                    # "run#":[1299, 1300, 1301, 1302, 1308, 1309, 1310, 1311, 1312], 
                    "run#":[1312], # all others removed. Some are short run <100 spills; many dqc issues, missing CMB temp. etc.
                    "location": "/pnfs/mu2e/scratch/outstage/ehrlich/wideband24/"},
                "1315"           : {# 22 degC
                    "type":"cosmics", 
                    "config":"crvaging-007",
                    "run#":[1315],
                    "location": "/pnfs/mu2e/scratch/outstage/ehrlich/wideband25/"},
                "1327"           : {# 15 degC
                    "type":"cosmics", 
                    "config":"crvaging-007",
                    "run#":[1327],
                    "location": "/pnfs/mu2e/scratch/users/ehrlich/wideband28/"}
                }

# temperature scan data tags
ds_temp_scan_cosmics = ["1299-1312", "1315", "1327"]
ds_temp_scan_led = ["1314", "1317", "1316,1318-1321", "1328-1329"]

def getfilelist(taglist):
    filenamelist = []
    runnumlist = []
    configlist = []
    nFEBlist = []
    for tag in taglist:
        rootdir = data_location[tag]["location"]+file_reco_subdir
        tnamelsit = glob.glob(rootdir+"*.root")
        tnamelsit = list(filter(lambda x: int(x.split('/')[-1].split('.')[-2].split('_')[0]) in data_location[tag]["run#"], tnamelsit))
        tnamelsit.sort(key=lambda x:(int(x.split('/')[-1].split('.')[-2].split('_')[0]), int(x.split('/')[-1].split('.')[-2].split('_')[1])))
        filenamelist += tnamelsit
        runnumlist += data_location[tag]["run#"]
        for irun in data_location[tag]["run#"]:
            configlist.append(data_location[tag]["config"])
            nFEBlist.append(geometry_constants.setup_dict[data_location[tag]["config"]][0]["nFEB"])
    return filenamelist, runnumlist, configlist, nFEBlist

