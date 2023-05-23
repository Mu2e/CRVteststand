import os, subprocess
import geometry_constants

# preparation of environment 
os.system("setup mu2efiletools")
os.system("setup dhtools")

dataFileGroup_dict = {"recoROOT":"rec.mu2e.CRV_wideband_cosmics.CRVWB-000-001-000.root",
                      "rawROOT" :"ntd.mu2e.CRV_wideband_cosmics.CRVWB-000-001-000.root",
                      "caliPDF" :"etc.mu2e.CRV_wideband_cosmics_cali.CRVWB-000-001-000.pdf",
                      "caliTXT" :"etc.mu2e.CRV_wideband_cosmics_cali.CRVWB-000-001-000.txt",
                      "recoPDF" :"etc.mu2e.CRV_wideband_cosmics_reco.CRVWB-000-001-000.pdf",
                      "recoTXT" :"etc.mu2e.CRV_wideband_cosmics_reco.CRVWB-000-001-000.txt"
                     }

def filenameparser(fullfilename, item):
    localFilename = fullfilename.split('/')[-1]
    passVerStr = localFilename.split('.')[-3]
    passVer = int(passVerStr.split('-')[-2])
    runSubrunStr = localFilename.split('.')[-2]
    runNum = int(runSubrunStr.split('_')[0])
    subrunNum = int(runSubrunStr.split('_')[1])
    if item == 'verstr':
        return passVerStr
    elif item == 'ver':
        return passVer
    elif item == 'run':
        return runNum
    elif item == 'subrun':
        return subrunNum
    else:
        print(item, "is not a valid attribute name of filename.")
        return None

def getfullfilelist(label):
    if label not in dataFileGroup_dict:
        print(label, "is not a valid file group name.")
        print("available names are:")
        print(dataFileGroup_dict.keys())
        return []
    else:
        cmd = "mu2eDatasetFileList "+dataFileGroup_dict[label]
        buffer = subprocess.check_output(cmd, shell=True, text=True)
        file_list = [f for f in buffer.split('\n') if f.strip()]
        file_list.sort(key=lambda x:(filenameparser(x,'run'), filenameparser(x,'subrun')))
        return file_list

def findlinked(fullfilename, label):
    filelist = []
    for file in getfullfilelist(label):
        if filenameparser(fullfilename, 'run')==filenameparser(file, 'run') and filenameparser(fullfilename, 'subrun')==filenameparser(file, 'subrun'):
            filelist.append(file)
    if len(filelist)>1:
        print("WARNING: filepath.findlinked: multiple files associated with %s with type %s found."%(fullfilename,label))
        return filelist
    else:
        return filelist[0]

datatag = {
          "LED_temp_scan_bad": {# bad data due to long spills / external HD
              "type":"LED", 
              "config":"crvled-001",
              "run#":[1303, 1304, 1305, 1306]},
          "LED_temp_scan": {
              "type":"LED", 
              "config":"crvled-001",
              "run#":[1314]+ # start of 20s/220s time window
                     [1316, 1317, 1318, 1319, 1320, 1321]+ # start temp drop
                     [1328, 1329]},
          "LED_temp_scan_cosmics": {# cosmic data correstpond to temp scan
              "type":"cosmics", 
              "config":"crvaging-007",
              "run#": #[1299, 1300, 1301, 1302, 1308, 1309, 1310, 1311, 1312]+
                      [1312]+ # all others removed. Some are short run <100 spills; many dqc issues, missing CMB temp. etc.
                      [1315]+ # 22 degC
                      [1327]}, # 15 degC
          "LED_low_PE": {
              "type":"LED",
              "config":"crvled-001",
              "run#":[1361, 1362, 1363, 1364, 1365, 1366, 1367],
              "LEDbias":[0x600, 0x500, 0x580, 0x680, 0x700, 0x780, 0x800]},
          "cosmics_delay": {
              "type":"cosmics", 
              "config":"crvaging-007",
              "run#":[1370, 1354]} # 1370 with offset, 1354 as reference
}

# temperature scan data tags, exist for compatibility reasons
ds_temp_scan_cosmics = ["LED_temp_scan_cosmics"]
ds_temp_scan_led = ["LED_temp_scan"]
    
def getfilelist(taglist): # argument can be a list of keys from the above dictionary or run numbers
    filenamelist = []
    runnumlist = []
    configlist = []
    nFEBlist = []
    
    for tag in taglist:
        success = False
        if tag in datatag.keys():
            success = True
            for irun in datatag[tag]['run#']:
                if irun not in runnumlist:
                    runnumlist.append(irun)
                    configlist.append(datatag[tag]['config'])
                    nFEBlist.append(geometry_constants.setup_dict[datatag[tag]['config']][0]['nFEB'])
        else:
            # tag itself is a run number
            if isinstance(tag, str):
                if tag.isnumeric():
                    tag = int(tag)
            if isinstance(tag, int):
                for k, v in datatag.items():
                    if tag in v['run#']:
                        success = True
                        runnumlist.append(tag)
                        configlist.append(v['config'])
                        nFEBlist.append(geometry_constants.setup_dict[v['config']][0]['nFEB'])
                        break
        if not success:
            runnumlist.append(tag)
            configlist.append(None)
            nFEBlist.append(None)
            print("WARNING: filepath.getfilelist: unable to recognize data tag:", tag)
    
    fullfilelist = getfullfilelist("recoROOT")
    for runnumber in runnumlist:
        templist = []
        for file in fullfilelist:
            if filenameparser(file, 'run')==runnumber:
                templist.append(file)
        if templist:
            filenamelist += templist
        else:
            print("WARNING: filepath.getfilelist: no file with run number %i was found."%(runnumber))
    
    return filenamelist, runnumlist, configlist, nFEBlist
