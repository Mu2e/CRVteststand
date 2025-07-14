import os, subprocess
import geometry_constants
import numpy as np

# preparation of environment 
os.system("setup mu2efiletools")
os.system("setup dhtools")

dataFileGroup_dict = {"recoROOT":["rec.mu2e.CRV_wideband_cosmics","CRVWB-000-001-000","root"],
                      "rawROOT" :["ntd.mu2e.CRV_wideband_cosmics","CRVWB-000-001-000","root"],
                      "caliPDF" :["etc.mu2e.CRV_wideband_cosmics_cali","CRVWB-000-001-000","pdf"],
                      "caliTXT" :["etc.mu2e.CRV_wideband_cosmics_cali","CRVWB-000-001-000","txt"],
                      "recoPDF" :["etc.mu2e.CRV_wideband_cosmics_reco","CRVWB-000-001-000","pdf"],
                      "recoTXT" :["etc.mu2e.CRV_wideband_cosmics_reco","CRVWB-000-001-000","txt"]
                     }

def filenameparser(fullfilename, item):
    localFilename = fullfilename.split('/')[-1]
    # print(fullfilename)
    passVerStr = localFilename.split('.')[-3]
    passVer = int(passVerStr.split('-')[-2])
    runSubrunStr = localFilename.split('.')[-2]
    runNum = int(runSubrunStr.split('_')[0])
    runSubrunSplit = runSubrunStr.split('_')
    subrunNum = int(runSubrunSplit[1]) if len(runSubrunSplit)>1 else None
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

def getfullfilelist(label, passVer = 1):
    if label not in dataFileGroup_dict:
        print(label, "is not a valid file group name.")
        print("available names are:")
        print(dataFileGroup_dict.keys())
        return []
    else:
        namestrings = dataFileGroup_dict[label]
        if int(namestrings[1].split('-')[-2])!=passVer:
            namestrings[1] = "CRVWB-000-%03i-000"%(passVer)
        cmd = "mu2eDatasetFileList "+'.'.join(namestrings)
        buffer = subprocess.check_output(cmd, shell=True, text=True)
        file_list = [f for f in buffer.split('\n') if f.strip()]
        file_list.sort(key=lambda x:(filenameparser(x,'run'), (filenameparser(x,'subrun') if filenameparser(x,'subrun') else -1)))
        # print(file_list)
        return file_list

def findlinked(fullfilename, label, passVer = 1):
    filelist = []
    for file in getfullfilelist(label, passVer):
        if filenameparser(fullfilename, 'run')==filenameparser(file, 'run') and filenameparser(fullfilename, 'subrun')==filenameparser(file, 'subrun'):
            filelist.append(file)
    if len(filelist)>1:
        print("WARNING: filepath.findlinked: multiple files associated with %s with type %s found."%(fullfilename,label))
        return filelist
    else:
        return filelist[0]

datatag = { #updated to 1420. #FIXME: tag previous runs
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
                      [1327]+ # 15 degC
                      [1341, 1354, 1370]},
          "LED_low_PE": {
              "type":"LED",
              "config":"crvled-001",
              "run#":[1361, 1364, 1365, 1366, 1367],
              "LEDbias":[0x600, 0x680, 0x700, 0x780, 0x800]}, # 1362 1363 too low no light
              # "run#":[1361, 1362, 1363, 1364, 1365, 1366, 1367],
              # "LEDbias":[0x600, 0x500, 0x580, 0x680, 0x700, 0x780, 0x800]},
          "cosmics_delay": {
              "type":"cosmics", 
              "config":"crvaging-007",
              "run#":[1370, 1354]}, # 1370 with offset, 1354 as reference
          "cosmics_delay_phase_aligned": {
              "type":"cosmics", 
              "config":"crvaging-006",
              "run#":[1379, 1380]}, # 1380 with offset, 1379 as reference
          # stack 127
          "crvaging001": {
              "type":"cosmics", 
              "config":"crvaging-001",
              "run#":[66, 94, 105, 119, 1010, 1020, 1021, 1022, 
                      1031, 1033, 1034, 1035]}, 
          "crvaging007": {
              "type":"cosmics", 
              "config":"crvaging-007",
              "run#":[1243, 1244, 1245, 1246, 1247, 1248, 1251, # 1242, 1249 short
                      1262, 1264, 1276, 1279, 1280, 1281, # 1263, 1265, 1271, 1277 short
                      1284, 1285, 1286, 1287, 1288, #
                      1312, 1315, 1327, 1341, 1354, 1370]}, # 1340, 1352 too short
                     #[1299, 1300, 1301, 1302, 1308, 1309, 1310, 1311, 1312] many dqc issues, missing CMB temp. etc.
          # stack 168
          "crvaging003": {
              "type":"cosmics", 
              "config":"crvaging-003",
              "run#":[1059, 1091, 1116]}, # 1124 short
          "crvaging004": {
              "type":"cosmics", 
              "config":"crvaging-004",
              "run#":[1080, 1081, 1082]}, # 1079, 1133, 1134 short
          "crvaging005": {
              "type":"cosmics", 
              "config":"crvaging-005",
              "run#":[1137, 1138, 1146,   # 1141, 1142 short
                      1148, 1149, 1150, 1152, 1154]}, #  
          "crvaging006": {
              "type":"cosmics", 
              "config":"crvaging-006",
              "run#":[1159, 1160, 1161, 1167, 1168, 1169, 1170, 1171, #
                      1177, 1181, # 1172, 1173, 1176, 1178, 1182, 1208 short
                      1212, 1217, 1218, 1219, 1240, 1241, # 1209, 1220 short
                      1375, 1376, 1379, 1380, 1400, 1407]},  
          # stack 169
          "crvaging008": {
              "type":"cosmics", 
              "config":"crvaging-008",
              "run#":[1421, 1619, 1631, 1632]}, #1419, 1420 removed for being too short
          "crvaging009": {
              "type":"cosmics", 
              "config":"crvaging-009",
              "run#":[1634]}, #1633 removed for being too short
          "trim_scan": {
              "type":"cosmics", 
              "config":"crvaging-006",
              "run#":[1376, 1377, 1378],
              "trim":[   0, -500, +500]},
          "trim_scan_led_like": {
              "type":"led", 
              "config":"crvled-002",
              "run#":[1381, 1382, 1383, 1384, 1385, 1386, 1387, 1388, 1390,
                      1391, 1392, 1393, 1394, 1395, 1396],
              "trim":[   0, -500, +500, +375, -375, -250, +250, +125, -125,
                         0,    0,    0,    0,    0,    0],
              "bulk":[   0,    0,    0,    0,    0,    0,    0,    0,    0,
                       -50,  +50,  +30,  -30,  -10,  +10]},
          "trim_scan_stack_169": {
              "type":"led", 
              "config":"crvled-003",
              "run#":[1414, 1415, 1416, 1417, 1418],
              "trim":[   0, -500, +500, +250, -250]},
          "SiDetTempScan": {
              "type":"led",
              "config":"crvled-004", 
              "run#":[1433, 1434, 1435, 1436, 1437, 1438, 1439, 1440, 1441, 1442, 
                      1443, 1444, 1445, 1446, 1447, 1448, 1449, 1450, 1451, 1452, 1453, 1454, 1455,
                      1456, 1457, 1458, 1459, 1460, 1461, 1462, 1463, 1464, 1465, 1466, 1467, 1468,
                      1469, 1470, 1471, 1472, 1473, 1474, 1475, 1476, 1477, 1478, 1479, 1480, 1481,
                      1482, 1483, 1484, 1485, 1486, 1487, 1488, 1489, 1490, 1493, 1494, 1495, 1496,
                      1497, 1498, 1499, 1500, 1501, 1502, 1503, 1504, 1505, 1506,
                      # 1507, 1508, 1509, 1510, 1511, 1512, 1513, 1514, 1515, 
                      1516, 1517, 1518, 1519, 1520, 1521, 1522, 1523, 1524, 1525, 1526, 1527, 
                      # 1528, # corrupted ignore for now
                      1532, 1534, 1535, 1536, 1537, 1538, 1539, 1540, 1541, 1542, 1543, 1544, # 1546,
                      1547, 1548, 1549, 1550, 1551, 1552, 1553, 1554, 1555, 1556, 1557, 1558, 1559,
                      1560, 1561, 1562, 1563, 1564, 1565, 1566, 1567, 1568, 1569, 1570, 1571, 1572,
                      1573, 1574, 1575, 1576, 1577, 1578, 1579, 1580, 1581, 1582, 1583, 1584, 1585,
                      1586, 1587, 1588, 1589, 1590, 1591, 1592, 1593, 1594, 1595, 1596, 1597, 1598],
              "tempN":[25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 
                       30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 
                       35., 35., 35., 35., 35., 35., 35., 35., 35., 35., 35., 35., 35., 
                       40., 40., 40., 40., 40., 40., 40., 40., 40., 40., 40., 40., 40., 
                       20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 
                       40., 40., 40., 40., 40., 40., 40., 40., 40., 40., 
                       # 15., 15., 15., 15., 15., 15., 15., 15., 15., 
                       15., 15., 15., 15., 15., 15., 15., 15., 15., 15., 15., 15., # 15., 
                       10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., # 10., 
                        5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5., 
                        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 
                       -5., -5., -5., -5., -5., -5., -5., -5., -5., -5., -5., -5., -5., 
                       -10., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10., -10.],
              "FEBBias":[53.5, 53.0, 54.0, 54.5, 55.0, 55.5, 56.0, 56.5, 57.0, 57.5,
                         53.8, 54.3, 54.8, 55.3, 55.8, 56.3, 56.8, 57.3, 57.8, 53.3, 53.8, 54.3, 54.8, 
                         53.5, 54.0, 54.5, 55.0, 55.5, 56.0, 56.5, 57.0, 57.5, 58.0, 53.5, 54.0, 54.5, 
                         53.8, 54.3, 54.8, 55.3, 55.8, 56.3, 56.8, 57.3, 57.8, 58.3, 53.8, 54.3, 54.8, 
                         52.7, 53.2, 53.7, 54.2, 54.7, 55.2, 55.7, 56.2, 56.7, 57.2, 52.7, 53.2, 53.7, 
                         53.8, 54.3, 54.8, 55.3, 55.8, 56.3, 56.8, 57.3, 57.8, 58.3, 
                         # 52.5, 53.0, 53.5, 54.0, 54.5, 55.0, 55.5, 56.0, 56.5, 
                         52.5, 53.0, 53.5, 54.0, 54.5, 55.0, 55.5, 56.0, 56.5, 57.0, 52.5, 53.0, # 53.5, 
                         52.2, 52.7, 53.2, 53.7, 54.2, 54.7, 55.2, 55.7, 56.2, 56.7, 52.2, 52.7, # 53.2, 
                         51.9, 52.4, 52.9, 53.4, 53.9, 54.4, 54.9, 55.4, 55.9, 56.4, 51.9, 52.4, 52.9, 
                         51.7, 52.2, 52.7, 53.2, 53.7, 54.2, 54.7, 55.2, 55.7, 56.2, 51.7, 52.2, 52.7,
                         51.4, 51.9, 52.4, 52.9, 53.4, 53.9, 54.4, 54.9, 55.4, 55.9, 51.4, 51.9, 52.4, 
                         51.1, 51.6, 52.1, 52.6, 53.1, 53.6, 54.1, 54.6, 55.1, 55.6, 51.1, 51.6, 52.1], 
              "highGain":[False, False, False, False, False, False, False, False, False, False, 
                          False, False, False, False, False, False, False, False, False, False, True, True, True, 
                          False, False, False, False, False, False, False, False, False, False, True, True, True, 
                           True,  True,  True,  True,  True,  True,  True,  True,  True,  True, True, True, True,
                          False, False, False, False, False, False, False, False, False, False, True, True, True, 
                          False, False, False, False, False, False, False, False, False, False, 
                          # False, False, False, False, False, False, False, False, False,  
                          False, False, False, False, False, False, False, False, False, False, True, True, # True, 
                          False, False, False, False, False, False, False, False, False, False, True, True, # True, 
                          False, False, False, False, False, False, False, False, False, False, True, True, True, 
                          False, False, False, False, False, False, False, False, False, False, True, True, True, 
                          False, False, False, False, False, False, False, False, False, False, True, True, True, 
                          False, False, False, False, False, False, False, False, False, False, True, True, True]},
          "SiDetFEBScan": {
              "type":"led", 
              "config":"crvled-004",
              "run#":[1599, 1601, 1602], # run 1600 had 9/10 empty spills, ignore... 
              "tempN":[10., 30., 40.]},
              # "run#":[1599, 1600, 1601, 1602],
              # "tempN":[10., 20., 30., 40.]}
          "gainAFEScan": {
              "type":"led", 
              "config":"crvled-003",
              "run#":[1620, 1621, 1622, 1623, 1624, 1625, 
                      1626, 1627, 1628, 1629, 1630],
              "gainAFE":[0x000, 0x800, 0x400, 0xc00, 0x180, 0x600,
                         0x280, 0x100, 0x200, 0x080, 0x300]},
          "bulk_scan_crvaging011_0": {
              "type":"led", 
              "config":"crvled-005",
              "run#":[1694, 1695, 1696, 1697, 1698],
              "bulk":[-25, 0, +25, +50, +75]},
          "bulk_scan_crvaging011_1": {
              "type":"led", 
              "config":"crvled-005",
              "run#":[1700, 1701, 1702, 1703, 1704],
              "bulk":[-25, 0, +25, +50, +75]},
          "SiDet2GainCalib": {
              "type":"led",
              "config":"crvled-004",
              "run#":list(range(1755, 1824+1)),
              "tempN":[25.]*len(range(1755, 1824+1)),
              "VoverN":[1.5+i*0.5 for j in range(7) for i in range(10)],
              "gainAFE":[0x384]*10 + [0x354]*10 + [0x3b4]*10 + [0x3a4]*10 + [0x364]*10 + [0x374]*10 + [0x394]*10},
          "SiDet2CMBScan": {
              "type":"led",
              "config":"crvled-004",
              "run#":list(range(1827, 1936+1)),
              "VoverN":[1.5+i*0.5 for j in range(11) for i in range(10)],
              # "gainAFE": vairous gain, n/a
              "VppN":[53.0+i*0.5 for i in range(10)] + \
                     [53.3+i*0.5 for i in range(10)] + \
                     [53.5+i*0.5 for i in range(10)] + \
                     [53.8+i*0.5 for i in range(10)] + \
                     [52.7+i*0.5 for i in range(10)] + \
                     [52.5+i*0.5 for i in range(10)] + \
                     [52.2+i*0.5 for i in range(10)] + \
                     [51.9+i*0.5 for i in range(10)] + \
                     [51.7+i*0.5 for i in range(10)] + \
                     [51.4+i*0.5 for i in range(10)] + \
                     [51.1+i*0.5 for i in range(10)],
              "tempSetup": [25.0]*10 + [30.0]*10 + [35.0]*10 + [40.0]*10 + [20.0]*10 + [15.0]*10 + [10.0]*10 + [5.0]*10 + [0.0]*10 + [-5.0]*10 + [-10.0]*10,
              "tempSensor": [25.0]*10 + [29.9]*10 + [34.9]*10 + [39.8]*10 + [20.1]*10 + [15.2]*10 + [10.3]*10 + [5.3]*10 + [0.4]*10 + [-4.4]*10 + [-9.4]*10},
          "SiDet1CMBScan": {
              "type":"led",
              "config":"crvled-004",
              "run#":[1434, 1433, 1435, 1436, 1437, 1438, 1439, 1440, 1441, 1442, 
                      1452, 1443, 1444, 1445, 1446, 1447, 1448, 1449, 1450, 1451,
                      1456, 1457, 1458, 1459, 1460, 1461, 1462, 1463, 1464, 1465,
                      1497, 1498, 1499, 1500, 1501, 1502, 1503, 1504, 1505, 1506,
                      1482, 1483, 1484, 1485, 1486, 1487, 1488, 1489, 1490, 1493, 
                      1516, 1517, 1518, 1519, 1520, 1521, 1522, 1523, 1524, 1525, 
                      1532, 1534, 1535, 1536, 1537, 1538, 1539, 1540, 1541, 1542, 
                      1547, 1548, 1549, 1550, 1551, 1552, 1553, 1554, 1555, 1556,
                      1560, 1561, 1562, 1563, 1564, 1565, 1566, 1567, 1568, 1569,
                      1573, 1574, 1575, 1576, 1577, 1578, 1579, 1580, 1581, 1582,
                      1586, 1587, 1588, 1589, 1590, 1591, 1592, 1593, 1594, 1595],
              "VoverN":[1.5+i*0.5 for j in range(11) for i in range(10)],
              # "gainAFE": vairous gain, n/a
              "VppN":[53.0+i*0.5 for i in range(10)] + \
                     [53.3+i*0.5 for i in range(10)] + \
                     [53.5+i*0.5 for i in range(10)] + \
                     [53.8+i*0.5 for i in range(10)] + \
                     [52.7+i*0.5 for i in range(10)] + \
                     [52.5+i*0.5 for i in range(10)] + \
                     [52.2+i*0.5 for i in range(10)] + \
                     [51.9+i*0.5 for i in range(10)] + \
                     [51.7+i*0.5 for i in range(10)] + \
                     [51.4+i*0.5 for i in range(10)] + \
                     [51.1+i*0.5 for i in range(10)],
              "tempSetup": [25.0]*10 + [30.0]*10 + [35.0]*10 + [40.0]*10 + [20.0]*10 + [15.0]*10 + [10.0]*10 + [5.0]*10 + [0.0]*10 + [-5.0]*10 + [-10.0]*10,
              "tempSensor": [25.0]*10 + [30.0]*10 + [35.0]*10 + [40.0]*10 + [20.0]*10 + [15.3]*10 + [10.3]*10 + [5.3]*10 + [0.4]*10 + [-4.4]*10 + [-9.4]*10},
          "SiDet2DynamicVbias": {
              "type":"led",
              "config":"crvled-004",
              "run#":list(range(1939, 2002+1))
              # "VoverN": 54.5V at 25 degC ref temp, dynamic correction with CMB correction only. 
              #           AFE0 tracking CMB0 only.
              # "gainAFE": vairous gain, n/a
              # "VppN": n/a
          },
          "bulk_scan_crvaging016_0": {
              "type":"led", 
              "config":"crvaging-016",
              "run#":[2052, 2053, 2054, 2055, 2056],
              "bulk":[-25, 0, +25, +50, +75]},
          "bulk_scan_crvaging019": {
              "type":"led", 
              "config":"crvaging-019",
              "run#":[2093, 2094, 2095, 2096, 2097],
              "bulk":[-25, 0, +25, +50, +75]}
}

# temperature scan data tags, exist for compatibility reasons
ds_temp_scan_cosmics = ["LED_temp_scan_cosmics"]
ds_temp_scan_led = ["LED_temp_scan"]
    
def getfilelist(taglist, filetype = "recoROOT", passVer = 1): # argument can be a list of keys from the above dictionary or run numbers
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
    
    fullfilelist = getfullfilelist(filetype, passVer)
    for runnumber in runnumlist:
        templist = []
        for file in fullfilelist:
            if filenameparser(file, 'run')==runnumber:
                if filenameparser(file, 'subrun') is not None:
                    templist.append(file)
        if templist:
            filenamelist += templist
        else:
            print("WARNING: filepath.getfilelist: no file with run number %i was found."%(runnumber))
    
    return filenamelist, runnumlist, configlist, nFEBlist
