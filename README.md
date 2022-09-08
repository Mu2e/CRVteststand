# CRVteststand
Code to manage data from the CRV test stand.  

Build:  
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh  
muse setup  
mkdir workDir  
cd workDir  
git clone git@github.com:Mu2e/CRVteststand.git  
cd CRVteststand  
make  

Run:  
Modify config.txt  
./parserCrv Run_number  
Ex: ./parserCrv 001052_000  
./calibCrv Run_number  
./recoCrv Run_number  

Run on OSG:<br>
Requires the following things
<ul>
<li>Text file with a list of raw files to be processed.</li>
<li>A directory that is accessible by the OSG, where the executables (parserCrv, calibCrv, and recoCrv) 
  and the config.txt file are located. The paths mentioned in the config.txt are ignored by the OSG script</li>
<li>An output directory that is accessible by the OSG, where the parsed, calib, and reco files will arrive.</li>
<li>The submit script (processSubmit.sh) and script that runs at the OSG (processDo.sh). Both files are at github.</li>
</ul>
Submit the jobs in this way:
./processSubmit.sh -r listOfRawFiles.txt -e directoryOfExecutablesAndConfig -o outputDirectory<br>
<br>
Event Display<br>
Example: root -l "DisplayEventsNew.C(\"/pnfs/mu2e/scratch/outstage/ehrlich/wideband9/crvreco/rec.mu2e.CRV_wideband_cosmics.crvaging-005-junk.001137.root\",\"channelMapTest5.txt\")"
