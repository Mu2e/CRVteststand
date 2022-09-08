# CRVteststand
Code to manage data from the CRV test stand.  

<h3>Build</h3>  
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh<br>
muse setup<br>
mkdir workDir<wr> 
cd workDir<br>
git clone git@github.com:Mu2e/CRVteststand.git<br>
cd CRVteststand<br>
make<br>

<h3>Run locally (see below for OSG jobs)</h3>
<i>Make sure your files are pre-staged, if they are on tape.</i><br>
Modify config.txt<br>
./parserCrv Run_number<br>
Example: ./parserCrv 001052_000<br>  
./calibCrv Run_number<br>
./recoCrv Run_number<br>

<h3>Event Display</h3>
Example: root -l "DisplayEventsNew.C(\"/pnfs/mu2e/scratch/outstage/ehrlich/wideband9/crvreco/rec.mu2e.CRV_wideband_cosmics.crvaging-005-junk.001137.root\",\"channelMapTest5.txt\")"<br>

<h3>Run on OSG</h3>
<i>Make sure your files are pre-staged, if they are on tape.</i><br>
Requires the following things<br>
<ul>
<li>Text file with a list of raw files to be processed.</li>
<li>A directory that is accessible by the OSG, where the executables (parserCrv, calibCrv, and recoCrv) 
  and the config.txt file are located. The paths mentioned in the config.txt are ignored by the OSG script</li>
<li>An output directory that is accessible by the OSG, where the parsed, calib, and reco files will arrive.</li>
<li>The submit script (processSubmit.sh) and script that runs at the OSG (processDo.sh). Both files are at github.</li>
</ul>
Submit the jobs in this way:<br>
./processSubmit.sh -r listOfRawFiles.txt -e directoryOfExecutablesAndConfig -o outputDirectory<br>
