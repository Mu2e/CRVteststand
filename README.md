# CRVteststand
Code to manage data from the CRV test stand.

#Build:
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh 
muse setup
mkdir workDir
cd workDir
git clone git@github.com:Mu2e/CRVteststand.git
make

Run:
Modify config.txt 
./parserCrv Run_number
Ex: ./parserCrv 001052_000
./calibCrv Run_number
./recoCrv Run_number
