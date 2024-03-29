#!/bin/bash 

date
hostname
lsb_release -a

#setup environment
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh
setup mu2e
setup muse
muse setup
setup mu2egrid v6_09_00
ups active

setup ifdhc

# to prevent wasting time trying to copy non-existing stage 1 files
export IFDH_GRIDFTP_EXTRA="-st 100"
export IFDH_CP_MAXRETRIES=0

# copy and run executable
ifdh cp $RAWFILES rawfiles.txt
ifdh cp "$EXECDIR/parserCrv" parserCrv
ifdh cp "$EXECDIR/calibCrv" calibCrv
ifdh cp "$EXECDIR/recoCrv" recoCrv
ifdh cp "$EXECDIR/config.txt" config.txt

echo "===================================================="

chmod 777 parserCrv
chmod 777 calibCrv
chmod 777 recoCrv

inputfileRemote="$(sed "$((PROCESS+1))q;d" rawfiles.txt)"
inputfile="$(basename $inputfileRemote)"
IFS='.' read -ra filenameArray <<< "$inputfile"
runAndSubrun=${filenameArray[4]}
ifdh cp $inputfileRemote $inputfile
date

#remove directories from config.txt, if present
sed -i '/crvraw/d' config.txt
sed -i '/crvparsed/d' config.txt
sed -i '/crvcalib/d' config.txt
sed -i '/crvreco/d' config.txt

#add current directory to config.txt
mkdir crvparsed
mkdir crvcalib
mkdir crvreco
mkdir log
echo crvraw$'\t'$(pwd) >> config.txt
echo crvparsed$'\t'$(pwd)/crvparsed/ >> config.txt
echo crvcalib$'\t'$(pwd)/crvcalib/ >> config.txt
echo crvreco$'\t'$(pwd)/crvreco/ >> config.txt

echo "===================================================="

cat config.txt

echo "===================================================="
echo "This is process $PROCESS"
echo "Raw file $inputfileRemote"
echo "===================================================="
echo "Parsing run/subrub $runAndSubrun"
parserCrv $runAndSubrun 
ls -l crvparsed
date
echo "===================================================="
echo "Calibrating run/subrub $runAndSubrun"
calibCrv $runAndSubrun 
ls -l crvcalib
date
echo "===================================================="
echo "Reconstructing run/subrub $runAndSubrun"
if [ $POISSON -eq 0 ]; then
  recoCrv $runAndSubrun
else
  recoCrv $runAndSubrun -p
fi
ls -l crvreco
date
echo "===================================================="

ls -l
ls -l crvparsed
ls -l crvcalib
ls -l crvreco

p1=$(basename crvparsed/*.root)
c1=$(basename crvcalib/*.txt)
c2=$(basename crvcalib/*.pdf)
r1=$(basename crvreco/*.root)
r2=$(basename crvreco/*.txt)
r3=$(basename crvreco/*.pdf)
ifdh cp crvparsed/$p1 $OUTPUTDIR/crvparsed/$p1
ifdh cp crvcalib/$c1 $OUTPUTDIR/crvcalib/$c1
ifdh cp crvcalib/$c2 $OUTPUTDIR/crvcalib/$c2
ifdh cp crvreco/$r1 $OUTPUTDIR/crvreco/$r1
ifdh cp crvreco/$r2 $OUTPUTDIR/crvreco/$r2
ifdh cp crvreco/$r3 $OUTPUTDIR/crvreco/$r3
date

# create log file
filenameArray[0]="log"
filenameArray[5]="txt"
logfile=$(printf ".%s" "${filenameArray[@]}")
logfile=${logfile:1}
cat jsb_tmp/JOBSUB_LOG_FILE >> $logfile
echo "=========== error log file ==========" >> $logfile
cat jsb_tmp/JOBSUB_ERR_FILE >> $logfile
ifdh cp $logfile "$OUTPUTDIR/log/$logfile"
