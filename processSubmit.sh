#!/bin/bash

help()
{
  echo "Usage:"
  echo "./processSubmit.sh -r listOfRawFiles.txt -e directoryOfExecutablesAndConfig -o outputDirectory [-p]"
}

rawFiles=""
execDir=""
outputDir=""
poisson=0

OPTIND=1  #needed by getopts

while getopts "r:e:o:ph" opt; do
  case "$opt" in
    r)
       rawFiles=$OPTARG
       ;;
    e)
       execDir=$OPTARG
       ;;
    o)
       outputDir=$OPTARG
       ;;
    p)
       poisson=1
       ;;
    h)
       help
       return 0 2> /dev/null || exit 0
       ;;
  esac
done

if ! [[ "$*" == *"-r"* && "$*" == *"-e"* && "$*" == *"-o"* ]]; then
  echo "one or more options are missing"
  help
  return 1 2> /dev/null || exit 1
fi

if [ ! -f "$rawFiles" ]; then
  echo "File with list of raw files does not exist!"
  help
  return 1 2> /dev/null || exit 1
fi

if [ ! -d "$execDir" ]; then
  echo "Directory with executable does not exist!"
  help
  return 1 2> /dev/null || exit 1
fi

if [ ! -f "$execDir/parserCrv" ]; then
  echo "parserCrv does not exist in directory of executables!"
  help
  return 1 2> /dev/null || exit 1
fi

if [ ! -f "$execDir/calibCrv" ]; then
  echo "calibCrv does not exist in directory of executables!"
  help
  return 1 2> /dev/null || exit 1
fi

if [ ! -f "$execDir/recoCrv" ]; then
  echo "calibCrv does not exist in directory of executables!"
  help
  return 1 2> /dev/null || exit 1
fi

if [ ! -f "$execDir/config.txt" ]; then
  echo "config.txt does not exist in directory of executables!"
  help
  return 1 2> /dev/null || exit 1
fi

if [ ! -d "$outputDir" ]; then
  echo "Output directory does not exist!"
  help
  return 1 2> /dev/null || exit 1
fi

mkdir $outputDir/crvparsed
mkdir $outputDir/crvcalib
mkdir $outputDir/crvreco
mkdir $outputDir/log

processes=`< $rawFiles wc -l`
jobname="mu2eWideband"
currentDir=`pwd`

jobsub_submit \
-N $processes \
-e JOBNAME=$jobname \
-e RAWFILES=$rawFiles \
-e EXECDIR=$execDir \
-e OUTPUTDIR=$outputDir \
-e POISSON=$poisson \
--memory=1GB \
--disk=100GB \
--expected-lifetime=4h \
--resource-provides usage_model=OPPORTUNISTIC,DEDICATED,OFFSITE \
--group=mu2e \
--append_condor_requirements='(TARGET.HAS_SINGULARITY=?=true)' \
--lines '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' \
--lines '+FERMIHTC_AutoRelease=True' \
--lines '+FERMIHTC_GraceMemory=1024' \
--lines '+FERMIHTC_GraceLifetime=7200' \
file:"//$currentDir/processDo.sh"

