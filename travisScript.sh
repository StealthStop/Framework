#!/bin/bash
set -ex

cd $1
#source /cvmfs/sft.cern.ch/lcg/views/LCG_89/x86_64-slc6-gcc62-opt/setup.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc6_amd64_gcc630
scramv1 project CMSSW CMSSW_9_3_3
cd CMSSW_9_3_3/src
set +x
eval `scramv1 runtime -sh`
set -x

pwd

cd -

pwd

echo "Yeah Buddy"
