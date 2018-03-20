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

ls
pwd

cd -

ls
pwd

git clone --depth=50 --branch=master https://github.com/susy2015/TopTagger.git 
cd TopTagger/TopTagger/test
./configure
make -j4

cd $1
ls


echo "Yeah Buddy"
