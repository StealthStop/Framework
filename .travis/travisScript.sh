#!/bin/bash

#stop upon failed command
set -ex

TRAVIS_BUILD_DIR=$1
cd $TRAVIS_BUILD_DIR/..
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch                                        
export SCRAM_ARCH=slc6_amd64_gcc630
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git/                       
source $VO_CMS_SW_DIR/cmsset_default.sh 
scramv1 project CMSSW CMSSW_9_3_3
cd CMSSW_9_3_3/src/
#suppress huge printout from "cmsenv"
set +x
echo "============================"
echo "SETTING UP CMSSW ENVIRONMENT"
echo "============================"
eval `scramv1 runtime -sh`
set -x
cp -r $TRAVIS_BUILD_DIR .
ls -l
pwd
git clone https://github.com/susy2015/TopTagger.git
cd TopTagger/TopTagger/test
echo "========================================================================="
./configure TENSORFLOWDIR=
make -j4
echo "========================================================================="
cd ../../../
cd Framework/Framework/test
source runTest.sh
#getTaggerCfg.sh -t Tensorflow_Simple_Example_v1.0.0 -o
#make -j4
#echo "========================================================================="
#./MyAnalysis -s -H myoutputfile.root -D TT -E 2000
echo "Yeah Buddy"
