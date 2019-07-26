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
echo "============================"
echo "Testing the top tagger"
echo "============================"
git clone https://github.com/susy2015/TopTaggerTools.git
git clone https://github.com/susy2015/SusyAnaTools.git
git clone -b adding2018 https://github.com/StealthStop/Analyzer.git
git clone https://github.com/susy2015/TopTagger.git
cd TopTagger/TopTagger/test
#./configure TENSORFLOWDIR=
./configure
make -j4
echo "============================"
echo "Testing Analyzer"
echo "============================"
cd ../../../
cd Analyzer/Analyzer/test
source ${CMSSW_BASE}/src/TopTagger/TopTagger/test/taggerSetup.sh
make -j4
#source runTest.sh
#getTaggerCfg.sh -t Tensorflow_Simple_Example_v1.0.0 -o
#make -j4
#echo "========================================================================="
#./MyAnalysis -s -H myoutputfile.root -D TT -E 2000

echo "============================"
echo "Testing BackgroundMVA"
echo "============================"
cd ../../../
cd Framework/BackgroundMVA/test
make -j4

echo "Yeah Buddy"
