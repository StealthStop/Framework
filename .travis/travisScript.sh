#!/bin/bash

ls -a
ls -a /home/
#stop upon failed command
set -ex

TRAVIS_BUILD_DIR=$1
cd /home/cmsusr/
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch                                        
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git/                       
source $VO_CMS_SW_DIR/cmsset_default.sh 
scramv1 project CMSSW CMSSW_11_2_0_pre5
cd CMSSW_11_2_0_pre5/src/
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
git clone -b Stealth https://github.com/susy2015/TopTaggerTools.git
git clone https://github.com/susy2015/NTupleReader.git
git clone -b Run2_UL https://github.com/StealthStop/Analyzer.git
git clone https://github.com/susy2015/TopTagger.git
cd TopTagger/TopTagger/test
./configure
make -j4
echo "============================"
echo "Testing Analyzer"
echo "============================"
cd ${CMSSW_BASE}/src/Analyzer/Analyzer/test
source ${CMSSW_BASE}/src/TopTagger/TopTagger/test/taggerSetup.sh
export PATH=$PATH:${CMSSW_BASE}/src/Framework/Framework/scripts
./configure
make -j4
echo "============================"
echo "Testing BackgroundMVA"
echo "============================"
cd ${CMSSW_BASE}/src/Framework/BackgroundMVA/test
make -j4

echo "Yeah Buddy"
