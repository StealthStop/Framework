#!/bin/tcsh

# set up the top tagger libraries etc
echo "Set up top tagger"
source ${CMSSW_BASE}/src/TopTagger/TopTagger/test/taggerSetup.csh

# Copy over filelists if they are not present or have changed
if (! -d condor/filelists ) then
    echo "No filelists found, copying from eos"
    mkdir condor/filelists/
    xrdcp -r root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/filelists/ condor/filelists/
    ln -s condor/filelists filelists
else
    echo "You already have the filelists. To get the current version, delete condor/filelists and run this again"
endif

if (! -d condor/filelists_Kevin ) then
    echo "No filelists_Kevin found, copying from eos"
    mkdir condor/filelists_Kevin/
    xrdcp -r root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/filelists_Kevin/ condor/filelists_Kevin/
    ln -s condor/filelists_Kevin filelists_Kevin
else
    echo "You already have the filelists_Kevin. To get the current version, delete condor/filelists_Kevin and run this again"
endif
