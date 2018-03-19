#!/bin/tcsh

#--- You should copy this to a scratch disk where you will manage your job submission and output
#    and prepare the two tar files as described below.

##--- See this web page: http://uscms.org/uscms_at_work/computing/setup/condor_worker_node.shtml

##--- See also make_condor_jdl_files.c which makes jdl files that go with this executable script.

set dataset = $1
set nfiles = $2
set startfile = $3
set filelist = $4
set options = $5

set base_dir = `pwd`
printf "\n\n base dir is $base_dir\n\n"

source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH slc6_amd64_gcc630

printf "\n\n ls output\n"
ls -l


#--- You need to make a tar file of your CMSSW directory that contains the
#      TopTagger stuff and your compiled code from the Exploration directory.
#      You may need to change some symbolic links that reference full path
#      names to relative path names in order for the untarred file to work
#      on the batch node.  I forgot which specific changes I made.
#      Once you have done that, make the tar file with something like this
#
#        tar -cvf cmssw-toptagger.tar CMSSW_8_0_28
#

printf "\n\n unpacking CMSSW tar file.\n\n"
tar -xf CMSSW_9_3_3.tar.gz

printf "\n\n ls output\n"
ls -l

printf "\n\n changing to CMSSW_9_3_3/ dir\n"
cd CMSSW_9_3_3/
mkdir -p src
cd src
scram b ProjectRename
eval `scramv1 runtime -csh`

printf "\n\n ls output\n"
ls -l

printf "\n\n output of uname -s : "
uname -s
printf "\n\n"

cp ${base_dir}/exestuff.tar.gz .
tar xzvf exestuff.tar.gz

cd exestuff/

setenv LD_LIBRARY_PATH ${PWD}:${LD_LIBRARY_PATH}

printf "\n\n LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n\n"

printf "\n\n ls output\n"
ls -l

printf "Copy over the needed filelist"
mkdir filelists_Kevin
xrdcp root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/${filelist} filelists_Kevin/

printf "\n\n Attempting to run MyAnalysis executable.\n\n"
#./MyAnalysis root://cmseos.fnal.gov/${input_fpat} ${output_file} ${weight}
./MyAnalysis -${options} --condor -D ${dataset} -N ${nfiles} -M ${startfile}


printf "\n\n ls output\n"
ls -l

mv *.root ${base_dir}

cd ${base_dir}

printf "\n\n ls output\n"
ls -l

