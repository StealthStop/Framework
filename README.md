# Framework
Code used to read our nTuples and compute necessary variables for our analysis.


## Using the tensor-flow based top tagger

To have easy access to TensorFlow, we need to work in a CMSSW93 release:
```
cmsrel CMSSW_9_3_3
cd CMSSW_9_3_3/src/
cmsenv
```

Then, check out the latest tagged version of the top tagger repository. 

```
git clone git@github.com:susy2015/TopTagger.git
```

Then configure and compile the tagger:
```
cd TopTagger/TopTagger/test
./configure 
make -j4
```

Now also check out our repository if not done already:
```
cd $CMSSW_BASE/src
git clone git@github.com:StealthStop/Framework.git
git clone git@github.com:susy2015/TopTaggerTools.git
git clone -b TopTaggingLPC git@github.com:susy2015/SusyAnaTools.git
git clone git@github.com:StealthStop/Analyzer.git
cd Analyzer/Analyzer/test
source setup.csh
make -j4
```

Last step is to get the cfg and model files for the top tagger and deepESM.
```
cmsenv
getTaggerCfg.sh -t Tensorflow_Medium_Example_v1.0.2 -o
getDeepESMCfg.sh -t Keras_Tensorflow_v1.0.0 -o
```

No changes to the analysis code should be needed. 


## Condor submission

The condor directory contains some scripts to help submit jobs via condor on the cmslpc cluster. 
The requirements for condor submission are: 
 - A script to run on the worker node. This script should set up the area, copy any needed files, call your executable with the right options, and make sure the output gets copied to where you want. The example included here is [run_Exploration_condor.tcsh](condor/run_Exploration_condor.tcsh)
 - One or more tarballs to unpack on the worker node, these usually contain a slimmed down CMSSW area, and your executable with any needed libraries
 - A so-called jdl file that contains the condor setup and specifies the jobs to be submitted
The last two items are produced by a python script called [condorSubmit.py](condor/condorSubmit.py). 

```
[condor]$ python condorSubmit.py -h
Usage: condorSubmit.py [options]


Options:
  -h, --help         show this help message and exit
  -n NUMFILE         number of files per job
  -d DATASETS        List of datasets, comma separated
  -l                 List all datacollections
  -L                 List all datacollections and sub collections
  -c                 Do not submit jobs.  Only create condor_submit.txt.
  --explore=EXPLORE  ExploreTopTagger (t), ExploreBackground (b),
                     ExploreEventSelection (s)
```
As you can see from the above help menu, there are a number of options. 
With the `-n` option you can specify how many files to run over per job. The `--explore` option lets you pick which analyzer to use. 
The MyAnalysis program has been updated to have these same switches. 
MyAnalysis now also uses the samples code to keep track of datasets, their cross sections, and their names. 
To see a list of available datasets, you can call the submission script with the `-l` or `-L` options. Pass the list of datasets you want to run over to the script with the option `-d`. 
Before submitting jobs, make sure to have called `voms-proxy-init`. 

