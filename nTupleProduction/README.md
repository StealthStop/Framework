# Producing nTuples with TreeMaker
## (as of 8/2/2019)

This is a guide for producing nTuples from miniAOD files using [TreeMaker](https://github.com/TreeMaker/TreeMaker). It is meant to compliment the [README](https://github.com/TreeMaker/TreeMaker/blob/Run2_2017/README.md), and serves as a step-by-step guide rather than any kind of documentation of TreeMaker. It is not written, endorsed, nor maintained by the developers of TreeMaker. Please update this as necessary.

## Clone and install TreeMaker

Connect via `ssh` to either LPC or CMSConnect. Ensure that a grid proxy is alread set-up. Initiate a VOMS proxy with 

`voms-proxy-init -voms cms`

Now clone and install TreeMaker:

```
wget https://raw.githubusercontent.com/TreeMaker/TreeMaker/Run2_2017/
setup.sh
chmod +x setup.sh
./setup.sh
cd CMSSW_10_2_11_patch1/src/
cmsenv
```

If running `./setup.sh` gives an error indicating the repository does not exist, open `setup.sh` and change line 6 from `ACCESS=ssh` to `ACCESS=https`. Then save, close, and try again.


## For MC Samples

Optional: query DAS to ensure that the miniAOD sample names are correct. Example:
`dasgoclient --query="dataset dataset=/WW*/*Autumn18*/MINIAODSIM”`

Once the user has a list of miniAOD samples to nTupilize, proceed with the following.

### Overview:
* Generate `_cff.py` files
* Get neff values where appropriate
* Produce text files containing configuration lines and add them to `getWeightProducer_cff.py`
* Create new dictionaries for Condor submission
* Produce and submit jobs to Condor
* Manage jobs until complete

### Generating `_cff.py` files

Each sample needs a `_cff.py` that contains its corresponding list of root files.

```
cd Treemaker/Production/test`
source /cvmfs/cms.cern.ch/crab3/crab.csh
```

if in tcsh: `setenv SSL_CERT_DIR '/etc/pki/tls/certs:/etc/grid-security/certificates` 

if in bash: `export SSL_CERT_DIR='/etc/pki/tls/certs:/etc/grid-security/certificates` 

Create "dictionaries" (python files containing sample names, formatted depending on the script that will be using it) for each sample scenario (e.g. `Autumn18`, `Fall17`, `Summer16`, `Summer16v3`) in the style of [`dict.py`](https://github.com/TreeMaker/TreeMaker/blob/Run2_2017/Production/test/dict.py). Follow the formatting in the header assuming for now that there are no negative weights.

`cd ../test`

The script `get_py.py` run with the `-p` option is used to produce the `_cff.py` files. Since `get_py.py` is located in [`TreeMaker/Production/python`](https://github.com/TreeMaker/TreeMaker/tree/Run2_2017/Production/python), it is run in `/test` using `run_get_py.py`:

`python run_get_py.py -d [dictionaryname] -p -o ../python/[sample scenario]`

Run for each dictionary. For output destination, specify the relative path to the appropriate folder in [`TreeMaker/Production/python`](https://github.com/TreeMaker/TreeMaker/tree/Run2_2017/Production/python) for that dictionary's scenario. Alternatively, remove the `-o` option. This will produce the files in the current directory, and can be moved or copied to the appropriate folder.

### Getting neff values

Some MC samples have events with negative weights which must be handled appropriately. This is done by running a program on the samples in Condor that outputs histograms containing neff(= positive events - negative events) values. These neff values will be manually entered into the dictonaries created in the previous step. If the user is sure that there are no negative weighted events (it does not hurt to check) they may skip this part.

`cd condorSub`

A second type of dictionary must be created here in the style of [`dict_Fall17_neff.py`](https://github.com/TreeMaker/TreeMaker/blob/Run2_2017/Production/test/condorSub/dict_Fall17_neff.py). **The dictionary name must begin with** `dict_`. All of the samples, regardless of scenario, can be placed in the same dictionary.

```
cd ..
./lnbatch.sh myNeff
cd myNeff
```

This directory should contain a list of softlinks to files in `condorSub`.

Produce (or just check for) a tarball in a writeable EOS directory and create an output folder for the neff Condor jobs:

`./checkVomsTar.sh -i root://cmseos.fnal.gov//store/user/YOURUSERNAME/myNeff/`

Now produce and submit jobs. If the neff dictionary is named `dict_neff.py`,

`python submitJobsNeff.py -p -d neff -N 50 -o root://cmseos.fnal.gov//store/user/YOURUSERNAME/myNeff -s`

Note the prefix `dict_` is removed. Once jobs are done,

```
./haddEOS.sh -d /store/user/amercald/myNeff -g _part -r
python submitJobsNeff.py -g -d neff -N 50 -o root://cmseos.fnal.gov//store/user/amercald/myNeff
```

This should print each sample name in the neff dictionary with its neff value (and pos, neg, tot events). Now they need to be added to the dictionaries that were used to produce `_cff.py` files.

`cd ..`

Edit the dictionaries here used for `get_py.py` with the neff values for each sample, only entering the neff value if it's different from the number of tot events, i.e. `neg != 0`. 

### Produce configuration lines to add to `getWeightProducer_cff.py`

Now use `get_py.py` again with the `-w` option to generate text files containing the lines for configuring `getWeightProducer_cff.py`: 

`python run_get_py.py -d [dictionaryname] -w`

Repeat for each dictionary. This will produce a `.txt` file in the current directory (or a directory specified with the `-o` option) for that dictionary. Open a second terminal window or tab in LPC (or CMSConnect, whichever is being used) and navigate to `TreeMaker`.

`cd WeightProducer/python`

Here there is a list of files named `MCSamples_[scenario].py`. These files contain the configuration lines for that scenario's samples, which is appended into a sample list in `getWeightProducer_cff.py`.

In the first terminal window, open the `.txt` files and copy and paste each line into the respective `MCSamples` file. Save the files and close the second terminal window.

### Create dictionaries for Condor submission

`cd condorSub`

Another style of dictionary must be made in the style of e.g.[`dict_Autumn18_diboson.py`](https://github.com/TreeMaker/TreeMaker/blob/Run2_2017/Production/test/condorSub/dict_Autumn18_diboson.py) for each of your sample scenarios. **The dictionary names must begin with** `dict_`. Create them, and ensure to change the scenario name at the top where appropriate.

Once these dictionaries are made, the samples are ready to be submitted for nTuple production.

### Submit jobs to Condor

Condor submission here uses the [CondorProduction](https://github.com/kpedro88/CondorProduction) package.

```
cd ..
./lnbatch.sh myProduction
cd myProduction
```

There should again be a list of softlinks to the files in condorSub. The script `submitJobs.py` will be used to produce and submit jobs to Condor (see [here](https://github.com/TreeMaker/TreeMaker#submit-production-to-condor) for documentation). When running `submitJobs.py`, enter the dictionary names without the `dict_` prefix:

`python submitJobs.py -p -d [dictionary names separated by commas] -o root://cmseos.fnal.gov//store/user/YOURUSERNAME/NTUPLEDESTINAtION -s`

The jobs should now be submitted to Condor. The samples will be nTupilized and placed in the specified EOS directory.

### Manage jobs until completion

It is common for these jobs to be held for various reasons. They need to be monitored and resubmitted as necessary. To manage this, use either basic [HTCondor](https://research.cs.wisc.edu/htcondor/manual/v7.8/10_Command_Reference.html) commands or (better) the `manageJobs.py` script documented [here](https://github.com/kpedro88/CondorProduction#job-management).

#### Once all the jobs are complete, the nTuples have been produced. Add the samples to the appropriate place in dict_master.py to keep track of what is in the EOS space.


# getSampleList.py

This is a simple script that lists the samples in the master dictionary corresponding to a specified year and sample type (or all sample types). The output can be printed or saved to a `.txt` file.

### Summary of options:

* `-y`: Specify year (Summer16v3, Fall17, or Autumn18)
* `-s`: Specify one sample type (`QCD_HT`, `QCD_PT`, `QCD_PT_MuEnriched`, `diboson`, `dyjets`, `gjets`, `ttbar`, `tth`, `wjets`, `zjets`, or `singlet`)
* `-a`: Give list of all samples for that year
* `-p`: Print list in terminal
* `-f`: Save list to `.txt` file

The same information about options can also be obtained by running the script with the `-h` option for "help".

### Example usage:

Input: `python getSampleList.py -y Autumn18 -s diboson -pf` 

Output: 

```
Autumn18.WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8

Autumn18.WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8

Autumn18.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8

Autumn18.WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_ext1

Autumn18.WWZ_TuneCP5_13TeV-amcatnlo-pythia8_ext1

Autumn18.WZZ_TuneCP5_13TeV-amcatnlo-pythia8_ext1

Autumn18.ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_ext1

Saving all sample names for Autumn18 to sampleList_Autumn18_diboson.txt

```