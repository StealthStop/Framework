# BackgroundMVA
The code make_training_trees.C produces a small TTree that's used as input to the mva training that's done by train_tmva.C

## Setup

cmsrel CMSSW_9_3_3
cd CMSSW_9_3_3/src/
cmsenv
git clone git@github.com:StealthStop/Framework.git
git clone -b TopTaggingLPC git@github.com:susy2015/SusyAnaTools.git
cd Framework/BackgroundMVA/test/
$CMSSW_BASE/src/Framework/Framework/scripts/getSamplesCfg.sh 
make -j4 

## Example Code
Make trees to run the training on 

Make trees interactively
```
./make_training_trees -D rpv_stop_350 -H outputfiles/mva-trees-rpv_stop_350.root
./make_training_trees -D rpv_stop_450 -H outputfiles/mva-trees-rpv_stop_450.root
./make_training_trees -D rpv_stop_550 -H outputfiles/mva-trees-rpv_stop_550.root
./make_training_trees -D rpv_stop_650 -H outputfiles/mva-trees-rpv_stop_650.root
./make_training_trees -D rpv_stop_750 -H outputfiles/mva-trees-rpv_stop_750.root
./make_training_trees -D rpv_stop_850 -H outputfiles/mva-trees-rpv_stop_850.root

./make_training_trees -D TT -H outputfiles/mva-trees-ttbar.root
```

Make trees with condor (Same options as Analyzer)
```
cd condor
python condorSubmit.py -d AllSignal,TT,QCD -n 10
./combineTree.sh
```

Run the training
```
./train_tmva
```

## Output

  When you run train_tmva.C, the output will be in a couple of places.

    ==== outputfiles/tmva-train-example-output.root

      This root file contains TTrees of the train and test samples.
      For example, in TestTree, there will be branches for all of the
      input variables plus the spectator variables you asked for.
      The output of the fisher that was trained will be in FisherG.
      The branch classID is 0 for signal and 1 for background events.
      The events in TestTree were not used in the training.

      If you want to rerun the TMVA GUI, you can do it with this file like this
      from the root prompt

        TMVA::TMVAGui("outputfiles/tmva-train-example-output.root")

      This doesn't redo the training.  It just lets you draw the TMVA output
      plots from the output file.


    ==== weights/TMVAClassification_FisherG.class.C

      You probably recognize this as the code that does the Gaussian
      input variable transformations and computes the fisher.
