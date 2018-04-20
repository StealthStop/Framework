
  Wed Apr 18 13:56:43 PDT 2018

  The code make_mva_training_tree_example.c produces a small TTree that's used as input
  to the mva training that's done by tmva_train_example.c

  When you run tmva_train_example.c, the output will be in a couple of places.

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
