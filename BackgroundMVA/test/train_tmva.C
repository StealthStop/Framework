#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

void loadTree(TMVA::DataLoader* loader, const std::string& type, const double weight, const std::string& file, const std::string& treePath)
{
    TFile *f = TFile::Open( file.c_str() );
    TTree *t = (TTree*) f->Get( treePath.c_str() );
    if (type=="Signal")
    {
        loader->AddSignalTree( t, weight );
    }
    else if (type=="Background")
    {
        loader->AddBackgroundTree( t, weight );
    }
}

int train_tmva()
{
    TString myMethodList = "" ;
    
    //---------------------------------------------------------------
    // This loads the library
    TMVA::Tools::Instance();
    
    // Default MVA methods to be trained + tested
    std::map<std::string,int> Use;
    // ---------------------------------------------------------------

    std::cout << std::endl;
    std::cout << "==> Start TMVAClassification" << std::endl;

    // Select methods (don't look at this code - not of interest)
    if (myMethodList != "") 
    {
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

        std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
        for (UInt_t i=0; i<mlist.size(); i++) 
        {
            std::string regMethod(mlist[i]);

            if (Use.find(regMethod) == Use.end()) 
            {
                std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
                for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
                std::cout << std::endl;
                return 1;
            }
            Use[regMethod] = 1;
        }
    }

    // --------------------------------------------------------------------------------------------------
    // --- Here the preparation phase begins

    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    gSystem -> Exec( "mkdir -p outputfiles" ) ;
    TString outfileName( "outputfiles/tmva-train-example-output.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    // Create the factory object. Later you can choose the methods
    // whose performance you'd like to investigate. The factory is 
    // the only TMVA object you have to interact with
    //
    // The first argument is the base of the name of all the
    // weightfiles in the directory weight/
    //
    // The second argument is the output file for the training results
    // All TMVA output can be suppressed by removing the "!" (not) in
    // front of the "Silent" argument in the option string
    TMVA::Factory* factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                                "V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
    TMVA::DataLoader* loader = new TMVA::DataLoader("fisherLoader");

    loader->AddVariable( "fwm2_top6", 'D' ) ;
    loader->AddVariable( "fwm3_top6", 'D' ) ;
    loader->AddVariable( "fwm4_top6", 'D' ) ;
    loader->AddVariable( "fwm5_top6", 'D' ) ;
    loader->AddVariable( "fwm6_top6", 'D' ) ;
    loader->AddVariable( "jmt_ev0_top6", 'D' ) ;
    loader->AddVariable( "jmt_ev1_top6", 'D' ) ;
    loader->AddVariable( "jmt_ev2_top6", 'D' ) ;

    //-------------
    // You can add so-called "Spectator variables", which are not used in the MVA training,
    // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
    // input variables, the response values of all trained MVAs, and the spectator variables
    //-------------

    //loader -> AddSpectator( "ds_index", "ds_index" ) ;

    // Read training and test data
    //loadTree(TMVA::DataLoader* loader, const std::string& type, const double weight, const std::string& file)

    loadTree(loader, "Signal", 1.0, "condor/output-files/AllSignal/make_training_trees_stealth_stop_350_SHuHd_0.root", "mvatraintt");
    loadTree(loader, "Signal", 1.0, "condor/output-files/AllSignal/make_training_trees_stealth_stop_450_SHuHd_0.root", "mvatraintt");
    loadTree(loader, "Signal", 1.0, "condor/output-files/AllSignal/make_training_trees_stealth_stop_550_SHuHd_0.root", "mvatraintt");
    loadTree(loader, "Signal", 1.0, "condor/output-files/AllSignal/make_training_trees_stealth_stop_650_SHuHd_0.root", "mvatraintt");
    loadTree(loader, "Signal", 1.0, "condor/output-files/AllSignal/make_training_trees_stealth_stop_750_SHuHd_0.root", "mvatraintt");
    loadTree(loader, "Signal", 1.0, "condor/output-files/AllSignal/make_training_trees_stealth_stop_850_SHuHd_0.root", "mvatraintt");

    //loadTree(loader, "Background", 0.388, "condor/output-files/TT/make_training_trees_TT.root"  , "mvatraintt");
    //loadTree(loader, "Background", 0.388, "condor/output-files/QCD/mva-trees-QCD.root", "mvatraintt");
    
    loadTree(loader, "Background", 1.0, "condor/output-files/TT/make_training_trees_TT.root"  , "mvatraintt");
    loadTree(loader, "Background", 1.0, "condor/output-files/QCD/make_training_trees_QCD.root", "mvatraintt");

    //Define Event wise weight
    loader->SetBackgroundWeightExpression("Weight");
    loader->SetSignalWeightExpression    ("Weight");

    // 0 Lepton Selection
    TCut mycuts = "passBaseline0l";
    TCut mycutb = "passBaseline0l";

    // 1 Lepton Selection
    //TCut mycuts = "passBaseline1l && Mbl>30 && Mbl<180";
    //TCut mycutb = "passBaseline1l && Mbl>30 && Mbl<180";

    // Tell the factory how to use the training and testing events
    loader->PrepareTrainingAndTestTree( mycuts, mycutb, "SplitMode=Random:NormMode=None:!V" );

    // ---- Book MVA methods
    //
    // Please lookup the various method configuration options in the corresponding cxx files, eg:
    // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
    // it is possible to preset ranges in the option string in which the cut optimisation should be done:
    // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

    // Fisher with Gauss-transformed input variables
    factory->BookMethod(loader,  TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss_Background" );

    // For an example of the category classifier usage, see: TMVAClassificationCategory

    // --------------------------------------------------------------------------------------------------
    // ---- Now you can tell the factory to train, test, and evaluate the MVAs

    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // ---- Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // ----- Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // --------------------------------------------------------------

    // Save the output
    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
    delete loader;

    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

    return 0;
}

int main()
{
    return train_tmva(); 
}
