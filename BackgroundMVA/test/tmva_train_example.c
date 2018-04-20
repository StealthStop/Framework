
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

int tmva_train_example()
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
    if (myMethodList != "") {
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

        std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
        for (UInt_t i=0; i<mlist.size(); i++) {
            std::string regMethod(mlist[i]);

            if (Use.find(regMethod) == Use.end()) {
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



    loader -> AddVariable( "fwm2_top6", 'D' ) ;
    loader -> AddVariable( "fwm3_top6", 'D' ) ;
    loader -> AddVariable( "fwm4_top6", 'D' ) ;
    loader -> AddVariable( "fwm5_top6", 'D' ) ;
    loader -> AddVariable( "fwm6_top6", 'D' ) ;

    loader -> AddVariable( "jmt_ev0_top6", 'D' ) ;
    loader -> AddVariable( "jmt_ev1_top6", 'D' ) ;
    loader -> AddVariable( "jmt_ev2_top6", 'D' ) ;



    //-------------

    // You can add so-called "Spectator variables", which are not used in the MVA training,
    // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
    // input variables, the response values of all trained MVAs, and the spectator variables
    //////////factory->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
    //////////factory->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );

    //--- owen : What will happen if I add spectators for vars used in MVA?
    //           Duplicated variables will be duplicated in output tree, so not good.


    loader -> AddSpectator( "ds_index", "ds_index" ) ;
    loader -> AddSpectator( "mva_train_weight", "mva_train_weight" ) ;

    loader -> AddSpectator( "njets_pt45_eta24", "njets_pt45_eta24" ) ;
    loader -> AddSpectator( "njets_pt30_eta24", "njets_pt30_eta24" ) ;
    loader -> AddSpectator( "njets_pt20_eta50", "njets_pt20_eta50" ) ;
    loader -> AddSpectator( "nbtag_csv85_pt30_eta24", "nbtag_csv85_pt30_eta24" ) ;
    loader -> AddSpectator( "pfht_pt40_eta24", "pfht_pt40_eta24" ) ;
    loader -> AddSpectator( "pfht_pt45_eta24", "pfht_pt45_eta24" ) ;
    loader -> AddSpectator( "nleptons", "nleptons" ) ;
    loader -> AddSpectator( "leppt1", "leppt1" ) ;
    loader -> AddSpectator( "m_lep1_b", "m_lep1_b" ) ;
    loader -> AddSpectator( "leppt2", "leppt2" ) ;
    loader -> AddSpectator( "m_lep2_b", "m_lep2_b" ) ;

    loader -> AddSpectator( "evt_count", "evt_count" ) ;
    loader -> AddSpectator( "run", "run" ) ;
    loader -> AddSpectator( "lumi", "lumi" ) ;
    loader -> AddSpectator( "event", "event" ) ;





    // Read training and test data

    TFile *input_signal350  = TFile::Open( "outputfiles/mva-train-example-signal-rpv_stop_350.root" ) ;
    TFile *input_signal450  = TFile::Open( "outputfiles/mva-train-example-signal-rpv_stop_450.root" ) ;
    TFile *input_signal550  = TFile::Open( "outputfiles/mva-train-example-signal-rpv_stop_550.root" ) ;
    TFile *input_signal650  = TFile::Open( "outputfiles/mva-train-example-signal-rpv_stop_650.root" ) ;

    TFile *input_ttbar   = TFile::Open( "outputfiles/mva-train-example-ttbar.root" ) ;




    // --- Register the training and test trees

    TTree *tt_signal350  = (TTree*) input_signal350  -> Get( "mvatraintt" ) ;
    TTree *tt_signal450  = (TTree*) input_signal450  -> Get( "mvatraintt" ) ;
    TTree *tt_signal550  = (TTree*) input_signal550  -> Get( "mvatraintt" ) ;
    TTree *tt_signal650  = (TTree*) input_signal650  -> Get( "mvatraintt" ) ;
    TTree *tt_ttbar      = (TTree*) input_ttbar      -> Get( "mvatraintt" ) ;



    // You can add an arbitrary number of signal or background trees
    loader->AddSignalTree    ( tt_signal350 , 1.    );
    loader->AddSignalTree    ( tt_signal450 , 1.    );
    loader->AddSignalTree    ( tt_signal550 , 1.    );
    loader->AddSignalTree    ( tt_signal650 , 1.    );
    loader->AddBackgroundTree( tt_ttbar, 0.388 );

    TCut mycuts = "njets_pt30_eta24>=6 && nleptons>=1 && ( (leppt1>30 && m_lep1_b > 30 && m_lep1_b < 180) || (leppt2>30 && m_lep2_b > 30 && m_lep2_b < 180) )";
    TCut mycutb = "njets_pt30_eta24>=6 && nleptons>=1 && ( (leppt1>30 && m_lep1_b > 30 && m_lep1_b < 180) || (leppt2>30 && m_lep2_b > 30 && m_lep2_b < 180) )";


    // Tell the factory how to use the training and testing events
    //
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

    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

    return 0;
}

int main( int argc, char** argv )
{
    return tmva_train_example(); 
}
