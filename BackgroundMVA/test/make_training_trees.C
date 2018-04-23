#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TClassTable.h"
#include "TMath.h"
#include "TSystem.h"
#include "TH1F.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Framework/BackgroundMVA/include/MakeMVAVariables.h"

void make_mva_training_tree_example( const char* ntuple_dir = "prod-hadlep-skim-v2/",
                                     const char* sample_string = "rpv_stop_350",
                                     const char* outfile = "outputfiles/mva-train-rpv_stop_350.root",
                                     int arg_ds_index = 11,
                                     float lumi_times_xsec = ( 3.79 * 35.9 * 1000.0 ),
                                     bool verb = false) 
{
    char fpat[1000] ;
    TChain* tt_in = new TChain( "TreeMaker2/PreSelection", "" ) ;
    sprintf( fpat, "%s/*%s*.root", ntuple_dir, sample_string ) ;
    printf("\n\n Loading files that match %s\n", fpat ) ;

    int n_added = tt_in->Add( fpat ) ;
    if ( n_added <= 0 ) { printf("\n\n *** No files match %s\n\n", fpat ) ; gSystem->Exit(-1) ; }
    printf("  Added %d files to chain.\n", n_added ) ;
    int n_entries = tt_in->GetEntries() ;
    if ( n_entries <= 0 ) { printf("\n\n *** No entries in ntuple chain.\n\n" ) ; gSystem->Exit(-1) ; }
    printf("  Number of entries: %d\n\n", n_entries ) ;
    NTupleReader tr(tt_in);

    //**** Need to use a pre-skim histogram in the file(s) in this chain to get the number of generated
    //     events run on in order to correctly compute the dataset weight factor.

    double n_entries_pre_skim(0) ;
    TObjArray* files = tt_in->GetListOfFiles() ;
    if ( files == nullptr ) { printf("\n\n *** List of files returned null pointer.\n\n") ; gSystem->Exit(-1) ; } 
    printf("\n\n Number of files in chain %d\n", files->GetEntries() ) ;
    for ( int fi=0; fi<files->GetEntries(); fi++ ) 
    {
        TObject* op = files->At(fi) ;
        if ( op == nullptr ) { printf("\n\n *** null pointer!\n\n") ; gSystem->Exit(-1) ; }
        printf("  file %2d : name %s  title %s\n", fi, op->GetName(), op->GetTitle() ) ;
        TFile* fp = new TFile( op->GetTitle(), "READ" ) ;
        if ( fp == nullptr ) { printf("\n\n *** pointer to TFile is null.\n\n" ) ; gSystem->Exit(-1) ; }
        if ( ! fp->IsOpen() ) { printf("\n\n *** file is not open.\n\n" ) ; gSystem->Exit(-1) ; }
        TH1F* hp = (TH1F*) fp->Get( "h_njets_pt45" ) ;
        if ( hp == nullptr ) 
        { 
            printf("\n\n *** can't find histogram h_njets_pt45 in file %d\n", fi ) ; 
            n_entries_pre_skim += n_entries;
            break;
        }
        printf("  hist pointer : %p\n", hp ) ;
        printf("  hist name : %s \n", hp->GetName() ) ;
        double entries = hp->GetEntries() ;
        printf("  file %2d : entries = %9.1f\n", fi, entries ) ;
        n_entries_pre_skim += entries ;
        delete fp;
    } // fi

    if ( n_entries_pre_skim <= 0 ) 
    {
        printf("\n\n *** n_entries_pre_skim <= 0 ??? %f\n\n", n_entries_pre_skim ) ; 
        gSystem->Exit(-1) ; 
    }

    float ds_weight = lumi_times_xsec / n_entries_pre_skim ;
    printf("  Dataset weight : %.8f\n\n", ds_weight ) ;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    gSystem->Exec( "mkdir -p outputfiles" ) ;

    TFile* tf_output = new TFile( outfile, "RECREATE" ) ;
    if ( tf_output == nullptr || !tf_output->IsOpen() ) 
    {
        printf("\n\n *** bad output file: %s\n\n", outfile ) ; 
        gSystem->Exit(-1) ;
    }
    
    TTree* tt_out = new TTree( "mvatraintt", "MVA training ttree" ) ;

    //--- Extra output histograms
    TH1F* h_costheta_ppweight = new TH1F( "h_costheta_ppweight", "cos(theta_ij) in CM frame, pipj/Esq weight", 110, -1.05, 1.05 ) ;
    TH1F* h_costheta_ppweight_noieqj = new TH1F( "h_costheta_ppweight_noieqj", "cos(theta_ij) in CM frame, pipj/Esq weight, excluding i=j", 110, -1.05, 1.05 ) ;

    //--- New branches for output.
    int ds_index = arg_ds_index ;
    float mva_train_weight = ds_weight ;
    double fwm2_top6 ;
    double fwm3_top6 ;
    double fwm4_top6 ;
    double fwm5_top6 ;
    double fwm6_top6 ;
    double jmt_ev0_top6 ;
    double jmt_ev1_top6 ;
    double jmt_ev2_top6 ;
    double fwm2_top6_tr_v3pt30 ;
    double fwm3_top6_tr_v3pt30 ;
    double fwm4_top6_tr_v3pt30 ;
    double fwm5_top6_tr_v3pt30 ;
    double fwm6_top6_tr_v3pt30 ;
    double jmt_ev0_top6_tr_v3pt30 ;
    double jmt_ev1_top6_tr_v3pt30 ;
    double jmt_ev2_top6_tr_v3pt30 ;
    double event_beta_z ;
    int    njets_pt45_eta24 ;
    int    njets_pt30_eta24 ;
    int    njets_pt20_eta24 ;
    int    njets_pt45_eta50 ;
    int    njets_pt30_eta50 ;
    int    njets_pt20_eta50 ;
    double pfht_pt40_eta24 ;
    double pfht_pt45_eta24 ;
    int    nleptons ;
    int    nbtag_csv85_pt30_eta24 ;
    double leppt1 ;
    double m_lep1_b ;
    double leppt2 ;
    double m_lep2_b ;
    int    evt_count ;
    int    run ;
    int    lumi ;
    ULong64_t  event ;

    tt_out->Branch( "mva_train_weight", &mva_train_weight, "mva_train_weight/F" ) ;
    tt_out->Branch( "ds_index", &ds_index, "ds_index/I" ) ;
    tt_out->Branch( "njets_pt45_eta24", &njets_pt45_eta24, "njets_pt45_eta24/I" ) ;
    tt_out->Branch( "njets_pt30_eta24", &njets_pt30_eta24, "njets_pt30_eta24/I" ) ;
    tt_out->Branch( "njets_pt20_eta24", &njets_pt20_eta24, "njets_pt20_eta24/I" ) ;
    tt_out->Branch( "njets_pt45_eta50", &njets_pt45_eta50, "njets_pt45_eta50/I" ) ;
    tt_out->Branch( "njets_pt30_eta50", &njets_pt30_eta50, "njets_pt30_eta50/I" ) ;
    tt_out->Branch( "njets_pt20_eta50", &njets_pt20_eta50, "njets_pt20_eta50/I" ) ;
    tt_out->Branch( "nbtag_csv85_pt30_eta24", &nbtag_csv85_pt30_eta24, "nbtag_csv85_pt30_eta24/I" ) ;
    tt_out->Branch( "pfht_pt40_eta24"       , &pfht_pt40_eta24        , "pfht_pt40_eta24/D" ) ;
    tt_out->Branch( "pfht_pt45_eta24"       , &pfht_pt45_eta24        , "pfht_pt45_eta24/D" ) ;
    tt_out->Branch( "nleptons"              , &nleptons               , "nleptons/I" ) ;
    tt_out->Branch( "leppt1"                , &leppt1                 , "leppt1/D" ) ;
    tt_out->Branch( "m_lep1_b"              , &m_lep1_b               , "m_lep1_b/D" ) ;
    tt_out->Branch( "leppt2"                , &leppt2                 , "leppt2/D" ) ;
    tt_out->Branch( "m_lep2_b"              , &m_lep2_b               , "m_lep2_b/D" ) ;
    tt_out->Branch( "fwm2_top6", &fwm2_top6, "fwm2_top6/D" ) ;
    tt_out->Branch( "fwm3_top6", &fwm3_top6, "fwm3_top6/D" ) ;
    tt_out->Branch( "fwm4_top6", &fwm4_top6, "fwm4_top6/D" ) ;
    tt_out->Branch( "fwm5_top6", &fwm5_top6, "fwm5_top6/D" ) ;
    tt_out->Branch( "fwm6_top6", &fwm6_top6, "fwm6_top6/D" ) ;
    tt_out->Branch( "jmt_ev0_top6", &jmt_ev0_top6, "jmt_ev0_top6/D" ) ;
    tt_out->Branch( "jmt_ev1_top6", &jmt_ev1_top6, "jmt_ev1_top6/D" ) ;
    tt_out->Branch( "jmt_ev2_top6", &jmt_ev2_top6, "jmt_ev2_top6/D" ) ;
    tt_out->Branch( "evt_count", &evt_count, "evt_count/I" ) ;
    tt_out->Branch( "run", &run, "run/I" ) ;
    tt_out->Branch( "lumi", &lumi, "lumi/I" ) ;
    tt_out->Branch( "event", &event, "event/l" ) ;
    tt_out->Branch( "event_beta_z", &event_beta_z, "event_beta_z/D" ) ;

    //--- Loop over events

    Long64_t nevts_ttree = tt_in->GetEntries() ;
    tr.registerDerivedVar("nevts_ttree", nevts_ttree);
    printf("\n\n Number of events in input tree: %lld\n\n", nevts_ttree ) ;

    int modnum(1) ;
    if ( nevts_ttree > 0 ) modnum = nevts_ttree / 100 ;
    if ( modnum <= 0 ) modnum = 1 ;
    tr.registerDerivedVar("modnum", modnum);

    int nsave(0) ;

    MakeMVAVariables makeMVAVariables(false);
    tr.registerFunction( std::move(makeMVAVariables) );

    while( tr.getNextEvent() )
    {
        const auto& cm_jets = tr.getVec<math::RThetaPhiVector>("cm_jets");
        fwm2_top6 = tr.getVar<double>("fwm2_top6");
        fwm3_top6 = tr.getVar<double>("fwm3_top6");
        fwm4_top6 = tr.getVar<double>("fwm4_top6");
        fwm5_top6 = tr.getVar<double>("fwm5_top6");
        fwm6_top6 = tr.getVar<double>("fwm6_top6");
        jmt_ev0_top6 = tr.getVar<double>("jmt_ev0_top6");
        jmt_ev1_top6 = tr.getVar<double>("jmt_ev1_top6");
        jmt_ev2_top6 = tr.getVar<double>("jmt_ev2_top6");
        fwm2_top6_tr_v3pt30 = tr.getVar<double>("fwm2_top6_tr_v3pt30");
        fwm3_top6_tr_v3pt30 = tr.getVar<double>("fwm3_top6_tr_v3pt30");
        fwm4_top6_tr_v3pt30 = tr.getVar<double>("fwm4_top6_tr_v3pt30");
        fwm5_top6_tr_v3pt30 = tr.getVar<double>("fwm5_top6_tr_v3pt30");
        fwm6_top6_tr_v3pt30 = tr.getVar<double>("fwm6_top6_tr_v3pt30");
        jmt_ev0_top6_tr_v3pt30 = tr.getVar<double>("jmt_ev0_top6_tr_v3pt30");
        jmt_ev1_top6_tr_v3pt30 = tr.getVar<double>("jmt_ev1_top6_tr_v3pt30");
        jmt_ev2_top6_tr_v3pt30 = tr.getVar<double>("jmt_ev2_top6_tr_v3pt30");
        event_beta_z = tr.getVar<double>("event_beta_z");
        njets_pt45_eta24 = tr.getVar<int>("njets_pt45_eta24");
        njets_pt30_eta24 = tr.getVar<int>("njets_pt30_eta24");
        njets_pt20_eta24 = tr.getVar<int>("njets_pt20_eta24");
        njets_pt45_eta50 = tr.getVar<int>("njets_pt45_eta50");
        njets_pt30_eta50 = tr.getVar<int>("njets_pt30_eta50");
        njets_pt20_eta50 = tr.getVar<int>("njets_pt20_eta50");
        pfht_pt40_eta24  = tr.getVar<double>("pfht_pt40_eta24");
        pfht_pt45_eta24  = tr.getVar<double>("pfht_pt45_eta24");
        nleptons         = tr.getVar<int>("nleptons");
        nbtag_csv85_pt30_eta24 = tr.getVar<int>("nbtag_csv85_pt30_eta24");
        leppt1    = tr.getVar<double>("leppt1");
        m_lep1_b  = tr.getVar<double>("m_lep1_b");
        leppt2    = tr.getVar<double>("leppt2");
        m_lep2_b  = tr.getVar<double>("m_lep2_b");
        evt_count = tr.getVar<int>("evt_count");
        run       = tr.getVar<int>("run");
        lumi      = tr.getVar<int>("lumi");
        event     = tr.getVar<ULong64_t>("event");
        
        {
            double esum_total(0.) ;
            for ( unsigned int i = 0 ; i < cm_jets.size() ; i ++ ) 
            {
                esum_total += cm_jets[i].R() ;
            } // i
            double esum_total_sq = esum_total * esum_total ;

            for ( unsigned int i = 0 ; i < cm_jets.size() ; i ++ ) 
            {
                double p_i = cm_jets[i].R() ;
                if ( p_i <= 0 ) continue ;

                for ( unsigned int j = 0 ; j < cm_jets.size() ; j ++ ) 
                {
                    double p_j = cm_jets[j].R() ;
                    if ( p_j <= 0 ) continue ;

                    double cosTheta = cm_jets[i].Dot( cm_jets[j] ) / (p_i * p_j) ;
                    double pi_pj_over_etot2 = p_i * p_j / esum_total_sq ;

                    h_costheta_ppweight->Fill( cosTheta, pi_pj_over_etot2 ) ;
                    if ( i != j ) h_costheta_ppweight_noieqj->Fill( cosTheta, pi_pj_over_etot2 ) ;

                } // j
            } // i
        }

        nsave++ ;
        tt_out->Fill() ;

    } // ei

    tt_out->AutoSave() ;

    h_costheta_ppweight->Write() ;
    h_costheta_ppweight_noieqj->Write() ;

    tf_output->Close() ;

    printf("\n\n Done.\n") ;
    if ( nevts_ttree > 0 ) printf("  Saved %9d / %9lld  (%.4f)\n", nsave, nevts_ttree, (1.*nsave)/(1.*nevts_ttree) ) ;
    printf("  ds_weight = %.3f\n", ds_weight ) ;
    printf("  Created %s\n\n", outfile ) ;

    // Cleaning up dynamic memory
    delete tt_in;
    delete tf_output;

} // skim_subset

int main()
{
    //make_mva_training_tree_example("temp/", "rpv_stop_350", "outputfiles/mva-train-rpv_stop_350.root");
    make_mva_training_tree_example("temp/", "TT", "outputfiles/mva-train-example-ttbar.root");
    return 0;
}
