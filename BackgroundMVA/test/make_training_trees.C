#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Framework/BackgroundMVA/include/MakeMVAVariables.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TClassTable.h"
#include "TMath.h"
#include "TSystem.h"
#include "TH1F.h"

#include <iostream>
#include <getopt.h>
#include <string>

std::set<AnaSamples::FileSummary> setFS(const std::string& dataSets, const bool& isCondor)
{
    AnaSamples::SampleSet        ss("sampleSets.cfg", isCondor);
    AnaSamples::SampleCollection sc("sampleCollections.cfg", ss);

    std::map<std::string, std::vector<AnaSamples::FileSummary>> fileMap;
    if(ss[dataSets] != ss.null())
    {
        fileMap[dataSets] = {ss[dataSets]};
        for(const auto& colls : ss[dataSets].getCollections())
        {
            fileMap[colls] = {ss[dataSets]};
        }
    }
    else if(sc[dataSets] != sc.null())
    {
        fileMap[dataSets] = {sc[dataSets]};
        int i = 0;
        for(const auto& fs : sc[dataSets])
        {
            fileMap[sc.getSampleLabels(dataSets)[i++]].push_back(fs);
        }
    }
    std::set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    return vvf;
}

void make_mva_training_tree_example( NTupleReader& tr, TFile* tf_output, const int maxEvts = -1, 
                                     const int totalEvts = -1, int arg_ds_index = 11,
                                     float lumi_times_xsec = ( 3.79 * 35.9 * 1000.0 ),
                                     bool verb = false) 
{
    Long64_t nevts_ttree = tr.getNEntries() ;
    //Long64_t nevts_ttree = totalEvts;
    tr.registerDerivedVar("nevts_ttree", nevts_ttree);
    printf("\n\n Number of events in input tree: %lld\n\n", nevts_ttree ) ;

    float ds_weight = lumi_times_xsec / nevts_ttree ;
    printf("  Dataset weight : %.8f\n\n", ds_weight ) ;

    gSystem->Exec( "mkdir -p outputfiles" ) ;

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
    tt_out->Branch( "fwm2_top6", &fwm2_top6, "fwm2_top6/D" ) ;
    tt_out->Branch( "fwm3_top6", &fwm3_top6, "fwm3_top6/D" ) ;
    tt_out->Branch( "fwm4_top6", &fwm4_top6, "fwm4_top6/D" ) ;
    tt_out->Branch( "fwm5_top6", &fwm5_top6, "fwm5_top6/D" ) ;
    tt_out->Branch( "fwm6_top6", &fwm6_top6, "fwm6_top6/D" ) ;
    tt_out->Branch( "jmt_ev0_top6", &jmt_ev0_top6, "jmt_ev0_top6/D" ) ;
    tt_out->Branch( "jmt_ev1_top6", &jmt_ev1_top6, "jmt_ev1_top6/D" ) ;
    tt_out->Branch( "jmt_ev2_top6", &jmt_ev2_top6, "jmt_ev2_top6/D" ) ;
    tt_out->Branch( "event_beta_z", &event_beta_z, "event_beta_z/D" ) ;
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
    tt_out->Branch( "evt_count", &evt_count, "evt_count/I" ) ;
    tt_out->Branch( "run", &run, "run/I" ) ;
    tt_out->Branch( "lumi", &lumi, "lumi/I" ) ;
    tt_out->Branch( "event", &event, "event/l" ) ;

    //--- Loop over events

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

        // ------------------------
        // -- Print event number
        // -----------------------        

        if(maxEvts != -1 && tr.getEvtNum() >= maxEvts) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // Make cos theta ppweight histos
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

        nsave++ ;
        tt_out->Fill() ;

    } // event loop

    tt_out->AutoSave() ;

    h_costheta_ppweight->Write() ;
    h_costheta_ppweight_noieqj->Write() ;

    tf_output->Close() ;

    printf("\n\n Done.\n") ;
    if ( nevts_ttree > 0 ) printf("  Saved %9d / %9lld  (%.4f)\n", nsave, nevts_ttree, (1.*nsave)/(1.*nevts_ttree) ) ;
    printf("  ds_weight = %.3f\n", ds_weight ) ;

    // Cleaning up dynamic memory
    delete tf_output;

}

int main(int argc, char *argv[])
{
    int opt, option_index = 0;
    bool runOnCondor = false;
    std::string histFile = "", dataSets = "";
    int nFiles = -1, startFile = 0, maxEvts = -1;

    static struct option long_options[] = {
        {"condor",             no_argument, 0, 'c'},
        {"histFile",     required_argument, 0, 'H'},
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
    };

    while((opt = getopt_long(argc, argv, "cH:D:N:M:E:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
            case 'c': runOnCondor      = true;              break;
            case 'H': histFile         = optarg;            break;
            case 'D': dataSets         = optarg;            break;
            case 'N': nFiles           = int(atoi(optarg)); break;
            case 'M': startFile        = int(atoi(optarg)); break;
            case 'E': maxEvts          = int(atoi(optarg)); break;
        }
    }

    if(runOnCondor)
    {
        char thistFile[128];
        sprintf(thistFile, "make_training_trees_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
    }

    std::set<AnaSamples::FileSummary> vvf = setFS(dataSets, runOnCondor); 
    TFile* outfile = TFile::Open(histFile.c_str(), "RECREATE");
    printf("  Created %s\n\n", histFile.c_str() ) ;

    for(const AnaSamples::FileSummary& file : vvf)
    {
        TChain* ch = new TChain( (file.treePath).c_str() );
        file.addFilesToChain(ch, startFile, nFiles);
        NTupleReader tr(ch);
        
        // Loop over all of the events and fill trees
        make_mva_training_tree_example(tr, outfile, maxEvts, file.nEvts);

        // Cleaning up dynamic memory
        delete ch;
            
    }

    //make_mva_training_tree_example("temp/", "rpv_stop_350", "outputfiles/mva-train-rpv_stop_350.root");
    //make_mva_training_tree_example("temp/", "TT", "outputfiles/mva-train-example-ttbar.root");
    return 0;
}
