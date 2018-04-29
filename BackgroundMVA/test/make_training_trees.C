#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include "Framework/BackgroundMVA/include/MakeMVAVariables.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/Baseline.h"

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
    //Long64_t nevts_ttree = tr.getNEntries() ;
    //Long64_t nevts_ttree = totalEvts;
    //tr.registerDerivedVar("nevts_ttree", nevts_ttree);
    //printf("\n\n Number of events in input tree: %lld\n\n", nevts_ttree ) ;

    //float ds_weight = lumi_times_xsec / nevts_ttree ;
    //printf("  Dataset weight : %.8f\n\n", ds_weight ) ;

    gSystem->Exec( "mkdir -p outputfiles" ) ;

    TTree* tt_out = new TTree( "mvatraintt", "MVA training ttree" ) ;

    //--- Extra output histograms
    TH1F* h_costheta_ppweight = new TH1F( "h_costheta_ppweight", "cos(theta_ij) in CM frame, pipj/Esq weight", 110, -1.05, 1.05 ) ;
    TH1F* h_costheta_ppweight_noieqj = new TH1F( "h_costheta_ppweight_noieqj", "cos(theta_ij) in CM frame, pipj/Esq weight, excluding i=j", 110, -1.05, 1.05 ) ;

    //--- New branches for output.
    bool passBaseline1l;
    double Mbl;          
    double fwm2_top6;
    double fwm3_top6;
    double fwm4_top6;
    double fwm5_top6;
    double fwm6_top6;
    double jmt_ev0_top6;
    double jmt_ev1_top6;
    double jmt_ev2_top6;
    double event_beta_z;

    tt_out->Branch( "passBaseline1l", &passBaseline1l, "passBaseline1l/B" ) ;
    tt_out->Branch( "Mbl",            &Mbl,            "Mbl/D" ) ;
    tt_out->Branch( "fwm2_top6", &fwm2_top6, "fwm2_top6/D" ) ;
    tt_out->Branch( "fwm3_top6", &fwm3_top6, "fwm3_top6/D" ) ;
    tt_out->Branch( "fwm4_top6", &fwm4_top6, "fwm4_top6/D" ) ;
    tt_out->Branch( "fwm5_top6", &fwm5_top6, "fwm5_top6/D" ) ;
    tt_out->Branch( "fwm6_top6", &fwm6_top6, "fwm6_top6/D" ) ;
    tt_out->Branch( "jmt_ev0_top6", &jmt_ev0_top6, "jmt_ev0_top6/D" ) ;
    tt_out->Branch( "jmt_ev1_top6", &jmt_ev1_top6, "jmt_ev1_top6/D" ) ;
    tt_out->Branch( "jmt_ev2_top6", &jmt_ev2_top6, "jmt_ev2_top6/D" ) ;
    tt_out->Branch( "event_beta_z", &event_beta_z, "event_beta_z/D" ) ;
    //tt_out->Branch( "njets_pt45_eta24", &njets_pt45_eta24, "njets_pt45_eta24/I" ) ;
    //tt_out->Branch( "njets_pt30_eta24", &njets_pt30_eta24, "njets_pt30_eta24/I" ) ;
    //tt_out->Branch( "njets_pt20_eta24", &njets_pt20_eta24, "njets_pt20_eta24/I" ) ;
    //tt_out->Branch( "njets_pt45_eta50", &njets_pt45_eta50, "njets_pt45_eta50/I" ) ;
    //tt_out->Branch( "njets_pt30_eta50", &njets_pt30_eta50, "njets_pt30_eta50/I" ) ;
    //tt_out->Branch( "njets_pt20_eta50", &njets_pt20_eta50, "njets_pt20_eta50/I" ) ;
    //tt_out->Branch( "nbtag_csv85_pt30_eta24", &nbtag_csv85_pt30_eta24, "nbtag_csv85_pt30_eta24/I" ) ;
    //tt_out->Branch( "pfht_pt40_eta24"       , &pfht_pt40_eta24        , "pfht_pt40_eta24/D" ) ;
    //tt_out->Branch( "pfht_pt45_eta24"       , &pfht_pt45_eta24        , "pfht_pt45_eta24/D" ) ;
    //tt_out->Branch( "nleptons"              , &nleptons               , "nleptons/I" ) ;
    //tt_out->Branch( "leppt1"                , &leppt1                 , "leppt1/D" ) ;
    //tt_out->Branch( "m_lep1_b"              , &m_lep1_b               , "m_lep1_b/D" ) ;
    //tt_out->Branch( "leppt2"                , &leppt2                 , "leppt2/D" ) ;
    //tt_out->Branch( "m_lep2_b"              , &m_lep2_b               , "m_lep2_b/D" ) ;
    //tt_out->Branch( "evt_count", &evt_count, "evt_count/I" ) ;
    //tt_out->Branch( "run", &run, "run/I" ) ;
    //tt_out->Branch( "lumi", &lumi, "lumi/I" ) ;
    //tt_out->Branch( "event", &event, "event/l" ) ;

    //--- Loop over events
    MakeMVAVariables makeMVAVariables(false);
    Jet jet;
    Muon muon;
    Electron electron;
    BJet bjet;
    CommonVariables commonVariables;
    Baseline baseline;
    tr.registerFunction( std::move(makeMVAVariables) );
    tr.registerFunction( std::move(jet) );
    tr.registerFunction( std::move(muon) );
    tr.registerFunction( std::move(electron) );
    tr.registerFunction( std::move(bjet) );
    tr.registerFunction( std::move(commonVariables) );
    tr.registerFunction( std::move(baseline) );

    int nsave = 0;
    while( tr.getNextEvent() )
    {
        const auto& cm_jets = tr.getVec<math::RThetaPhiVector>("cm_jets");

        passBaseline1l = tr.getVar<bool>("passBaseline1l");
        Mbl            = tr.getVar<double>("Mbl");
        fwm2_top6 = tr.getVar<double>("fwm2_top6");
        fwm3_top6 = tr.getVar<double>("fwm3_top6");
        fwm4_top6 = tr.getVar<double>("fwm4_top6");
        fwm5_top6 = tr.getVar<double>("fwm5_top6");
        fwm6_top6 = tr.getVar<double>("fwm6_top6");
        jmt_ev0_top6 = tr.getVar<double>("jmt_ev0_top6");
        jmt_ev1_top6 = tr.getVar<double>("jmt_ev1_top6");
        jmt_ev2_top6 = tr.getVar<double>("jmt_ev2_top6");
        event_beta_z = tr.getVar<double>("event_beta_z");

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
        double weight = file.getWeight(); // not used currently
        std::string runtype = (file.tag.find("Data") != std::string::npos) ? "Data" : "MC";
        std::cout << "Starting loop (in run)" << std::endl;
        printf( "runtype: %s fileWeight: %f nFiles: %i startFile: %i maxEvts: %i \n",runtype.c_str(),weight,nFiles,startFile,maxEvts ); fflush( stdout );
        tr.registerDerivedVar<std::string>("runtype",runtype);
        tr.registerDerivedVar<std::string>("filetag",file.tag);
        tr.registerDerivedVar<double>("etaCut",2.4);
        tr.registerDerivedVar<bool>("blind",true);
        
        // Loop over all of the events and fill trees
        make_mva_training_tree_example(tr, outfile, maxEvts, file.nEvts);

        // Cleaning up dynamic memory
        delete ch;
            
    }

    //make_mva_training_tree_example("temp/", "rpv_stop_350", "outputfiles/mva-train-rpv_stop_350.root");
    //make_mva_training_tree_example("temp/", "TT", "outputfiles/mva-train-example-ttbar.root");
    return 0;
}
