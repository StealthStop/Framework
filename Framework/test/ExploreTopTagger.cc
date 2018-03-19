#define ExploreTopTagger_cxx
#include "ExploreTopTagger.h"

#include "Utility.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

//mandatory includes to use top tagger
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/include/TTException.h"
#include "SetUpTopTagger.h"

void ExploreTopTagger::InitHistos()
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("HT",new TH1D("HT","HT",60,0,3000));

    // -----------------------
    // make some histograms
    // -----------------------
    my_histos.emplace("myHisto", new TH1D("njets","njets", 20, 0, 20));
    my_histos.emplace("h_met", new TH1D("h_met","h_met", 20, 0, 200));
    my_histos.emplace("h_ht", new TH1D("h_ht","h_ht", 60, 0, 3000));
    my_histos.emplace("h_ntops", new TH1D("h_ntops","h_ntops", 5, 0, 5));
    my_histos.emplace("h_ntops_presel", new TH1D("h_ntops_presel","h_ntops_presel;Ntop;Events", 5, 0, 5));
    my_histos.emplace("h_ntops_3jet", new TH1D("h_ntops_3jet","h_ntops_3jet", 5, 0, 5));
    my_histos.emplace("h_ntops_3jet_presel", new TH1D("h_ntops_3jet_presel","h_ntops_3jet_presel;Ntop(trijet);Events", 5, 0, 5));
    my_histos.emplace("h_ntops_2jet", new TH1D("h_ntops_2jet","h_ntops_2jet", 5, 0, 5));
    my_histos.emplace("h_ntops_1jet", new TH1D("h_ntops_1jet","h_ntops_1jet", 5, 0, 5));
    my_histos.emplace("h_baseline_ntops", new TH1D("h_baseline_ntops","h_baseline_ntops", 5, 0, 5));
    my_histos.emplace("h_baseline_ntops_3jet", new TH1D("h_baseline_ntops_3jet","h_baseline_ntops_3jet", 5, 0, 5));
    my_histos.emplace("h_baseline_ntops_2jet", new TH1D("h_baseline_ntops_2jet","h_baseline_ntops_2jet", 5, 0, 5));
    my_histos.emplace("h_baseline_ntops_1jet", new TH1D("h_baseline_ntops_1jet","h_baseline_ntops_1jet", 5, 0, 5));

    my_histos.emplace("h_dphi_2tops", new TH1D("dphi_2tops","dphi_2tops", 40, -4, 4));

    // -----------------------
    // Histograms without event selection (only presence of 2 hadronic gen tops
    // -----------------------
    my_2d_histos.emplace("h_gentop_pT_daughterDR", new TH2D("h_gentop_pT_daughterDR","h_gentop_pT_daughterDR",50,0,1000,60,0,3));
    my_histos.emplace("h_gentop_pT", new TH1D("h_gentop_pT","h_gentop_pT",50,0,1000));
    my_histos.emplace("h_gentop_pT_type1", new TH1D("h_gentop_pT_type1","h_gentop_pT_type1",50,0,1000));
    my_histos.emplace("h_gentop_pT_type2", new TH1D("h_gentop_pT_type2","h_gentop_pT_type2",50,0,1000));
    my_histos.emplace("h_gentop_pT_type3", new TH1D("h_gentop_pT_type3","h_gentop_pT_type3",50,0,1000));

    my_histos.emplace("h_top_gentop_minDR", new TH1D("h_top_gentop_minDR","h_top_gentop_minDR", 60, 0, 3));
    my_histos.emplace("h_top_gentop_Dpt", new TH1D("h_top_gentop_Dpt","h_top_gentop_Dpt", 50, 0, 5));
    my_2d_histos.emplace("h_top_gentop_minDR_Dpt", new TH2D("h_top_gentop_minDR_Dpt", "h_top_gentop_minDR_Dpt", 60, 0, 3, 50, 0, 5));

    my_histos.emplace("h_top_3jet_gentop_minDR", new TH1D("h_top_3jet_gentop_minDR","h_top_3jet_gentop_minDR", 60, 0, 3));
    my_histos.emplace("h_top_3jet_gentop_Dpt", new TH1D("h_top_3jet_gentop_Dpt","h_top_3jet_gentop_Dpt", 50, 0, 5));
    my_2d_histos.emplace("h_top_3jet_gentop_minDR_Dpt", new TH2D("h_top_3jet_gentop_minDR_Dpt", "h_top_3jet_gentop_minDR_Dpt", 60, 0, 3, 50, 0, 5));
    my_histos.emplace("h_top_2jet_gentop_minDR", new TH1D("h_top_2jet_gentop_minDR","h_top_2jet_gentop_minDR", 60, 0, 3));
    my_histos.emplace("h_top_2jet_gentop_Dpt", new TH1D("h_top_2jet_gentop_Dpt","h_top_2jet_gentop_Dpt", 50, 0, 5));
    my_2d_histos.emplace("h_top_2jet_gentop_minDR_Dpt", new TH2D("h_top_2jet_gentop_minDR_Dpt", "h_top_2jet_gentop_minDR_Dpt", 60, 0, 3, 50, 0, 5));
    my_histos.emplace("h_top_1jet_gentop_minDR", new TH1D("h_top_1jet_gentop_minDR","h_top_1jet_gentop_minDR", 60, 0, 3));
    my_histos.emplace("h_top_1jet_gentop_Dpt", new TH1D("h_top_1jet_gentop_Dpt","h_top_1jet_gentop_Dpt", 50, 0, 5));
    my_2d_histos.emplace("h_top_1jet_gentop_minDR_Dpt", new TH2D("h_top_1jet_gentop_minDR_Dpt", "h_top_1jet_gentop_minDR_Dpt", 60, 0, 3, 50, 0, 5));

    my_histos.emplace("h_top_gentop_minDR_3jet_daughters", new TH1D("h_top_gentop_minDR_3jet_daughters", "h_top_gentop_minDR_3jet_daughters", 60, 0, 3)); 
    my_histos.emplace("h_top_gentop_Dpt_3jet_daughters", new TH1D("h_top_gentop_Dpt_3jet_daughters", "h_top_gentop_Dpt_3jet_daughters", 50, 0, 5)); 
    my_2d_histos.emplace("h_top_gentop_minDR_Dpt_3jet_daughters", new TH2D("h_top_gentop_minDR_Dpt_3jet_daughters", "h_top_gentop_minDR_Dpt_3jet_daughters", 60, 0, 3, 50, 0, 5)); 
    my_histos.emplace("h_top_gentop_topmatch_minDR_3jet_daughters", new TH1D("h_top_gentop_topmatch_minDR_3jet_daughters", "h_top_gentop_topmatch_minDR_3jet_daughters", 60, 0, 3)); 
    my_histos.emplace("h_top_gentop_topmatch_Dpt_3jet_daughters", new TH1D("h_top_gentop_topmatch_Dpt_3jet_daughters", "h_top_gentop_topmatch_Dpt_3jet_daughters", 50, 0, 5)); 
    my_2d_histos.emplace("h_top_gentop_topmatch_minDR_Dpt_3jet_daughters", new TH2D("h_top_gentop_topmatch_minDR_Dpt_3jet_daughters", "h_top_gentop_topmatch_minDR_Dpt_3jet_daughters", 60, 0, 3, 50, 0, 5)); 

    my_histos.emplace("h_top_trijet_n_matched_constituents", new TH1D("h_top_trijet_n_matched_constituents", "h_top_trijet_n_matched_constituents", 4, -0.5, 3.5));
    my_histos.emplace("h_top_trijet_match_n_matched_constituents", new TH1D("h_top_trijet_match_n_matched_constituents", "h_top_trijet_match_n_matched_constituents", 4, -0.5, 3.5));
    my_2d_histos.emplace("h_top_gentop_minDR_Dpt_anymatch", new TH2D("h_top_gentop_minDR_Dpt_anymatch", "h_top_gentop_minDR_Dpt_anymatch", 60, 0, 3, 50, 0, 5));
    my_2d_histos.emplace("h_top_gentop_minDR_Dpt_3match", new TH2D("h_top_gentop_minDR_Dpt_3match", "h_top_gentop_minDR_Dpt_3match", 60, 0, 3, 50, 0, 5));
    my_2d_histos.emplace("h_top_gentop_minDR_Dpt_2match", new TH2D("h_top_gentop_minDR_Dpt_2match", "h_top_gentop_minDR_Dpt_2match", 60, 0, 3, 50, 0, 5));
    my_2d_histos.emplace("h_top_gentop_minDR_Dpt_1match", new TH2D("h_top_gentop_minDR_Dpt_1match", "h_top_gentop_minDR_Dpt_1match", 60, 0, 3, 50, 0, 5));
    my_2d_histos.emplace("h_top_gentop_minDR_Dpt_0match", new TH2D("h_top_gentop_minDR_Dpt_0match", "h_top_gentop_minDR_Dpt_0match", 60, 0, 3, 50, 0, 5));
    
    my_histos.emplace("h_top_gentop_discr_anymatch", new TH1D("h_top_gentop_discr_anymatch", "h_top_gentop_discr_anymatch", 40, 0.8, 1));
    my_histos.emplace("h_top_gentop_discr_3match", new TH1D("h_top_gentop_discr_3match", "h_top_gentop_discr_3match", 40, 0.8, 1));
    my_histos.emplace("h_top_gentop_discr_2match", new TH1D("h_top_gentop_discr_2match", "h_top_gentop_discr_2match", 40, 0.8, 1));
    my_histos.emplace("h_top_gentop_discr_1match", new TH1D("h_top_gentop_discr_1match", "h_top_gentop_discr_1match", 40, 0.8, 1));
    my_histos.emplace("h_top_gentop_discr_0match", new TH1D("h_top_gentop_discr_0match", "h_top_gentop_discr_0match", 40, 0.8, 1));

    my_histos.emplace("h_top_gentop_discr_anymatch_topmatch", new TH1D("h_top_gentop_discr_anymatch_topmatch", "h_top_gentop_discr_anymatch_topmatch", 40, 0.8, 1));
    my_histos.emplace("h_top_gentop_discr_3match_topmatch", new TH1D("h_top_gentop_discr_3match_topmatch", "h_top_gentop_discr_3match_topmatch", 40, 0.8, 1));
    my_histos.emplace("h_top_gentop_discr_2match_topmatch", new TH1D("h_top_gentop_discr_2match_topmatch", "h_top_gentop_discr_2match_topmatch", 40, 0.8, 1));
    my_histos.emplace("h_top_gentop_discr_1match_topmatch", new TH1D("h_top_gentop_discr_1match_topmatch", "h_top_gentop_discr_1match_topmatch", 40, 0.8, 1));
    my_histos.emplace("h_top_gentop_discr_0match_topmatch", new TH1D("h_top_gentop_discr_0match_topmatch", "h_top_gentop_discr_0match_topmatch", 40, 0.8, 1));

    my_histos.emplace("h_gentop_top_minDR", new TH1D("h_gentop_top_minDR",     "h_gentop_top_minDR", 60, 0, 3));
    my_histos.emplace("h_gentop_top_Dpt", new TH1D("h_gentop_top_Dpt",       "h_gentop_top_Dpt", 50, 0, 5));
    my_2d_histos.emplace("h_gentop_top_minDR_Dpt", new TH2D("h_gentop_top_minDR_Dpt", "h_gentop_top_minDR_Dpt", 60, 0, 3, 50, 0, 5));

    my_histos.emplace("h_top_type1_matched_nsub", new TH1D("h_top_type1_matched_nsub","h_top_type1_matched_nsub",50,0,1));
    my_histos.emplace("h_top_type1_unmatched_nsub", new TH1D("h_top_type1_unmatched_nsub","h_top_type1_unmatched_nsub",50,0,1));
    my_histos.emplace("h_top_type1_matched_softdrop", new TH1D("h_top_type1_matched_softdrop","h_top_type1_matched_softdrop",60,0,300));
    my_histos.emplace("h_top_type1_unmatched_softdrop", new TH1D("h_top_type1_unmatched_softdrop","h_top_type1_unmatched_softdrop",60,0,300));

    my_efficiencies.emplace("toptag_eff", new TEfficiency("toptag_eff","Top tagging efficiency;gentop p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_eff_type1", new TEfficiency("toptag_eff_type1","Top tagging efficiency;type1 gentop p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_eff_type2", new TEfficiency("toptag_eff_type2","Top tagging efficiency;type2 gentop p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_eff_type3", new TEfficiency("toptag_eff_type3","Top tagging efficiency;type3 gentop p_T;#epsilon",10,0,1000));

    // -----------------------
    // Histograms after baseline event selection
    // -----------------------
    my_2d_histos.emplace("h_baseline_gentop_pT_daughterDR", new TH2D("h_baseline_gentop_pT_daughterDR","h_baseline_gentop_pT_daughterDR",50,0,1000,60,0,3));
    my_histos.emplace("h_baseline_gentop_pT", new TH1D("h_baseline_gentop_pT","h_baseline_gentop_pT",50,0,1000));
    my_histos.emplace("h_baseline_gentop_pT_type1", new TH1D("h_baseline_gentop_pT_type1","h_baseline_gentop_pT_type1",50,0,1000));
    my_histos.emplace("h_baseline_gentop_pT_type2", new TH1D("h_baseline_gentop_pT_type2","h_baseline_gentop_pT_type2",50,0,1000));
    my_histos.emplace("h_baseline_gentop_pT_type3", new TH1D("h_baseline_gentop_pT_type3","h_baseline_gentop_pT_type3",50,0,1000));

    my_histos.emplace("h_baseline_top_gentop_minDR", new TH1D("h_baseline_top_gentop_minDR","h_baseline_top_gentop_minDR", 60, 0, 3));
    my_histos.emplace("h_baseline_top_gentop_Dpt", new TH1D("h_baseline_top_gentop_Dpt","h_baseline_top_gentop_Dpt", 50, 0, 5));
    my_2d_histos.emplace("h_baseline_top_gentop_minDR_Dpt", new TH2D("h_baseline_top_gentop_minDR_Dpt", "h_baseline_top_gentop_minDR_Dpt", 60, 0, 3, 50, 0, 5));

    my_histos.emplace("h_baseline_top_3jet_gentop_minDR", new TH1D("h_baseline_top_3jet_gentop_minDR","h_baseline_top_3jet_gentop_minDR", 60, 0, 3));
    my_histos.emplace("h_baseline_top_3jet_gentop_Dpt", new TH1D("h_baseline_top_3jet_gentop_Dpt","h_baseline_top_3jet_gentop_Dpt", 50, 0, 5));
    my_2d_histos.emplace("h_baseline_top_3jet_gentop_minDR_Dpt", new TH2D("h_baseline_top_3jet_gentop_minDR_Dpt", "h_baseline_top_3jet_gentop_minDR_Dpt", 60, 0, 3, 50, 0, 5));
    my_histos.emplace("h_baseline_top_2jet_gentop_minDR", new TH1D("h_baseline_top_2jet_gentop_minDR","h_baseline_top_2jet_gentop_minDR", 60, 0, 3));
    my_histos.emplace("h_baseline_top_2jet_gentop_Dpt", new TH1D("h_baseline_top_2jet_gentop_Dpt","h_baseline_top_2jet_gentop_Dpt", 50, 0, 5));
    my_2d_histos.emplace("h_baseline_top_2jet_gentop_minDR_Dpt", new TH2D("h_baseline_top_2jet_gentop_minDR_Dpt", "h_baseline_top_2jet_gentop_minDR_Dpt", 60, 0, 3, 50, 0, 5));
    my_histos.emplace("h_baseline_top_1jet_gentop_minDR", new TH1D("h_baseline_top_1jet_gentop_minDR","h_baseline_top_1jet_gentop_minDR", 60, 0, 3));
    my_histos.emplace("h_baseline_top_1jet_gentop_Dpt", new TH1D("h_baseline_top_1jet_gentop_Dpt","h_baseline_top_1jet_gentop_Dpt", 50, 0, 5));
    my_2d_histos.emplace("h_baseline_top_1jet_gentop_minDR_Dpt", new TH2D("h_baseline_top_1jet_gentop_minDR_Dpt", "h_baseline_top_1jet_gentop_minDR_Dpt", 60, 0, 3, 50, 0, 5));

    my_histos.emplace("h_baseline_top_gentop_minDR_3jet_daughters", new TH1D("h_baseline_top_gentop_minDR_3jet_daughters", "h_baseline_top_gentop_minDR_3jet_daughters", 60, 0, 3)); 
    my_histos.emplace("h_baseline_top_gentop_Dpt_3jet_daughters", new TH1D("h_baseline_top_gentop_Dpt_3jet_daughters", "h_baseline_top_gentop_Dpt_3jet_daughters", 50, 0, 5)); 
    my_2d_histos.emplace("h_baseline_top_gentop_minDR_Dpt_3jet_daughters", new TH2D("h_baseline_top_gentop_minDR_Dpt_3jet_daughters", "h_baseline_top_gentop_minDR_Dpt_3jet_daughters", 60, 0, 3, 50, 0, 5)); 
    my_histos.emplace("h_baseline_top_gentop_topmatch_minDR_3jet_daughters", new TH1D("h_baseline_top_gentop_topmatch_minDR_3jet_daughters", "h_baseline_top_gentop_topmatch_minDR_3jet_daughters", 60, 0, 3)); 
    my_histos.emplace("h_baseline_top_gentop_topmatch_Dpt_3jet_daughters", new TH1D("h_baseline_top_gentop_topmatch_Dpt_3jet_daughters", "h_baseline_top_gentop_topmatch_Dpt_3jet_daughters", 50, 0, 5)); 
    my_2d_histos.emplace("h_baseline_top_gentop_topmatch_minDR_Dpt_3jet_daughters", new TH2D("h_baseline_top_gentop_topmatch_minDR_Dpt_3jet_daughters", "h_baseline_top_gentop_topmatch_minDR_Dpt_3jet_daughters", 60, 0, 3, 50, 0, 5)); 

    my_histos.emplace("h_baseline_top_trijet_n_matched_constituents", new TH1D("h_baseline_top_trijet_n_matched_constituents", "h_baseline_top_trijet_n_matched_constituents", 4, -0.5, 3.5));
    my_histos.emplace("h_baseline_top_trijet_match_n_matched_constituents", new TH1D("h_baseline_top_trijet_match_n_matched_constituents", "h_baseline_top_trijet_match_n_matched_constituents", 4, -0.5, 3.5));
    my_2d_histos.emplace("h_baseline_top_gentop_minDR_Dpt_anymatch", new TH2D("h_baseline_top_gentop_minDR_Dpt_anymatch", "h_baseline_top_gentop_minDR_Dpt_anymatch", 60, 0, 3, 50, 0, 5));
    my_2d_histos.emplace("h_baseline_top_gentop_minDR_Dpt_3match", new TH2D("h_baseline_top_gentop_minDR_Dpt_3match", "h_baseline_top_gentop_minDR_Dpt_3match", 60, 0, 3, 50, 0, 5));
    my_2d_histos.emplace("h_baseline_top_gentop_minDR_Dpt_2match", new TH2D("h_baseline_top_gentop_minDR_Dpt_2match", "h_baseline_top_gentop_minDR_Dpt_2match", 60, 0, 3, 50, 0, 5));
    my_2d_histos.emplace("h_baseline_top_gentop_minDR_Dpt_1match", new TH2D("h_baseline_top_gentop_minDR_Dpt_1match", "h_baseline_top_gentop_minDR_Dpt_1match", 60, 0, 3, 50, 0, 5));
    my_2d_histos.emplace("h_baseline_top_gentop_minDR_Dpt_0match", new TH2D("h_baseline_top_gentop_minDR_Dpt_0match", "h_baseline_top_gentop_minDR_Dpt_0match", 60, 0, 3, 50, 0, 5));

    my_histos.emplace("h_baseline_top_gentop_discr_anymatch", new TH1D("h_baseline_top_gentop_discr_anymatch", "h_baseline_top_gentop_discr_anymatch", 40, 0.8, 1));
    my_histos.emplace("h_baseline_top_gentop_discr_3match", new TH1D("h_baseline_top_gentop_discr_3match", "h_baseline_top_gentop_discr_3match", 40, 0.8, 1));
    my_histos.emplace("h_baseline_top_gentop_discr_2match", new TH1D("h_baseline_top_gentop_discr_2match", "h_baseline_top_gentop_discr_2match", 40, 0.8, 1));
    my_histos.emplace("h_baseline_top_gentop_discr_1match", new TH1D("h_baseline_top_gentop_discr_1match", "h_baseline_top_gentop_discr_1match", 40, 0.8, 1));
    my_histos.emplace("h_baseline_top_gentop_discr_0match", new TH1D("h_baseline_top_gentop_discr_0match", "h_baseline_top_gentop_discr_0match", 40, 0.8, 1));

    my_histos.emplace("h_baseline_top_gentop_discr_anymatch_topmatch", new TH1D("h_baseline_top_gentop_discr_anymatch_topmatch", "h_baseline_top_gentop_discr_anymatch_topmatch", 40, 0.8, 1));
    my_histos.emplace("h_baseline_top_gentop_discr_3match_topmatch", new TH1D("h_baseline_top_gentop_discr_3match_topmatch", "h_baseline_top_gentop_discr_3match_topmatch", 40, 0.8, 1));
    my_histos.emplace("h_baseline_top_gentop_discr_2match_topmatch", new TH1D("h_baseline_top_gentop_discr_2match_topmatch", "h_baseline_top_gentop_discr_2match_topmatch", 40, 0.8, 1));
    my_histos.emplace("h_baseline_top_gentop_discr_1match_topmatch", new TH1D("h_baseline_top_gentop_discr_1match_topmatch", "h_baseline_top_gentop_discr_1match_topmatch", 40, 0.8, 1));
    my_histos.emplace("h_baseline_top_gentop_discr_0match_topmatch", new TH1D("h_baseline_top_gentop_discr_0match_topmatch", "h_baseline_top_gentop_discr_0match_topmatch", 40, 0.8, 1));

    my_histos.emplace("h_baseline_gentop_top_minDR", new TH1D("h_baseline_gentop_top_minDR",     "h_baseline_gentop_top_minDR", 60, 0, 3));
    my_histos.emplace("h_baseline_gentop_top_Dpt", new TH1D("h_baseline_gentop_top_Dpt",       "h_baseline_gentop_top_Dpt", 50, 0, 5));
    my_2d_histos.emplace("h_baseline_gentop_top_minDR_Dpt", new TH2D("h_baseline_gentop_top_minDR_Dpt", "h_baseline_gentop_top_minDR_Dpt", 60, 0, 3, 50, 0, 5));

    my_histos.emplace("h_baseline_top_type1_matched_nsub", new TH1D("h_baseline_top_type1_matched_nsub","h_baseline_top_type1_matched_nsub",50,0,1));
    my_histos.emplace("h_baseline_top_type1_unmatched_nsub", new TH1D("h_baseline_top_type1_unmatched_nsub","h_baseline_top_type1_unmatched_nsub",50,0,1));
    my_histos.emplace("h_baseline_top_type1_matched_softdrop", new TH1D("h_baseline_top_type1_matched_softdrop","h_baseline_top_type1_matched_softdrop",60,0,300));
    my_histos.emplace("h_baseline_top_type1_unmatched_softdrop", new TH1D("h_baseline_top_type1_unmatched_softdrop","h_baseline_top_type1_unmatched_softdrop",60,0,300));

    my_efficiencies.emplace("toptag_eff_baseline", new TEfficiency("toptag_eff_baseline","Top tagging efficiency;gentop p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_eff_type1_baseline", new TEfficiency("toptag_eff_type1_baseline","Top tagging efficiency;type1 gentop p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_eff_type2_baseline", new TEfficiency("toptag_eff_type2_baseline","Top tagging efficiency;type2 gentop p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_eff_type3_baseline", new TEfficiency("toptag_eff_type3_baseline","Top tagging efficiency;type3 gentop p_T;#epsilon",10,0,1000));

    my_efficiencies.emplace("toptag_fakerate_baseline", new TEfficiency("toptag_fakerate_baseline","Top tagging fake rate;reco top p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_fakerate", new TEfficiency("toptag_fakerate","Top tagging fake rate;reco top p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_fakerate_excl_baseline", new TEfficiency("toptag_fakerate_excl_baseline","Top tagging fake rate;reco top p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_fakerate_excl", new TEfficiency("toptag_fakerate_excl","Top tagging fake rate;reco top p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_chitagrate", new TEfficiency("toptag_chitagrate","Tagged top fraction that matches a neutralino;reco top p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_chitagrate_baseline", new TEfficiency("toptag_chitagrate_baseline","Tagged top fraction that matches a neutralino;reco top p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_chitagrate_excl", new TEfficiency("toptag_chitagrate_excl","Tagged top fraction that matches a neutralino and no tops;reco top p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_chitagrate_excl_baseline", new TEfficiency("toptag_chitagrate_excl_baseline","Tagged top fraction that matches a neutralino and no tops;reco top p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_singletrate", new TEfficiency("toptag_singletrate","Tagged top fraction that matches a singlet or singlino;reco top p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_singletrate_baseline", new TEfficiency("toptag_singletrate_baseline","Tagged top fraction that matches a singlet or singlino;reco top p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_singletrate_excl", new TEfficiency("toptag_singletrate_excl","Tagged top fraction that matches a singlet or singlino and no tops;reco top p_T;#epsilon",10,0,1000));
    my_efficiencies.emplace("toptag_singletrate_excl_baseline", new TEfficiency("toptag_singletrate_excl_baseline","Tagged top fraction that matches a singlet or singlino and no tops;reco top p_T;#epsilon",10,0,1000));

    my_efficiencies.emplace("toptag_fully_matched", new TEfficiency("toptag_fully_matched","Both tops matched;top[0] category;top[1] category;#epsilon",3,0.5,3.5,3,0.5,3.5));
    my_efficiencies.emplace("toptag_partially_matched", new TEfficiency("toptag_partially_matched","One top matched;top[0] category;top[1] category;#epsilon",3,0.5,3.5,3,0.5,3.5));
    my_efficiencies.emplace("toptag_unmatched", new TEfficiency("toptag_unmatched","Both tops not matched;top[0] category;top[1] category;#epsilon",3,0.5,3.5,3,0.5,3.5));

    my_efficiencies.emplace("real_fakerate", new TEfficiency("real_fakerate", "Fake rate;reco cand pT[GeV];fake rate",10,0,1000 ));
    my_efficiencies.emplace("real_fakerate_weighted", new TEfficiency("real_fakerate_weighted", "Fake rate;reco cand pT[GeV];fake rate",10,0,1000 ));

    my_2d_histos.emplace("toptag_breakdown", new TH2D("toptag_breakdown","Number of events/category;top[0] category;top[1] category;#Events",3,0.5,3.5,3,0.5,3.5));

    // Histograms of tagger inputs
    my_histos.emplace("h_cand_m", new TH1D("h_cand_m", "h_cand_m", 50, 50, 300));
    my_histos.emplace("h_cand_p", new TH1D("h_cand_p", "h_cand_p", 50, 0, 2000));
    my_histos.emplace("h_j12_m", new TH1D("h_j12_m", "h_j12_m", 50, 0, 250));
    my_histos.emplace("h_j23_m", new TH1D("h_j23_m", "h_j23_m", 50, 0, 200));
    my_histos.emplace("h_j13_m", new TH1D("h_j13_m", "h_j13_m", 50, 0, 200));
    my_histos.emplace("h_dTheta12", new TH1D("h_dTheta12", "h_dTheta12", 60, 1, 4));
    my_histos.emplace("h_dTheta23", new TH1D("h_dTheta23", "h_dTheta23", 60, 0, 3));
    my_histos.emplace("h_dTheta13", new TH1D("h_dTheta13", "h_dTheta13", 60, 1, 4));
    my_histos.emplace("h_j1_m", new TH1D("h_j1_m", "h_j1_m", 50, 0, 100));
    my_histos.emplace("h_j1_p", new TH1D("h_j1_p", "h_j1_p", 50, 0, 250));
    my_histos.emplace("h_j1_QGL", new TH1D("h_j1_QGL", "h_j1_QGL", 50, 0, 1));
    my_histos.emplace("h_j1_CSV", new TH1D("h_j1_CSV", "h_j1_CSV", 50, 0, 1));
    my_histos.emplace("h_j2_m", new TH1D("h_j2_m", "h_j2_m", 50, 0, 100));
    my_histos.emplace("h_j2_p", new TH1D("h_j2_p", "h_j2_p", 50, 0, 250));
    my_histos.emplace("h_j2_QGL", new TH1D("h_j2_QGL", "h_j2_QGL", 50, 0, 1));
    my_histos.emplace("h_j2_CSV", new TH1D("h_j2_CSV", "h_j2_CSV", 50, 0, 1));
    my_histos.emplace("h_j3_m", new TH1D("h_j3_m", "h_j3_m", 50, 0, 100));
    my_histos.emplace("h_j3_p", new TH1D("h_j3_p", "h_j3_p", 50, 0, 250));
    my_histos.emplace("h_j3_QGL", new TH1D("h_j3_QGL", "h_j3_QGL", 50, 0, 1));
    my_histos.emplace("h_j3_CSV", new TH1D("h_j3_CSV", "h_j3_CSV", 50, 0, 1));
    
    my_histos.emplace("h_cand_m_topmatch", new TH1D("h_cand_m_topmatch", "h_cand_m_topmatch", 50, 50, 300));
    my_histos.emplace("h_cand_p_topmatch", new TH1D("h_cand_p_topmatch", "h_cand_p_topmatch", 50, 0, 2000));
    my_histos.emplace("h_j12_m_topmatch", new TH1D("h_j12_m_topmatch", "h_j12_m_topmatch", 50, 0, 250));
    my_histos.emplace("h_j23_m_topmatch", new TH1D("h_j23_m_topmatch", "h_j23_m_topmatch", 50, 0, 200));
    my_histos.emplace("h_j13_m_topmatch", new TH1D("h_j13_m_topmatch", "h_j13_m_topmatch", 50, 0, 200));
    my_histos.emplace("h_dTheta12_topmatch", new TH1D("h_dTheta12_topmatch", "h_dTheta12_topmatch", 60, 1, 4));
    my_histos.emplace("h_dTheta23_topmatch", new TH1D("h_dTheta23_topmatch", "h_dTheta23_topmatch", 60, 0, 3));
    my_histos.emplace("h_dTheta13_topmatch", new TH1D("h_dTheta13_topmatch", "h_dTheta13_topmatch", 60, 1, 4));
    my_histos.emplace("h_j1_m_topmatch", new TH1D("h_j1_m_topmatch", "h_j1_m_topmatch", 50, 0, 100));
    my_histos.emplace("h_j1_p_topmatch", new TH1D("h_j1_p_topmatch", "h_j1_p_topmatch", 50, 0, 250));
    my_histos.emplace("h_j1_QGL_topmatch", new TH1D("h_j1_QGL_topmatch", "h_j1_QGL_topmatch", 50, 0, 1));
    my_histos.emplace("h_j1_CSV_topmatch", new TH1D("h_j1_CSV_topmatch", "h_j1_CSV_topmatch", 50, 0, 1));
    my_histos.emplace("h_j2_m_topmatch", new TH1D("h_j2_m_topmatch", "h_j2_m_topmatch", 50, 0, 100));
    my_histos.emplace("h_j2_p_topmatch", new TH1D("h_j2_p_topmatch", "h_j2_p_topmatch", 50, 0, 250));
    my_histos.emplace("h_j2_QGL_topmatch", new TH1D("h_j2_QGL_topmatch", "h_j2_QGL_topmatch", 50, 0, 1));
    my_histos.emplace("h_j2_CSV_topmatch", new TH1D("h_j2_CSV_topmatch", "h_j2_CSV_topmatch", 50, 0, 1));
    my_histos.emplace("h_j3_m_topmatch", new TH1D("h_j3_m_topmatch", "h_j3_m_topmatch", 50, 0, 100));
    my_histos.emplace("h_j3_p_topmatch", new TH1D("h_j3_p_topmatch", "h_j3_p_topmatch", 50, 0, 250));
    my_histos.emplace("h_j3_QGL_topmatch", new TH1D("h_j3_QGL_topmatch", "h_j3_QGL_topmatch", 50, 0, 1));
    my_histos.emplace("h_j3_CSV_topmatch", new TH1D("h_j3_CSV_topmatch", "h_j3_CSV_topmatch", 50, 0, 1));
    
    my_histos.emplace("h_cand_m_susymatch", new TH1D("h_cand_m_susymatch", "h_cand_m_susymatch", 50, 50, 300));
    my_histos.emplace("h_cand_p_susymatch", new TH1D("h_cand_p_susymatch", "h_cand_p_susymatch", 50, 0, 2000));
    my_histos.emplace("h_j12_m_susymatch", new TH1D("h_j12_m_susymatch", "h_j12_m_susymatch", 50, 0, 250));
    my_histos.emplace("h_j23_m_susymatch", new TH1D("h_j23_m_susymatch", "h_j23_m_susymatch", 50, 0, 200));
    my_histos.emplace("h_j13_m_susymatch", new TH1D("h_j13_m_susymatch", "h_j13_m_susymatch", 50, 0, 200));
    my_histos.emplace("h_dTheta12_susymatch", new TH1D("h_dTheta12_susymatch", "h_dTheta12_susymatch", 60, 1, 4));
    my_histos.emplace("h_dTheta23_susymatch", new TH1D("h_dTheta23_susymatch", "h_dTheta23_susymatch", 60, 0, 3));
    my_histos.emplace("h_dTheta13_susymatch", new TH1D("h_dTheta13_susymatch", "h_dTheta13_susymatch", 60, 1, 4));
    my_histos.emplace("h_j1_m_susymatch", new TH1D("h_j1_m_susymatch", "h_j1_m_susymatch", 50, 0, 100));
    my_histos.emplace("h_j1_p_susymatch", new TH1D("h_j1_p_susymatch", "h_j1_p_susymatch", 50, 0, 250));
    my_histos.emplace("h_j1_QGL_susymatch", new TH1D("h_j1_QGL_susymatch", "h_j1_QGL_susymatch", 50, 0, 1));
    my_histos.emplace("h_j1_CSV_susymatch", new TH1D("h_j1_CSV_susymatch", "h_j1_CSV_susymatch", 50, 0, 1));
    my_histos.emplace("h_j2_m_susymatch", new TH1D("h_j2_m_susymatch", "h_j2_m_susymatch", 50, 0, 100));
    my_histos.emplace("h_j2_p_susymatch", new TH1D("h_j2_p_susymatch", "h_j2_p_susymatch", 50, 0, 250));
    my_histos.emplace("h_j2_QGL_susymatch", new TH1D("h_j2_QGL_susymatch", "h_j2_QGL_susymatch", 50, 0, 1));
    my_histos.emplace("h_j2_CSV_susymatch", new TH1D("h_j2_CSV_susymatch", "h_j2_CSV_susymatch", 50, 0, 1));
    my_histos.emplace("h_j3_m_susymatch", new TH1D("h_j3_m_susymatch", "h_j3_m_susymatch", 50, 0, 100));
    my_histos.emplace("h_j3_p_susymatch", new TH1D("h_j3_p_susymatch", "h_j3_p_susymatch", 50, 0, 250));
    my_histos.emplace("h_j3_QGL_susymatch", new TH1D("h_j3_QGL_susymatch", "h_j3_QGL_susymatch", 50, 0, 1));
    my_histos.emplace("h_j3_CSV_susymatch", new TH1D("h_j3_CSV_susymatch", "h_j3_CSV_susymatch", 50, 0, 1));
    
    my_histos.emplace("h_cand_m_nomatch", new TH1D("h_cand_m_nomatch", "h_cand_m_nomatch", 50, 50, 300));
    my_histos.emplace("h_cand_p_nomatch", new TH1D("h_cand_p_nomatch", "h_cand_p_nomatch", 50, 0, 2000));
    my_histos.emplace("h_j12_m_nomatch", new TH1D("h_j12_m_nomatch", "h_j12_m_nomatch", 50, 0, 250));
    my_histos.emplace("h_j23_m_nomatch", new TH1D("h_j23_m_nomatch", "h_j23_m_nomatch", 50, 0, 200));
    my_histos.emplace("h_j13_m_nomatch", new TH1D("h_j13_m_nomatch", "h_j13_m_nomatch", 50, 0, 200));
    my_histos.emplace("h_dTheta12_nomatch", new TH1D("h_dTheta12_nomatch", "h_dTheta12_nomatch", 60, 1, 4));
    my_histos.emplace("h_dTheta23_nomatch", new TH1D("h_dTheta23_nomatch", "h_dTheta23_nomatch", 60, 0, 3));
    my_histos.emplace("h_dTheta13_nomatch", new TH1D("h_dTheta13_nomatch", "h_dTheta13_nomatch", 60, 1, 4));
    my_histos.emplace("h_j1_m_nomatch", new TH1D("h_j1_m_nomatch", "h_j1_m_nomatch", 50, 0, 100));
    my_histos.emplace("h_j1_p_nomatch", new TH1D("h_j1_p_nomatch", "h_j1_p_nomatch", 50, 0, 250));
    my_histos.emplace("h_j1_QGL_nomatch", new TH1D("h_j1_QGL_nomatch", "h_j1_QGL_nomatch", 50, 0, 1));
    my_histos.emplace("h_j1_CSV_nomatch", new TH1D("h_j1_CSV_nomatch", "h_j1_CSV_nomatch", 50, 0, 1));
    my_histos.emplace("h_j2_m_nomatch", new TH1D("h_j2_m_nomatch", "h_j2_m_nomatch", 50, 0, 100));
    my_histos.emplace("h_j2_p_nomatch", new TH1D("h_j2_p_nomatch", "h_j2_p_nomatch", 50, 0, 250));
    my_histos.emplace("h_j2_QGL_nomatch", new TH1D("h_j2_QGL_nomatch", "h_j2_QGL_nomatch", 50, 0, 1));
    my_histos.emplace("h_j2_CSV_nomatch", new TH1D("h_j2_CSV_nomatch", "h_j2_CSV_nomatch", 50, 0, 1));
    my_histos.emplace("h_j3_m_nomatch", new TH1D("h_j3_m_nomatch", "h_j3_m_nomatch", 50, 0, 100));
    my_histos.emplace("h_j3_p_nomatch", new TH1D("h_j3_p_nomatch", "h_j3_p_nomatch", 50, 0, 250));
    my_histos.emplace("h_j3_QGL_nomatch", new TH1D("h_j3_QGL_nomatch", "h_j3_QGL_nomatch", 50, 0, 1));
    my_histos.emplace("h_j3_CSV_nomatch", new TH1D("h_j3_CSV_nomatch", "h_j3_CSV_nomatch", 50, 0, 1));
    
    // Cut flows
    my_efficiencies.emplace("event_sel", new TEfficiency("event_sel","Event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total", new TEfficiency("event_sel_total","Total event selection efficiency;Cut;#epsilon",8,0,8));

}

void ExploreTopTagger::Loop(std::string runtype, double weight, int maxevents, bool isQuiet)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes_total = 0, nbytes = 0;

   TopTagger tt;
   tt.setCfgFile("TopTagger.cfg");

   TRandom3 rand = TRandom3(123);

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if(maxevents != -1 && jentry >= maxevents) break;

      nbytes = fChain->GetEntry(jentry);   
      nbytes_total += nbytes;

      if ( jentry % (nentries/10) == 0 ) printf("  Event %9llu / %9llu  (%2.0f%%)\n", jentry, nentries, 100*(jentry*1.)/(nentries*1.) ) ;

      // -----------------
      // check for number of hadronic tops at gen level
      // -----------------
      int nhadWs = 0;
      std::vector<TLorentzVector> hadtops;
      std::vector<TLorentzVector> hadWs;
      std::vector<int> hadtops_idx;
      std::vector<std::vector<const TLorentzVector*> > hadtopdaughters;
      std::vector<TLorentzVector> neutralinos;
      std::vector<TLorentzVector> singlets;
      std::vector<TLorentzVector> singlinos;
      for ( unsigned int gpi=0; gpi < GenParticles->size() ; gpi++ ) 
      {
          int pdgid = abs( GenParticles_PdgId->at(gpi) ) ;
          int momid = abs( GenParticles_ParentId->at(gpi) ) ;
          int momidx = GenParticles_ParentIdx->at(gpi);
          int status = GenParticles_Status->at(gpi);
          if(pdgid == 1000022 && (status==22 || status == 52))
          {
              neutralinos.push_back(GenParticles->at(gpi));
          }
          if(pdgid == 5000001 && (status == 22 || status == 52))
          {
              singlinos.push_back(GenParticles->at(gpi));
          }
          if(pdgid == 5000002 && (status == 22 || status == 52))
          {
              singlets.push_back(GenParticles->at(gpi));
          }
          if(status == 23 && momid == 24 && pdgid < 6)
          {
              // Should be the quarks from W decay
              nhadWs++;
              // find the top
              int Wmotherid = GenParticles_ParentId->at(momidx);
              if (abs(Wmotherid) == 6){
                  int Wmotheridx = GenParticles_ParentIdx->at(momidx);
                  std::vector<int>::iterator found = std::find(hadtops_idx.begin(), hadtops_idx.end(), Wmotheridx);
                  if (found != hadtops_idx.end())
                  {
                      // already found before
                      // std::cout << "Found this top before: " << *found << std::endl;
                      int position = distance(hadtops_idx.begin(),found);
                      // add the daughter to the list
                      hadtopdaughters[position].push_back(&(GenParticles->at(gpi)));
                  } else
                  {
                      // not yet found
                      hadtops_idx.push_back(Wmotheridx);
                      hadtops.push_back(GenParticles->at(Wmotheridx));
                      hadWs.push_back(GenParticles->at(momidx));
                      std::vector<const TLorentzVector*> daughters;
                      daughters.push_back(&(GenParticles->at(gpi)));
                      hadtopdaughters.push_back(daughters);
                      //std::cout << "Found a new top at idx " << Wmotheridx << std::endl;
                  }
              }
          } 
      }
      if(neutralinos.size() != 2 && neutralinos.size() != 0)
          std::cout << "Found " << neutralinos.size() << " neutralinos!" << std::endl;
      if(singlinos.size() != 2 && singlinos.size() != 0)
          std::cout << "Found " << singlinos.size() << " singlinos!" << std::endl;
      if(singlets.size() != 2 && singlets.size() != 0)
          std::cout << "Found " << singlets.size() << " singlets!" << std::endl;

      // Now check the b quarks (we only want the ones associated with a hadronic W decay for now)
      for ( unsigned int gpi=0; gpi < GenParticles->size() ; gpi++ ) 
      {
          int pdgid = abs( GenParticles_PdgId->at(gpi) ) ;
          int momid = abs( GenParticles_ParentId->at(gpi) ) ;
          int momidx = GenParticles_ParentIdx->at(gpi);
          int status = GenParticles_Status->at(gpi);
          
          if(status == 23 && momid == 6 && pdgid == 5)
          {
              // found a b quark from top decay, need to add this to the list of daughters
              std::vector<int>::iterator found = std::find(hadtops_idx.begin(), hadtops_idx.end(), momidx);
              if (found != hadtops_idx.end())
              {
                  // already found
                  int position = distance(hadtops_idx.begin(),found);
                  hadtopdaughters[position].push_back(&(GenParticles->at(gpi)));
                  //std::cout << "(b) Found this top before: " << *found << std::endl;
              } 
              //else
              //{
                  // not yet found
                  //std::cout << "(b) Found a new leptonic top at idx " << momidx << std::endl;
              //}
          }
      }

      bool verbose = false;
      if (verbose)
      {
          for (int ht=0; ht<hadtops.size(); ++ht){
              std::cout << "Hadtop index = " << hadtops_idx[ht] << std::endl;
              std::cout << "       daughters: ";
              for (int htd=0; htd<hadtopdaughters[ht].size(); ++htd){
                  std::cout << (*hadtopdaughters[ht][htd]).Pt() << " " ;
              }
              std::cout << std::endl;
          }
      }

      // Only keep events with two hadronic top decays
      if (runtype.find("qcd") == std::string::npos && nhadWs != 4) continue;  

      // Figure out whether the gentop is more similar to a monojet, dijet or trijet reco top
      // Monojet criterion: pT>400, DR(daughter,top)<0.8
      std::vector<int> hadtoptype;
      for (int igentop=0; igentop<hadtops.size(); igentop++ )
      {
          double maxDR_gentop_daughter = 0;
          for(int idaughter=0; idaughter<hadtopdaughters[igentop].size(); idaughter++)
          {
              double DR_gentop_daughter = utility::calcDR(hadtops[igentop].Eta(), hadtopdaughters[igentop][idaughter]->Eta(), hadtops[igentop].Phi(), hadtopdaughters[igentop][idaughter]->Phi());
              if(DR_gentop_daughter>maxDR_gentop_daughter && hadtopdaughters[igentop][idaughter]->Pt()/hadtops[igentop].Pt() > 0.1)
                  maxDR_gentop_daughter = DR_gentop_daughter;
          }
          my_2d_histos["h_gentop_pT_daughterDR"]->Fill(hadtops[igentop].Pt(), maxDR_gentop_daughter);
          
          my_histos["h_gentop_pT"]->Fill(hadtops[igentop].Pt());
          // fully merged case
          if (hadtops[igentop].Pt() > 450 && maxDR_gentop_daughter < 0.8)
          {
              hadtoptype.push_back(1);
              my_histos["h_gentop_pT_type1"]->Fill(hadtops[igentop].Pt());
          }
          else if (hadWs[igentop].Pt() > 250) // merged W case
          {
              hadtoptype.push_back(2);
              my_histos["h_gentop_pT_type2"]->Fill(hadtops[igentop].Pt());
          }
          else // assume everything else would be a resolved top
          {
              hadtoptype.push_back(3);
              my_histos["h_gentop_pT_type3"]->Fill(hadtops[igentop].Pt());
          }
      }

      // ------------------
      // --- TOP TAGGER ---
      // ------------------
      
      // setup variables needed for top tagger
      SetUpTopTagger st(*static_cast<const NtupleClass*> (this) , hadtops, hadtopdaughters);
      std::vector<Constituent> constituents = st.getConstituents();
      
      // run the top tagger
      tt.runTagger(constituents);

      // retrieve the top tagger results object
      const TopTaggerResults& ttr = tt.getResults();

      // get reconstructed top
      const std::vector<TopObject*>& tops = ttr.getTops();
      my_histos["h_ntops"]->Fill(tops.size(), weight);

      // get set of all constituents (i.e. AK4 and AK8 jets) used in one of the tops
      std::set<Constituent const *> usedConstituents = ttr.getUsedConstituents();

      // count number of tops per type
      int ntops_3jet = 0;
      int ntops_2jet=0;
      int ntops_1jet=0;
      for (const TopObject* top : tops)
      {
          if(top->getNConstituents() == 3 )
          {
              ntops_3jet++;
          }
          else if(top->getNConstituents() == 2 )
          {
              ntops_2jet++;
          }
          else if(top->getNConstituents() == 1 )
          {
              ntops_1jet++;
          }
      }

      if (jentry < 10) 
      {
          printf("\tN tops: %ld\n", tops.size());

          // print top properties
          for(const TopObject* top : tops)
          {
              //print basic top properties (top->p() gives a TLorentzVector)
              //N constituents refers to the number of jets included in the top
              //3 for resolved tops 
              //2 for W+jet tops
              //1 for fully merged AK8 tops
              printf("\tTop properties: N constituents: %3d,   Pt: %6.1lf,   Eta: %7.3lf,   Phi: %7.3lf,   Mass: %6.1lf\n", top->getNConstituents(), top->p().Pt(), top->p().Eta(), top->p().Phi(), top->p().M());
              
              //get vector of top constituents 
              const std::vector<Constituent const *>& constituents = top->getConstituents();
              
              //Print properties of individual top constituent jets 
              for(const Constituent* constituent : constituents)
              {
                  printf("\t\tConstituent properties: Constituent type: %3d,   Pt: %6.1lf,   Eta: %7.3lf,   Phi: %7.3lf,   Mass: %6.1lf\n", constituent->getType(), constituent->p().Pt(), constituent->p().Eta(), constituent->p().Phi(), constituent->p().M());
              }        
          }        

          std::cout << "Properties of all used constituents" << std::endl;
          // Print properties of individual top constituent jets 
          for(const Constituent* constituent : usedConstituents)
          {
              printf("\t\tConstituent properties: Constituent type: %3d,   Pt: %6.1lf,   Eta: %7.3lf,   Phi: %7.3lf,   Mass: %6.1lf\n", constituent->getType(), constituent->p().Pt(), constituent->p().Eta(), constituent->p().Phi(), constituent->p().M());
          }        
          
      }


      // -------------------------------
      // -- Basic event selection stuff
      // -------------------------------

      // Check whether event would pass the trigger requirement
      bool passTrigger = true;
      int rec_njet_pt45(0) ;
      int rec_njet_pt20(0) ;
      int rec_njet_pt45_btag(0) ;
      double HT_pt40 = 0.0;
      for ( unsigned int rji=0; rji < Jets->size() ; rji++ ) {
          TLorentzVector jlv( Jets->at(rji) ) ;
          if (abs(jlv.Eta()) > 2.4) continue;
          if ( jlv.Pt() > 20 ) 
              rec_njet_pt20++;
          if (jlv.Pt() > 40)
              HT_pt40 += jlv.Pt();
          if ( jlv.Pt() > 45 ) 
          {
              rec_njet_pt45++ ;
              if ( Jets_bDiscriminatorCSV->at(rji) > 0.8484) 
                  rec_njet_pt45_btag++;
          }
      } 
      if ( !( HT_pt40>500 && rec_njet_pt45>=6 ) ) 
          passTrigger = false;

      bool passLoose = passTrigger && rec_njet_pt45_btag>1;
      bool passBaseline = HT_pt40>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>1 && tops.size()>1;

      if(passLoose)
      {
          my_histos["h_ntops_presel"]->Fill(tops.size(), weight);
      }

      if(passBaseline)
      {
          for (int igentop=0; igentop<hadtops.size(); igentop++ )
          {
              double maxDR_gentop_daughter = 0;
              for(int idaughter=0; idaughter<hadtopdaughters[igentop].size(); idaughter++)
              {
                  double DR_gentop_daughter = utility::calcDR(hadtops[igentop].Eta(), hadtopdaughters[igentop][idaughter]->Eta(), hadtops[igentop].Phi(), hadtopdaughters[igentop][idaughter]->Phi());
                  if(DR_gentop_daughter>maxDR_gentop_daughter && hadtopdaughters[igentop][idaughter]->Pt()/hadtops[igentop].Pt() > 0.1)
                      maxDR_gentop_daughter = DR_gentop_daughter;
              }
              my_2d_histos["h_baseline_gentop_pT_daughterDR"]->Fill(hadtops[igentop].Pt(), maxDR_gentop_daughter);
              
              my_histos["h_baseline_gentop_pT"]->Fill(hadtops[igentop].Pt());
              // fully merged case
              if (hadtops[igentop].Pt() > 450 && maxDR_gentop_daughter < 0.8)
              {
                  my_histos["h_baseline_gentop_pT_type1"]->Fill(hadtops[igentop].Pt());
              }
              else if (hadWs[igentop].Pt() > 250) // merged W case
              {
                  my_histos["h_baseline_gentop_pT_type2"]->Fill(hadtops[igentop].Pt());
              }
              else // assume everything else would be a resolved top
              {
                  my_histos["h_baseline_gentop_pT_type3"]->Fill(hadtops[igentop].Pt());
              }
          }

      }



      // Fill event selection efficiencies
      my_efficiencies["event_sel_total"]->Fill(true,0);
      my_efficiencies["event_sel_total"]->Fill(HT_pt40>500,1);
      my_efficiencies["event_sel_total"]->Fill(HT_pt40>500 && rec_njet_pt45>=6 ,2);
      my_efficiencies["event_sel_total"]->Fill(HT_pt40>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>0 ,3);
      my_efficiencies["event_sel_total"]->Fill(HT_pt40>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>0 && tops.size()>0 ,4);
      my_efficiencies["event_sel_total"]->Fill(HT_pt40>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>0 && tops.size()>0 && rec_njet_pt45_btag>1 ,5);
      my_efficiencies["event_sel_total"]->Fill(HT_pt40>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>0 && tops.size()>0 && rec_njet_pt45_btag>1 && tops.size()>1 ,6);
      my_efficiencies["event_sel_total"]->Fill(HT_pt40>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>0 && tops.size()>0 && rec_njet_pt45_btag>1 && tops.size()>1 && rec_njet_pt20>=8 ,7);
      
      my_efficiencies["event_sel"]->Fill(true,0);
      my_efficiencies["event_sel"]->Fill(HT_pt40>500,1);
      if(HT_pt40>500)
      {
          my_efficiencies["event_sel"]->Fill(rec_njet_pt45>=6,2);
          if (rec_njet_pt45>=6)
          {
              my_efficiencies["event_sel"]->Fill(rec_njet_pt45_btag>0,3);
              if (rec_njet_pt45_btag>0)
              {
                  my_efficiencies["event_sel"]->Fill(tops.size()>0,4);
                  if (tops.size()>0)
                  {
                      my_efficiencies["event_sel"]->Fill(rec_njet_pt45_btag>1,5);
                      if (rec_njet_pt45_btag>1)
                      {
                          my_efficiencies["event_sel"]->Fill(tops.size()>1,6);
                          if (tops.size()>1)
                          {
                              my_efficiencies["event_sel"]->Fill(rec_njet_pt20>=8,7);
                          }
                      }
                  }
              }
          }
      }
 
      my_histos["myHisto"]->Fill(NJets, weight);
      my_histos["h_met"]->Fill(MET, weight);
      my_histos["h_ht"]->Fill(HT, weight);

      // ----------------------------------
      // -- Study top tagger performance --
      // ----------------------------------

      ttUtility::TrijetInputCalculator TIC = ttUtility::TrijetInputCalculator();
      std::vector<std::string> myvars;
      myvars.push_back("cand_m");
      myvars.push_back("cand_p");
      myvars.push_back("j12_m");
      myvars.push_back("j13_m");
      myvars.push_back("j23_m");
      myvars.push_back("dTheta12");
      myvars.push_back("dTheta23");
      myvars.push_back("dTheta13");
      myvars.push_back("j1_m");
      myvars.push_back("j1_p");
      myvars.push_back("j1_QGL");
      myvars.push_back("j1_CSV");
      myvars.push_back("j2_m");
      myvars.push_back("j2_p");
      myvars.push_back("j2_QGL");
      myvars.push_back("j2_CSV");
      myvars.push_back("j3_m");
      myvars.push_back("j3_p");
      myvars.push_back("j3_QGL");
      myvars.push_back("j3_CSV");
      std::vector<float> mydata;
      mydata.resize(myvars.size());
      TIC.mapVars(myvars);
      TIC.setPtr(mydata.data());

      // --- Check input variables for resolved tagger ---
      const std::vector<TopObject> topcandidates = ttr.getTopCandidates();
      for(const TopObject top : topcandidates)
      {
          TIC.calculateVars(top, 0);
          my_histos["h_cand_m"]->Fill(mydata[0], weight);
          my_histos["h_cand_p"]->Fill(mydata[1], weight);
          my_histos["h_j12_m"]->Fill(mydata[2], weight);
          my_histos["h_j23_m"]->Fill(mydata[3], weight);
          my_histos["h_j13_m"]->Fill(mydata[4], weight);
          my_histos["h_dTheta12"]->Fill(mydata[5], weight);
          my_histos["h_dTheta23"]->Fill(mydata[6], weight);
          my_histos["h_dTheta13"]->Fill(mydata[7], weight);
          my_histos["h_j1_m"]->Fill(mydata[8], weight);
          my_histos["h_j1_p"]->Fill(mydata[9], weight);
          my_histos["h_j1_QGL"]->Fill(mydata[10], weight);
          my_histos["h_j1_CSV"]->Fill(mydata[11], weight);
          my_histos["h_j2_m"]->Fill(mydata[12], weight);
          my_histos["h_j2_p"]->Fill(mydata[13], weight);
          my_histos["h_j2_QGL"]->Fill(mydata[14], weight);
          my_histos["h_j2_CSV"]->Fill(mydata[15], weight);
          my_histos["h_j3_m"]->Fill(mydata[16], weight);
          my_histos["h_j3_p"]->Fill(mydata[17], weight);
          my_histos["h_j3_QGL"]->Fill(mydata[18], weight);
          my_histos["h_j3_CSV"]->Fill(mydata[19], weight);

          // Also divide this into categories
          bool matches_top = false;
          bool matches_susy = false;

          for (TLorentzVector hadtop : hadtops)
          {
              double DR_top_gentop = utility::calcDR(top.p().Eta(), hadtop.Eta(), top.p().Phi(), hadtop.Phi());
              double Dpt_top_gentop = abs(top.p().Pt() - hadtop.Pt())/top.p().Pt();

              if(DR_top_gentop<0.4 && Dpt_top_gentop < 0.5)
              {
                  matches_top = true;
                  break;
              }
          }
          for (TLorentzVector hadtop : neutralinos)
          {
              double DR_top_gentop = utility::calcDR(top.p().Eta(), hadtop.Eta(), top.p().Phi(), hadtop.Phi());
              double Dpt_top_gentop = abs(top.p().Pt() - hadtop.Pt())/top.p().Pt();

              if(DR_top_gentop<0.4 && Dpt_top_gentop < 0.5)
              {
                  matches_susy = true;
                  break;
              }
          }
          for (TLorentzVector hadtop : singlets)
          {
              double DR_top_gentop = utility::calcDR(top.p().Eta(), hadtop.Eta(), top.p().Phi(), hadtop.Phi());
              double Dpt_top_gentop = abs(top.p().Pt() - hadtop.Pt())/top.p().Pt();

              if(DR_top_gentop<0.4 && Dpt_top_gentop < 0.5)
              {
                  matches_susy = true;
                  break;
              }
          }
          for (TLorentzVector hadtop : singlinos)
          {
              double DR_top_gentop = utility::calcDR(top.p().Eta(), hadtop.Eta(), top.p().Phi(), hadtop.Phi());
              double Dpt_top_gentop = abs(top.p().Pt() - hadtop.Pt())/top.p().Pt();

              if(DR_top_gentop<0.4 && Dpt_top_gentop < 0.5)
              {
                  matches_susy = true;
                  break;
              }
          }

          if(matches_top)
          {
              my_histos["h_cand_m_topmatch"]->Fill(mydata[0], weight);
              my_histos["h_cand_p_topmatch"]->Fill(mydata[1], weight);
              my_histos["h_j12_m_topmatch"]->Fill(mydata[2], weight);
              my_histos["h_j23_m_topmatch"]->Fill(mydata[3], weight);
              my_histos["h_j13_m_topmatch"]->Fill(mydata[4], weight);
              my_histos["h_dTheta12_topmatch"]->Fill(mydata[5], weight);
              my_histos["h_dTheta23_topmatch"]->Fill(mydata[6], weight);
              my_histos["h_dTheta13_topmatch"]->Fill(mydata[7], weight);
              my_histos["h_j1_m_topmatch"]->Fill(mydata[8], weight);
              my_histos["h_j1_p_topmatch"]->Fill(mydata[9], weight);
              my_histos["h_j1_QGL_topmatch"]->Fill(mydata[10], weight);
              my_histos["h_j1_CSV_topmatch"]->Fill(mydata[11], weight);
              my_histos["h_j2_m_topmatch"]->Fill(mydata[12], weight);
              my_histos["h_j2_p_topmatch"]->Fill(mydata[13], weight);
              my_histos["h_j2_QGL_topmatch"]->Fill(mydata[14], weight);
              my_histos["h_j2_CSV_topmatch"]->Fill(mydata[15], weight);
              my_histos["h_j3_m_topmatch"]->Fill(mydata[16], weight);
              my_histos["h_j3_p_topmatch"]->Fill(mydata[17], weight);
              my_histos["h_j3_QGL_topmatch"]->Fill(mydata[18], weight);
              my_histos["h_j3_CSV_topmatch"]->Fill(mydata[19], weight);
          }
          else if(matches_susy)
          {
              my_histos["h_cand_m_susymatch"]->Fill(mydata[0], weight);
              my_histos["h_cand_p_susymatch"]->Fill(mydata[1], weight);
              my_histos["h_j12_m_susymatch"]->Fill(mydata[2], weight);
              my_histos["h_j23_m_susymatch"]->Fill(mydata[3], weight);
              my_histos["h_j13_m_susymatch"]->Fill(mydata[4], weight);
              my_histos["h_dTheta12_susymatch"]->Fill(mydata[5], weight);
              my_histos["h_dTheta23_susymatch"]->Fill(mydata[6], weight);
              my_histos["h_dTheta13_susymatch"]->Fill(mydata[7], weight);
              my_histos["h_j1_m_susymatch"]->Fill(mydata[8], weight);
              my_histos["h_j1_p_susymatch"]->Fill(mydata[9], weight);
              my_histos["h_j1_QGL_susymatch"]->Fill(mydata[10], weight);
              my_histos["h_j1_CSV_susymatch"]->Fill(mydata[11], weight);
              my_histos["h_j2_m_susymatch"]->Fill(mydata[12], weight);
              my_histos["h_j2_p_susymatch"]->Fill(mydata[13], weight);
              my_histos["h_j2_QGL_susymatch"]->Fill(mydata[14], weight);
              my_histos["h_j2_CSV_susymatch"]->Fill(mydata[15], weight);
              my_histos["h_j3_m_susymatch"]->Fill(mydata[16], weight);
              my_histos["h_j3_p_susymatch"]->Fill(mydata[17], weight);
              my_histos["h_j3_QGL_susymatch"]->Fill(mydata[18], weight);
              my_histos["h_j3_CSV_susymatch"]->Fill(mydata[19], weight);
          }
          else 
          {
              my_histos["h_cand_m_nomatch"]->Fill(mydata[0], weight);
              my_histos["h_cand_p_nomatch"]->Fill(mydata[1], weight);
              my_histos["h_j12_m_nomatch"]->Fill(mydata[2], weight);
              my_histos["h_j23_m_nomatch"]->Fill(mydata[3], weight);
              my_histos["h_j13_m_nomatch"]->Fill(mydata[4], weight);
              my_histos["h_dTheta12_nomatch"]->Fill(mydata[5], weight);
              my_histos["h_dTheta23_nomatch"]->Fill(mydata[6], weight);
              my_histos["h_dTheta13_nomatch"]->Fill(mydata[7], weight);
              my_histos["h_j1_m_nomatch"]->Fill(mydata[8], weight);
              my_histos["h_j1_p_nomatch"]->Fill(mydata[9], weight);
              my_histos["h_j1_QGL_nomatch"]->Fill(mydata[10], weight);
              my_histos["h_j1_CSV_nomatch"]->Fill(mydata[11], weight);
              my_histos["h_j2_m_nomatch"]->Fill(mydata[12], weight);
              my_histos["h_j2_p_nomatch"]->Fill(mydata[13], weight);
              my_histos["h_j2_QGL_nomatch"]->Fill(mydata[14], weight);
              my_histos["h_j2_CSV_nomatch"]->Fill(mydata[15], weight);
              my_histos["h_j3_m_nomatch"]->Fill(mydata[16], weight);
              my_histos["h_j3_p_nomatch"]->Fill(mydata[17], weight);
              my_histos["h_j3_QGL_nomatch"]->Fill(mydata[18], weight);
              my_histos["h_j3_CSV_nomatch"]->Fill(mydata[19], weight);
          }

      }



      // --- Gen matching ---
      int n_matched_recotops = 0;
      int n_matched_other = 0;
      int n_matched_recotops_auto = 0;
      // How often does a tagged top match with genlevel?
      std::vector<bool> top_matches;
      if(hadtops.size() > 0)
      {
          for (const TopObject* top : tops)
          {
              bool matched_top = false;
              bool matched_neutralino = false;
              bool matched_singlet = false;
              bool matched_singlino = false;

              TLorentzVector matched_top_LV;
              std::vector<const TLorentzVector*> matched_top_constituents;
              double minDR = 999;
              for (int i_gentop=0; i_gentop<hadtops.size(); ++i_gentop)
              {
                  double DR_top_gentop = utility::calcDR(top->p().Eta(), hadtops[i_gentop].Eta(), top->p().Phi(), hadtops[i_gentop].Phi());
                  if (DR_top_gentop < minDR)
                  {
                      minDR = DR_top_gentop;
                      matched_top_LV = hadtops[i_gentop];
                      matched_top_constituents = hadtopdaughters[i_gentop];
                  }
              }
              // std::cout << "Top Pt, Eta, Phi: " << top->p().Pt() << " " << top->p().Eta() << " " << top->p().Phi() << std::endl;
              //std::cout << "Gen top Pt, Eta, Phi, DR: " << matched_top.Pt() << " " << matched_top.Eta() << " " << matched_top.Phi() << " " << minDR << std::endl;

              double Dpt_top_gentop = abs(top->p().Pt() - matched_top_LV.Pt())/top->p().Pt();
              my_histos["h_top_gentop_minDR"]->Fill(minDR);
              my_histos["h_top_gentop_Dpt"]->Fill(Dpt_top_gentop);
              my_2d_histos["h_top_gentop_minDR_Dpt"]->Fill(minDR, Dpt_top_gentop);
              matched_top = minDR<0.4 && Dpt_top_gentop < 0.5;
              top_matches.push_back(matched_top);

              double minDR_top_neutralino = 999;
              double Dpt_top_neutralino = -1;
              for (int i_chi=0; i_chi<neutralinos.size(); ++i_chi)
              {
                  double DR_top_neutralino = utility::calcDR(top->p().Eta(), neutralinos[i_chi].Eta(), top->p().Phi(), neutralinos[i_chi].Phi());
                  if(DR_top_neutralino < minDR_top_neutralino)
                  {
                      minDR_top_neutralino = DR_top_neutralino;
                      Dpt_top_neutralino = abs(top->p().Pt() - neutralinos[i_chi].Pt())/top->p().Pt();
                  }
              }
              matched_neutralino = minDR_top_neutralino < 0.4 && Dpt_top_neutralino < 0.5;
              if (verbose && matched_neutralino && matched_top)
              {
                  std::cout << "Matched to a neutralino and to a top quark..." << std::endl;
                  std::cout << "DR and Dpt for neutralino: " << minDR_top_neutralino << ", " << Dpt_top_neutralino << std::endl;
                  std::cout << "DR and Dpt for top:        " << minDR << ", " << Dpt_top_gentop << std::endl;
              }
              my_efficiencies["toptag_chitagrate"]->Fill( matched_neutralino, top->p().Pt());
              my_efficiencies["toptag_chitagrate_excl"]->Fill( matched_neutralino && (!matched_top), top->p().Pt());

              double minDR_top_singlet = 999;
              double Dpt_top_singlet = -1;
              for (int i_chi=0; i_chi<singlets.size(); ++i_chi)
              {
                  double DR_top_singlet = utility::calcDR(top->p().Eta(), singlets[i_chi].Eta(), top->p().Phi(), singlets[i_chi].Phi());
                  if(DR_top_singlet < minDR_top_singlet)
                  {
                      minDR_top_singlet = DR_top_singlet;
                      Dpt_top_singlet = abs(top->p().Pt() - singlets[i_chi].Pt())/top->p().Pt();
                  }
              }
              matched_singlet = minDR_top_singlet < 0.4 && Dpt_top_singlet < 0.5;
              if (verbose && matched_singlet && matched_top)
              {
                  std::cout << "Matched to a singlet and to a top quark..." << std::endl;
                  std::cout << "DR and Dpt for singlet: " << minDR_top_singlet << ", " << Dpt_top_singlet << std::endl;
                  std::cout << "DR and Dpt for top:        " << minDR << ", " << Dpt_top_gentop << std::endl;
              }
              double minDR_top_singlino = 999;
              double Dpt_top_singlino = -1;
              for (int i_chi=0; i_chi<singlinos.size(); ++i_chi)
              {
                  double DR_top_singlino = utility::calcDR(top->p().Eta(), singlinos[i_chi].Eta(), top->p().Phi(), singlinos[i_chi].Phi());
                  if(DR_top_singlino < minDR_top_singlino)
                  {
                      minDR_top_singlino = DR_top_singlino;
                      Dpt_top_singlino = abs(top->p().Pt() - singlinos[i_chi].Pt())/top->p().Pt();
                  }
              }
              matched_singlino = minDR_top_singlino < 0.4 && Dpt_top_singlino < 0.5;
              if (verbose && matched_singlino && matched_top)
              {
                  std::cout << "Matched to a singlino and to a top quark..." << std::endl;
                  std::cout << "DR and Dpt for singlino: " << minDR_top_singlino << ", " << Dpt_top_singlino << std::endl;
                  std::cout << "DR and Dpt for top:        " << minDR << ", " << Dpt_top_gentop << std::endl;
              }
              my_efficiencies["toptag_singletrate"]->Fill( matched_singlet || matched_singlino, top->p().Pt());
              my_efficiencies["toptag_singletrate_excl"]->Fill( (matched_singlet || matched_singlino) && (!matched_top), top->p().Pt());

              // Fake rate defined as not matching a top, neutralino, singlet, or singlino
              my_efficiencies["toptag_fakerate"]->Fill( !matched_top, top->p().Pt());
              my_efficiencies["toptag_fakerate_excl"]->Fill( !( matched_top || matched_neutralino || matched_singlet || matched_singlino ), top->p().Pt());
              if(passBaseline)
              {
                  my_histos["h_baseline_top_gentop_minDR"]->Fill(minDR);
                  my_histos["h_baseline_top_gentop_Dpt"]->Fill(Dpt_top_gentop);
                  my_2d_histos["h_baseline_top_gentop_minDR_Dpt"]->Fill(minDR, Dpt_top_gentop);
                  my_efficiencies["toptag_fakerate_baseline"]->Fill( !matched_top , top->p().Pt());
                  my_efficiencies["toptag_fakerate_excl_baseline"]->Fill( !( matched_top || matched_neutralino || matched_singlet || matched_singlino ) , top->p().Pt());
                  my_efficiencies["toptag_chitagrate_baseline"]->Fill( matched_neutralino, top->p().Pt());
                  my_efficiencies["toptag_chitagrate_excl_baseline"]->Fill( matched_neutralino && (!matched_top), top->p().Pt());
                  my_efficiencies["toptag_singletrate_baseline"]->Fill(matched_singlet || matched_singlino, top->p().Pt());
                  my_efficiencies["toptag_singletrate_excl_baseline"]->Fill( (matched_singlet || matched_singlino) && (!matched_top), top->p().Pt());
              }
              if(matched_top) n_matched_recotops++;
              if((matched_neutralino || matched_singlet || matched_singlino) && (!matched_top))
                  n_matched_other++;
              // Compare with what comes out of the top tagger itself: 
              const TLorentzVector* bestgentop = top->getBestGenTopMatch(0.4);
              if(bestgentop != nullptr) n_matched_recotops_auto++;

              if(top->getNConstituents() == 3 )
              {
                  // do stuff for trijet
                  my_histos["h_top_3jet_gentop_minDR"]->Fill(minDR);
                  my_histos["h_top_3jet_gentop_Dpt"]->Fill(Dpt_top_gentop);
                  my_2d_histos["h_top_3jet_gentop_minDR_Dpt"]->Fill(minDR, Dpt_top_gentop);

                  int n_matched_constituents = 0;
                  for(const Constituent* constituent : top->getConstituents())
                  {
                      double minDR_AK4_daughter = 999;
                      const TLorentzVector *matched_daughter;

                      for(const TLorentzVector* daughter : matched_top_constituents)
                      {
                          // for each AK4, check whether we find a matched daughter
                          double DR_daughter_constituent = utility::calcDR(daughter->Eta(), constituent->p().Eta(), daughter->Phi(), constituent->p().Phi());
                          if(DR_daughter_constituent < minDR_AK4_daughter)
                          {
                              minDR_AK4_daughter = DR_daughter_constituent;
                              matched_daughter = daughter;
                          }
                      }
                      double Dpt_3jet_daughter = abs(constituent->p().Pt() - matched_daughter->Pt())/constituent->p().Pt();
                      my_histos["h_top_gentop_minDR_3jet_daughters"]->Fill(minDR_AK4_daughter);
                      my_histos["h_top_gentop_Dpt_3jet_daughters"]->Fill( Dpt_3jet_daughter );
                      my_2d_histos["h_top_gentop_minDR_Dpt_3jet_daughters"]->Fill(minDR_AK4_daughter, Dpt_3jet_daughter );
                      if(passBaseline)
                      {
                          my_histos["h_baseline_top_gentop_minDR_3jet_daughters"]->Fill(minDR_AK4_daughter);
                          my_histos["h_baseline_top_gentop_Dpt_3jet_daughters"]->Fill( Dpt_3jet_daughter );
                          my_2d_histos["h_baseline_top_gentop_minDR_Dpt_3jet_daughters"]->Fill(minDR_AK4_daughter, Dpt_3jet_daughter );                      
                      }
                      if(minDR<0.4)
                      {
                          my_histos["h_top_gentop_topmatch_minDR_3jet_daughters"]->Fill(minDR_AK4_daughter);
                          my_histos["h_top_gentop_topmatch_Dpt_3jet_daughters"]->Fill( Dpt_3jet_daughter );
                          my_2d_histos["h_top_gentop_topmatch_minDR_Dpt_3jet_daughters"]->Fill(minDR_AK4_daughter, Dpt_3jet_daughter );
                          if(passBaseline)
                          {
                              my_histos["h_baseline_top_gentop_topmatch_minDR_3jet_daughters"]->Fill(minDR_AK4_daughter);
                              my_histos["h_baseline_top_gentop_topmatch_Dpt_3jet_daughters"]->Fill( Dpt_3jet_daughter );
                              my_2d_histos["h_baseline_top_gentop_topmatch_minDR_Dpt_3jet_daughters"]->Fill(minDR_AK4_daughter, Dpt_3jet_daughter );
                          }
                      }

                      if(minDR_AK4_daughter < 0.3 && Dpt_3jet_daughter < 0.5)
                      {
                          n_matched_constituents++;
                          //std::cout << "Found a match for this AK4 jet" << std::endl;
                          //std::cout << "\t matched daughter pt, eta, phi: " << matched_daughter->Pt() << ", " << matched_daughter->Eta() << ", " << matched_daughter->Phi() << std::endl;
                          //std::cout << "\t AK4 jet pt, eta, phi: " << constituent->p().Pt() << ", " << constituent->p().Eta() << ", " << constituent->p().Phi() << std::endl;
                      }
                      else 
                      {
                          //std::cout << "Did not find a match for this AK4 jet" << std::endl;
                          //std::cout << "\t AK4 jet pt, eta, phi: " << constituent->p().Pt() << ", " << constituent->p().Eta() << ", " << constituent->p().Phi() << std::endl;
                      }
                  }
                  //std::cout << "Was able to match " << n_matched_constituents << " AK4 constituents to a genlevel top daughter" << std::endl;
                  my_histos["h_top_trijet_n_matched_constituents"]->Fill(n_matched_constituents);
                  if(matched_top)
                      my_histos["h_top_trijet_match_n_matched_constituents"]->Fill(n_matched_constituents);
                  my_2d_histos["h_top_gentop_minDR_Dpt_anymatch"]->Fill(minDR, Dpt_top_gentop);
                  my_histos["h_top_gentop_discr_anymatch"]->Fill(top->getDiscriminator());
                  if(passBaseline)
                  {
                      my_histos["h_baseline_top_3jet_gentop_minDR"]->Fill(minDR);
                      my_histos["h_baseline_top_3jet_gentop_Dpt"]->Fill(Dpt_top_gentop);
                      my_2d_histos["h_baseline_top_3jet_gentop_minDR_Dpt"]->Fill(minDR, Dpt_top_gentop);

                      my_histos["h_baseline_top_trijet_n_matched_constituents"]->Fill(n_matched_constituents);
                      my_2d_histos["h_baseline_top_gentop_minDR_Dpt_anymatch"]->Fill(minDR, Dpt_top_gentop);
                      my_histos["h_baseline_top_gentop_discr_anymatch"]->Fill(top->getDiscriminator());
                      if(matched_top)
                          my_histos["h_baseline_top_trijet_match_n_matched_constituents"]->Fill(n_matched_constituents);
                  }
                  if(matched_top)
                  {
                      my_histos["h_top_gentop_discr_anymatch_topmatch"]->Fill(top->getDiscriminator());
                      if(passBaseline)
                          my_histos["h_baseline_top_gentop_discr_anymatch_topmatch"]->Fill(top->getDiscriminator());
                  }
                  if(n_matched_constituents == 3)
                  {
                      my_2d_histos["h_top_gentop_minDR_Dpt_3match"]->Fill(minDR, Dpt_top_gentop);
                      my_histos["h_top_gentop_discr_3match"]->Fill(top->getDiscriminator());
                      if(minDR<0.4)
                          my_histos["h_top_gentop_discr_3match_topmatch"]->Fill(top->getDiscriminator());
                  }
                  else if(n_matched_constituents == 2)
                  {
                      my_2d_histos["h_top_gentop_minDR_Dpt_2match"]->Fill(minDR, Dpt_top_gentop);
                      my_histos["h_top_gentop_discr_2match"]->Fill(top->getDiscriminator());
                      if(minDR<0.4)
                          my_histos["h_top_gentop_discr_2match_topmatch"]->Fill(top->getDiscriminator());
                  }
                  else if(n_matched_constituents == 1)
                  {
                      my_2d_histos["h_top_gentop_minDR_Dpt_1match"]->Fill(minDR, Dpt_top_gentop);
                      my_histos["h_top_gentop_discr_1match"]->Fill(top->getDiscriminator());
                      if(minDR<0.4)
                          my_histos["h_top_gentop_discr_1match_topmatch"]->Fill(top->getDiscriminator());
                  }
                  else if(n_matched_constituents == 0)
                  {
                      my_2d_histos["h_top_gentop_minDR_Dpt_0match"]->Fill(minDR, Dpt_top_gentop);
                      my_histos["h_top_gentop_discr_0match"]->Fill(top->getDiscriminator());
                      if(minDR<0.4)
                          my_histos["h_top_gentop_discr_0match_topmatch"]->Fill(top->getDiscriminator());
                  }
                  if(passBaseline)
                  {
                      if(n_matched_constituents == 3)
                      {
                          my_2d_histos["h_baseline_top_gentop_minDR_Dpt_3match"]->Fill(minDR, Dpt_top_gentop);
                          my_histos["h_baseline_top_gentop_discr_3match"]->Fill(top->getDiscriminator());
                          if(minDR<0.4)
                              my_histos["h_baseline_top_gentop_discr_3match_topmatch"]->Fill(top->getDiscriminator());
                      }
                      else if(n_matched_constituents == 2)
                      {
                          my_2d_histos["h_baseline_top_gentop_minDR_Dpt_2match"]->Fill(minDR, Dpt_top_gentop);
                          my_histos["h_baseline_top_gentop_discr_2match"]->Fill(top->getDiscriminator());
                          if(minDR<0.4)
                              my_histos["h_baseline_top_gentop_discr_2match_topmatch"]->Fill(top->getDiscriminator());
                      }
                      else if(n_matched_constituents == 1)
                      {
                          my_2d_histos["h_baseline_top_gentop_minDR_Dpt_1match"]->Fill(minDR, Dpt_top_gentop);
                          my_histos["h_baseline_top_gentop_discr_1match"]->Fill(top->getDiscriminator());
                          if(minDR<0.4)
                              my_histos["h_baseline_top_gentop_discr_1match_topmatch"]->Fill(top->getDiscriminator());
                      }
                      else if(n_matched_constituents == 0)
                      {
                          my_2d_histos["h_baseline_top_gentop_minDR_Dpt_0match"]->Fill(minDR, Dpt_top_gentop);
                          my_histos["h_baseline_top_gentop_discr_0match"]->Fill(top->getDiscriminator());
                          if(minDR<0.4)
                              my_histos["h_baseline_top_gentop_discr_0match_topmatch"]->Fill(top->getDiscriminator());
                      }
                  
                  }
              }
              else if(top->getNConstituents() == 2 )
              {
                  // do stuff for dijet
                  my_histos["h_top_2jet_gentop_minDR"]->Fill(minDR);
                  my_histos["h_top_2jet_gentop_Dpt"]->Fill(Dpt_top_gentop);
                  my_2d_histos["h_top_2jet_gentop_minDR_Dpt"]->Fill(minDR, Dpt_top_gentop);
                  if(passBaseline)
                  {
                      my_histos["h_baseline_top_2jet_gentop_minDR"]->Fill(minDR);
                      my_histos["h_baseline_top_2jet_gentop_Dpt"]->Fill(Dpt_top_gentop);
                      my_2d_histos["h_baseline_top_2jet_gentop_minDR_Dpt"]->Fill(minDR, Dpt_top_gentop);
                  }
              }
              else if(top->getNConstituents() == 1 )
              {
                  my_histos["h_top_1jet_gentop_minDR"]->Fill(minDR);
                  my_histos["h_top_1jet_gentop_Dpt"]->Fill(Dpt_top_gentop);
                  my_2d_histos["h_top_1jet_gentop_minDR_Dpt"]->Fill(minDR, Dpt_top_gentop);
                  if(passBaseline)
                  {
                      my_histos["h_baseline_top_1jet_gentop_minDR"]->Fill(minDR);
                      my_histos["h_baseline_top_1jet_gentop_Dpt"]->Fill(Dpt_top_gentop);
                      my_2d_histos["h_baseline_top_1jet_gentop_minDR_Dpt"]->Fill(minDR, Dpt_top_gentop);
                  }

                  Constituent const *thistop = top->getConstituents()[0];
                  // do stuff for monojet
                  // check mindr, dpt/pt, nsubjettiness
                  //TLorentzVector matched_top;
                  //std::vector<const TLorentzVector*> matched_top_constituents;
                  //double minDR = 999;
                  // need to look at the AK4s that were removed as well. 
                  if (minDR<0.4 && Dpt_top_gentop < 0.5)
                  {
                      my_histos["h_top_type1_matched_nsub"]->Fill(thistop->getTau3()/thistop->getTau2());
                      my_histos["h_top_type1_matched_softdrop"]->Fill(thistop->getSoftDropMass());
                  } 
                  else 
                  {
                      my_histos["h_top_type1_unmatched_nsub"]->Fill(thistop->getTau3()/thistop->getTau2());
                      my_histos["h_top_type1_unmatched_softdrop"]->Fill(thistop->getSoftDropMass());
                  }
                  if (passBaseline)
                  {
                      if (minDR<0.4 && Dpt_top_gentop < 0.5)
                      {
                          my_histos["h_baseline_top_type1_matched_nsub"]->Fill(thistop->getTau3()/thistop->getTau2());
                          my_histos["h_baseline_top_type1_matched_softdrop"]->Fill(thistop->getSoftDropMass());
                      } 
                      else 
                      {
                          my_histos["h_baseline_top_type1_unmatched_nsub"]->Fill(thistop->getTau3()/thistop->getTau2());
                          my_histos["h_baseline_top_type1_unmatched_softdrop"]->Fill(thistop->getSoftDropMass());
                      }                  
                  }
              }
          }
      }
      my_histos["h_ntops_3jet"]->Fill(ntops_3jet, weight);
      if(passLoose)
          my_histos["h_ntops_3jet_presel"]->Fill(ntops_3jet, weight);
      my_histos["h_ntops_2jet"]->Fill(ntops_2jet, weight);
      my_histos["h_ntops_1jet"]->Fill(ntops_1jet, weight);
      if(n_matched_recotops>2)
          std::cout << "Double matched a gentop!" << std::endl;

      //if (n_matched_recotops_auto != n_matched_recotops)
      //    std::cout << "top tagger code found different number of matches" << std::endl;


      // How often is a gentop actually reconstructed as a final reco top
      for (int i_gentop=0; i_gentop<hadtops.size(); ++i_gentop)
      {
          //std::cout << "Gen top Pt, Eta, Phi: " << hadtops[i_gentop].Pt() << " " << hadtops[i_gentop].Eta() << " " << hadtops[i_gentop].Phi()  << std::endl;

          const TopObject* matched_top;
          double minDR = 999;
          for (const TopObject* top : tops)
          {
              double DR_top_gentop = utility::calcDR(top->p().Eta(), hadtops[i_gentop].Eta(), top->p().Phi(), hadtops[i_gentop].Phi());
              if (DR_top_gentop < minDR)
              {
                  minDR = DR_top_gentop;
                  matched_top = top;
              }          
          }
          if (minDR < 999)
          {
              //std::cout << "For this gentop, the closest reco top is at DR= " << minDR << " and has pT, eta, phi of: " << matched_top->p().Pt() << " " << matched_top->p().Eta() << " " << matched_top->p().Phi() << std::endl;
              double Dpt_gentop_top = abs(hadtops[i_gentop].Pt() - matched_top->p().Pt())/hadtops[i_gentop].Pt();
              my_histos["h_gentop_top_minDR"]->Fill(minDR);
              my_histos["h_gentop_top_Dpt"]->Fill(Dpt_gentop_top);
              my_2d_histos["h_gentop_top_minDR_Dpt"]->Fill(minDR,Dpt_gentop_top);

              // Call it a match if dR<0.4 and dPT/PT < 0.5
              my_efficiencies["toptag_eff"]->Fill( (minDR < 0.4 && Dpt_gentop_top < 0.5) , hadtops[i_gentop].Pt());
              if (hadtoptype[i_gentop] == 1)
                  my_efficiencies["toptag_eff_type1"]->Fill( (minDR < 0.4 && Dpt_gentop_top < 0.5) , hadtops[i_gentop].Pt());
              else if (hadtoptype[i_gentop] == 2)
                  my_efficiencies["toptag_eff_type2"]->Fill( (minDR < 0.4 && Dpt_gentop_top < 0.5) , hadtops[i_gentop].Pt());
              else if (hadtoptype[i_gentop] == 3)
                  my_efficiencies["toptag_eff_type3"]->Fill( (minDR < 0.4 && Dpt_gentop_top < 0.5) , hadtops[i_gentop].Pt());

              if(passBaseline)
              {
                  my_histos["h_baseline_gentop_top_minDR"]->Fill(minDR);
                  my_histos["h_baseline_gentop_top_Dpt"]->Fill(Dpt_gentop_top);
                  my_2d_histos["h_baseline_gentop_top_minDR_Dpt"]->Fill(minDR,Dpt_gentop_top);
                  
                  // Call it a match if dR<0.4 and dPT/PT < 0.5
                  my_efficiencies["toptag_eff_baseline"]->Fill( (minDR < 0.4 && Dpt_gentop_top < 0.5) , hadtops[i_gentop].Pt());
                  if (hadtoptype[i_gentop] == 1)
                      my_efficiencies["toptag_eff_type1_baseline"]->Fill( (minDR < 0.4 && Dpt_gentop_top < 0.5) , hadtops[i_gentop].Pt());
                  else if (hadtoptype[i_gentop] == 2)
                      my_efficiencies["toptag_eff_type2_baseline"]->Fill( (minDR < 0.4 && Dpt_gentop_top < 0.5) , hadtops[i_gentop].Pt());
                  else if (hadtoptype[i_gentop] == 3)
                      my_efficiencies["toptag_eff_type3_baseline"]->Fill( (minDR < 0.4 && Dpt_gentop_top < 0.5) , hadtops[i_gentop].Pt());
              }
          }
          else
          {
              //std::cout << "No reco top found" << std::endl;
          }

      }


      if(!passLoose) continue;

      // Now check the actual fake rate (not the impurity)
      // Take a random candidate that is not genmatched and check how often it passes the tagger      
      // Only do this for 3-jet category
      std::vector<TopObject> trijetcandidates;
      for(TopObject mycand : topcandidates)
      {
          if(mycand.getNConstituents() == 3)
              trijetcandidates.push_back(mycand);
      }
      int ncands = trijetcandidates.size();
      //std::cout << "Start fake rate check " << std::endl;
      if(ncands > 0)
      {
          int myrandint = rand.Integer(ncands);
          TopObject mycand = trijetcandidates[myrandint];
          bool found_not_matched = false;
          int ntries = 0;
          //std::cout << "Cand top pT, eta, phi: " << mycand.p().Pt() << " " << mycand.p().Eta() << " " << mycand.p().Phi() << std::endl; 
          while (!found_not_matched && (ntries < 100) && (hadtops.size() > 0))
          {
              ntries++;
              double minDR = 999;
              double Dpt = 0;
              for (int gt=0; gt<hadtops.size(); gt++)
              {
                  //std::cout << "had top pT, eta, phi: " << hadtops[gt].Pt() << " " << hadtops[gt].Eta() << " " << hadtops[gt].Phi() << std::endl; 
                  
                  double DR_top_gentop = utility::calcDR(mycand.p().Eta(), hadtops[gt].Eta(), mycand.p().Phi(), hadtops[gt].Phi());
                  double Dpt_top_gentop = abs(mycand.p().Pt() - hadtops[gt].Pt())/mycand.p().Pt();
                  if (DR_top_gentop < minDR)
                  {
                      minDR = DR_top_gentop;
                      Dpt = Dpt_top_gentop;
                  }
              }
              if(minDR > 0.4 || Dpt > 0.5)
              {
                  found_not_matched = true;
                  //std::cout << "Found a candidate not matching any gentops" << std::endl;
                  }
              else
              {
                  mycand = topcandidates[rand.Integer(ncands)];
                  //std::cout << "updating candidate: " << mycand.p().Pt() << " " << mycand.p().Eta() << " " << mycand.p().Phi() << std::endl;
              }
          }
          
          bool found_match_reco = false;
          for (const TopObject* mytop : tops)
          {
              //std::cout << "reco top: " << mytop->p().Pt() << " " << mytop->p().Eta() << " " << mytop->p().Phi() << std::endl;
              // check whether mytop and mycand match
              double DR_top_recotop = utility::calcDR(mycand.p().Eta(), mytop->p().Eta(), mycand.p().Phi(), mytop->p().Phi());
              double Dpt_top_recotop = abs(mycand.p().Pt() - mytop->p().Pt())/mycand.p().Pt();
              //std::cout << "DR, Dpt: " << DR_top_recotop << " " << Dpt_top_recotop << std::endl;
              if(DR_top_recotop < 0.1 && Dpt_top_recotop < 0.1)
              {
                  found_match_reco = true;
                  //std::cout << "Found matching reco top" << std::endl;
                  break;
              }
          }
          my_efficiencies["real_fakerate"]->Fill(found_match_reco, mycand.p().Pt());
          my_efficiencies["real_fakerate_weighted"]->FillWeighted(found_match_reco, weight, mycand.p().Pt());
      }
      
      if(!passBaseline) continue;
      
      my_histos["h_baseline_ntops"]->Fill(tops.size(), weight);
      my_histos["h_baseline_ntops_3jet"]->Fill(ntops_3jet, weight);
      my_histos["h_baseline_ntops_2jet"]->Fill(ntops_2jet, weight);
      my_histos["h_baseline_ntops_1jet"]->Fill(ntops_1jet, weight);
      
      my_histos["h_dphi_2tops"]->Fill( utility::calcDPhi(tops[0]->p().Phi(), tops[1]->p().Phi()) );

      if(n_matched_other == 0)
      {
          // fully matched or not
          bool fully_matched = false;
          bool partially_matched = false;
          bool unmatched = false;
          if (n_matched_recotops >= 2)
              fully_matched = true;
          else if(n_matched_recotops == 1)
              partially_matched = true;
          else if(n_matched_recotops == 0)
              unmatched = true;
          
          my_efficiencies["toptag_fully_matched"]->Fill(fully_matched, tops[0]->getNConstituents(), tops[1]->getNConstituents());
          my_efficiencies["toptag_partially_matched"]->Fill(partially_matched, tops[0]->getNConstituents(), tops[1]->getNConstituents());
          my_efficiencies["toptag_unmatched"]->Fill(unmatched, tops[0]->getNConstituents(), tops[1]->getNConstituents());
      }
      my_2d_histos["toptag_breakdown"]->Fill(tops[0]->getNConstituents(), tops[1]->getNConstituents(), weight);



   } // end of event loop

}

void ExploreTopTagger::WriteHistos()
{
    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }
    
}
