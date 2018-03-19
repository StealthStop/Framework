#define ExploreEventSelection_cxx
#include "ExploreEventSelection.h"

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

// includes for the event shapes
#include "bdt_350to650_fwm10_jmtev_top6.h"
#include "EventShapeVariables.h"
#include "get_cmframe_jets.c"
//#include "fisher_350to650_fwm10_jmtev_top6.h"

void ExploreEventSelection::InitHistos()
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("h_met", new TH1D("h_met","h_met", 20, 0, 200));
    my_histos.emplace("h_ht", new TH1D("h_ht","h_ht", 60, 0, 3000));
    my_histos.emplace("h_ntops", new TH1D("h_ntops","h_ntops", 5, 0, 5));
    my_histos.emplace("h_mbl_2l_test", new TH1D("h_mbl_2l_test","h_mbl_2l_test", 50, 0, 200));
    
    // 0 lepton plots
    // "6j" is the control region only. Only look at data for that region
    // Add cuts on Owen's BDT 
    std::vector<std::string> mycuts_0l {"g6j_HT500_g1b", "g6j_HT500_g1b_0t", "g6j_HT500_g1b_1t", "g6j_HT500_g1b_2t",
            "6j_HT500_g1b_0t", "6j_HT500_g1b_1t1", "6j_HT500_g1b_1t2","6j_HT500_g1b_1t3","6j_HT500_g1b_2t", 
            "g7j_HT500_g1b_0t", "g7j_HT500_g1b_1t1", "g7j_HT500_g1b_1t2", "g7j_HT500_g1b_1t3", "g7j_HT500_g1b_2t"};
    for(std::string mycut : mycuts_0l)
    {
        my_histos.emplace("h_njets_0l_"+mycut, new TH1D(("h_njets_0l_"+mycut).c_str(),("h_njets_0l_"+mycut).c_str(), 19, 0, 19));
        my_histos.emplace("h_ntops_0l_"+mycut, new TH1D(("h_ntops_0l_"+mycut).c_str(),("h_ntops_0l_"+mycut).c_str(), 5, 0, 5));
        my_histos.emplace("h_nb_0l_"+mycut, new TH1D(("h_nb_0l_"+mycut).c_str(),("h_nb_0l_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_0l_"+mycut, new TH1D(("h_HT_0l_"+mycut).c_str(),("h_HT_0l_"+mycut).c_str(), 60, 0, 3000));
        my_histos.emplace("h_bdt_0l_"+mycut, new TH1D(("h_bdt_0l_"+mycut).c_str(),("h_bdt_0l_"+mycut).c_str(), 40, -0.5, 0.5));

        my_2d_histos.emplace("h_njets_bdt_0l_"+mycut, new TH2D(("h_njets_bdt_0l_"+mycut).c_str(),("h_njets_bdt_0l_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));
    }

    // 1 lepton plots
    // attempt to have "6j" be the control region, also include 6j in case that works better
    // TODO: add BDT cuts
    std::vector<std::string> mycuts_1l {"g6j","g6j_g1b", "g6j_g1b_mbl", "g6j_g1b_mbl_0t", "g6j_g1b_mbl_1t1", "g6j_g1b_mbl_1t2", "g6j_g1b_mbl_1t3",
            "6j","6j_g1b", "6j_g1b_mbl", "6j_g1b_mbl_0t", "6j_g1b_mbl_1t1", "6j_g1b_mbl_1t2", "6j_g1b_mbl_1t3",
            "g7j","g7j_g1b", "g7j_g1b_mbl", "g7j_g1b_mbl_0t", "g7j_g1b_mbl_1t1", "g7j_g1b_mbl_1t2", "g7j_g1b_mbl_1t3",
            };
    std::vector<std::string> mycuts_1mu {
            "6j","6j_g1b", "6j_g1b_mbl", "6j_g1b_mbl_0t", "6j_g1b_mbl_1t1", "6j_g1b_mbl_1t2", "6j_g1b_mbl_1t3",
            };
    std::vector<std::string> mycuts_1el {
            "6j","6j_g1b", "6j_g1b_mbl", "6j_g1b_mbl_0t", "6j_g1b_mbl_1t1", "6j_g1b_mbl_1t2", "6j_g1b_mbl_1t3",
            };
    for(std::string mycut : mycuts_1l)
    {
        my_histos.emplace("h_njets_1l_"+mycut, new TH1D(("h_njets_1l_"+mycut).c_str(),("h_njets_1l_"+mycut).c_str(), 15, 0, 15));
        my_histos.emplace("h_ntops_1l_"+mycut, new TH1D(("h_ntops_1l_"+mycut).c_str(),("h_ntops_1l_"+mycut).c_str(), 5, 0, 5));
        my_histos.emplace("h_nb_1l_"+mycut, new TH1D(("h_nb_1l_"+mycut).c_str(),("h_nb_1l_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_1l_"+mycut, new TH1D(("h_HT_1l_"+mycut).c_str(),("h_HT_1l_"+mycut).c_str(), 60, 0, 3000));
        my_histos.emplace("h_mbl_1l_"+mycut, new TH1D(("h_mbl_1l_"+mycut).c_str(),("h_mbl_1l_"+mycut).c_str(), 30, 0, 300));
        my_histos.emplace("h_bdt_1l_"+mycut, new TH1D(("h_bdt_1l_"+mycut).c_str(),("h_bdt_1l_"+mycut).c_str(), 40, -0.5, 0.5));

        my_2d_histos.emplace("h_njets_bdt_1l_"+mycut, new TH2D(("h_njets_bdt_1l_"+mycut).c_str(),("h_njets_bdt_1l_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));
    }
    for(std::string mycut : mycuts_1mu)
    {
        my_histos.emplace("h_mupt_1mu_"+mycut, new TH1D(("h_mupt_1mu_"+mycut).c_str(),("h_mupt_1mu_"+mycut).c_str(), 50, 0, 500));
        my_histos.emplace("h_njets_1mu_"+mycut, new TH1D(("h_njets_1mu_"+mycut).c_str(),("h_njets_1mu_"+mycut).c_str(), 15, 0, 15));
        my_histos.emplace("h_ntops_1mu_"+mycut, new TH1D(("h_ntops_1mu_"+mycut).c_str(),("h_ntops_1mu_"+mycut).c_str(), 5, 0, 5));
        my_histos.emplace("h_nb_1mu_"+mycut, new TH1D(("h_nb_1mu_"+mycut).c_str(),("h_nb_1mu_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_1mu_"+mycut, new TH1D(("h_HT_1mu_"+mycut).c_str(),("h_HT_1mu_"+mycut).c_str(), 60, 0, 3000));
        my_histos.emplace("h_mbl_1mu_"+mycut, new TH1D(("h_mbl_1mu_"+mycut).c_str(),("h_mbl_1mu_"+mycut).c_str(), 30, 0, 300));
        my_histos.emplace("h_bdt_1mu_"+mycut, new TH1D(("h_bdt_1mu_"+mycut).c_str(),("h_bdt_1mu_"+mycut).c_str(), 40, -0.5, 0.5));

        my_2d_histos.emplace("h_njets_bdt_1mu_"+mycut, new TH2D(("h_njets_bdt_1mu_"+mycut).c_str(),("h_njets_bdt_1mu_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));
    }
    for(std::string mycut : mycuts_1el)
    {
        my_histos.emplace("h_elpt_1el_"+mycut, new TH1D(("h_elpt_1el_"+mycut).c_str(),("h_elpt_1el_"+mycut).c_str(), 50, 0, 500));
        my_histos.emplace("h_njets_1el_"+mycut, new TH1D(("h_njets_1el_"+mycut).c_str(),("h_njets_1el_"+mycut).c_str(), 15, 0, 15));
        my_histos.emplace("h_ntops_1el_"+mycut, new TH1D(("h_ntops_1el_"+mycut).c_str(),("h_ntops_1el_"+mycut).c_str(), 5, 0, 5));
        my_histos.emplace("h_nb_1el_"+mycut, new TH1D(("h_nb_1el_"+mycut).c_str(),("h_nb_1el_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_1el_"+mycut, new TH1D(("h_HT_1el_"+mycut).c_str(),("h_HT_1el_"+mycut).c_str(), 60, 0, 3000));
        my_histos.emplace("h_mbl_1el_"+mycut, new TH1D(("h_mbl_1el_"+mycut).c_str(),("h_mbl_1el_"+mycut).c_str(), 30, 0, 300));
        my_histos.emplace("h_bdt_1el_"+mycut, new TH1D(("h_bdt_1el_"+mycut).c_str(),("h_bdt_1el_"+mycut).c_str(), 40, -0.5, 0.5));

        my_2d_histos.emplace("h_njets_bdt_1el_"+mycut, new TH2D(("h_njets_bdt_1el_"+mycut).c_str(),("h_njets_bdt_1el_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));
    }

    
    // 2 lepton plots
    // onZ control region only for now, also check for 1top to see if a fake top changes any behavior
    std::vector<std::string> mycuts_2l {"onZ", "onZ_g1b", "onZ_g1b_nombl", "onZ_g1b_g1t", "onZ_g1b_nombl_g1t", "2b",
            "onZ_g1b_nombl_bdt1","onZ_g1b_nombl_bdt2","onZ_g1b_nombl_bdt3","onZ_g1b_nombl_bdt4"};
    for(std::string mycut : mycuts_2l)
    {
        my_histos.emplace("h_njets_2l_"+mycut, new TH1D(("h_njets_2l_"+mycut).c_str(),("h_njets_2l_"+mycut).c_str(), 15, 0, 15));
        my_histos.emplace("h_ntops_2l_"+mycut, new TH1D(("h_ntops_2l_"+mycut).c_str(),("h_ntops_2l_"+mycut).c_str(), 5, 0, 5));
        my_histos.emplace("h_nb_2l_"+mycut, new TH1D(("h_nb_2l_"+mycut).c_str(),("h_nb_2l_"+mycut).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_2l_"+mycut, new TH1D(("h_HT_2l_"+mycut).c_str(),("h_HT_2l_"+mycut).c_str(), 60, 0, 3000));
        my_histos.emplace("h_bdt_2l_"+mycut, new TH1D(("h_bdt_2l_"+mycut).c_str(),("h_bdt_2l_"+mycut).c_str(), 40, -0.5, 0.5));

        my_2d_histos.emplace("h_njets_bdt_2l_"+mycut, new TH2D(("h_njets_bdt_2l_"+mycut).c_str(),("h_njets_bdt_2l_"+mycut).c_str(), 15, 0, 15, 40, -0.5, 0.5));

    }

    // Cut flows
    my_efficiencies.emplace("event_sel", new TEfficiency("event_sel","Event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total", new TEfficiency("event_sel_total","Total event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_0l", new TEfficiency("event_sel_0l","0 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_0l", new TEfficiency("event_sel_total_0l","Total 0 lepton event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_1l", new TEfficiency("event_sel_1l","1 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_1l", new TEfficiency("event_sel_total_1l","Total 1 lepton event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_2l", new TEfficiency("event_sel_2l","2 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_2l", new TEfficiency("event_sel_total_2l","Total 2 lepton event selection efficiency;Cut;#epsilon",8,0,8));

}

void ExploreEventSelection::Loop(double weight, int maxevents, std::string type, std::string filetag, bool isQuiet)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes_total = 0, nbytes = 0;

   // make a toptagger object
   TopTagger tt;
   tt.setCfgFile("TopTagger.cfg");

   // Set up Event shape BDT
   std::vector<std::string> inputVarNames_top6 ;
   std::vector<double> bdtInputVals_top6 ;

   {
       std::string vname ;
       vname = "fwm2_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm3_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm4_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm5_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm6_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm7_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm8_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm9_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "fwm10_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "jmt_ev0_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "jmt_ev1_top6" ; inputVarNames_top6.push_back( vname ) ;
       vname = "jmt_ev2_top6" ; inputVarNames_top6.push_back( vname ) ;
       
       for ( unsigned int i=0; i < inputVarNames_top6.size() ; i++ ) {
           bdtInputVals_top6.push_back( 0.5 ) ; //--- load vector with dummy values.
       } // i
       
   }
   ReadBDT_350to650_fwm10_jmtev_top6 eventshapeBDT( inputVarNames_top6 ) ;
   //ReadFisher_350to650_fwm10_jmtev_top6 read_fisher_350to650_fwm10_jmtev_top6( inputVarNames_top6 ) ;


   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if(maxevents != -1 && jentry >= maxevents) break;

      nbytes = fChain->GetEntry(jentry);   
      nbytes_total += nbytes;

      if ( jentry % 10000 == 0 ) printf("  Event %9llu\n", jentry ) ;

      // Exclude events with MadGraph HT > 100 from the DY inclusive sample
      if(filetag == "DYJetsToLL_M-50_Incl" && madHT > 100) continue;

      // Make sure event weight is not 0 for data
      double eventweight = 1.;
      if(type != "Data")
          eventweight = Weight;
      
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
      if(type != "Data")
      {
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
      }

      // ------------------------------
      // -- Trigger for data
      // ------------------------------
      
      bool passTriggerAllHad = PassTriggerAllHad();
      bool passTriggerMuon = PassTriggerMuon();
      bool passTriggerElectron = PassTriggerElectron();
      if (type == "Data")
      {
          if (filetag == "Data_JetHT" && !passTriggerAllHad) continue;
          if (filetag == "Data_SingleMuon" && !passTriggerMuon) continue;
          if (filetag == "Data_SingleElectron" && !passTriggerElectron) continue;
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
      my_histos["h_ntops"]->Fill(tops.size(), eventweight);

      // get set of all constituents (i.e. AK4 and AK8 jets) used in one of the tops
      std::set<Constituent const *> usedConstituents = ttr.getUsedConstituents();

      // count number of tops per type
      int ntops_3jet=0;
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


      // -------------------------------
      // -- Event shape BDT
      // -------------------------------

      std::vector<math::RThetaPhiVector> cm_frame_jets ;
      get_cmframe_jets( Jets, cm_frame_jets, 6 ) ;
      EventShapeVariables esv_top6( cm_frame_jets ) ;
      TVectorD eigen_vals_norm_top6 = esv_top6.getEigenValues() ;

      {
          int vi(0) ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(2) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(3) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(4) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(5) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(6) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(7) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(8) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(9) ; vi++ ;
          bdtInputVals_top6.at(vi) = esv_top6.getFWmoment(10) ; vi++ ;
          bdtInputVals_top6.at(vi) = eigen_vals_norm_top6[0] ; vi++ ;
          bdtInputVals_top6.at(vi) = eigen_vals_norm_top6[1] ; vi++ ;
          bdtInputVals_top6.at(vi) = eigen_vals_norm_top6[2] ; vi++ ;
      }

      double eventshape_bdt_val = eventshapeBDT.GetMvaValue( bdtInputVals_top6 ) ;
      //double fisher_val = read_fisher_350to650_fwm10_jmtev_top6.GetMvaValue( bdtInputVals_top6 ) ;



      // -------------------------------
      // -- Basic event selection stuff
      // -------------------------------

      // Count jets & bjets
      int rec_njet_pt45(0) ;
      int rec_njet_pt30(0) ;
      int rec_njet_pt30_btag(0) ;
      int rec_njet_pt45_btag(0) ;
      double HT_trigger = 0.0;
      std::vector<TLorentzVector> rec_bjets_pt30;
      for ( unsigned int rji=0; rji < Jets->size() ; rji++ ) {
          TLorentzVector jlv( Jets->at(rji) ) ;
          if (abs(jlv.Eta()) > 2.4) continue;
          if ( jlv.Pt() > 30 )
          { 
              rec_njet_pt30++;
              if ( Jets_bDiscriminatorCSV->at(rji) > 0.8484) 
              {
                  rec_njet_pt30_btag++;
                  rec_bjets_pt30.push_back(jlv);
              }
          }
          if (jlv.Pt() > 40)
              HT_trigger += jlv.Pt();
          if ( jlv.Pt() > 45 ) 
          {
              rec_njet_pt45++ ;
              if ( Jets_bDiscriminatorCSV->at(rji) > 0.8484) 
                  rec_njet_pt45_btag++;
          }
      } 

      // Count leptons > 30 GeV
      std::vector<TLorentzVector> rec_muon_pt30;
      std::vector<int> rec_charge_muon_pt30;
      for (unsigned int imu = 0; imu < Muons->size(); ++imu)
      {
          TLorentzVector lvmu(Muons->at(imu));
          if( abs(lvmu.Eta()) < 2.4 && lvmu.Pt() > 30 && Muons_passIso->at(imu))
          {
              rec_muon_pt30.push_back(lvmu);
              rec_charge_muon_pt30.push_back(Muons_charge->at(imu));
          }
      }
      std::vector<TLorentzVector> rec_electron_pt30;
      std::vector<int> rec_charge_electron_pt30;
      for (unsigned int iel = 0; iel < Electrons->size(); ++iel)
      {
          TLorentzVector lvel(Electrons->at(iel));
          if( abs(lvel.Eta()) < 2.4 && lvel.Pt() > 30 && Electrons_tightID->at(iel) && Electrons_passIso->at(iel))
          {
              rec_electron_pt30.push_back(lvel);
              rec_charge_electron_pt30.push_back(Electrons_charge->at(iel));
          }
      }
      int nleptons = rec_muon_pt30.size() + rec_electron_pt30.size();
      bool onZ = false;
      bool passMbl_2l = false;
      if ( nleptons == 2 )
      {
          if ( (rec_muon_pt30.size() == 2) && (rec_charge_muon_pt30[0] != rec_charge_muon_pt30[1]) )
          {
              double mll = (rec_muon_pt30[0] + rec_muon_pt30[1]).M();
              if( mll > 81 && mll < 101)
                  onZ = true; 
              // check whether a bl pair passes the M(b,l) cut
              for (TLorentzVector myb : rec_bjets_pt30)
              {
                  double mass_bl_1 = (rec_muon_pt30[0] + myb).M();
                  if(mass_bl_1 < 180 && mass_bl_1 > 30)
                      passMbl_2l = true;
                  double mass_bl_2 = (rec_muon_pt30[1] + myb).M();
                  if(mass_bl_2 < 180 && mass_bl_2 > 30)
                      passMbl_2l = true;
                  my_histos["h_mbl_2l_test"]->Fill(mass_bl_1,eventweight);
                  my_histos["h_mbl_2l_test"]->Fill(mass_bl_2,eventweight);
              }
          } 
          else if ( (rec_electron_pt30.size() == 2) && (rec_charge_electron_pt30[0] != rec_charge_electron_pt30[1]) )
          {
              double mll = (rec_electron_pt30[0] + rec_electron_pt30[1]).M();
              if( mll > 81 && mll < 101)
                  onZ = true;  
              // check whether a bl pair passes the M(b,l) cut
              for (TLorentzVector myb : rec_bjets_pt30)
              {
                  double mass_bl_1 = (rec_electron_pt30[0] + myb).M();
                  if(mass_bl_1 < 180 && mass_bl_1 > 30)
                      passMbl_2l = true;
                  double mass_bl_2 = (rec_electron_pt30[1] + myb).M();
                  if(mass_bl_2 < 180 && mass_bl_2 > 30)
                      passMbl_2l = true;
                  my_histos["h_mbl_2l_test"]->Fill(mass_bl_1,eventweight);
                  my_histos["h_mbl_2l_test"]->Fill(mass_bl_2,eventweight);
              }
          }
      }


      bool passBaseline0l = nleptons==0 && rec_njet_pt45>=6 && HT_trigger > 500 && rec_njet_pt45_btag >= 1;
      bool passBaseline1l = nleptons==1 && rec_njet_pt30>=6 ;
      bool passBaseline1mu = rec_muon_pt30.size()==1 && rec_njet_pt30>=6 ;
      bool passBaseline1el = rec_electron_pt30.size()==1 && rec_njet_pt30>=6 ;
      bool passBaseline2l = nleptons==2;
      if(type == "Data")
      {
          passBaseline0l = passBaseline0l && passTriggerAllHad && (filetag == "Data_JetHT");
          if (rec_muon_pt30.size() > 0)
          {
              passBaseline1l = passBaseline1l && passTriggerMuon && (filetag == "Data_SingleMuon");
              passBaseline2l = passBaseline2l && passTriggerMuon && (filetag == "Data_SingleMuon");
          } 
          else if (rec_electron_pt30.size() > 0)
          {
              passBaseline1l = passBaseline1l && passTriggerElectron && (filetag == "Data_SingleElectron");
              passBaseline2l = passBaseline2l && passTriggerElectron && (filetag == "Data_SingleElectron");
          }
      }
      bool pass_g1b = rec_njet_pt30_btag >= 1;
      bool pass_0t = tops.size()==0, pass_1t = tops.size()==1, pass_2t = tops.size()==2;
      bool pass_1t1 = tops.size()==1 && ntops_1jet==1, pass_1t2 = tops.size()==1 && ntops_2jet==1, pass_1t3 = tops.size()==1 && ntops_3jet==1;
      double mbl = -1;
      TLorentzVector used_bjet;
      double mblmet = -1;
      TLorentzVector metlv;
      metlv.SetPtEtaPhiM(MET, 0, METPhi, 0);
      if(nleptons == 1 && pass_g1b)
      {
          TLorentzVector mylepton = (rec_electron_pt30.size() == 1) ? rec_electron_pt30[0] : rec_muon_pt30[0];
          bool passMtop = false;
          //std::cout << "found lepton and " << rec_njet_pt30_btag << " bjets" << std::endl;
          for (TLorentzVector myb : rec_bjets_pt30)
          {
              double mass_bl = (mylepton + myb).M();
              double mass_blmet = (mylepton + myb + metlv).M();
              //std::cout << "mbl and mblmet are " << mass_bl << " and " << mass_blmet << std::endl;
              if (mbl == -1)
              {
                  mbl = mass_bl;
                  mblmet = mass_blmet;
                  used_bjet = myb;
              }
              else if( abs(mass_bl-172.5) < abs(mbl-172.5))
              {
                  mbl = mass_bl;
                  mblmet = mass_blmet;
                  used_bjet = myb;
              }
          }
      }
      bool pass_mbl = mbl > 30 && mbl < 180;
      // Now check that used_bjet isn't also used for the top tagger
      if (pass_mbl)
      {
          int top_type1_to_remove = 0;
          int top_type2_to_remove = 0;
          int top_type3_to_remove = 0;
          for(const TopObject* mytop : tops)
          {
              const std::vector<Constituent const *> mytop_constituents = mytop->getConstituents();
              bool usedup = false;
              for(const Constituent* c: mytop_constituents)
              {
                  if(c->p() == used_bjet)
                  {
                      usedup = true;
                      //std::cout << "Already used this b for the leptonic top" << std::endl;
                  }
              }
              if (usedup)
              {
                  if (mytop_constituents.size() == 1)
                      top_type1_to_remove++;
                  else if(mytop_constituents.size() == 2)
                      top_type2_to_remove++;
                  else if(mytop_constituents.size() == 3)
                      top_type3_to_remove++;
              }
          }
          int ntops_to_remove = top_type1_to_remove + top_type2_to_remove + top_type3_to_remove;
          //std::cout << "Old top counting: " << pass_0t << " " << pass_1t << " " << pass_2t << " " << pass_1t1 << " " << pass_1t2 << " " << pass_1t3 << std::endl;
          if (ntops_to_remove > 0)
          {
              pass_0t = (tops.size() - ntops_to_remove) == 0;
              pass_1t = (tops.size() - ntops_to_remove) == 1;
              pass_2t = (tops.size() - ntops_to_remove) == 2;
              pass_1t1 = (ntops_1jet - top_type1_to_remove) == 1;
              pass_1t2 = (ntops_2jet - top_type2_to_remove) == 1;
              pass_1t3 = (ntops_3jet - top_type3_to_remove) == 1;
          }
          //std::cout << "New top counting: " << pass_0t << " " << pass_1t << " " << pass_2t << " " << pass_1t1 << " " << pass_1t2 << " " << pass_1t3 << std::endl;
      }

      bool bdt_bin1 = eventshape_bdt_val > -1.   && eventshape_bdt_val <= -0.04;
      bool bdt_bin2 = eventshape_bdt_val > -0.04 && eventshape_bdt_val <= 0;
      bool bdt_bin3 = eventshape_bdt_val > 0     && eventshape_bdt_val <= 0.04;
      bool bdt_bin4 = eventshape_bdt_val > 0.04  && eventshape_bdt_val <= 1;

      
      const std::map<std::string, bool> cut_map_0l {
          {"g6j_HT500_g1b", passBaseline0l},
          {"g6j_HT500_g1b_0t", passBaseline0l && pass_0t},
          {"g6j_HT500_g1b_1t", passBaseline0l && pass_1t},
          {"g6j_HT500_g1b_2t", passBaseline0l && pass_2t},
          {"6j_HT500_g1b_0t",  passBaseline0l && rec_njet_pt30==6 && pass_0t},
          {"6j_HT500_g1b_1t1", passBaseline0l && rec_njet_pt30==6 && pass_1t1},
          {"6j_HT500_g1b_1t2", passBaseline0l && rec_njet_pt30==6 && pass_1t2},
          {"6j_HT500_g1b_1t3", passBaseline0l && rec_njet_pt30==6 && pass_1t3},
          {"6j_HT500_g1b_2t",  passBaseline0l && rec_njet_pt30==6 && pass_2t},
          {"g7j_HT500_g1b_0t",  passBaseline0l && rec_njet_pt30>=7 && pass_0t},
          {"g7j_HT500_g1b_1t1", passBaseline0l && rec_njet_pt30>=7 && pass_1t1},
          {"g7j_HT500_g1b_1t2", passBaseline0l && rec_njet_pt30>=7 && pass_1t2},
          {"g7j_HT500_g1b_1t3", passBaseline0l && rec_njet_pt30>=7 && pass_1t3},
          {"g7j_HT500_g1b_2t",  passBaseline0l && rec_njet_pt30>=7 && pass_2t}
      };

      for(auto& kv : cut_map_0l)
      {
          if(kv.second)
          {
              my_histos["h_njets_0l_"+kv.first]->Fill(rec_njet_pt30, eventweight);
              my_histos["h_ntops_0l_"+kv.first]->Fill(tops.size(), eventweight);
              my_histos["h_nb_0l_"+kv.first]->Fill(rec_njet_pt30_btag, eventweight);
              my_histos["h_HT_0l_"+kv.first]->Fill(HT_trigger, eventweight);
              my_histos["h_bdt_0l_"+kv.first]->Fill(eventshape_bdt_val, eventweight);
              my_2d_histos["h_njets_bdt_0l_"+kv.first]->Fill(rec_njet_pt30, eventshape_bdt_val, eventweight);
          }
      }

      const std::map<std::string, bool> cut_map_1l {
          {"g6j", passBaseline1l},
          {"g6j_g1b", passBaseline1l && pass_g1b},
          {"g6j_g1b_mbl", passBaseline1l && pass_g1b && pass_mbl},
          {"g6j_g1b_mbl_0t",  passBaseline1l && pass_g1b && pass_mbl && pass_0t},
          {"g6j_g1b_mbl_1t1", passBaseline1l && pass_g1b && pass_mbl && pass_1t1},
          {"g6j_g1b_mbl_1t2", passBaseline1l && pass_g1b && pass_mbl && pass_1t2},
          {"g6j_g1b_mbl_1t3", passBaseline1l && pass_g1b && pass_mbl && pass_1t3},
          {"6j", passBaseline1l && rec_njet_pt30==6},
          {"6j_g1b", passBaseline1l && rec_njet_pt30==6 && pass_g1b},
          {"6j_g1b_mbl", passBaseline1l && rec_njet_pt30==6 && pass_g1b && pass_mbl},
          {"6j_g1b_mbl_0t",  passBaseline1l && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_0t},
          {"6j_g1b_mbl_1t1", passBaseline1l && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_1t1},
          {"6j_g1b_mbl_1t2", passBaseline1l && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_1t2},
          {"6j_g1b_mbl_1t3", passBaseline1l && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_1t3},
          {"g7j", passBaseline1l && rec_njet_pt30>=7},
          {"g7j_g1b", passBaseline1l && rec_njet_pt30>=7 && pass_g1b},
          {"g7j_g1b_mbl", passBaseline1l && rec_njet_pt30>=7 && pass_g1b && pass_mbl},
          {"g7j_g1b_mbl_0t",  passBaseline1l && rec_njet_pt30>=7 && pass_g1b && pass_mbl && pass_0t},
          {"g7j_g1b_mbl_1t1", passBaseline1l && rec_njet_pt30>=7 && pass_g1b && pass_mbl && pass_1t1},
          {"g7j_g1b_mbl_1t2", passBaseline1l && rec_njet_pt30>=7 && pass_g1b && pass_mbl && pass_1t2},
          {"g7j_g1b_mbl_1t3", passBaseline1l && rec_njet_pt30>=7 && pass_g1b && pass_mbl && pass_1t3},
      };
      const std::map<std::string, bool> cut_map_1mu {
          {"6j", passBaseline1mu && rec_njet_pt30==6},
          {"6j_g1b", passBaseline1mu && rec_njet_pt30==6 && pass_g1b},
          {"6j_g1b_mbl", passBaseline1mu && rec_njet_pt30==6 && pass_g1b && pass_mbl},
          {"6j_g1b_mbl_0t",  passBaseline1mu && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_0t},
          {"6j_g1b_mbl_1t1", passBaseline1mu && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_1t1},
          {"6j_g1b_mbl_1t2", passBaseline1mu && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_1t2},
          {"6j_g1b_mbl_1t3", passBaseline1mu && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_1t3},
      };
      const std::map<std::string, bool> cut_map_1el {
          {"6j", passBaseline1el && rec_njet_pt30==6},
          {"6j_g1b", passBaseline1el && rec_njet_pt30==6 && pass_g1b},
          {"6j_g1b_mbl", passBaseline1el && rec_njet_pt30==6 && pass_g1b && pass_mbl},
          {"6j_g1b_mbl_0t",  passBaseline1el && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_0t},
          {"6j_g1b_mbl_1t1", passBaseline1el && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_1t1},
          {"6j_g1b_mbl_1t2", passBaseline1el && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_1t2},
          {"6j_g1b_mbl_1t3", passBaseline1el && rec_njet_pt30==6 && pass_g1b && pass_mbl && pass_1t3},
      };


      for(auto& kv : cut_map_1l)
      {
          if(kv.second)
          {
              my_histos["h_njets_1l_"+kv.first]->Fill(rec_njet_pt30, eventweight);
              my_histos["h_ntops_1l_"+kv.first]->Fill(tops.size(), eventweight);
              my_histos["h_nb_1l_"+kv.first]->Fill(rec_njet_pt30_btag, eventweight);
              my_histos["h_HT_1l_"+kv.first]->Fill(HT_trigger, eventweight);
              my_histos["h_bdt_1l_"+kv.first]->Fill(eventshape_bdt_val, eventweight);
              my_histos["h_mbl_1l_"+kv.first]->Fill(mbl, eventweight);
              my_2d_histos["h_njets_bdt_1l_"+kv.first]->Fill(rec_njet_pt30, eventshape_bdt_val, eventweight);
          }
      }
      for(auto& kv : cut_map_1mu)
      {
          if(kv.second)
          {
              my_histos["h_mupt_1mu_"+kv.first]->Fill(rec_muon_pt30[0].Pt(), eventweight);
              my_histos["h_njets_1mu_"+kv.first]->Fill(rec_njet_pt30, eventweight);
              my_histos["h_ntops_1mu_"+kv.first]->Fill(tops.size(), eventweight);
              my_histos["h_nb_1mu_"+kv.first]->Fill(rec_njet_pt30_btag, eventweight);
              my_histos["h_HT_1mu_"+kv.first]->Fill(HT_trigger, eventweight);
              my_histos["h_bdt_1mu_"+kv.first]->Fill(eventshape_bdt_val, eventweight);
              my_histos["h_mbl_1mu_"+kv.first]->Fill(mbl, eventweight);
              my_2d_histos["h_njets_bdt_1mu_"+kv.first]->Fill(rec_njet_pt30, eventshape_bdt_val, eventweight);
          }
      }
      for(auto& kv : cut_map_1el)
      {
          if(kv.second)
          {
              my_histos["h_elpt_1el_"+kv.first]->Fill(rec_electron_pt30[0].Pt(), eventweight);
              my_histos["h_njets_1el_"+kv.first]->Fill(rec_njet_pt30, eventweight);
              my_histos["h_ntops_1el_"+kv.first]->Fill(tops.size(), eventweight);
              my_histos["h_nb_1el_"+kv.first]->Fill(rec_njet_pt30_btag, eventweight);
              my_histos["h_HT_1el_"+kv.first]->Fill(HT_trigger, eventweight);
              my_histos["h_bdt_1el_"+kv.first]->Fill(eventshape_bdt_val, eventweight);
              my_histos["h_mbl_1el_"+kv.first]->Fill(mbl, eventweight);
              my_2d_histos["h_njets_bdt_1el_"+kv.first]->Fill(rec_njet_pt30, eventshape_bdt_val, eventweight);
          }
      }

      const std::map<std::string, bool> cut_map_2l {
          {"onZ", passBaseline2l && onZ},
          {"onZ_g1b", passBaseline2l && onZ && pass_g1b},
          {"onZ_g1b_nombl", passBaseline2l && onZ && pass_g1b && !passMbl_2l},
          {"onZ_g1b_nombl_bdt1", passBaseline2l && onZ && pass_g1b && !passMbl_2l && bdt_bin1},
          {"onZ_g1b_nombl_bdt2", passBaseline2l && onZ && pass_g1b && !passMbl_2l && bdt_bin2},
          {"onZ_g1b_nombl_bdt3", passBaseline2l && onZ && pass_g1b && !passMbl_2l && bdt_bin3},
          {"onZ_g1b_nombl_bdt4", passBaseline2l && onZ && pass_g1b && !passMbl_2l && bdt_bin4},
          {"onZ_g1b_g1t", passBaseline2l && onZ && pass_g1b && pass_1t}, 
          {"onZ_g1b_nombl_g1t", passBaseline2l && onZ && pass_g1b && !passMbl_2l && pass_1t}, 
          {"2b", passBaseline2l && rec_njet_pt30_btag == 2} 
      };

      for(auto& kv : cut_map_2l)
      {
          if(kv.second)
          {
              my_histos["h_njets_2l_"+kv.first]->Fill(rec_njet_pt30, eventweight);
              my_histos["h_ntops_2l_"+kv.first]->Fill(tops.size(), eventweight);
              my_histos["h_nb_2l_"+kv.first]->Fill(rec_njet_pt30_btag, eventweight);
              my_histos["h_HT_2l_"+kv.first]->Fill(HT_trigger, eventweight);
              my_histos["h_bdt_2l_"+kv.first]->Fill(eventshape_bdt_val, eventweight);
              my_2d_histos["h_njets_bdt_2l_"+kv.first]->Fill(rec_njet_pt30, eventshape_bdt_val, eventweight);
          }
      }

      // Fill event selection efficiencies
      my_efficiencies["event_sel_total"]->Fill(true,0);
      my_efficiencies["event_sel_total"]->Fill(HT_trigger>500,1);
      my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 ,2);
      my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>0 ,3);
      my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>0 && tops.size()>0 ,4);
      my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>0 && tops.size()>0 && rec_njet_pt45_btag>1 ,5);
      my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>0 && tops.size()>0 && rec_njet_pt45_btag>1 && tops.size()>1 ,6);
      my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>0 && tops.size()>0 && rec_njet_pt45_btag>1 && tops.size()>1 && rec_njet_pt30>=8 ,7);
      
      my_efficiencies["event_sel"]->Fill(true,0);
      my_efficiencies["event_sel"]->Fill(HT_trigger>500,1);
      if(HT_trigger>500)
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
                              my_efficiencies["event_sel"]->Fill(rec_njet_pt30>=8,7);
                          }
                      }
                  }
              }
          }
      }
 
      my_histos["h_met"]->Fill(MET, eventweight);
      my_histos["h_ht"]->Fill(HT, eventweight);


   } // end of event loop

}

bool ExploreEventSelection::PassTriggerGeneral(std::vector<std::string> &mytriggers)
{
    bool passTrigger = false;
    for(unsigned int i=0; i<TriggerNames->size(); ++i)
    {
        if(TriggerPass->at(i) != 1)
            continue;
        std::string trigname = TriggerNames->at(i);
        if( std::any_of(mytriggers.begin(), mytriggers.end(), [&] (std::string s) { return trigname.find(s)!=std::string::npos; }) )
        {
            passTrigger = true;
            break;
        }
    }
    return passTrigger;

}


bool ExploreEventSelection::PassTriggerAllHad()
{
    std::vector<std::string> mytriggers {
        //"HLT_PFHT1050", // 2017 trigger
        //"HLT_PFHT900"
            //"HLT_PFHT380_SixPFJet32_DoublePFBTagCSV", // 2017 trigger
            //"HLT_PFHT430_SixPFJet40_PFBTagCSV", // 2017 trigger
            "HLT_PFHT450_SixJet40_BTagCSV",
            "HLT_PFHT400_SixJet30_DoubleBTagCSV",            
            };
    return PassTriggerGeneral(mytriggers);
}

bool ExploreEventSelection::PassTriggerMuon()
{
    std::vector<std::string> mytriggers {"HLT_IsoMu24","HLT_IsoTkMu24_v"};
    return PassTriggerGeneral(mytriggers);
}

bool ExploreEventSelection::PassTriggerElectron()
{
    std::vector<std::string> mytriggers {"HLT_Ele27_WPTight_Gsf"};
    return PassTriggerGeneral(mytriggers);
}

void ExploreEventSelection::WriteHistos()
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
