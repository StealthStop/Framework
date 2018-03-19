#define NtupleClass_cxx
#include "NtupleClass.h"

#include "Utility.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <iostream>

//mandatory includes to use top tagger
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/include/TTException.h"
#include "SetUpTopTagger.h"

void NtupleClass::Loop(std::string runtype)
{
//   In a ROOT session, you can do:
//      root> .L NtupleClass.C
//      root> NtupleClass t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   // -----------------------
   // make some histograms
   // -----------------------
   TH1D *myHisto  = new TH1D("njets","njets", 20, 0, 20);
   TH1D *h_met  = new TH1D("h_met","h_met", 20, 0, 200);
   TH1D *h_ht  = new TH1D("h_ht","h_ht", 60, 0, 3000);
   TH1D *h_ntops  = new TH1D("h_ntops","h_ntops", 5, 0, 5);

   TopTagger tt;
   tt.setCfgFile("TopTagger.cfg");

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
   
      nb = fChain->GetEntry(jentry);   
      nbytes += nb;

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
      h_ntops->Fill(tops.size());

      // get set of all constituents (i.e. AK4 and AK8 jets) used in one of the tops
      std::set<Constituent const *> usedConstituents = ttr.getUsedConstituents();

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

      bool passBaseline = HT_pt40>500 && rec_njet_pt45>=6 && rec_njet_pt45_btag>1 && tops.size()>1;
 
      if (passBaseline)
      {
          myHisto->Fill(NJets);
          h_met->Fill(MET);
          h_ht->Fill(HT);
      }
   }

   myHisto->Write();
   h_met->Write();
   h_ht->Write();
   h_ntops->Write();
}
