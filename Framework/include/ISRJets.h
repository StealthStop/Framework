#ifndef ISRJets_h
#define ISRJets_h

#include "Framework/Framework/include/Utility.h"

#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopObject.h"
#include "TopTagger/TopTagger/interface/Constituent.h"

#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/interface/lester_mt2_bisect.h"

#include "TLorentzVector.h"
#include <iostream> 
#include <vector>
#include <cmath>

class ISRJets
{
private:
    std::string myVarSuffix_;

    // ----------------------------------------
    // function to find the particle parents
    // ----------------------------------------
    inline int findParent(const int p, const int idx, const std::vector<int>& GenParticles_ParentId, const std::vector<int>& GenParticles_ParentIdx) const
    {
        //std::cout << "GenParticles_ParentIdx size: " << GenParticles_ParentIdx.size() << " --- " << "Idx: " << idx << std::endl;
     
        if (idx == -1) 
        {
            return -1; 
        }   
        else if (GenParticles_ParentId[idx] == p)
        {
            return GenParticles_ParentId[idx];   
        }
        else
        {
            return findParent(p, GenParticles_ParentIdx[idx], GenParticles_ParentId, GenParticles_ParentIdx);
        }
    }   

    // -----------------------------------------------------------------------------------------------------
    // function to generate all possible matches between gen particles and gen jets if pass DR and pt cut
    // -----------------------------------------------------------------------------------------------------
    inline std::vector<std::tuple<int, int, double>> findAllDR(const std::vector<TLorentzVector>& GenParticles, 
                                                               const std::vector<TLorentzVector>& Jets,
                                                               const std::vector<bool>& GoodGenParticles,
                                                               const int resPartID,
                                                               const std::vector<int>& GenParticles_ParentId,
                                                               const std::vector<int>& GenParticles_ParentIdx,
                                                               const int& ISR_Idx, 
                                                               const double maxDR, const double maxPtRatio) const
    {
        bool check_ISR      = true;
        int check_resPartID = resPartID;

        std::vector<std::tuple<int, int, double>> AllDR;
        std::tuple<int, int, double> DRtup;
       
        for (unsigned int g = 0; g < GenParticles.size(); g++)
        {
            if (resPartID == 2212 )
            {
                check_ISR = GenParticles_ParentIdx.at(g) == ISR_Idx;
            }
            
            for (unsigned int j = 0; j < Jets.size(); j++)
            {
                bool passDR = GenParticles.at(g).DeltaR(Jets.at(j)) < maxDR;
                bool passPt = abs( 1 - Jets.at(j).Pt() / GenParticles.at(g).Pt() ) < maxPtRatio;
                if (findParent(check_resPartID, g, GenParticles_ParentId, GenParticles_ParentIdx) == check_resPartID && GoodGenParticles.at(g) && check_ISR && passDR && passPt)
                {
                    std::get<0>(DRtup) = g;
                    std::get<1>(DRtup) = j;
                    std::get<2>(DRtup) = GenParticles.at(g).DeltaR(Jets.at(j));
                    AllDR.push_back(DRtup);
                }
            }
        }
        return AllDR;      
    }

    // ------------------------------------------------------------
    // fuction to sort the best GenParticles and Jets matches
    // ------------------------------------------------------------
    void getMatches(const std::vector<std::tuple<int, int, double>>& AllDR, 
                    std::vector<std::pair<int, int>>& Matches, 
                    std::vector<bool> availableDR) const
    {
        std::tuple<int, int, double> bestDR;
        double minDR = 999;
        for ( unsigned int d = 0; d < AllDR.size(); d++ )
        {
            if ( std::get<2>(AllDR.at(d)) < minDR && availableDR.at(d) )
            {
                bestDR = AllDR.at(d);
                minDR  = std::get<2>(AllDR.at(d));
            }
        }

        bool allgone = true;
        for ( const auto& u : availableDR )
        {
            if (u) allgone = false;
        }

        if (!allgone)
        {
            Matches.push_back( std::make_pair( std::get<0>(bestDR), std::get<1>(bestDR) ) );
        
            for ( unsigned int d = 0; d < AllDR.size(); d++ )
            {
                if ( std::get<0>(AllDR.at(d)) == std::get<0>(bestDR) || std::get<1>(AllDR.at(d)) == std::get<1>(bestDR) )
                {
                    availableDR.at(d) = false;
                }
            } 
            getMatches(AllDR, Matches, availableDR);
        }
     }


    // ----------------------
    // Remove the ISR jets
    // ----------------------
    void removeISRJets(NTupleReader& tr) const
    {
        const auto& runtype = tr.getVar<std::string>("runtype");

        if(runtype != "Data")
        {
            const auto& Jets                   = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
            const auto& GoodJets_pt20          = tr.getVec<bool>("GoodJets_pt20"+myVarSuffix_);
            const auto& GenParticles           = tr.getVec<TLorentzVector>("GenParticles"+myVarSuffix_);
            const auto& GenParticles_PdgId     = tr.getVec<int>("GenParticles_PdgId"+myVarSuffix_);
            const auto& GenParticles_ParentId  = tr.getVec<int>("GenParticles_ParentId"+myVarSuffix_);
            const auto& GenParticles_ParentIdx = tr.getVec<int>("GenParticles_ParentIdx"+myVarSuffix_);
            const auto& GenParticles_Status    = tr.getVec<int>("GenParticles_Status"+myVarSuffix_); 
        
            // ---------------------------
            // loop over the RecoParticles
            // ---------------------------
            std::vector<TLorentzVector> RecoISR;
            int nRecoISR = 0;
            for (unsigned int j = 0; j < Jets.size(); j++)
            {
                //if (!GoodJets_pt20[j]) continue;
                RecoISR.push_back(Jets.at(j));
                nRecoISR++;
            }
 
            // -----------------------------------------
            // loop over the GenParticles to get Gen ISR 
            // ----------------------------------------- 
            std::vector<bool> GenISR(GenParticles.size(), false);
            int nGenISR = 0;
            for (unsigned int g = 0; g < GenParticles.size(); g++)
            {
                int pdgId  = GenParticles_PdgId.at(g);
                int momId  = GenParticles_ParentId.at(g);
                int status = GenParticles_Status.at(g);
                
                //std::cout << "pdgId: " << pdgId  << " --- " << "momId: " << momId << " --- " << "status code: " << status << std::endl;                           

                bool is_jet   = (abs(pdgId) <= 6 || abs(pdgId) == 21);
                bool pass_isr = is_jet ? status == 23 && (abs(momId) == 2212) : false;  
                //bool filter   = GenParticles.at(g).Pt() > 20 && abs(GenParticles.at(g).Eta()) < 2.4;           

                if (pass_isr)
                {
                    GenISR.at(g) = true;
                    nGenISR++;
                }
            }  

            // ----------------
            // ISR Gen matching
            // ----------------
            auto& GM_ISRmatching_DR = tr.createDerivedVec<double>("GM_ISRmatching_DR"+myVarSuffix_);              
            auto& GM_ISRmatching_Pt = tr.createDerivedVec<double>("GM_ISRmatching_Pt"+myVarSuffix_);
            auto& ISRmatched        = tr.createDerivedVec<bool>("ISRmatched"+myVarSuffix_, RecoISR.size());

            std::vector<int> ListParticles{1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6, 2212};
            std::vector<TLorentzVector> isrGenList;
            std::vector<TLorentzVector> isrRecoList;
            int ISR_Idx  = -1;
            int nISRJets = 0;

            // matching
            for (unsigned int m = 0; m < ListParticles.size(); m++)            
            {
                std::vector<std::pair<int, int>> Matches;
                double maxDR      = 0.1; // set max DR allowed for matching
                double maxPtRatio = 0.5; // set max pT allowed for matching
                std::vector<std::tuple<int, int, double>> AllDR = findAllDR(GenParticles, RecoISR, GenISR, ListParticles[m], GenParticles_ParentId, GenParticles_ParentIdx, ISR_Idx, maxDR, maxPtRatio);
                std::vector<bool> availableDR(AllDR.size(), true);

                getMatches(AllDR, Matches, availableDR);   
                
                std::vector<double> DRvec;
                std::vector<double> PtVec;
                TLorentzVector isrGenMatched;
                TLorentzVector isrRecoMatched;
                for (unsigned int match = 0; match < Matches.size(); match++)
                {
                    isrGenMatched  = GenParticles.at(Matches.at(match).first);
                    isrRecoMatched = RecoISR.at(Matches.at(match).second);
                    
                    GM_ISRmatching_DR.push_back( GenParticles.at(Matches.at(match).first).DeltaR(RecoISR.at(Matches.at(match).second)) );
                    GM_ISRmatching_Pt.push_back( abs( 1 - GenParticles.at(Matches.at(match).first).Pt() / RecoISR.at(Matches.at(match).second).Pt() ) );                    

                    // getting ISRmatched jets
                    ISRmatched.at(Matches.at(match).second) = true;
                    nISRJets++;
                    
                }
                //isrGenList.push_back(isrGenMatched);
                //isrRecoList.push_back(isrRecoMatched);
            }
          
             
            tr.registerDerivedVar("nGenISR"+myVarSuffix_, nGenISR); 
            //tr.registerDerivedVar("GM_genISR"+myVarSuffix_, isrGenList.at(0));
            tr.registerDerivedVar("nRecoISR"+myVarSuffix_, nRecoISR);
            //tr.registerDerivedVar("GM_recoISR"+myVarSuffix_, isrRecoList.at(0));
            tr.registerDerivedVar("nISRJets"+myVarSuffix_, nISRJets);       
            //tr.registerDerivedVec("ISRmatched"+myVarSuffix_, ISRmatched);
 
        }
    }

public:    
    ISRJets(const std::string& myVarSuffix = "")
        :myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up ISRJets"<<std::endl;;
    }

    void operator()(NTupleReader& tr)
    {
        removeISRJets(tr);
    }
};

#endif
