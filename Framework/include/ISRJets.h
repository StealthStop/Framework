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
    inline void findBestDR(const std::vector<TLorentzVector>& GenParticles, 
                           const std::vector<TLorentzVector>& Jets,
                           const std::vector<bool>& GoodGenParticles,
                           const int resPartID,
                           const std::vector<int>& GenParticles_ParentId,
                           const std::vector<int>& GenParticles_ParentIdx,
                           const int& ISR_Idx, 
                           const double maxDR, const double maxPtRatio,
                           std::vector<std::tuple<int, int, double>>& allDR,
                           std::vector<std::tuple<int, int, double>>& bestAllDR,
                           std::vector<std::tuple<int, int, double>>& justDR,
                           std::vector<std::tuple<int, int, double>>& justPtRatio) const
    {
        bool check_ISR      = true;
        int check_resPartID = resPartID;
 
        std::tuple<int, int, double> DRtup;
       
        for (unsigned int g = 0; g < GenParticles.size(); g++)
        {
            if (resPartID <= 6 || resPartID == 21) 
            {
                check_ISR = GenParticles_ParentIdx.at(g) == ISR_Idx;
            }
            
            for (unsigned int j = 0; j < Jets.size(); j++)
            {
                bool passBestDR  = GenParticles.at(g).DeltaR(Jets.at(j)) < maxDR;
                bool passBestPt  = ( Jets.at(j).Pt() / GenParticles.at(g).Pt() ) > ( 1 - maxPtRatio ) && ( Jets.at(j).Pt() / GenParticles.at(g).Pt() ) < ( 1 + maxPtRatio );
                bool findParents = (findParent(check_resPartID, g, GenParticles_ParentId, GenParticles_ParentIdx) == check_resPartID && GoodGenParticles.at(g) && check_ISR);

                if (findParents)
                {
                    // all possible matching
                    std::get<0>(DRtup) = g;
                    std::get<1>(DRtup) = j;
                    std::get<2>(DRtup) = GenParticles.at(g).DeltaR(Jets.at(j));
                    allDR.push_back(DRtup);

                    // best matching
                    if (passBestDR && passBestPt)
                    {
                        std::get<0>(DRtup) = g;
                        std::get<1>(DRtup) = j;
                        std::get<2>(DRtup) = GenParticles.at(g).DeltaR(Jets.at(j));
                        bestAllDR.push_back(DRtup);
                    }
                    
                    // matching with cutting on deltaR
                    if (passBestDR) 
                    {
                        std::get<0>(DRtup) = g;
                        std::get<1>(DRtup) = j;
                        std::get<2>(DRtup) = GenParticles.at(g).DeltaR(Jets.at(j));
                        justDR.push_back(DRtup);                  
                    }

                    // matching with cutting on pt ratio
                    if (passBestPt)
                    {
                        std::get<0>(DRtup) = g;
                        std::get<1>(DRtup) = j;
                        std::get<2>(DRtup) = GenParticles.at(g).DeltaR(Jets.at(j));
                        justPtRatio.push_back(DRtup);
                    }    
                }
            }
        }
    }

    // ----------------------------------------------------------
    // function to sort the best GenParticles and Jets matches
    // ----------------------------------------------------------
    void getMatches(const std::vector<std::tuple<int, int, double>>& bestAllDR, 
                    std::vector<std::pair<int, int>>& bestMathches, 
                    std::vector<bool> availableDR) const
    {
        std::tuple<int, int, double> bestDR;
        double minDR = 999;
        for ( unsigned int d = 0; d < bestAllDR.size(); d++ )
        {
            if ( std::get<2>(bestAllDR.at(d)) < minDR && availableDR.at(d) )
            {
                bestDR = bestAllDR.at(d);
                minDR  = std::get<2>(bestAllDR.at(d));
            }
        }

        bool allgone = true;
        for ( const auto& u : availableDR )
        {
            if (u) allgone = false;
        }

        if (!allgone)
        {
            bestMathches.push_back( std::make_pair( std::get<0>(bestDR), std::get<1>(bestDR) ) );
        
            for ( unsigned int d = 0; d < bestAllDR.size(); d++ )
            {
                if ( std::get<0>(bestAllDR.at(d)) == std::get<0>(bestDR) || std::get<1>(bestAllDR.at(d)) == std::get<1>(bestDR) )
                {
                    availableDR.at(d) = false;
                }
            } 
            getMatches(bestAllDR, bestMathches, availableDR);
        }
     }

    // -------------------------
    // -- Remove the ISR jets
    // -------------------------
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
        
            // ---------------------
            // loop over the Jets
            // ---------------------
            std::vector<TLorentzVector> RecoISR;
            int nRecoISR = 0;
            for (unsigned int j = 0; j < Jets.size(); j++)
            {
                if (!GoodJets_pt20[j]) continue;
                RecoISR.push_back(Jets.at(j));
                nRecoISR++;
            }
 
            // --------------------------------------------
            // loop over the GenParticles to get Gen ISR 
            // -------------------------------------------- 
            auto& GenISR  = tr.createDerivedVec<bool>("GenISR"+myVarSuffix_, GenParticles.size(), false);
            int   nGenISR = 0;
            
            for (unsigned int g = 0; g < GenParticles.size(); g++)
            {
                int pdgId    = GenParticles_PdgId.at(g);
                int momPdgId = GenParticles_ParentId.at(g);
                int status   = GenParticles_Status.at(g);

                //if ( ( abs(pdgId) <= 6 || abs(pdgId) == 21 )
                //       && ( abs(momPdgId) <= 6 || abs(momPdgId) == 21 )
                //       && ( status == 23 ))
                //{
                //    std::cout << "momPdgId: " <<  momPdgId << " || " << "pdgId : " << pdgId  << " || " << "status: " << status << std::endl; 
                //}

                // mom-dau: qg, gq, gg
                bool pass_ISR = ( abs(pdgId) <= 6  && abs(momPdgId) == 21 && status == 23 ) ||
                                ( abs(pdgId) == 21 && abs(momPdgId) <= 6  && status == 23 ) || 
                                ( abs(pdgId) == 21 && abs(momPdgId) == 21 && status == 23 ) ; 
                bool filter   = GenParticles.at(g).Pt() > 20 && abs(GenParticles.at(g).Eta()) < 2.4;

                if (pass_ISR && filter)
                {
                    GenISR.at(g) = true;
                    nGenISR++;
                }
            }  

            // ---------------------------------------
            // ISR filter by truth definition of it  
            // ---------------------------------------
            auto& GM_ISRmatching_allDR                    = tr.createDerivedVec<double>("GM_ISRmatching_allDR"+myVarSuffix_);
            auto& GM_ISRmatching_bestDR                   = tr.createDerivedVec<double>("GM_ISRmatching_bestDR"+myVarSuffix_);
            auto& GM_ISRmatching_justCutOnDR_DR           = tr.createDerivedVec<double>("GM_ISRmatching_justCutOnDR_DR"+myVarSuffix_);
            auto& GM_ISRmatching_justCutOnPtRatio_DR      = tr.createDerivedVec<double>("GM_ISRmatching_justCutOnPtRatio_DR"+myVarSuffix_);
            auto& GM_ISRmatching_allPtRatio               = tr.createDerivedVec<double>("GM_ISRmatching_allPtRatio"+myVarSuffix_);
            auto& GM_ISRmatching_bestPtRatio              = tr.createDerivedVec<double>("GM_ISRmatching_bestPtRatio"+myVarSuffix_);
            auto& GM_ISRmatching_justCutOnDR_PtRatio      = tr.createDerivedVec<double>("GM_ISRmatching_justCutOnDR_PtRatio"+myVarSuffix_);
            auto& GM_ISRmatching_justCutOnPtRatio_PtRatio = tr.createDerivedVec<double>("GM_ISRmatching_justCutOnPtRatio_PtRatio"+myVarSuffix_);
            auto& ISRmatched_dr_ptr                       = tr.createDerivedVec<bool>("ISRmatched_dr_ptr"+myVarSuffix_, RecoISR.size(), false);
            auto& ISRmatched_dr                           = tr.createDerivedVec<bool>("ISRmatched_dr"+myVarSuffix_, RecoISR.size(), false);

            std::vector<int> ListParticles{1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6, 21};
            int ISR_Idx         = -1;
            int nISRJets_dr_ptr = 0;
            int nISRJets_dr     = 0;

            for (unsigned int m = 0; m < ListParticles.size(); m++)            
            {
                std::vector<std::tuple<int, int, double>> allMatches;               
                std::vector<std::tuple<int, int, double>> bestMatches;
                std::vector<std::tuple<int, int, double>> justDRMatches;
                std::vector<std::tuple<int, int, double>> justPtRatioMatches;
                double maxDR      = 0.3; // set max DR allowed for matching / = 0.1
                double maxPtRatio = 0.5; // set max pT allowed for matching
                
                findBestDR(GenParticles, RecoISR, GenISR, ListParticles[m], GenParticles_ParentId, GenParticles_ParentIdx, ISR_Idx, maxDR, maxPtRatio, allMatches, bestMatches, justDRMatches, justPtRatioMatches);

                // --------------------
                // all possible matches
                // --------------------
                for (unsigned int match = 0; match < allMatches.size(); match++)
                {
                    GM_ISRmatching_allDR.push_back( GenParticles.at(std::get<0>(allMatches.at(match))).DeltaR(RecoISR.at(std::get<1>(allMatches.at(match)))) );
                    GM_ISRmatching_allPtRatio.push_back( ( RecoISR.at(std::get<1>(allMatches.at(match))).Pt() / GenParticles.at(std::get<0>(allMatches.at(match))).Pt() ) );
                }       

                // -------------------------------------------
                // matches with cutting on deltaR and pt ratio
                // -------------------------------------------
                for (unsigned int match = 0; match < bestMatches.size(); match++)
                {
                    GM_ISRmatching_bestDR.push_back( GenParticles.at(std::get<0>(bestMatches.at(match))).DeltaR(RecoISR.at(std::get<1>(bestMatches.at(match)))) );
                    GM_ISRmatching_bestPtRatio.push_back( ( RecoISR.at(std::get<1>(bestMatches.at(match))).Pt() / GenParticles.at(std::get<0>(bestMatches.at(match))).Pt() ) );

                    // getting ISRmatched_dr_ptr jets
                    ISRmatched_dr_ptr.at(std::get<1>(bestMatches.at(match))) = true;
                    nISRJets_dr_ptr++;
                }

                // ------------------------------
                // matches with cutting on deltaR
                // ------------------------------
                for (unsigned int match = 0; match < justDRMatches.size(); match++)
                {
                    GM_ISRmatching_justCutOnDR_DR.push_back( GenParticles.at(std::get<0>(justDRMatches.at(match))).DeltaR(RecoISR.at(std::get<1>(justDRMatches.at(match)))) );
                    GM_ISRmatching_justCutOnDR_PtRatio.push_back( ( RecoISR.at(std::get<1>(justDRMatches.at(match))).Pt() / GenParticles.at(std::get<0>(justDRMatches.at(match))).Pt() ) );                

                    // getting ISRmatched_dr jets
                    ISRmatched_dr.at(std::get<1>(justDRMatches.at(match))) = true;
                    nISRJets_dr++;
                }   

                // --------------------------------
                // matches with cutting on pt ratio
                // --------------------------------
                for (unsigned int match = 0; match < justPtRatioMatches.size(); match++)
                {
                    GM_ISRmatching_justCutOnPtRatio_DR.push_back( GenParticles.at(std::get<0>(justPtRatioMatches.at(match))).DeltaR(RecoISR.at(std::get<1>(justPtRatioMatches.at(match)))) );
                    GM_ISRmatching_justCutOnPtRatio_PtRatio.push_back( ( RecoISR.at(std::get<1>(justPtRatioMatches.at(match))).Pt() / GenParticles.at(std::get<0>(justPtRatioMatches.at(match))).Pt() ) );
                }       

                tr.registerDerivedVar("nGenISR"+myVarSuffix_, nGenISR);
                tr.registerDerivedVar("nRecoISR"+myVarSuffix_, nRecoISR);
                tr.registerDerivedVar("nISRJets_dr_ptr"+myVarSuffix_, nISRJets_dr_ptr);      
                tr.registerDerivedVar("nISRJets_dr"+myVarSuffix_, nISRJets_dr); 
            
            } // ISR truth level matching loop 

            // ------------------------------------------
            // NonISR filter by using TreeMaker method
            // ------------------------------------------
            auto& NonISRmatched = tr.createDerivedVec<bool>("NonISRmatched"+myVarSuffix_, Jets.size(), false);
            int nISR = 0;

            for (unsigned int j = 0; j < Jets.size(); j++)
            {
                if (!GoodJets_pt20[j]) continue;
        
                bool matched = false;

                // gen matching
                for (unsigned int g = 0; g < GenParticles.size(); g++)
                {
                    if (matched) break;

                    int pdgId    = GenParticles_PdgId.at(g);
                    int momPdgId = GenParticles_ParentId.at(g);
                    int status   = GenParticles_Status.at(g);
 
                    if (abs(pdgId) > 5 || abs(pdgId) != 21 || status != 23) continue;    
                    if(!(momPdgId == 6 || momPdgId == 23   || momPdgId == 24 || momPdgId == 25 || momPdgId > 1e6)) continue;            
                    
                    // check against daughter in case of hard initial splitting
                    double dR = GenParticles.at(g).DeltaR(Jets.at(j)); 
                    if (dR < 0.3)
                    {
                        matched = true;
                        break;
                    }        
                } 
                  
                // NonISR Jets filter      
                if(!matched)
                {
                    NonISRmatched.at(j) = true;
                    nISR++;
                } 
                
                tr.registerDerivedVar("nISR"+myVarSuffix_, nISR);

            } // TreeMaker matching loop       
                

        } // data/MC loop
    } // event loop

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
