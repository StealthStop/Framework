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

    // ----------------
    // -- Event Loop
    // ----------------
    void removeISRJets(NTupleReader& tr) const
    {
        const auto& runtype = tr.getVar<std::string>("runtype");

        if(runtype != "Data")
        {
            const auto& Jets                  = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
            const auto& GoodJets_pt20         = tr.getVec<bool>("GoodJets_pt20"+myVarSuffix_);
            const auto& GenParticles          = tr.getVec<TLorentzVector>("GenParticles"+myVarSuffix_);
            const auto& GenParticles_PdgId    = tr.getVec<int>("GenParticles_PdgId"+myVarSuffix_);
            const auto& GenParticles_ParentId = tr.getVec<int>("GenParticles_ParentId"+myVarSuffix_);
            const auto& GenParticles_Status   = tr.getVec<int>("GenParticles_Status"+myVarSuffix_); 
       
            // ---------------------------------------
            // ISR filter by truth definition of it  
            //     -- ISR gen matching
            // ---------------------------------------
            auto& GM_ISRmatching_allDR                    = tr.createDerivedVec<double>("GM_ISRmatching_allDR"+myVarSuffix_);
            auto& GM_ISRmatching_bestDR                   = tr.createDerivedVec<double>("GM_ISRmatching_bestDR"+myVarSuffix_);
            auto& GM_ISRmatching_justCutOnDR_DR           = tr.createDerivedVec<double>("GM_ISRmatching_justCutOnDR_DR"+myVarSuffix_);
            auto& GM_ISRmatching_justCutOnPtRatio_DR      = tr.createDerivedVec<double>("GM_ISRmatching_justCutOnPtRatio_DR"+myVarSuffix_);
            auto& GM_ISRmatching_allPtRatio               = tr.createDerivedVec<double>("GM_ISRmatching_allPtRatio"+myVarSuffix_);
            auto& GM_ISRmatching_bestPtRatio              = tr.createDerivedVec<double>("GM_ISRmatching_bestPtRatio"+myVarSuffix_);
            auto& GM_ISRmatching_justCutOnDR_PtRatio      = tr.createDerivedVec<double>("GM_ISRmatching_justCutOnDR_PtRatio"+myVarSuffix_);
            auto& GM_ISRmatching_justCutOnPtRatio_PtRatio = tr.createDerivedVec<double>("GM_ISRmatching_justCutOnPtRatio_PtRatio"+myVarSuffix_);
            auto& ISRmatched_dr_ptr                       = tr.createDerivedVec<bool>("ISRmatched_dr_ptr"+myVarSuffix_, Jets.size(), false);
            auto& ISRmatched_dr                           = tr.createDerivedVec<bool>("ISRmatched_dr"+myVarSuffix_, Jets.size(), false);

            // ---------------------------
            // Loop over the gen particles
            // --------------------------- 
            for (unsigned int g = 0; g < GenParticles.size(); g++)
            {
                int pdgId    = GenParticles_PdgId.at(g);
                int momPdgId = GenParticles_ParentId.at(g);
                int status   = GenParticles_Status.at(g);

                // ISR definition - mom-dau: qg, gq, gg
                bool pass_ISR = ( abs(pdgId) <= 6  && abs(momPdgId) == 21 && status == 23 ) ||
                                ( abs(pdgId) == 21 && abs(momPdgId) <= 6  && status == 23 ) || 
                                ( abs(pdgId) == 21 && abs(momPdgId) == 21 && status == 23 ) ; 
       
                // ----------------------------
                // Loop over the reco particles
                // ----------------------------
                for (unsigned int j = 0; j < Jets.size(); j++)
                {
                    double maxDR      = 0.3; // set max DR allowed for matching / = 0.1
                    double maxPtRatio = 0.5; // set max pT allowed for matching

                    bool passBestDR = GenParticles.at(g).DeltaR(Jets.at(j)) < maxDR;
                    bool passBestPt = ( Jets.at(j).Pt() / GenParticles.at(g).Pt() ) > ( 1 - maxPtRatio ) && 
                                      ( Jets.at(j).Pt() / GenParticles.at(g).Pt() ) < ( 1 + maxPtRatio );

                    // ----------------------------
                    // matching without any cutting
                    // ----------------------------           
                    if (pass_ISR)
                    {
                        GM_ISRmatching_allDR.push_back( GenParticles.at(g).DeltaR(Jets.at(j)) );                         
                        GM_ISRmatching_allPtRatio.push_back( Jets.at(j).Pt() / GenParticles.at(g).Pt() );
                    }

                    // --------------------------------------------
                    // matching with cutting on deltaR and pt ratio
                    // --------------------------------------------           
                    if (pass_ISR && passBestDR && passBestPt)
                    {
                        GM_ISRmatching_bestDR.push_back( GenParticles.at(g).DeltaR(Jets.at(j)) );
                        GM_ISRmatching_bestPtRatio.push_back( ( Jets.at(j).Pt() / GenParticles.at(g).Pt() ) );

                        // ISR jet filter
                        ISRmatched_dr_ptr.at(j) = true;
                    }

                    // -------------------------------
                    // matching with cutting on deltaR
                    // -------------------------------               
                    if (pass_ISR && passBestDR)
                    {
                        GM_ISRmatching_justCutOnDR_DR.push_back( GenParticles.at(g).DeltaR(Jets.at(j)) );
                        GM_ISRmatching_justCutOnDR_PtRatio.push_back( ( Jets.at(j).Pt() / GenParticles.at(g).Pt() ) );

                        // ISR jet filter
                        ISRmatched_dr.at(j) = true;
                    } 

                    // ---------------------------------
                    // matching with cutting on pt ratio
                    // ---------------------------------
                    if (pass_ISR && passBestPt)
                    {
                        GM_ISRmatching_justCutOnPtRatio_DR.push_back( GenParticles.at(g).DeltaR(Jets.at(j)) );
                        GM_ISRmatching_justCutOnPtRatio_PtRatio.push_back( ( Jets.at(j).Pt() / GenParticles.at(g).Pt() ) );                        
                    }     
 
                } // reco particles
            } // Gen particles 

            // -------------------------------------------------------
            // NonISR & Other Jets filter by using TreeMaker method
            //     -- NonISR    : signal jets
            //     -- OtherJets : ISR + FSR + underlying
            // -------------------------------------------------------
            auto& NonISRmatched_dr     = tr.createDerivedVec<bool>("NonISRmatched_dr"+myVarSuffix_, Jets.size(), false);
            auto& NonISRmatched_dr_ptr = tr.createDerivedVec<bool>("NonISRmatched_dr_ptr"+myVarSuffix_, Jets.size(), false);
            auto& OtherJets            = tr.createDerivedVec<bool>("OtherJets"+myVarSuffix_, Jets.size(), false);
            int nNonISR_dr     = 0;
            int nNonISR_dr_ptr = 0;
            int nOtherJets     = 0;

            // ----------------------------
            // Loop over the reco particles
            // ----------------------------
            for (unsigned int j = 0; j < Jets.size(); j++)
            {
                if (!GoodJets_pt20[j]) continue;
        
                bool matched = false;

                // ----------------------------
                // Loop over the gen pareticles
                // ----------------------------
                for (unsigned int g = 0; g < GenParticles.size(); g++)
                {
                    if (matched) break;

                    int pdgId    = GenParticles_PdgId.at(g);
                    int momPdgId = GenParticles_ParentId.at(g);
                    int status   = GenParticles_Status.at(g);
 
                    if (abs(pdgId) > 5      || status != 23) continue;  
                    if(!(abs(momPdgId) == 6 || abs(momPdgId) == 24 || abs(momPdgId) == 1000022 || abs(momPdgId) == 1000006)) continue;                     

                    double dR    = GenParticles.at(g).DeltaR(Jets.at(j));
                    bool ptRatio = ( Jets.at(j).Pt() / GenParticles.at(g).Pt() ) > ( 0.5 ) &&
                                   ( Jets.at(j).Pt() / GenParticles.at(g).Pt() ) < ( 1.5 ); 
                    
                    // -------------------------------
                    // matching with cutting on deltaR
                    // -------------------------------
                    if (dR < 0.3)
                    {
                        matched = true;
                        NonISRmatched_dr.at(j) = true;
                        nNonISR_dr++;
                    }

                    // --------------------------------------------
                    // matching with cutting on deltaR and pt ratio            
                    // --------------------------------------------
                    if (dR < 0.3 && ptRatio)
                    {
                        matched = true;
                        NonISRmatched_dr_ptr.at(j) = true;
                        nNonISR_dr_ptr++;
                    }
                } 
                  
                // -------------------------------------
                // OtherJets (ISR+FSR+underlying) filter
                // -------------------------------------      
                if(!matched)
                {
                    OtherJets.at(j) = true;
                    nOtherJets++;
                } 
                
                tr.registerDerivedVar("nNonISR_dr"+myVarSuffix_, nNonISR_dr);
                tr.registerDerivedVar("nNonISR_dr_ptr"+myVarSuffix_, nNonISR_dr_ptr);
                tr.registerDerivedVar("nOtherJets"+myVarSuffix_, nOtherJets);

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