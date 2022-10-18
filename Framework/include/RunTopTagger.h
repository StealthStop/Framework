#ifndef RUNTOPTAGGER_H
#define RUNTOPTAGGER_H

#include "Framework/Framework/include/SetUpTopTagger.h"

class RunTopTagger
{
private:
    std::string taggerCfg_;
    std::string myVarSuffix_;

    std::unique_ptr<TopTagger> tt_;
    std::vector<TLorentzVector>* hadtops_;
    std::vector<std::vector<const TLorentzVector*>>* hadtopdaughters_;
    std::vector<std::vector<int>>* hadtopdaughters_id_;
    std::vector<TLorentzVector>* neutralinos_;
    std::vector<TLorentzVector>* singlinos_;
    std::vector<TLorentzVector>* singlets_;
    std::vector<int>* hadtops_idx_;
    int ntops_;
    int ntops_3jet_;
    int ntops_2jet_;
    int ntops_1jet_;
    
    inline int findParent(const int p, const int idx, const std::vector<int>& GenParticles_ParentId, const std::vector<int>& GenParticles_ParentIdx) const
    {
        if (idx == -1)
        {
            return -1;
        }
        else if(abs(GenParticles_ParentId[idx]) == p)
        {
            return GenParticles_ParentIdx[idx];
        }
        else
        {
            return findParent(p, GenParticles_ParentIdx[idx], GenParticles_ParentId, GenParticles_ParentIdx);
        }
    }

    void genMatch(NTupleReader& tr)
    {
        // ----------------------------------------------
        // check for number of hadronic tops at gen level
        // ----------------------------------------------
        const auto& runtype = tr.getVar<std::string>("runtype");

        if(runtype != "Data")
        {
            hadtops_            = &tr.createDerivedVec<TLorentzVector>("hadtops"+myVarSuffix_);
            hadtopdaughters_    = &tr.createDerivedVec<std::vector<const TLorentzVector*>>("hadtopdaughters"+myVarSuffix_);
            hadtopdaughters_id_ = &tr.createDerivedVec<std::vector<int>>("hadtopdaughters_id"+myVarSuffix_);
            neutralinos_        = &tr.createDerivedVec<TLorentzVector>("neutralinos"+myVarSuffix_);
            singlinos_          = &tr.createDerivedVec<TLorentzVector>("singlinos"+myVarSuffix_);
            singlets_           = &tr.createDerivedVec<TLorentzVector>("singlets"+myVarSuffix_);
            hadtops_idx_        = &tr.createDerivedVec<int>("hadtops_idx"+myVarSuffix_);
        
            const auto& GenParticles            = utility::convertVectorOfLV<TLorentzVector, utility::LorentzVector>(tr.getVec<utility::LorentzVector>("GenParticles"));
            const auto& GenParticles_PdgId      = tr.getVec<int>("GenParticles_PdgId");
            const auto& GenParticles_ParentId   = tr.getVec<int>("GenParticles_ParentId");
            const auto& GenParticles_ParentIdx  = tr.getVec<int>("GenParticles_ParentIdx");
            const auto& GenParticles_Status     = tr.getVec<int>("GenParticles_Status");            

            for(unsigned int gpi=0; gpi < GenParticles.size(); gpi++ ) 
            {
                int pdgid = GenParticles_PdgId.at(gpi);
                int status = GenParticles_Status.at(gpi);
                int momid = GenParticles_ParentId.at(gpi) ;
                int momidx = GenParticles_ParentIdx.at(gpi);
                int momstatus = (momidx == -1) ? -1 : GenParticles_Status.at(momidx);
                int topIdx = findParent(6, gpi, GenParticles_ParentId, GenParticles_ParentIdx);
                bool passWMomStatus = false;

                if((abs(momid) == 24) && (momstatus != 1 || status == 2 || momstatus != 22 || momstatus != 23 || momstatus != 52) )
                {
                    passWMomStatus = true;
                }
                if(abs(pdgid) == 1000022 && (status==22 || status == 52))
                {
                    neutralinos_->push_back(GenParticles.at(gpi));
                }
                if(abs(pdgid) == 5000001 && (status == 22 || status == 52))
                {
                    singlinos_->push_back(GenParticles.at(gpi));
                }
                if(abs(pdgid) == 5000002 && (status == 22 || status == 52))
                {
                    singlets_->push_back(GenParticles.at(gpi));
                }
                if( topIdx >= 0 && (abs(pdgid) != 24) && (passWMomStatus || abs(pdgid) == 5))
                {
                    
                    unsigned int position = 0;
                    for(;position < hadtops_idx_->size() && (*hadtops_idx_)[position] != topIdx; ++position);
                    if( position < hadtops_idx_->size() )
                    {
                        (*hadtopdaughters_)[position].push_back(&(GenParticles.at(gpi)));
                        (*hadtopdaughters_id_)[position].push_back(gpi);
                    } 
                    else
                    {
                        hadtops_idx_->push_back(topIdx);
                        hadtops_->push_back(GenParticles.at(topIdx));
                        hadtopdaughters_->push_back( {&(GenParticles.at(gpi))} );
                        hadtopdaughters_id_->push_back( {static_cast<int>(gpi)} );
                    }
                }
            }            
        }
    }

    void countTops(const std::vector<TopObject*>& tops)
    {
        ntops_ = 0;
        ntops_3jet_ = 0;
        ntops_2jet_ = 0;
        ntops_1jet_ = 0;
       
        for (const TopObject* top : tops)
        {
            ntops_++;

            if(top->getNConstituents() == 3 )
            {
                ntops_3jet_++;
            }
            else if(top->getNConstituents() == 2 )
            {
                ntops_2jet_++;
            }
            else if(top->getNConstituents() == 1 )
            {
                ntops_1jet_++;
            }
        }
    }

    void runTopTagger(NTupleReader& tr)
    {
        // Prepare class to run Top Tagger
        genMatch(tr);

        // Setup variables needed for top tagger
        SetUpTopTagger st( tr, *hadtops_, *hadtopdaughters_, myVarSuffix_ );
        const std::vector<Constituent>& constituents = st.getConstituents();

        // Run the top tagger             
        tt_->runTagger(constituents);

        // Get the top tagger results object     
        const TopTaggerResults& ttr = tt_->getResults();
       
        // Get tagged objects the new top tagger returns more than just tops now (MERGED_TOP, SEMIMERGEDWB_TOP, RESOLVED_TOP, MERGED_W, SEMIMERGEDQB_TOP)
        // For now we will only use merged and resolved tops
        const std::vector<TopObject*>& taggedObjects = ttr.getTops();
        std::vector<TopObject*> mergedTops;
        std::vector<TopObject*> resolvedTops;
        for(auto* o : taggedObjects)
        {
            if     (o->getType()==TopObject::MERGED_TOP)   mergedTops.push_back(o);
            else if(o->getType()==TopObject::RESOLVED_TOP) resolvedTops.push_back(o);
        }

        // --------------------------------------------------
        // Get reconstructed tops and derive needed variables
        // --------------------------------------------------                            
        // Regardless of how we are running the top tagger, we will always include resolved tops
        // Thus, initialize the general tops vector with resolved tops
        std::vector<TopObject*> tops(resolvedTops);                          

        // If _not_ doing the resolved-only tagger i.e. doing whole merged+resolved tagger, add in the merged tops to the tops vector
        if (taggerCfg_.find("Resolved") == std::string::npos)
            tops.insert(tops.end(), mergedTops.begin(), mergedTops.end());
        
        countTops(tops);
        
        // -------------------------------------
        // -- Calculate DeltaR between 2 tops
        // -------------------------------------
        double dR_top1_top2 = -1;
        if(ntops_ == 2)
        {
            dR_top1_top2 = tops.at(0)->p().DeltaR(tops.at(1)->p()); 
        }

        // -------------------------------------------------------------
        // -- get vectors of the top mass & eta & pT & LV in the event
        // -------------------------------------------------------------
        auto& topsMass = tr.createDerivedVec<double>("topsMass"+myVarSuffix_);
        auto& topsEta  = tr.createDerivedVec<double>("topsEta"+myVarSuffix_);
        auto& topsPhi  = tr.createDerivedVec<double>("topsPhi"+myVarSuffix_);
        auto& topsPt   = tr.createDerivedVec<double>("topsPt"+myVarSuffix_);
        auto& topsLV   = tr.createDerivedVec<utility::LorentzVector>("topsLV"+myVarSuffix_);

        for(const auto* t : tops)
        {
            topsMass.push_back(t->p().M());
            topsEta.push_back(t->p().Eta());
            topsPhi.push_back(t->p().Phi());
            topsPt.push_back(t->p().Pt());
            topsLV.push_back(utility::convertLV<utility::LorentzVector, TLorentzVector>(t->p()));
        }

        // ----------------------------------
        // -- get best top mass & eta & pt 
        // ----------------------------------
        double bestTopMass             = -9999.9;
        double bestTopEta              = -9999.9;
        double bestTopPhi              = -9999.9;
        double bestTopPt               = -9999.9;
        const TopObject* bestTopMassLV = nullptr;
        bool bestTopMassGenMatch       = false;
        bool bestTopMassTopTag         = false;
        for(unsigned int iTop = 0; iTop < tops.size(); ++iTop)
        {
            auto* top = tops[iTop];

            if(fabs(top->p().M() - 173.21) < fabs(bestTopMass - 173.21))
            {
                bestTopMass   = top->p().M();
                bestTopEta    = top->p().Eta();
                bestTopPhi    = top->p().Phi();
                bestTopPt     = top->p().Pt();
                bestTopMassLV = top;
            }     
        }
        bestTopMassGenMatch = (bestTopMassLV)?(bestTopMassLV->getBestGenTopMatch(0.6) != nullptr):(false);
        
        for(const auto& topPtr : tops)
        {
            if(topPtr == bestTopMassLV)
            {
                bestTopMassTopTag = true;
                break;
            }
        }

        // -----------------------------------------------
        // -- get variables for fake rate & efficiency 
        // -----------------------------------------------
        const auto& candidateTops          = ttr.getTopCandidates();
        float highestDisc                  = -9999.9;
        float bestTopMassCand              = -9999.9;
        float bestTopEtaCand               = -9999.9;
        const TopObject* bestTopMassLVCand = nullptr;
        float bestTopMassTopTagDisc        = -999.9;
        bool bestTopMassGenMatchCand       = false;
        bool bestTopMassTopTagCand         = false;
        for(unsigned int iTop = 0; iTop < candidateTops.size(); ++iTop)
        {
            auto& top = candidateTops[iTop];

            highestDisc = (top.getDiscriminator() > highestDisc ? top.getDiscriminator() : highestDisc);

            if(fabs(top.p().M() - 173.5) < fabs(bestTopMassCand - 173.5) && top.getNConstituents() == 3)
            {
                bestTopMassCand   = top.p().M();
                bestTopEtaCand    = top.p().Eta();
                bestTopMassLVCand = &top;
                bestTopMassTopTagDisc = top.getDiscriminator();
            }
        }
        bestTopMassGenMatchCand = (bestTopMassLVCand)?(bestTopMassLVCand->getBestGenTopMatch(0.6) != nullptr):(false);

        for(const auto& topPtr : tops) 
        {
            if(topPtr == bestTopMassLVCand) 
            {
                bestTopMassTopTagCand = true;
                break;
            }
        }
       
        // --------------------------------------------------------
        // -- get the variables for Efficiency & Fake Rate & ROC
        // --------------------------------------------------------
        double genTopMatchMass_R = -9999.9, genTopMatchEta_R = -9999.9, genTopMatchPhi_R = -9999.9, genTopMatchPt_R = -9999.9;
        double genTopMatchMass_M = -9999.9, genTopMatchEta_M = -9999.9, genTopMatchPhi_M = -9999.9, genTopMatchPt_M = -9999.9;
        double genTopMatchMass_C = -9999.9, genTopMatchEta_C = -9999.9, genTopMatchPhi_C = -9999.9, genTopMatchPt_C = -9999.9;
        bool isFakeRate_R = false, isFakeRate_M = false, isFakeRate_C = false;
        double topDiscGenMatch = -9999.9, topDiscNotGenMatch = -9999.9; 
 
        for(const auto* top : tops)
        {
            // ----------
            // Efficiency
            // ---------- 
            const auto* genMatchTop = top->getBestGenTopMatch();

            if (genMatchTop)
            {
                // gen match top variables for Resolved
                if (top->getNConstituents() == 3)
                {
                    genTopMatchMass_R = top->p().M();
                    genTopMatchEta_R  = top->p().Eta();
                    genTopMatchPhi_R  = top->p().Phi();
                    genTopMatchPt_R   = top->p().Pt();
                    
                }

                // gen match top variables for Merged
                else if (top->getNConstituents() == 1)
                {
                    genTopMatchMass_M = top->p().M();
                    genTopMatchEta_M  = top->p().Eta();
                    genTopMatchPhi_M  = top->p().Phi();
                    genTopMatchPt_M   = top->p().Pt();
                
                }
                 
                // gen match top variables for Resolved + Merged
                else if (top->getNConstituents() == 3 && top->getNConstituents() == 1)
                {
                    genTopMatchMass_C = top->p().M();
                    genTopMatchEta_C  = top->p().Eta();
                    genTopMatchPhi_C  = top->p().Phi();
                    genTopMatchPt_C   = top->p().Pt();

                }
            }   
        
            // ---------
            // Fake Rate
            // --------- 
            if (top->getBestGenTopMatch() == nullptr)
            {

                if (top->getNConstituents() == 3)
                {
                    isFakeRate_R = true;
                }

                else if (top->getNConstituents() == 1)
                {
                    isFakeRate_M = true;
                }

                else if (top->getNConstituents() == 3 && top->getNConstituents() == 1)
                {
                    isFakeRate_C = true;
                }
            }
            
            // ---------------
            // MVA Score - ROC
            // ---------------
            if (top->getBestGenTopMatch() != nullptr)
            {
                topDiscGenMatch = top->getDiscriminator();
            }

            else
            {
                topDiscNotGenMatch = top->getDiscriminator();;
            }

        } 


        // Making tight photon lv (should live somewhere else: is needed for HistoContainer.h)
        const auto& Photons        = tr.getVec<utility::LorentzVector>("Photons");
        const auto& Photons_fullID = tr.getVec<bool>("Photons_fullID");

        auto& tightPhotons = tr.createDerivedVec<TLorentzVector>("tightPhotons"+myVarSuffix_);
        for(unsigned int i = 0; i < Photons.size(); ++i)
        {
            if(Photons_fullID[i])
            {
                tightPhotons.push_back(utility::convertLV<TLorentzVector, utility::LorentzVector>(Photons[i]));
            }
        }
        
        // Register Variables
        tr.registerDerivedVar("ttr"+myVarSuffix_, &ttr);
        tr.registerDerivedVar("ntops"+myVarSuffix_, ntops_);
        tr.registerDerivedVar("ntops_3jet"+myVarSuffix_, ntops_3jet_);
        tr.registerDerivedVar("ntops_2jet"+myVarSuffix_, ntops_2jet_);
        tr.registerDerivedVar("ntops_1jet"+myVarSuffix_, ntops_1jet_);
        tr.registerDerivedVar("dR_top1_top2"+myVarSuffix_,dR_top1_top2);
        tr.registerDerivedVar("bestTopMass"+myVarSuffix_, bestTopMass);
        tr.registerDerivedVar("bestTopEta"+myVarSuffix_, bestTopEta);
        tr.registerDerivedVar("bestTopPt"+myVarSuffix_, bestTopPt);
        tr.registerDerivedVar("bestTopPhi"+myVarSuffix_, bestTopPhi);
        tr.registerDerivedVar("bestTopMassLV"+myVarSuffix_, bestTopMassLV?(bestTopMassLV->p()):(TLorentzVector()));
        tr.registerDerivedVar("bestTopMassGenMatch"+myVarSuffix_, bestTopMassGenMatch);
        tr.registerDerivedVar("bestTopMassTopTag"+myVarSuffix_, bestTopMassTopTag);
        tr.registerDerivedVar("highestDisc"+myVarSuffix_, highestDisc);
        tr.registerDerivedVar("bestTopMassCand"+myVarSuffix_, bestTopMassCand);
        tr.registerDerivedVar("bestTopEtaCand"+myVarSuffix_, bestTopEtaCand);
        tr.registerDerivedVar("bestTopMassLVCand"+myVarSuffix_, bestTopMassLVCand?(bestTopMassLVCand->p()):(TLorentzVector()));
        tr.registerDerivedVar("bestTopMassTopTagDisc"+myVarSuffix_, bestTopMassTopTagDisc);
        tr.registerDerivedVar("bestTopMassGenMatchCand"+myVarSuffix_, bestTopMassGenMatchCand);
        tr.registerDerivedVar("bestTopMassTopTagCand"+myVarSuffix_, bestTopMassTopTagCand);
        // variables for efficiency and fake rate for Resolved, Merged, Resolved + Merged
        tr.registerDerivedVar("genTopMatchMass_R"+myVarSuffix_,  genTopMatchMass_R );
        tr.registerDerivedVar("genTopMatchEta_R"+myVarSuffix_,   genTopMatchEta_R  );
        tr.registerDerivedVar("genTopMatchPhi_R"+myVarSuffix_,   genTopMatchPhi_R  );
        tr.registerDerivedVar("genTopMatchPt_R"+myVarSuffix_,    genTopMatchPt_R   );
        tr.registerDerivedVar("genTopMatchMass_M"+myVarSuffix_,  genTopMatchMass_M );
        tr.registerDerivedVar("genTopMatchEta_M"+myVarSuffix_,   genTopMatchEta_M  );
        tr.registerDerivedVar("genTopMatchPhi_M"+myVarSuffix_,   genTopMatchPhi_M  );
        tr.registerDerivedVar("genTopMatchPt_M"+myVarSuffix_,    genTopMatchPt_M   );
        tr.registerDerivedVar("genTopMatchMass_C"+myVarSuffix_,  genTopMatchMass_C );
        tr.registerDerivedVar("genTopMatchEta_C"+myVarSuffix_,   genTopMatchEta_C  );
        tr.registerDerivedVar("genTopMatchPhi_C"+myVarSuffix_,   genTopMatchPhi_C  );
        tr.registerDerivedVar("genTopMatchPt_C"+myVarSuffix_,    genTopMatchPt_C   );
        tr.registerDerivedVar("isFakeRate_R"+myVarSuffix_,       isFakeRate_R      );
        tr.registerDerivedVar("isFakeRate_M"+myVarSuffix_,       isFakeRate_M      );
        tr.registerDerivedVar("isFakeRate_C"+myVarSuffix_,       isFakeRate_C      );
        tr.registerDerivedVar("topDiscGenMatch"+myVarSuffix_,    topDiscGenMatch   );
        tr.registerDerivedVar("topDiscNotGenMatch"+myVarSuffix_, topDiscNotGenMatch);

    }

public:
    RunTopTagger(std::string taggerCfg = "TopTagger.cfg", std::string myVarSuffix = "") 
        : taggerCfg_          (taggerCfg)
        , myVarSuffix_        (myVarSuffix)
        , tt_                 (new TopTagger())
        , hadtops_            (nullptr)
        , hadtopdaughters_    (nullptr)
        , hadtopdaughters_id_ (nullptr) 
        , neutralinos_        (nullptr)
        , singlinos_          (nullptr) 
        , singlets_           (nullptr)
        , hadtops_idx_        (nullptr)
    {                
        std::cout<<"Setting up RunTopTagger"<<std::endl;
        tt_->setCfgFile(taggerCfg_);
        std::cout<<"Using "+taggerCfg+" as the TopTagger config file"<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        runTopTagger(tr);
    }
};

#endif
