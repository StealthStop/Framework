#ifndef RUNTOPTAGGER_H
#define RUNTOPTAGGER_H

#include "Framework/Framework/include/SetUpTopTagger.h"

class RunTopTagger
{
private:
    std::shared_ptr<TopTagger> tt_;
    std::vector<TLorentzVector>* hadtops_;
    std::vector<std::vector<const TLorentzVector*>>* hadtopdaughters_;
    std::vector<TLorentzVector>* hadWs_;
    std::vector<TLorentzVector>* neutralinos_;
    std::vector<TLorentzVector>* singlinos_;
    std::vector<TLorentzVector>* singlets_;
    std::vector<int>* hadtops_idx_;
    std::string taggerCfg_;
    int nhadWs_ ;
    int ntops_;
    int ntops_3jet_;
    int ntops_2jet_;
    int ntops_1jet_;

    void genMatch(NTupleReader& tr)
    {
        // ----------------------------------------------
        // check for number of hadronic tops at gen level
        // ----------------------------------------------
        hadtops_         = new std::vector<TLorentzVector>();
        hadtopdaughters_ = new std::vector<std::vector<const TLorentzVector*>>();
        hadWs_           = new std::vector<TLorentzVector>();
        neutralinos_     = new std::vector<TLorentzVector>();
        singlinos_       = new std::vector<TLorentzVector>();
        singlets_        = new std::vector<TLorentzVector>();
        hadtops_idx_     = new std::vector<int>();
        nhadWs_          = 0;

        const std::string& runtype = tr.getVar<std::string>("runtype");

        if(runtype != "Data")
        {
            const auto& GenParticles            = tr.getVec<TLorentzVector>("GenParticles");
            const auto& GenParticles_PdgId      = tr.getVec<int>("GenParticles_PdgId");
            const auto& GenParticles_ParentId   = tr.getVec<int>("GenParticles_ParentId");
            const auto& GenParticles_ParentIdx  = tr.getVec<int>("GenParticles_ParentIdx");
            const auto& GenParticles_Status     = tr.getVec<int>("GenParticles_Status");            

            for ( unsigned int gpi=0; gpi < GenParticles.size() ; gpi++ ) 
            {
                int pdgid = abs( GenParticles_PdgId.at(gpi) ) ;
                int momid = abs( GenParticles_ParentId.at(gpi) ) ;
                int momidx = GenParticles_ParentIdx.at(gpi);
                int status = GenParticles_Status.at(gpi);
                if(pdgid == 1000022 && (status==22 || status == 52))
                {
                    neutralinos_->push_back(GenParticles.at(gpi));
                }
                if(pdgid == 5000001 && (status == 22 || status == 52))
                {
                    singlinos_->push_back(GenParticles.at(gpi));
                }
                if(pdgid == 5000002 && (status == 22 || status == 52))
                {
                    singlets_->push_back(GenParticles.at(gpi));
                }
                if(status == 23 && momid == 24 && pdgid < 6)
                {
                    // Should be the quarks from W decay
                    nhadWs_++;
                    // find the top
                    int Wmotherid = GenParticles_ParentId.at(momidx);
                    if (abs(Wmotherid) == 6){
                        int Wmotheridx = GenParticles_ParentIdx.at(momidx);
                        std::vector<int>::iterator found = std::find(hadtops_idx_->begin(), hadtops_idx_->end(), Wmotheridx);
                        if (found != hadtops_idx_->end())
                        {
                            // already found before
                            // std::cout << "Found this top before: " << *found << std::endl;
                            int position = distance(hadtops_idx_->begin(),found);
                            // add the daughter to the list
                            (*hadtopdaughters_)[position].push_back(&(GenParticles.at(gpi)));
                        } else
                        {
                            // not yet found
                            hadtops_idx_->push_back(Wmotheridx);
                            hadtops_->push_back(GenParticles.at(Wmotheridx));
                            hadWs_->push_back(GenParticles.at(momidx));
                            std::vector<const TLorentzVector*> daughters;
                            daughters.push_back(&(GenParticles.at(gpi)));
                            hadtopdaughters_->push_back(daughters);
                            //std::cout << "Found a new top at idx " << Wmotheridx << std::endl;
                        }
                    }
                } 
            }
            // Now check the b quarks (we only want the ones associated with a hadronic W decay for now)
            for ( unsigned int gpi=0; gpi < GenParticles.size() ; gpi++ ) 
            {
                int pdgid = abs( GenParticles_PdgId.at(gpi) ) ;
                int momid = abs( GenParticles_ParentId.at(gpi) ) ;
                int momidx = GenParticles_ParentIdx.at(gpi);
                int status = GenParticles_Status.at(gpi);
              
                if(status == 23 && momid == 6 && pdgid == 5)
                {
                    // found a b quark from top decay, need to add this to the list of daughters
                    std::vector<int>::iterator found = std::find(hadtops_idx_->begin(), hadtops_idx_->end(), momidx);
                    if (found != hadtops_idx_->end())
                    {
                        // already found
                        int position = distance(hadtops_idx_->begin(),found);
                        (*hadtopdaughters_)[position].push_back(&(GenParticles.at(gpi)));
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
    }

    void countTops(std::vector<TopObject*>* tops)
    {
        ntops_ = 0;
        ntops_3jet_ = 0;
        ntops_2jet_ = 0;
        ntops_1jet_ = 0;
        for (const TopObject* top : *tops)
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
        SetUpTopTagger st( tr, (*hadtops_), (*hadtopdaughters_) );
        std::vector<Constituent> constituents = st.getConstituents();

        // Run the top tagger             
        tt_->runTagger(constituents);

        // Get the top tagger results object     
        const TopTaggerResults& ttr = tt_->getResults();

        // Get reconstructed tops and derive needed variables                            
        std::vector<TopObject*>* tops = new std::vector<TopObject*>(ttr.getTops());
        countTops(tops);

        const auto& candidateTops = ttr.getTopCandidates();
        double bestTopMass = -9999.9;
        double bestTopEta = -9999.9;
        const TopObject* bestTopMassLV = nullptr;
        bool bestTopMassGenMatch = false;
        bool bestTopMassTopTag = false;
        for(int iTop = 0; iTop < candidateTops.size(); ++iTop)
        {
            auto& top = candidateTops[iTop];

            if(fabs(top.p().M() - 173.5) < fabs(bestTopMass - 173.5) && top.getNConstituents() == 3)
            {
                bestTopMass = top.p().M();
                bestTopEta = top.p().Eta();
                bestTopMassLV = &top;
            }
        }

        bestTopMassGenMatch = (bestTopMassLV)?(bestTopMassLV->getBestGenTopMatch(0.6) != nullptr):(false);
        for(const auto& topPtr : (*tops)) 
        {
            if(topPtr == bestTopMassLV) 
            {
                bestTopMassTopTag = true;
                break;
            }
        }

        
        // Making tight photon lv (should live somewhere else: is needed for HistoContainer.h)
        const auto& Photons        = tr.getVec<TLorentzVector>("Photons");
        const auto& Photons_fullID = tr.getVec<bool>("Photons_fullID");

        std::vector<TLorentzVector>* tightPhotons = new std::vector<TLorentzVector>();
        for(int i = 0; i < Photons.size(); ++i)
        {
            if(Photons_fullID[i])
            {
                tightPhotons->push_back(Photons[i]);
            }
        }


        // Register Variables
        tr.registerDerivedVar("ttr", &ttr);
        tr.registerDerivedVar("ntops", ntops_);
        tr.registerDerivedVar("ntops_3jet", ntops_3jet_);
        tr.registerDerivedVar("ntops_2jet", ntops_2jet_);
        tr.registerDerivedVar("ntops_1jet", ntops_1jet_);
        tr.registerDerivedVar("nhadWs", nhadWs_);
        tr.registerDerivedVec("hadtops", hadtops_);
        tr.registerDerivedVec("hadtopdaughters", hadtopdaughters_);
        tr.registerDerivedVec("hadWs", hadWs_);
        tr.registerDerivedVec("neutralinos", neutralinos_);
        tr.registerDerivedVec("singlinos", singlinos_);
        tr.registerDerivedVec("singlets", singlets_);
        tr.registerDerivedVec("hadtops_idx", hadtops_idx_);
        tr.registerDerivedVar("bestTopMassLV", bestTopMassLV?(bestTopMassLV->p()):(TLorentzVector()));
        tr.registerDerivedVar("bestTopMass", bestTopMass);
        tr.registerDerivedVar("bestTopMassTopTag", bestTopMassTopTag);
        tr.registerDerivedVar("bestTopMassGenMatch", bestTopMassGenMatch);
        tr.registerDerivedVar("bestTopEta", bestTopEta);
        tr.registerDerivedVec("tightPhotons",tightPhotons);
    }

public:
    RunTopTagger(std::string taggerCfg = "TopTagger.cfg") :  
        taggerCfg_      (taggerCfg), 
        tt_             (new TopTagger()),
        hadtops_        (nullptr),
        hadtopdaughters_(nullptr),
        hadWs_          (nullptr),
        neutralinos_    (nullptr),
        singlinos_      (nullptr),
        singlets_       (nullptr),
        hadtops_idx_    (nullptr)
    {                
        tt_->setCfgFile(taggerCfg_);
    }

    ~RunTopTagger(){}

    void operator()(NTupleReader& tr)
    {
        runTopTagger(tr);
    }
};

#endif