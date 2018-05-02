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

    int findParent(const int idx, const std::vector<int>& GenParticles_ParentId, const std::vector<int>& GenParticles_ParentIdx)
    {
        if (idx == -1)
        {
            return -1;
        }
        else if(abs(GenParticles_ParentId[idx]) == 6)
        {
            return GenParticles_ParentIdx[idx];
        }
        else
        {
            return findParent(GenParticles_ParentIdx[idx], GenParticles_ParentId, GenParticles_ParentIdx);
        }
    }

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
                int pdgid = GenParticles_PdgId.at(gpi);
                int momid = GenParticles_ParentId.at(gpi) ;
                int momidx = GenParticles_ParentIdx.at(gpi);
                int status = GenParticles_Status.at(gpi);
                int topIdx = findParent(gpi, GenParticles_ParentId, GenParticles_ParentIdx);
                //printf(" %6i: status: %6i pdg: %6i motherID: %6i motherIDX: %6i ", gpi,  GenParticles_Status[gpi], GenParticles_PdgId[gpi], GenParticles_ParentId[gpi], GenParticles_ParentIdx[gpi]); fflush(stdout);
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
                if( topIdx >= 0 && (abs(pdgid) != 24) )
                {
                    //printf(" topIdx: %i particle: %i\n", topIdx, pdgid); fflush(stdout);
                    
                    int position = 0;
                    for(;position < hadtops_idx_->size() && (*hadtops_idx_)[position] != topIdx; ++position);
                    if( position < hadtops_idx_->size() )
                    {
                        (*hadtopdaughters_)[position].push_back(&(GenParticles.at(gpi)));
                    } 
                    else
                    {
                        hadtops_idx_->push_back(topIdx);
                        hadtops_->push_back(GenParticles.at(topIdx));
                        std::vector<const TLorentzVector*> daughters;
                        daughters.push_back(&(GenParticles.at(gpi)));
                        hadtopdaughters_->push_back(daughters);
                    }
                }
                else
                {
                    //printf("\n");
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
        const std::vector<Constituent>& constituents = st.getConstituents();

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

        ////std::cout<<" Size Yo "<<hadtops_->size()<<"  "<<hadtopdaughters_->size()<<std::endl;
        //for (int i = 0; i < hadtops_->size(); ++i)
        //{
        //    TLorentzVector dSum;
        //    for (int j = 0; j < ((*hadtopdaughters_)[i]).size(); j++)
        //    {
        //        dSum += *(((*hadtopdaughters_)[i])[j]);
        //    }
        //    printf("nTops: %i ndaughters %i   top: (pt %4.5lf , eta %4.5lf, phi %4.5lf, mass %4.5lf) dSum: (pt %4.5lf , eta %4.5lf, phi %4.5lf, mass %4.5lf)\n", hadtops_->size(), (*hadtopdaughters_)[i].size(), 
        //           (*hadtops_)[i].Pt(), (*hadtops_)[i].Eta(), (*hadtops_)[i].Phi(), (*hadtops_)[i].M(), 
        //           dSum.Pt(), dSum.Eta(), dSum.Phi(), dSum.M()
        //          );
        //
        //}
        //
        //printf("=========================================================================================\n");

        // Register Variables
        tr.registerDerivedVar("ttr", &ttr);
        tr.registerDerivedVar("ntops", ntops_);
        tr.registerDerivedVar("ntops_3jet", ntops_3jet_);
        tr.registerDerivedVar("ntops_2jet", ntops_2jet_);
        tr.registerDerivedVar("ntops_1jet", ntops_1jet_);
        //tr.registerDerivedVar("nhadWs", nhadWs_);
        tr.registerDerivedVec("hadtops", hadtops_);
        tr.registerDerivedVec("hadtopdaughters", hadtopdaughters_);
        //tr.registerDerivedVec("hadWs", hadWs_);
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
        std::cout<<"Using "+taggerCfg+" as the TopTagger config file"<<std::endl;
    }

    ~RunTopTagger(){}

    void operator()(NTupleReader& tr)
    {
        runTopTagger(tr);
    }
};

#endif
