#ifndef StopJets_h
#define StopJets_h

// for top-tagged jets
#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopObject.h"
#include "TopTagger/TopTagger/interface/Constituent.h"

// for gen level study
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/interface/lester_mt2_bisect.h"

#include "TLorentzVector.h"
#include <iostream> 
#include <vector>
#include <cmath>

class StopJets
{
private:
    std::string myVarSuffix_;

    // -------------------------------------
    // -- Top-Tagged Jets for hemispheres
    // -------------------------------------
    void getStopJets(NTupleReader& tr) const
    {
        const auto* ttr           = tr.getVar<TopTaggerResults*>("ttr");
        const auto& Jets          = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt20 = tr.getVec<bool>("GoodJets_pt20");

        // create an index for resolved tops
        std::vector<TLorentzVector> StopJets_;
        std::set<unsigned int> usedIndex;
        const std::vector<TopObject*>& taggedObjects = ttr->getTops();
        for(auto* t : taggedObjects)
        {
            if(t->getType()==TopObject::RESOLVED_TOP) 
            {
                TLorentzVector top;
                const std::vector<const Constituent*>& constituents = t->getConstituents();
                for(const auto& c : constituents)
                {
                    unsigned int index = c->getIndex(); 
                    usedIndex.insert(index);
                    top += c->P();
                }
            
                StopJets_.push_back(top);
            }
        }

        // get the resolved jets' Mass, Eta, Phi, Pt
        double resolvedMass = -9999.9, resolvedEta = -9999.9, resolvedPhi = -9999.9, resolvedPt = -9999.9;
        for(auto& i : usedIndex)
        {
            TLorentzVector resolved = Jets[i];
            resolvedMass = resolved.M();
            resolvedEta  = resolved.Eta();
            resolvedPhi  = resolved.Phi();
            resolvedPt   = resolved.Pt();
        }

        // get the notTopJets by using 'usedIndex' 
        std::vector<TLorentzVector> notTopJets;
        for(unsigned int i = 0; i < Jets.size(); ++i)
        {
            if(!GoodJets_pt20[i]) continue;
            if ( std::find(usedIndex.begin(), usedIndex.end(), i) == usedIndex.end() ) 
            {
                StopJets_.push_back(Jets[i]);
            }
        }
    
        // get the StopJets : add tops and not tops jets to each other
        auto& StopJets = tr.createDerivedVec<TLorentzVector>("StopJets"+myVarSuffix_, StopJets_);
        auto& GoodStopJets = tr.createDerivedVec<bool>("GoodStopJets"+myVarSuffix_, StopJets.size(), true);      
        tr.createDerivedVar<int>("NGoodStopJets"+myVarSuffix_, GoodStopJets.size());
        
        tr.registerDerivedVar("resolvedMass"+myVarSuffix_, resolvedMass);
        tr.registerDerivedVar("resolvedEta"+myVarSuffix_, resolvedEta);
        tr.registerDerivedVar("resolvedPhi"+myVarSuffix_, resolvedPhi);
        tr.registerDerivedVar("resolvedPt"+myVarSuffix_, resolvedPt);   
    }

public:    
    StopJets(const std::string& myVarSuffix = "")
        :myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up StopJets"<<std::endl;;
    }

    void operator()(NTupleReader& tr)
    {
        getStopJets(tr);
    }
};

#endif
