#ifndef MT2Jets_h
#define MT2Jets_h

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

class MT2Jets
{
private:
    std::string myVarSuffix_;

    // -------------------------------------
    // -- Top-Tagged Jets for hemispheres
    // -------------------------------------
    void getMT2Jets(NTupleReader& tr) const
    {
        const auto* ttr           = tr.getVar<TopTaggerResults*>("ttr");
        const auto& Jets          = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt45 = tr.getVec<bool>("GoodJets_pt45");

        // create an index for resolved tops
        std::vector<TLorentzVector> topJets;
        std::set<unsigned int> usedIndex;
        //auto& usedIndex = tr.createDerivedVec<unsigned int>("usedIndex"+myVarSuffix_);
        const std::vector<TopObject*>& taggedObjects = ttr->getTops();
        for(auto* t : taggedObjects)
        {
            if(t->getType()==TopObject::RESOLVED_TOP) 
            {
                //topJets.push_back(t->P());
                TLorentzVector top;
                const std::vector<const Constituent*>& constituents = t->getConstituents();
                for(const auto& c : constituents)
                {
                    unsigned int index = c->getIndex(); 
                    usedIndex.insert(index);
                    //usedIndex.push_back(index);
                    //if(GoodJets_pt45[index]) 
                    top += c->P();
                }
            
                topJets.push_back(top);
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
            if(!GoodJets_pt45[i]) continue;
            if ( std::find(usedIndex.begin(), usedIndex.end(), i) == usedIndex.end() ) 
            {
                notTopJets.push_back(Jets[i]);
            }
        }
    
        // get the MT2Jets : add tops and not tops jets to each other
        auto& MT2Jets = tr.createDerivedVec<TLorentzVector>("MT2Jets"+myVarSuffix_, topJets);
        MT2Jets.insert(MT2Jets.end(), notTopJets.begin(), notTopJets.end());
        auto& GoodMT2Jets = tr.createDerivedVec<bool>("GoodMT2Jets"+myVarSuffix_, MT2Jets.size(), true);      
        tr.createDerivedVar<int>("NGoodMT2Jets"+myVarSuffix_, GoodMT2Jets.size());
        
        tr.registerDerivedVar("resolvedMass"+myVarSuffix_, resolvedMass);
        tr.registerDerivedVar("resolvedEta"+myVarSuffix_, resolvedEta);
        tr.registerDerivedVar("resolvedPhi"+myVarSuffix_, resolvedPhi);
        tr.registerDerivedVar("resolvedPt"+myVarSuffix_, resolvedPt);   
    }

public:    
    MT2Jets(const std::string& myVarSuffix = "")
        :myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up MT2Jets"<<std::endl;;
    }

    void operator()(NTupleReader& tr)
    {
        getMT2Jets(tr);
    }
};

#endif
