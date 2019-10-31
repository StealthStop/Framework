#ifndef MT2Jets_h
#define MT2Jets_h

#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopObject.h"
#include "TopTagger/TopTagger/interface/Constituent.h"
#include "TLorentzVector.h"
#include <iostream> 
#include <vector>
#include <cmath>

class MT2Jets
{
private:
    std::string myVarSuffix_;

    void getMT2Jets(NTupleReader& tr) const
    {
        const auto* ttr       = tr.getVar<TopTaggerResults*>("ttr");
        const auto& Jets      = tr.getVec<TLorentzVector>("Jets");

        // create an index for resolved tops
        std::vector<TLorentzVector> topJets;
        std::set<unsigned int> usedIndex;
        const std::vector<TopObject*>& taggedObjects = ttr->getTops();
        for(auto* t : taggedObjects)
        {
            if(t->getType()==TopObject::RESOLVED_TOP) 
            {
                topJets.push_back(t->P());
                const std::vector<const Constituent*>& constituents = t->getConstituents();
                for (const auto& c : constituents)
                {
                    usedIndex.insert(c->getIndex());
                }
            }
        }
        
        // get the notTopJets by using 'usedIndex' 
        std::vector<TLorentzVector> notTopJets;
        for(unsigned int i = 0; i < Jets.size(); ++i)
        {
            if( std::find(usedIndex.begin(), usedIndex.end(), i) == usedIndex.end() ) 
            {
                notTopJets.push_back(Jets[i]);
            }
        }
    
        // Get the MT2Jets : add tops and not tops jets to each other
        auto& MT2Jets = tr.createDerivedVec<TLorentzVector>("MT2Jets"+myVarSuffix_, notTopJets);
        MT2Jets.insert(MT2Jets.end(), topJets.begin(), topJets.end());
        auto& GoodMT2Jets = tr.createDerivedVec<bool>("GoodMT2Jets"+myVarSuffix_, MT2Jets.size(), true);        
        tr.createDerivedVar<int>("NGoodMT2Jets"+myVarSuffix_, GoodMT2Jets.size());
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
