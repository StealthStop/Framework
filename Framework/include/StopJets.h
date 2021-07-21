#ifndef StopJets_h
#define StopJets_h

#include "Framework/Framework/include/Utility.h"

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
        const auto* ttr                   = tr.getVar<TopTaggerResults*>("ttr"+myVarSuffix_);
        const auto& Jets                  = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
        const auto& GoodJets_pt20         = tr.getVec<bool>("GoodJets_pt20"+myVarSuffix_);
        auto& StopJets                    = tr.createDerivedVec<TLorentzVector>("StopJets"+myVarSuffix_);
 
        // --------------------------------- 
        // create an index for resolved tops
        // --------------------------------- 
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
         
                // get top jets   
                StopJets.push_back(top);
                
            }
        }

        // -----------------------
        // sort the StopJets by pt
        // -----------------------
        std::sort(StopJets.begin(), StopJets.end(), utility::compare_pt_TLV);
        if( StopJets.size() > 2 ) 
        {
            std::sort(StopJets.begin()+1, StopJets.end(), [StopJets](const TLorentzVector& v1, const TLorentzVector& v2) 
                {return v1.Pt()*v1.DeltaR(StopJets[0]) > v2.Pt()*v2.DeltaR(StopJets[0]);}
            );
        }
        
        // -------------------------------------
        // get the notTopJets by using usedIndex 
        // -------------------------------------
        for(unsigned int i = 0; i < Jets.size(); ++i)
        {
            if(!GoodJets_pt20[i]) continue;
            if ( std::find(usedIndex.begin(), usedIndex.end(), i) == usedIndex.end() ) 
            {
                StopJets.push_back(Jets[i]);
            }
        }
        auto& GoodStopJets = tr.createDerivedVec<bool>("GoodStopJets"+myVarSuffix_, StopJets.size(), true);
        tr.createDerivedVar<int>("NGoodStopJets"+myVarSuffix_, GoodStopJets.size());   

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
