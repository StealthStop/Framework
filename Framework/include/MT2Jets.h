#ifndef MT2Jets_h
#define MT2Jets_h

#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TLorentzVector.h" 
#include <vector>
#include <iostream>
#include <cmath>

class MT2Jets
{
private:
    std::string jetMaskName_, nJetName_, myVarSuffix_;

    std::vector<TLorentzVector> vector_subtract(std::vector<TLorentzVector> vector1, std::vector<TLorentzVector> vector2)
    {
        std::vector<TLorentzVector> subtraction;
        for (unsigned int i = 0; i< vector1.size(); i++)
        {
                bool add = true;
                for (unsigned int j = 0; j < vector2.size(); j++)
                {
                        if (vector1.at(i) == vector2.at(j))
                        {
                                add = false;
                        }
                }
                if (add)
                {
                    subtraction.push_back(vector1.at(i));
                }
        }
        return subtraction;
    }

    void getHemisphereJets(NTupleReader& tr) const
    {
        const auto& Jets          = tr.getVec<TLorentzVector>("Jets");
        const auto& topsLV        = tr.getVec<TLorentzVector>("topsLV");
        const auto& NGoodJets     = tr.getVar<int>(nJetName_);
        const auto& GoodJets_pt45 = tr.getVec<bool>("GoodJets_pt45");
 
        std::vector<TLorentzVector> topJets, notTopJets;
        if(NGoodJets >= 2)
        {
            for(unsigned int i = 0; i < Jets.size(); ++i)
            {
                if(GoodJets_pt45[i])
                {
                    // to get the non-top jets 
                    bool add = true; 
                    for (unsigned int j = 0; j < topsLV.size(); j++)
                    {
                        if (Jets.at(i) == topsLV.at(j))
                        {
                                add = false;
                        }
                    }
                    if (add)
                    {
                        notTopJets.push_back(Jets.at(i));
                    }                                  
            
                    // sum over of the top jets
                    for(unsigned int t = 0; t < topsLV.size(); t++)
                    {
                        topJets.push_back(topsLV.at(t)); 
                    }
                }
            }
            //std::vector<TLorentzVector> MT2Jets;
            //MT2Jets = notTopJets.push_back(topJets);

        } 
        //tr.registerDerivedVar("MT2Jets"+myVarSuffix_,MT2Jets);
            
    }

public:    
    MT2Jets(const std::string& myVarSuffix = "")
        :myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up MT2Jets"<<std::endl;;
    }

    void operator()(NTupleReader& tr)
    {
        getHemisphereJets(tr);
    }
};

#endif
