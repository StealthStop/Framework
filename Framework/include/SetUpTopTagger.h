#ifndef SETUPTOPTAGGER_H
#define SETUPTOPTAGGER_H

//mandatory includes to use top tagger
#include "TopTagger/TopTagger/interface/TopTagger.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include <vector>

#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

class SetUpTopTagger
{
private:
    NTupleReader& tr_;
    std::string myVarSuffix_;
    ttUtility::ConstAK4Inputs<double>* AK4Inputs_;
    ttUtility::ConstAK8Inputs<double>* AK8Inputs_;
    std::vector<uint8_t>* ak4Filter_;
    const std::vector<TLorentzVector>& Jets_;                      
    const std::vector<float>& Jets_bJetTagDeepCSVtotb_;
    const std::vector<float>& Jets_qgLikelihood_;        
    const std::vector<bool>& GoodJets_;
    const std::vector<bool>& GoodJets_pt20_; 
    const std::vector<TLorentzVector>& JetsAK8_;                   
    const std::vector<float>& JetsAK8_DeepTagTvsQCD_;
    const std::vector<float>& JetsAK8_DeepTagWvsQCD_;
    const std::vector<float>& JetsAK8_softDropMass_;      
    const std::vector<std::vector<TLorentzVector>>& JetsAK8_subjets_;           
    const std::vector<TLorentzVector>& hadtops_;                   
    const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters_;
               
        
    std::vector<double>* intVecTodoubleVec(NTupleReader& tr, const std::string& vType);
    std::vector<double>* floatVecTodoubleVec(const std::vector<float>& vF);
    std::vector<std::vector<double>>* VecVecintToVecVecdouble(NTupleReader& tr, const std::string& name);

    template<typename I, typename J> std::vector<I>* add2Vec(NTupleReader& tr, const std::string& name1, const std::string& name2)
    {
        const auto& vec1 = tr.getVec<J>(name1);
        const auto& vec2 = tr.getVec<J>(name2);
        std::vector<I>* sumVec = new std::vector<I>(vec1.size());
        for(unsigned int i = 0; i < vec1.size(); i++)
        {
            (*sumVec)[i] = vec1[i] + vec2[i];
        }
        tr.registerDerivedVec(name1+"AddedTo"+name2, sumVec);
        return sumVec;
    }

    void addVariables();
    
public:  
    SetUpTopTagger(NTupleReader& tr, const std::vector<TLorentzVector>& hadtops, const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters, const std::string& myVarSuffix);
    ~SetUpTopTagger();
    std::vector<Constituent> getConstituents() const;
};

#endif
