#ifndef SETUPTOPTAGGER_H
#define SETUPTOPTAGGER_H

//mandatory includes to use top tagger
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include <vector>

#include "SusyAnaTools/Tools/NTupleReader.h"

class SetUpTopTagger
{
private:
    NTupleReader& tr_;
    const std::vector<TLorentzVector>& Jets_;                      
    const std::vector<double>& Jets_bDiscriminatorCSV_;    
    const std::vector<double>& Jets_qgLikelihood_;         
    const std::vector<TLorentzVector>& JetsAK8_;                   
    const std::vector<double>& JetsAK8_NsubjettinessTau1_; 
    const std::vector<double>& JetsAK8_NsubjettinessTau2_; 
    const std::vector<double>& JetsAK8_NsubjettinessTau3_; 
    const std::vector<double>& JetsAK8_softDropMass_;      
    const std::vector<std::vector<TLorentzVector>>& JetsAK8_subjets_;           
    const std::vector<TLorentzVector>& hadtops_;                   
    const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters_;           
    
    ttUtility::ConstAK4Inputs* AK4Inputs_;
    ttUtility::ConstAK8Inputs* AK8Inputs_;
    
    std::vector<double>* intVecTodoubleVec(NTupleReader& tr, const std::string& vType);

    template<typename I> std::vector<I>* add2Vec(NTupleReader& tr, const std::string& name1, const std::string& name2)
    {
        const auto& vec1 = tr.getVec<I>(name1);
        const auto& vec2 = tr.getVec<I>(name2);
        std::vector<I>* sumVec = new std::vector<I>(vec1.size());
        for(int i = 0; i < vec1.size(); i++)
        {
            (*sumVec)[i] = vec1[i] + vec2[i];
        }
        tr.registerDerivedVec(name1+"AddedTo"+name2, sumVec);
        return sumVec;
    }

    void addVariables();
    
public:  
    SetUpTopTagger(NTupleReader& tr, const std::vector<TLorentzVector>& hadtops, const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters);
    ~SetUpTopTagger();
    std::vector<Constituent> getConstituents() const;
};

#endif
