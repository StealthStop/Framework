#include "Framework/Framework/include/SetUpTopTagger.h"
#include "Framework/Framework/include/NtupleClass.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <string>

SetUpTopTagger::SetUpTopTagger(const NtupleClass& nTupleClass, 
			       const std::vector<TLorentzVector>& hadtops, 
			       const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters) : 
    nTupleClass_               (nTupleClass),
    tr_                        (nullptr),
    AK4Inputs_                 (nullptr),
    AK8Inputs_                 (nullptr),
    Jets_                      (*nTupleClass_.Jets),
    Jets_bDiscriminatorCSV_    (*nTupleClass_.Jets_bDiscriminatorCSV),
    Jets_qgLikelihood_         (*nTupleClass_.Jets_qgLikelihood),
    JetsAK8_                   (*nTupleClass_.JetsAK8),
    JetsAK8_NsubjettinessTau1_ (*nTupleClass_.JetsAK8_NsubjettinessTau1),
    JetsAK8_NsubjettinessTau2_ (*nTupleClass_.JetsAK8_NsubjettinessTau2),
    JetsAK8_NsubjettinessTau3_ (*nTupleClass_.JetsAK8_NsubjettinessTau3),
    JetsAK8_softDropMass_      (*nTupleClass_.JetsAK8_softDropMass),
    JetsAK8_subjets_           (*nTupleClass_.JetsAK8_subjets),
    Jets_ptD_                  (*nTupleClass_.Jets_ptD),
    Jets_axismajor_            (*nTupleClass_.Jets_axismajor),
    Jets_axisminor_            (*nTupleClass_.Jets_axisminor),
    Jets_multiplicity_         (*nTupleClass_.Jets_multiplicity),
    hadtops_                   (hadtops),
    hadtopdaughters_           (hadtopdaughters)
{
}

SetUpTopTagger::SetUpTopTagger(const NTupleReader& tr, 
			       const std::vector<TLorentzVector>& hadtops, 
			       const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters) : 
    nTupleClass_               (nullptr),
    tr_                        (tr),
    AK4Inputs_                 (nullptr),
    AK8Inputs_                 (nullptr),
    Jets_                      (tr_.getVec<TLorentzVector>("Jets")),
    Jets_bDiscriminatorCSV_    (tr_.getVec<double>("Jets_bDiscriminatorCSV")),
    Jets_qgLikelihood_         (tr_.getVec<double>("Jets_qgLikelihood")),
    JetsAK8_                   (tr_.getVec<TLorentzVector>("JetsAK8")),
    JetsAK8_NsubjettinessTau1_ (tr_.getVec<double>("JetsAK8_NsubjettinessTau1")),
    JetsAK8_NsubjettinessTau2_ (tr_.getVec<double>("JetsAK8_NsubjettinessTau2")),
    JetsAK8_NsubjettinessTau3_ (tr_.getVec<double>("JetsAK8_NsubjettinessTau3")),
    JetsAK8_softDropMass_      (tr_.getVec<double>("JetsAK8_softDropMass")),
    JetsAK8_subjets_           (tr_.getVec<std::vector<TLorentzVector>>("JetsAK8_subjets")),
    Jets_ptD_                  (tr_.getVec<double>("Jets_ptD")),
    Jets_axismajor_            (tr_.getVec<double>("Jets_axismajor")),
    Jets_axisminor_            (tr_.getVec<double>("Jets_axisminor")),
    Jets_multiplicity_         (tr_.getVec<int>("Jets_multiplicity")),
    hadtops_                   (hadtops),
    hadtopdaughters_           (hadtopdaughters)
{
}

SetUpTopTagger::~SetUpTopTagger()
{
    delete AK4Inputs_;
    delete AK8Inputs_;
}

std::vector<double> SetUpTopTagger::intVecTodoubleVec(const std::vector<int>& vecInt)
{
    vecDouble_.clear();
    for(int i = 0; i < vecInt.size(); i++)
    {
        vecDouble_.push_back(vecInt[i]);
    }
    return vecDouble_;
}

void SetUpTopTagger::addVariables()
{
    AK4Inputs_->addSupplamentalVector("qgPtD",   Jets_ptD_);
    AK4Inputs_->addSupplamentalVector("qgAxis1", Jets_axismajor_);
    AK4Inputs_->addSupplamentalVector("qgAxis2", Jets_axisminor_);
    AK4Inputs_->addSupplamentalVector("qgMult",  intVecTodoubleVec(Jets_multiplicity_));
}

// Create jets constituents list combining AK4 and AK8 jets, these are used to construct top candiates
// The vector of input constituents can also be constructed "by hand"    
std::vector<Constituent> SetUpTopTagger::getConstituents()
{
    // Use helper function to create input list 
    // Create AK4 inputs object
    AK4Inputs_ = new ttUtility::ConstAK4Inputs(
                                               Jets_, 
                                               Jets_bDiscriminatorCSV_,
                                               Jets_qgLikelihood_, 
                                               hadtops_, 
                                               hadtopdaughters_);  
    // Create AK8 inputs object
    AK8Inputs_ = new ttUtility::ConstAK8Inputs(
                                               JetsAK8_,
                                               JetsAK8_NsubjettinessTau1_,
                                               JetsAK8_NsubjettinessTau2_,
                                               JetsAK8_NsubjettinessTau3_,
                                               JetsAK8_softDropMass_,
                                               JetsAK8_subjets_,
                                               hadtops_,
                                               hadtopdaughters_);  
    //Add variables that are not passed to the constructor
    addVariables();  

    //Return Constituents
    return ttUtility::packageConstituents(*AK4Inputs_, *AK8Inputs_);
}
