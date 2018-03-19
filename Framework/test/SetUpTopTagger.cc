#include "SetUpTopTagger.h"
#include "NtupleClass.h"

#include <string>

SetUpTopTagger::SetUpTopTagger(const NtupleClass& nTupleClass, 
			       const std::vector<TLorentzVector>& hadtops, 
			       const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters) : 
  nTupleClass_(nTupleClass), 
  AK4Inputs_(nullptr), 
  AK8Inputs_(nullptr)
{
  
  // Use helper function to create input list 
  // Create AK4 inputs object
  AK4Inputs_ = new ttUtility::ConstAK4Inputs(
					 *nTupleClass_.Jets, 
					 *nTupleClass_.Jets_bDiscriminatorCSV,
					 *nTupleClass_.Jets_qgLikelihood, 
					 hadtops, 
					 hadtopdaughters);
  
  // Create AK8 inputs object
  AK8Inputs_ = new ttUtility::ConstAK8Inputs(
					 *nTupleClass_.JetsAK8,
					 *nTupleClass_.JetsAK8_NsubjettinessTau1,
					 *nTupleClass_.JetsAK8_NsubjettinessTau2,
					 *nTupleClass_.JetsAK8_NsubjettinessTau3,
					 *nTupleClass_.JetsAK8_softDropMass,
					 *nTupleClass_.JetsAK8_subjets,    // These should be the subjets!
					 hadtops,
					 hadtopdaughters);  
  //Add variables that are not passed to the constructor
  addVariables();
  
}

SetUpTopTagger::~SetUpTopTagger(){
  delete AK4Inputs_;
  delete AK8Inputs_;
}

std::vector<double> SetUpTopTagger::intVecTodoubleVec(const std::vector<int>& vecInt){
  vecDouble_.clear();
  for(int i = 0; i < vecInt.size(); i++){
    vecDouble_.push_back(vecInt[i]);
  }
  return vecDouble_;
}

void SetUpTopTagger::addVariables(){
  AK4Inputs_->addSupplamentalVector("qgPtD",   *nTupleClass_.Jets_ptD);
  AK4Inputs_->addSupplamentalVector("qgAxis1", *nTupleClass_.Jets_axismajor);
  AK4Inputs_->addSupplamentalVector("qgAxis2", *nTupleClass_.Jets_axisminor);
  AK4Inputs_->addSupplamentalVector("qgMult",  intVecTodoubleVec(*nTupleClass_.Jets_multiplicity));
}
  
// Create jets constituents list combining AK4 and AK8 jets, these are used to construct top candiates
// The vector of input constituents can also be constructed "by hand"    
std::vector<Constituent> SetUpTopTagger::getConstituents(){
  return ttUtility::packageConstituents(*AK4Inputs_, *AK8Inputs_);
}
