#ifndef SETUPTOPTAGGER_H
#define SETUPTOPTAGGER_H

//mandatory includes to use top tagger
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include <vector>

class NtupleClass;
class NTupleReader;

class SetUpTopTagger{
 private:
  const NtupleClass& nTupleClass_;
  const NTupleReader& tr_;
  const std::vector<TLorentzVector>& Jets_;                      
  const std::vector<double>& Jets_bDiscriminatorCSV_;    
  const std::vector<double>& Jets_qgLikelihood_;         
  const std::vector<TLorentzVector>& JetsAK8_;                   
  const std::vector<double>& JetsAK8_NsubjettinessTau1_; 
  const std::vector<double>& JetsAK8_NsubjettinessTau2_; 
  const std::vector<double>& JetsAK8_NsubjettinessTau3_; 
  const std::vector<double>& JetsAK8_softDropMass_;      
  const std::vector<std::vector<TLorentzVector>>& JetsAK8_subjets_;           
  const std::vector<double>& Jets_ptD_;                  
  const std::vector<double>& Jets_axismajor_;            
  const std::vector<double>& Jets_axisminor_;            
  const std::vector<int>& Jets_multiplicity_;         
  const std::vector<TLorentzVector>& hadtops_;                   
  const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters_;           

  ttUtility::ConstAK4Inputs* AK4Inputs_;
  ttUtility::ConstAK8Inputs* AK8Inputs_;
  std::vector<double> vecDouble_;

  std::vector<double> intVecTodoubleVec(const std::vector<int>& vecInt);
  void addVariables();

 public:  
  SetUpTopTagger(const NtupleClass& nTupleClass, const std::vector<TLorentzVector>& hadtops, const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters);
  SetUpTopTagger(const NTupleReader& tr, const std::vector<TLorentzVector>& hadtops, const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters);
  ~SetUpTopTagger();
  std::vector<Constituent> getConstituents();
};

#endif
