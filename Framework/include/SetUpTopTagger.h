#ifndef SETUPTOPTAGGER_H
#define SETUPTOPTAGGER_H

//mandatory includes to use top tagger
#include "TopTagger/TopTagger/include/TopTagger.h"
#include "TopTagger/TopTagger/include/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/include/TopTaggerResults.h"
#include <vector>

class NtupleClass;

class SetUpTopTagger{
 private:
  const NtupleClass& nTupleClass_;
  ttUtility::ConstAK4Inputs* AK4Inputs_;
  ttUtility::ConstAK8Inputs* AK8Inputs_;
  std::vector<double> vecDouble_;

 public:  
  SetUpTopTagger(const NtupleClass& nTupleClass, const std::vector<TLorentzVector>& hadtops, const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters);
  ~SetUpTopTagger();
  std::vector<double> intVecTodoubleVec(const std::vector<int>& vecInt);
  void addVariables();
  std::vector<Constituent> getConstituents();
};

#endif
