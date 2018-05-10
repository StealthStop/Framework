#include "Framework/Framework/include/SetUpTopTagger.h"

#include <string>
#include <iostream>

SetUpTopTagger::SetUpTopTagger(NTupleReader& tr, 
			       const std::vector<TLorentzVector>& hadtops, 
			       const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters) : 
    tr_                        (tr),
    AK4Inputs_                 (nullptr),
    AK8Inputs_                 (nullptr),
    Jets_                      (tr.getVec<TLorentzVector>("Jets")),
    Jets_bDiscriminatorCSV_    (tr.getVec<double>("Jets_bDiscriminatorCSV")),
    Jets_qgLikelihood_         (tr.getVec<double>("Jets_qgLikelihood")),
    JetsAK8_                   (tr.getVec<TLorentzVector>("JetsAK8")),
    JetsAK8_NsubjettinessTau1_ (tr.getVec<double>("JetsAK8_NsubjettinessTau1")),
    JetsAK8_NsubjettinessTau2_ (tr.getVec<double>("JetsAK8_NsubjettinessTau2")),
    JetsAK8_NsubjettinessTau3_ (tr.getVec<double>("JetsAK8_NsubjettinessTau3")),
    JetsAK8_softDropMass_      (tr.getVec<double>("JetsAK8_softDropMass")),
    JetsAK8_subjets_           (tr.getVec<std::vector<TLorentzVector>>("JetsAK8_subjets")),
    //JetsAk8_subjets_bDiscriminatorCSV_(tr.getVec<std::vector<double>>("JetsAk8_subjets_bDiscriminatorCSV")),
    //JetsAk8_subjets_multiplicity_     (tr.getVec<std::vector<double>>("JetsAk8_subjets_multiplicity")),
    //JetsAk8_subjets_ptD_              (tr.getVec<std::vector<double>>("JetsAk8_subjets_ptD")),
    //JetsAk8_subjets_axismajor_        (tr.getVec<std::vector<double>>("JetsAk8_subjets_axismajor")),
    //JetsAk8_subjets_axisminor_        (tr.getVec<std::vector<double>>("JetsAk8_subjets_axisminor")),
    hadtops_                   (hadtops),
    hadtopdaughters_           (hadtopdaughters)
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
                                               //JetsAk8_subjets_bDiscriminatorCSV_,
                                               //JetsAk8_subjets_multiplicity_,
                                               //JetsAk8_subjets_ptD_,
                                               //JetsAk8_subjets_axismajor_,
                                               //JetsAk8_subjets_axisminor_,
                                               hadtops_,
                                               hadtopdaughters_);  
    //Add variables that are not passed to the constructor
    addVariables();      
}

SetUpTopTagger::~SetUpTopTagger()
{
    delete AK4Inputs_;
    delete AK8Inputs_;
}

std::vector<double>* SetUpTopTagger::intVecTodoubleVec(NTupleReader& tr, const std::string& name)
{
    const auto& vI = tr.getVec<int>(name);
    std::vector<double>* vD = new std::vector<double>(vI.size());
    for(int i = 0; i < vI.size(); i++)
    {
        (*vD)[i] = vI[i];
    }
    tr.registerDerivedVec(name+"ConvertedToDouble", vD);
    return vD;
}

void SetUpTopTagger::addVariables()
{
    //AK4Inputs_->addSupplamentalVector("qgLikelihood",                        tr_.getVec<double>("Jets_qgLikelihood"));
    AK4Inputs_->addSupplamentalVector("qgPtD",                               tr_.getVec<double>("Jets_ptD"));
    AK4Inputs_->addSupplamentalVector("qgAxis1",                             tr_.getVec<double>("Jets_axismajor"));
    AK4Inputs_->addSupplamentalVector("qgAxis2",                             tr_.getVec<double>("Jets_axisminor"));
    //AK4Inputs_->addSupplamentalVector("recoJetschargedHadronEnergyFraction", tr_.getVec<double>("Jets_chargedHadronEnergyFraction"));
    //AK4Inputs_->addSupplamentalVector("recoJetschargedEmEnergyFraction",     tr_.getVec<double>("Jets_chargedEmEnergyFraction"));
    //AK4Inputs_->addSupplamentalVector("recoJetsneutralEmEnergyFraction",     tr_.getVec<double>("Jets_neutralEmEnergyFraction"));
    //AK4Inputs_->addSupplamentalVector("recoJetsmuonEnergyFraction",          tr_.getVec<double>("Jets_muonEnergyFraction"));
    //AK4Inputs_->addSupplamentalVector("recoJetsHFHadronEnergyFraction",      tr_.getVec<double>("Jets_hfHadronEnergyFraction"));
    //AK4Inputs_->addSupplamentalVector("recoJetsHFEMEnergyFraction",          tr_.getVec<double>("Jets_hfEMEnergyFraction"));
    //AK4Inputs_->addSupplamentalVector("recoJetsneutralEnergyFraction",       *add2Vec<double>(tr_,"Jets_neutralEmEnergyFraction","Jets_neutralHadronEnergyFraction"));
    //AK4Inputs_->addSupplamentalVector("PhotonEnergyFraction",                tr_.getVec<double>("Jets_photonEnergyFraction"));
    //AK4Inputs_->addSupplamentalVector("ElectronEnergyFraction",              tr_.getVec<double>("Jets_electronEnergyFraction"));
    //AK4Inputs_->addSupplamentalVector("ChargedHadronMultiplicity",           *intVecTodoubleVec(tr_,"Jets_chargedHadronMultiplicity"));
    //AK4Inputs_->addSupplamentalVector("NeutralHadronMultiplicity",           *intVecTodoubleVec(tr_,"Jets_neutralMultiplicity"));
    //AK4Inputs_->addSupplamentalVector("PhotonMultiplicity",                  *intVecTodoubleVec(tr_,"Jets_photonMultiplicity"));
    //AK4Inputs_->addSupplamentalVector("ElectronMultiplicity",                *intVecTodoubleVec(tr_,"Jets_electronMultiplicity"));
    //AK4Inputs_->addSupplamentalVector("MuonMultiplicity",                    *intVecTodoubleVec(tr_,"Jets_muonMultiplicity"));
    //AK4Inputs_->addSupplamentalVector("DeepCSVb",                            tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("DeepCSVc",                            tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("DeepCSVl",                            tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("DeepCSVbb",                           tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("DeepCSVcc",                           tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("DeepFlavorb",                         tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("DeepFlavorbb",                        tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("DeepFlavorlepb",                      tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("DeepFlavorc",                         tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("DeepFlavoruds",                       tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("DeepFlavorg",                         tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("CvsL",                                tr.getVec<double>(""));
    //AK4Inputs_->addSupplamentalVector("CvsB",                                tr.getVec<double>(""));
    AK4Inputs_->addSupplamentalVector("qgMult",  *intVecTodoubleVec(tr_,"Jets_multiplicity"));
}

// Create jets constituents list combining AK4 and AK8 jets, these are used to construct top candiates
// The vector of input constituents can also be constructed "by hand"    
std::vector<Constituent> SetUpTopTagger::getConstituents() const
{
    return ttUtility::packageConstituents(*AK4Inputs_, *AK8Inputs_);
}
