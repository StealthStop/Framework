#include "Framework/Framework/include/SetUpTopTagger.h"

#include <string>
#include <iostream>

SetUpTopTagger::SetUpTopTagger(NTupleReader& tr, 
			       const std::vector<TLorentzVector>& hadtops, 
			       const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters,
                               const std::string& myVarSuffix) : 
    tr_                        (tr),
    myVarSuffix_               (myVarSuffix),
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
    //JetsAK8_subjets_bDiscriminatorCSV_(tr.getVec<std::vector<double>>("JetsAK8_subjets_bDiscriminatorCSV")),
    //JetsAK8_subjets_ptD_              (tr.getVec<std::vector<double>>("JetsAK8_subjets_ptD")),
    //JetsAK8_subjets_axismajor_        (tr.getVec<std::vector<double>>("JetsAK8_subjets_axismajor")),
    //JetsAK8_subjets_axisminor_        (tr.getVec<std::vector<double>>("JetsAK8_subjets_axisminor")),
    hadtops_                   (hadtops),
    hadtopdaughters_           (hadtopdaughters)
{
    // Use helper function to create input list 
    // Create AK4 inputs object
    AK4Inputs_ = new ttUtility::ConstAK4Inputs<double>(
        Jets_, 
        Jets_bDiscriminatorCSV_,
        Jets_qgLikelihood_, 
        hadtops_, 
        hadtopdaughters_);  
    // Create AK8 inputs object
    AK8Inputs_ = new ttUtility::ConstAK8Inputs<double>(
        JetsAK8_,
        JetsAK8_NsubjettinessTau1_,
        JetsAK8_NsubjettinessTau2_,
        JetsAK8_NsubjettinessTau3_,
        JetsAK8_softDropMass_,
        JetsAK8_subjets_,
        //JetsAK8_subjets_bDiscriminatorCSV_,
        //*VecVecintToVecVecdouble(tr_, "JetsAK8_subjets_multiplicity"),
        //JetsAK8_subjets_ptD_,
        //JetsAK8_subjets_axismajor_,
        //JetsAK8_subjets_axisminor_,
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

std::vector<std::vector<double>>* SetUpTopTagger::VecVecintToVecVecdouble(NTupleReader& tr, const std::string& name)
{
    const auto& vvI = tr.getVec<std::vector<int>>(name);
    std::vector<std::vector<double>>* vvD = new std::vector<std::vector<double>>();
    for(int i = 0; i < vvI.size(); i++)
    {
        std::vector<double> vD;
        for(int j = 0; j < vvI[i].size(); j++)
        {
            vD.push_back(vvI[i][j]);
        }
        vvD->push_back(vD);
    }
    tr.registerDerivedVec(name+"ConvertedToDouble", vvD);
    return vvD;
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
    //AK4Inputs_->addSupplamentalVector("DeepCSVb",                            tr_.getVec<double>("Jets_bJetTagDeepCSVprobb"));
    //AK4Inputs_->addSupplamentalVector("DeepCSVc",                            tr_.getVec<double>("Jets_bJetTagDeepCSVprobc"));
    //AK4Inputs_->addSupplamentalVector("DeepCSVl",                            tr_.getVec<double>("Jets_bJetTagDeepCSVprobudsg"));
    //AK4Inputs_->addSupplamentalVector("DeepCSVbb",                           tr_.getVec<double>("Jets_bJetTagDeepCSVprobbb"));
    //AK4Inputs_->addSupplamentalVector("DeepCSVcc",                           std::vector<double>(tr_.getVec<double>("Jets_bJetTagDeepCSVprobbb").size(), 0.0));
    //AK4Inputs_->addSupplamentalVector("DeepFlavorb",                         tr_.getVec<double>("Jets_bJetTagDeepFlavourprobb"));
    //AK4Inputs_->addSupplamentalVector("DeepFlavorbb",                        tr_.getVec<double>("Jets_bJetTagDeepFlavourprobbb"));
    //AK4Inputs_->addSupplamentalVector("DeepFlavorlepb",                      tr_.getVec<double>("Jets_bJetTagDeepFlavourproblepb"));
    //AK4Inputs_->addSupplamentalVector("DeepFlavorc",                         tr_.getVec<double>("Jets_bJetTagDeepFlavourprobc"));
    //AK4Inputs_->addSupplamentalVector("DeepFlavoruds",                       tr_.getVec<double>("Jets_bJetTagDeepFlavourprobuds"));
    //AK4Inputs_->addSupplamentalVector("DeepFlavorg",                         tr_.getVec<double>("Jets_bJetTagDeepFlavourprobg"));
    //AK4Inputs_->addSupplamentalVector("CvsL",                                tr_.getVec<double>("Jets_bJetTagDeepCSVCvsL"));
    //AK4Inputs_->addSupplamentalVector("CvsB",                                tr_.getVec<double>("Jets_bJetTagDeepCSVCvsB"));
    AK4Inputs_->addSupplamentalVector("qgMult",  *intVecTodoubleVec(tr_,"Jets_multiplicity"));
}

// Create jets constituents list combining AK4 and AK8 jets, these are used to construct top candiates
// The vector of input constituents can also be constructed "by hand"    
std::vector<Constituent> SetUpTopTagger::getConstituents() const
{
    return ttUtility::packageConstituents(*AK4Inputs_, *AK8Inputs_);
}
