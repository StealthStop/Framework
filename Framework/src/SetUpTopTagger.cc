#include "Framework/Framework/include/SetUpTopTagger.h"

#include <string>
#include <iostream>

SetUpTopTagger::SetUpTopTagger(NTupleReader& tr, 
			       const std::vector<TLorentzVector>& hadtops, 
			       const std::vector<std::vector<const TLorentzVector*>>& hadtopdaughters,
                             const std::string& myVarSuffix): 
    tr_                      (tr),
    myVarSuffix_             (myVarSuffix),
    AK4Inputs_               (nullptr),
    AK8Inputs_               (nullptr),
    ak4Filter_               (nullptr),
    Jets_                    (tr.getVec<TLorentzVector>("Jets"+myVarSuffix_)),
    Jets_bJetTagDeepCSVtotb_ (tr.getVec<double>("Jets"+myVarSuffix_+"_bJetTagDeepCSVtotb")),
    Jets_qgLikelihood_       (tr.getVec<double>("Jets"+myVarSuffix_+"_qgLikelihood")),
    GoodJets_                (tr.getVec<bool>("GoodJets"+myVarSuffix_)),
    GoodJets_pt20_           (tr.getVec<bool>("GoodJets_pt20"+myVarSuffix_)), 
    JetsAK8_                 (tr.getVec<TLorentzVector>("JetsAK8"+myVarSuffix_)),
    JetsAK8_DeepTagTvsQCD_   (tr.getVec<double>("JetsAK8"+myVarSuffix_+"_DeepTagTvsQCD")),
    JetsAK8_DeepTagWvsQCD_   (tr.getVec<double>("JetsAK8"+myVarSuffix_+"_DeepTagWvsQCD")),
    JetsAK8_softDropMass_    (tr.getVec<double>("JetsAK8"+myVarSuffix_+"_softDropMass")),
    JetsAK8_subjets_         (tr.getVec<std::vector<TLorentzVector>>("JetsAK8"+myVarSuffix+"_subjets")),
    hadtops_                 (hadtops),
    hadtopdaughters_         (hadtopdaughters)
{
    // ------------------------------
    // -- Jet Filter
    // ------------------------------
    std::vector<uint8_t>* ak4Filter_ = new std::vector<uint8_t>(Jets_.size(), true);
    for( unsigned int i = 0; i < ak4Filter_->size(); ++i)
    {
        // Jet cleaning which is pT < 20, for resolved top candidate: pT > 40, 30, 20 GeV on the three jets respectively
        //if(Jets[i].Pt() < 20.0)
        //{
        //    ak4Filter_[i] = false;
        //}
    
        if (GoodJets_pt20_[i])
        {
            (*ak4Filter_)[i] = true;   
        } 
        else 
        {
            (*ak4Filter_)[i] = false;
        }
    }

    // Use helper function to create input list 
    // Create AK4 inputs object
    AK4Inputs_ = new ttUtility::ConstAK4Inputs<double>(
        Jets_,
        Jets_bJetTagDeepCSVtotb_,
        Jets_qgLikelihood_, 
        hadtops_, 
        hadtopdaughters_);  
    AK4Inputs_->setFilterVector(*ak4Filter_); // filter

    // Create AK8 inputs object
    AK8Inputs_ = new ttUtility::ConstAK8Inputs<double>(
        JetsAK8_,
        JetsAK8_DeepTagTvsQCD_,
        JetsAK8_DeepTagWvsQCD_,
        JetsAK8_softDropMass_,
        JetsAK8_subjets_,        
        hadtops_,
        hadtopdaughters_);  
    
    //Add variables that are not passed to the constructor
    addVariables();      
}

SetUpTopTagger::~SetUpTopTagger()
{
    delete AK4Inputs_;
    delete AK8Inputs_;
    delete ak4Filter_;
}

std::vector<double>* SetUpTopTagger::intVecTodoubleVec(NTupleReader& tr, const std::string& name)
{
    const auto& vI = tr.getVec<int>(name);
    std::vector<double>* vD = new std::vector<double>(vI.size());
    for(unsigned int i = 0; i < vI.size(); i++)
    {
        (*vD)[i] = vI[i];
    }
    tr.registerDerivedVec(name+"ConvertedToDouble"+myVarSuffix_, vD);
    return vD;
}

std::vector<std::vector<double>>* SetUpTopTagger::VecVecintToVecVecdouble(NTupleReader& tr, const std::string& name)
{
    const auto& vvI = tr.getVec<std::vector<int>>(name);
    std::vector<std::vector<double>>* vvD = new std::vector<std::vector<double>>();
    for(unsigned int i = 0; i < vvI.size(); i++)
    {
        std::vector<double> vD;
        for(unsigned int j = 0; j < vvI[i].size(); j++)
        {
            vD.push_back(vvI[i][j]);
        }
        vvD->push_back(vD);
    }
    tr.registerDerivedVec(name+"ConvertedToDouble"+myVarSuffix_, vvD);
    return vvD;
}

void SetUpTopTagger::addVariables()
{
    AK4Inputs_->addSupplamentalVector("qgPtD",                               tr_.getVec<double>("Jets"+myVarSuffix_+"_ptD"));
    AK4Inputs_->addSupplamentalVector("qgAxis1",                             tr_.getVec<double>("Jets"+myVarSuffix_+"_axismajor"));
    AK4Inputs_->addSupplamentalVector("qgAxis2",                             tr_.getVec<double>("Jets"+myVarSuffix_+"_axisminor"));
    AK4Inputs_->addSupplamentalVector("qgMult",                              *intVecTodoubleVec(tr_,"Jets"+myVarSuffix_+"_multiplicity"));
    AK4Inputs_->addSupplamentalVector("recoJetschargedHadronEnergyFraction", tr_.getVec<double>("Jets"+myVarSuffix_+"_chargedHadronEnergyFraction"));
    AK4Inputs_->addSupplamentalVector("recoJetschargedEmEnergyFraction",     tr_.getVec<double>("Jets"+myVarSuffix_+"_chargedEmEnergyFraction"));
    AK4Inputs_->addSupplamentalVector("recoJetsneutralEmEnergyFraction",     tr_.getVec<double>("Jets"+myVarSuffix_+"_neutralEmEnergyFraction"));
    AK4Inputs_->addSupplamentalVector("recoJetsmuonEnergyFraction",          tr_.getVec<double>("Jets"+myVarSuffix_+"_muonEnergyFraction"));
    AK4Inputs_->addSupplamentalVector("recoJetsHFHadronEnergyFraction",      tr_.getVec<double>("Jets"+myVarSuffix_+"_hfHadronEnergyFraction"));
    AK4Inputs_->addSupplamentalVector("recoJetsHFEMEnergyFraction",          tr_.getVec<double>("Jets"+myVarSuffix_+"_hfEMEnergyFraction"));
    AK4Inputs_->addSupplamentalVector("recoJetsneutralEnergyFraction",       *add2Vec<double>(tr_,"Jets"+myVarSuffix_+"_neutralEmEnergyFraction","Jets"+myVarSuffix_+"_neutralHadronEnergyFraction"));
    AK4Inputs_->addSupplamentalVector("PhotonEnergyFraction",                tr_.getVec<double>("Jets"+myVarSuffix_+"_photonEnergyFraction"));
    AK4Inputs_->addSupplamentalVector("ElectronEnergyFraction",              tr_.getVec<double>("Jets"+myVarSuffix_+"_electronEnergyFraction"));
    AK4Inputs_->addSupplamentalVector("ChargedHadronMultiplicity",           *intVecTodoubleVec(tr_,"Jets"+myVarSuffix_+"_chargedHadronMultiplicity"));
    AK4Inputs_->addSupplamentalVector("NeutralHadronMultiplicity",           *intVecTodoubleVec(tr_,"Jets"+myVarSuffix_+"_neutralHadronMultiplicity"));
    AK4Inputs_->addSupplamentalVector("PhotonMultiplicity",                  *intVecTodoubleVec(tr_,"Jets"+myVarSuffix_+"_photonMultiplicity"));
    AK4Inputs_->addSupplamentalVector("ElectronMultiplicity",                *intVecTodoubleVec(tr_,"Jets"+myVarSuffix_+"_electronMultiplicity"));
    AK4Inputs_->addSupplamentalVector("MuonMultiplicity",                    *intVecTodoubleVec(tr_,"Jets"+myVarSuffix_+"_muonMultiplicity"));
    AK4Inputs_->addSupplamentalVector("DeepCSVb",                            tr_.getVec<double>("Jets"+myVarSuffix_+"_bJetTagDeepCSVprobb"));
    AK4Inputs_->addSupplamentalVector("DeepCSVc",                            tr_.getVec<double>("Jets"+myVarSuffix_+"_bJetTagDeepCSVprobc"));
    AK4Inputs_->addSupplamentalVector("DeepCSVl",                            tr_.getVec<double>("Jets"+myVarSuffix_+"_bJetTagDeepCSVprobudsg"));
    AK4Inputs_->addSupplamentalVector("DeepCSVbb",                           tr_.getVec<double>("Jets"+myVarSuffix_+"_bJetTagDeepCSVprobbb"));
}

// Create jets constituents list combining AK4 and AK8 jets, these are used to construct top candiates
// The vector of input constituents can also be constructed "by hand"    
std::vector<Constituent> SetUpTopTagger::getConstituents() const
{
    return ttUtility::packageConstituents(*AK4Inputs_, *AK8Inputs_);
}
