#ifndef TRAININGNTUPLEVARS_H
#define TRAININGNTUPLEVARS_H

class TrainingNTupleVars
{
private:
    std::string myVarSuffix_;

    
    void moreVars(NTupleReader& tr)
    {
        const auto& Jets                    = tr.getVec<utility::LorentzVector>(("Jets"+myVarSuffix_));
        const auto& GoodJets_pt30           = tr.getVec<bool>("GoodJets_pt30"+myVarSuffix_);
        const auto& GoodLeptons             = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons"+myVarSuffix_);
        
        auto& GoodLeptonsMass               = tr.createDerivedVec<double>("GoodLeptonsMass"+myVarSuffix_);
        auto& GoodLeptonsPt                 = tr.createDerivedVec<double>("GoodLeptonsPt"+myVarSuffix_);
        auto& GoodLeptonsEta                = tr.createDerivedVec<double>("GoodLeptonsEta"+myVarSuffix_);
        auto& GoodLeptonsPhi                = tr.createDerivedVec<double>("GoodLeptonsPhi"+myVarSuffix_);
        auto& GoodLeptonsPx                 = tr.createDerivedVec<double>("GoodLeptonsPx"+myVarSuffix_);
        auto& GoodLeptonsPy                 = tr.createDerivedVec<double>("GoodLeptonsPy"+myVarSuffix_);
        auto& GoodLeptonsPz                 = tr.createDerivedVec<double>("GoodLeptonsPz"+myVarSuffix_);
        auto& GoodLeptonsE                  = tr.createDerivedVec<double>("GoodLeptonsE"+myVarSuffix_);

        for (const auto& lep : GoodLeptons)
        {
            GoodLeptonsMass.push_back(lep.second.M());
            GoodLeptonsPt.push_back(lep.second.Pt());
            GoodLeptonsEta.push_back(lep.second.Eta());
            GoodLeptonsPhi.push_back(lep.second.Phi());
            GoodLeptonsPx.push_back(lep.second.Px());
            GoodLeptonsPy.push_back(lep.second.Py());
            GoodLeptonsPz.push_back(lep.second.Pz());
            GoodLeptonsE.push_back(lep.second.E());
        }

        auto& GoodJetsMass                    = tr.createDerivedVec<double>("GoodJetsMass"+myVarSuffix_);
        auto& GoodJetsPt                      = tr.createDerivedVec<double>("GoodJetsPt"+myVarSuffix_);
        auto& GoodJetsEta                     = tr.createDerivedVec<double>("GoodJetsEta"+myVarSuffix_);
        auto& GoodJetsPhi                     = tr.createDerivedVec<double>("GoodJetsPhi"+myVarSuffix_);
        auto& GoodJetsPx                      = tr.createDerivedVec<double>("GoodJetsPx"+myVarSuffix_);
        auto& GoodJetsPy                      = tr.createDerivedVec<double>("GoodJetsPy"+myVarSuffix_);
        auto& GoodJetsPz                      = tr.createDerivedVec<double>("GoodJetsPz"+myVarSuffix_);
        auto& GoodJetsE                       = tr.createDerivedVec<double>("GoodJetsE"+myVarSuffix_);
        for (unsigned int j=0; j<Jets.size(); j++)
        {
            if (GoodJets_pt30[j])
            {
                GoodJetsMass.push_back(Jets[j].M());
                GoodJetsPt.push_back(Jets[j].Pt());
                GoodJetsEta.push_back(Jets[j].Eta());
                GoodJetsPhi.push_back(Jets[j].Phi());
                GoodJetsPx.push_back(Jets[j].Px());
                GoodJetsPy.push_back(Jets[j].Py());
                GoodJetsPz.push_back(Jets[j].Pz());
                GoodJetsE.push_back(Jets[j].E());
            }
        }
    }

public:
    TrainingNTupleVars(std::string myVarSuffix = "") 
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up TrainingNTupleVars"<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        moreVars(tr);
    }
};

#endif
