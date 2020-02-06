#ifndef TRAININGNTUPLEVARS_H
#define TRAININGNTUPLEVARS_H

class TrainingNTupleVars
{
private:
    std::string myVarSuffix_;

    
    void moreVars(NTupleReader& tr)
    {
        const auto& Jets                    = tr.getVec<TLorentzVector>(("Jets"+myVarSuffix_));
        const auto& GoodJets_pt30           = tr.getVec<bool>("GoodJets_pt30"+myVarSuffix_);
        const auto& GoodLeptons             = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons"+myVarSuffix_);

        






        auto& GoodLeptonsMass                 = tr.createDerivedVec<double>("GoodLeptonsMass"+myVarSuffix_);
        auto& GoodLeptonsPt                   = tr.createDerivedVec<double>("GoodLeptonsPt"+myVarSuffix_);
        auto& GoodLeptonsEta                  = tr.createDerivedVec<double>("GoodLeptonsEta"+myVarSuffix_);
        auto& GoodLeptonsPhi                  = tr.createDerivedVec<double>("GoodLeptonsPhi"+myVarSuffix_);

        for (const auto& lep : GoodLeptons)
        {
            GoodLeptonsMass.push_back(lep.second.M());
            GoodLeptonsPt.push_back(lep.second.Pt());
            GoodLeptonsEta.push_back(lep.second.Eta());
            GoodLeptonsPhi.push_back(lep.second.Phi());
        }

        auto& GoodJetsMass                    = tr.createDerivedVec<double>("GoodJetsMass"+myVarSuffix_);
        auto& GoodJetsPt                      = tr.createDerivedVec<double>("GoodJetsPt"+myVarSuffix_);
        auto& GoodJetsEta                  = tr.createDerivedVec<double>("GoodJetsEta"+myVarSuffix_);
        auto& GoodJetsPhi                  = tr.createDerivedVec<double>("GoodJetsPhi"+myVarSuffix_);

        for (unsigned int j=0; j<Jets.size(); j++)
        {
            if (GoodJets_pt30[j])
            {
                GoodJetsMass.push_back(Jets[j].M());
                GoodJetsPt.push_back(Jets[j].Pt());
                GoodJetsEta.push_back(Jets[j].Eta());
                GoodJetsPhi.push_back(Jets[j].Phi());
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
