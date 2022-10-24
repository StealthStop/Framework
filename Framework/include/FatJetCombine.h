#ifndef FATJETCOMBINE_H
#define FATJETCOMBINE_H

#include "Framework/Framework/include/Utility.h"

class FatJetCombine
{
private:
    std::string myVarSuffix_;

    void fatjetcombine(NTupleReader& tr)
    {
        const auto& JetsAK8               = tr.getVec<utility::LorentzVector>("JetsAK8"+myVarSuffix_);
        const auto& Muons                 = tr.getVec<utility::LorentzVector>("Muons");
        const auto& Electrons             = tr.getVec<utility::LorentzVector>("Electrons");

        //First clean out leptons from JetsAK8 collection
        auto& GoodJetsAK8                = tr.createDerivedVec<bool>("GoodJetsAK8"+myVarSuffix_, JetsAK8.size(), true);
        
        for(unsigned int mu=0; mu < Muons.size(); mu++)
        {
            double minDeltaR = 0.8;
            utility::LorentzVector myMuon = Muons.at(mu);
            int muonCand = -1;
            for(unsigned int j=0; j < JetsAK8.size(); j++)
            {
                utility::LorentzVector myJet = JetsAK8.at(j);
                if( std::fabs(myMuon.Pt() - myJet.Pt()) / myMuon.Pt() < 1 && utility::DeltaR(myMuon, myJet) < minDeltaR)
                {
                    minDeltaR = utility::DeltaR(myMuon, myJet);
                    muonCand = j;
                }
            }
            if(muonCand != -1) GoodJetsAK8.at(muonCand) = false;
        }

        for(unsigned int el=0; el < Electrons.size(); el++)
        {
            double minDeltaR = 0.8;
            utility::LorentzVector myElec = Electrons.at(el);
            int elecCand = -1;
            for(unsigned int j=0; j < JetsAK8.size(); j++)
            {
                utility::LorentzVector myJet = JetsAK8.at(j);
                if( std::fabs(myElec.Pt() - myJet.Pt()) / myElec.Pt() < 1 && utility::DeltaR(myElec, myJet) < minDeltaR)
                {
                    minDeltaR = utility::DeltaR(myElec, myJet);
                    elecCand = j;                  
                }
            }
            if(elecCand != -1) GoodJetsAK8.at(elecCand) = false;
        }

        int NGoodJetsAK8 = 0;
        for (unsigned int j = 0; j < JetsAK8.size(); j++)
        {
            if(abs(JetsAK8.at(j).Eta()) > 2.4 or JetsAK8.at(j).Pt() < 170) 
                GoodJetsAK8.at(j) = false;

            else if(GoodJetsAK8.at(j))
                NGoodJetsAK8++;
        }                

        tr.registerDerivedVar("NGoodJetsAK8"+myVarSuffix_, NGoodJetsAK8);
        }

public:
    FatJetCombine(std::string myVarSuffix = "") 
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up FatJetCombine"<<std::endl;   
    }
    
    void operator()(NTupleReader& tr)
    {
        fatjetcombine(tr);
    }
};

#endif
