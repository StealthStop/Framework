#ifndef ELECTRON_H
#define ELECTRON_H

#include "Framework/Framework/include/Utility.h"

class Electron
{
private:
    std::string myVarSuffix_;

    void electron(NTupleReader& tr)
    {
        const auto& allElectrons = tr.getVec<utility::LorentzVector>("Electrons");
        const auto& allElectrons_passIso = tr.getVec<bool>("Electrons_passIso");
        const auto& allElectrons_charge  = tr.getVec<int>("Electrons_charge");
        const auto& allElectrons_tightID = tr.getVec<bool>("Electrons_tightID");
        const auto& etaCut = tr.getVar<double>("etaCut");
        const auto& MET = tr.getVar<float>("MET"); 
        const auto& METPhi = tr.getVar<float>("METPhi");
        const auto& runYear = tr.getVar<std::string>("runYear");

        utility::LorentzVector lvMET;
        lvMET.SetPt(MET); lvMET.SetEta(0.0); lvMET.SetPhi(METPhi); lvMET.SetE(MET);

        auto* good_electrons_ = new std::vector<bool>();
        auto* good_electrons_pt20 = new std::vector<bool>();
        auto* electrons_mtw_ = new std::vector<double>();
        int NGoodElectrons = 0;
        int NGoodPlusElectrons = 0;
        int NGoodMinusElectrons = 0;
        int NGoodElectrons_pt20 = 0;
        int NGoodPlusElectrons_pt20 = 0;
        int NGoodMinusElectrons_pt20 = 0;
        double ptCut = 0.0;
        if (runYear.find("2016") != std::string::npos)
            ptCut = 30.0;
        else
            ptCut = 37.0; 
        for(unsigned int iel = 0; iel < allElectrons.size(); ++iel)
        {
            utility::LorentzVector lvel = allElectrons.at(iel);
            double mtw = sqrt( 2*( lvMET.Pt()*lvel.Pt() - (lvMET.Px()*lvel.Px() + lvMET.Py()*lvel.Py()) ) );
            electrons_mtw_->push_back(mtw);
            if( abs(lvel.Eta()) < etaCut && 
                lvel.Pt() > ptCut && 
                allElectrons_passIso.at(iel) &&
                allElectrons_tightID.at(iel) 
                )
            {
                good_electrons_->push_back(true);
                NGoodElectrons++;
                if( allElectrons_charge.at(iel) ==  1 ) NGoodPlusElectrons++;
                else if( allElectrons_charge.at(iel) == -1 ) NGoodMinusElectrons++;
                else std::cout<<utility::color("Charge values in nTuples are different", "red")<<std::endl;
            }
            else
            {
                good_electrons_->push_back(false);
            }
            if ( abs(lvel.Eta()) < etaCut && 
                lvel.Pt() > 20 && 
                allElectrons_passIso.at(iel) &&
                allElectrons_tightID.at(iel) 
                )
            {
                good_electrons_pt20->push_back(true);
                NGoodElectrons_pt20++;
                if( allElectrons_charge.at(iel) ==  1 ) NGoodPlusElectrons_pt20++;
                else if( allElectrons_charge.at(iel) == -1 ) NGoodMinusElectrons_pt20++;
                else std::cout<<utility::color("Charge values in nTuples are different", "red")<<std::endl;
            }
            else
            {
                good_electrons_pt20->push_back(false);
            }
        }

        tr.registerDerivedVec("GoodElectrons"+myVarSuffix_, good_electrons_);
        tr.registerDerivedVar("NGoodElectrons"+myVarSuffix_, NGoodElectrons);
        tr.registerDerivedVar("NGoodPlusElectrons"+myVarSuffix_, NGoodPlusElectrons);
        tr.registerDerivedVar("NGoodMinusElectrons"+myVarSuffix_, NGoodMinusElectrons);
        tr.registerDerivedVec("GoodElectrons_pt20"+myVarSuffix_, good_electrons_pt20);
        tr.registerDerivedVar("NGoodElectrons_pt20"+myVarSuffix_, NGoodElectrons_pt20);
        tr.registerDerivedVar("NGoodPlusElectrons_pt20"+myVarSuffix_, NGoodPlusElectrons_pt20);
        tr.registerDerivedVar("NGoodMinusElectrons_pt20"+myVarSuffix_, NGoodMinusElectrons_pt20);

        tr.registerDerivedVec("ElectronsMTW"+myVarSuffix_, electrons_mtw_);
    }

public:
    Electron(std::string myVarSuffix = "") 
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up Electron"<<std::endl;   
    }

    void operator()(NTupleReader& tr)
    {
        electron(tr);
    }
};

#endif
