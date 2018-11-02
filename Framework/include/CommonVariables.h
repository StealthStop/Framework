#ifndef COMMONVARIABLES_H
#define COMMONVARIABLES_H

class CommonVariables
{
private:
    std::string myVarSuffix_;

    void genMatch(NTupleReader& tr) const
    {
        const auto& runtype = tr.getVar<std::string>("runtype");
        int used_gen_bjet_for_mbl = -1;
        
        if(runtype != "Data")
        {
            const auto& GenParticles = tr.getVec<TLorentzVector>("GenParticles");
            const auto& Mbl_Index    = tr.getVar<int>("used_bjet_for_mbl"+myVarSuffix_);
            const auto& Jets         = tr.getVec<TLorentzVector>(("Jets"+myVarSuffix_));

            for(unsigned int gpi=0; gpi < GenParticles.size(); gpi++ ) 
            {
                if(Mbl_Index < 0) continue;
                double deltaR = Jets.at(Mbl_Index).DeltaR( GenParticles.at(gpi) );
                if(deltaR < 0.4) used_gen_bjet_for_mbl = gpi;
            }            
        }
        tr.registerDerivedVar("used_gen_bjet_for_mbl"+myVarSuffix_, used_gen_bjet_for_mbl);
    }

    void commonVariables(NTupleReader& tr)
    {
        // Get needed branches
        const auto& Jets = tr.getVec<TLorentzVector>(("Jets"+myVarSuffix_));
        const auto& GoodJets = tr.getVec<bool>("GoodJets"+myVarSuffix_);
        const auto& GoodBJets_pt30 = tr.getVec<bool>("GoodBJets_pt30"+myVarSuffix_);

        const auto& Muons = tr.getVec<TLorentzVector>("Muons");
        const auto& MuonsCharge = tr.getVec<int>("Muons_charge");
        const auto& MuonsMTW = tr.getVec<double>("MuonsMTW");
        const auto& GoodMuons = tr.getVec<bool>("GoodMuons");
        const auto& NGoodMuons = tr.getVar<int>("NGoodMuons");

        const auto& Electrons = tr.getVec<TLorentzVector>("Electrons");
        const auto& ElectronsCharge = tr.getVec<int>("Electrons_charge");
        const auto& ElectronsMTW = tr.getVec<double>("ElectronsMTW");
        const auto& GoodElectrons = tr.getVec<bool>("GoodElectrons");
        const auto& NGoodElectrons = tr.getVar<int>("NGoodElectrons");

        const auto& etaCut = tr.getVar<double>("etaCut");

        // HT of jets with pT>40
        double ht = 0;
        for(unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
        {            
            if(!GoodJets[ijet]) continue;
            TLorentzVector jet = Jets.at(ijet);

            if(jet.Pt() > 40 && abs(jet.Eta()) < etaCut)
                ht += jet.Pt();
        }
        tr.registerDerivedVar("HT_trigger"+myVarSuffix_, ht);

        // Put leptons together
        auto* GoodLeptons = new std::vector<std::pair<std::string, TLorentzVector>>();
        auto* GoodLeptonsCharge = new std::vector<int>();
        int NGoodLeptons = 0;
        for(unsigned int imu = 0; imu < Muons.size(); ++imu)
        {            
            if(!GoodMuons[imu]) continue;
            TLorentzVector muon = Muons.at(imu);
            GoodLeptons->push_back( std::make_pair("m", muon) );
            GoodLeptonsCharge->push_back( MuonsCharge.at(imu) );
            NGoodLeptons++;
        }
        for(unsigned int iel = 0; iel < Electrons.size(); ++iel)
        {
            if(!GoodElectrons[iel]) continue;
            TLorentzVector electron = Electrons.at(iel);
            GoodLeptons->push_back( std::make_pair("e", electron) );
            GoodLeptonsCharge->push_back( ElectronsCharge.at(iel) );
            NGoodLeptons++;
        }

        tr.registerDerivedVec("GoodLeptons"+myVarSuffix_, GoodLeptons);
        tr.registerDerivedVar("NGoodLeptons"+myVarSuffix_, NGoodLeptons);
        tr.registerDerivedVec("GoodLeptonsCharge"+myVarSuffix_, GoodLeptonsCharge);
        
        // M(l,b); closest to 105 GeV if multiple combinations (halfway between 30 and 180 GeV)
        double Mbl = 0;
        double Mbldiff = 999.;
        int used_bjet_for_mbl = -1;
        for(const auto& pair : *GoodLeptons)
        {
            TLorentzVector lepton = pair.second;
            for(unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
            {            
                if(!GoodBJets_pt30[ijet]) continue;
                TLorentzVector bjet = Jets.at(ijet);                
                double mbl = (lepton+bjet).M();
                if( abs(mbl - 105) < Mbldiff)
                {
                    Mbl = mbl;
                    Mbldiff = abs(mbl - 105);
                    used_bjet_for_mbl = ijet;
                }
            }
        }
        tr.registerDerivedVar("Mbl"+myVarSuffix_,Mbl);
        tr.registerDerivedVar("used_bjet_for_mbl"+myVarSuffix_,used_bjet_for_mbl);
        genMatch(tr);
        
        //Find single lepton for HistoContainer
        TLorentzVector singleLepton;
        double mTLep = 999.9;
        for(unsigned int i = 0; i < Muons.size(); ++i)
        {
            if(!GoodMuons[i]) continue;
            if(Muons[i].Pt() > 20)
            {
                singleLepton = Muons[i];
                mTLep = MuonsMTW[i];
                break;
            }
        }
        for(unsigned int i = 0; i < Electrons.size(); ++i)
        {
            if(!GoodElectrons[i]) continue;
            if(Electrons[i].Pt() > 20)
            {
                singleLepton = Electrons[i];
                mTLep = ElectronsMTW[i];
                break;
            }
        }
        tr.registerDerivedVar("singleLepton"+myVarSuffix_, singleLepton);            

        // 2 lepton onZ selection variables
        bool onZ = false;
        double mll = 0;
        if ( GoodLeptons->size() == 2 )
        {
            if( (NGoodMuons == 2 || NGoodElectrons == 2) && (GoodLeptonsCharge->at(0) != GoodLeptonsCharge->at(1)) )
            {
                mll = ( GoodLeptons->at(0).second + GoodLeptons->at(1).second ).M();
                if( mll > 81 && mll < 101)
                    onZ = true; 
            }
        }
        tr.registerDerivedVar("onZ"+myVarSuffix_, onZ);
        tr.registerDerivedVar("mll"+myVarSuffix_, mll);        
    }

public:
    CommonVariables(std::string myVarSuffix = "")
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up CommonVariables"<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        commonVariables(tr);
    }
};

#endif
