#ifndef COMMONVARIABLES_H
#define COMMONVARIABLES_H

class CommonVariables
{
private:
    std::string myVarSuffix_;

    bool objectInHEM(const std::vector<utility::LorentzVector>& objects, const double etalow, const double etahigh, const double philow, const double phihigh, const double ptcut) const
    {
        bool inHEM = false;
        for(const auto& o : objects)
        {
            if( (o.Eta() >= etalow) && (o.Eta() <= etahigh) && (o.Phi() >= philow) && (o.Phi() <= phihigh) && (o.Pt() > ptcut) )
            {
                inHEM = true;
            }
        }
        return inHEM;
    }

    bool objectInHEM(const std::vector<utility::LorentzVector>& objects, const std::vector<bool>& filter, const double etalow, const double etahigh, const double philow, const double phihigh, const double ptcut) const
    {
        bool inHEM = false;
        for(unsigned int i = 0; i < objects.size(); i++)
        {
            if( (objects[i].Eta() >= etalow) && (objects[i].Eta() <= etahigh) && (objects[i].Phi() >= philow) && (objects[i].Phi() <= phihigh) && (objects[i].Pt() > ptcut) && filter[i])
            {
                inHEM = true;
            }
        }
        return inHEM;
    }

    void commonVariables(NTupleReader& tr)
    {
        // Get needed branches
        const auto& runYear = tr.getVar<std::string>("runYear");
        const auto& filetag = tr.getVar<std::string>("filetag");
        const auto& runType = tr.getVar<std::string>("runtype");
        const auto& Jets = tr.getVec<utility::LorentzVector>(("Jets"+myVarSuffix_));
        const auto& GoodJets = tr.getVec<bool>("GoodJets"+myVarSuffix_);
        const auto& GoodBJets = tr.getVec<bool>("GoodBJets"+myVarSuffix_);
        const auto& GoodBJets_pt30 = tr.getVec<bool>("GoodBJets_pt30"+myVarSuffix_);
        const auto& GoodBJets_pt45 = tr.getVec<bool>("GoodBJets_pt45"+myVarSuffix_);
        const auto& NonIsoMuonJets_pt30 = tr.getVec<bool>("NonIsoMuonJets_pt30"+myVarSuffix_);
        const auto& NonIsoMuonJets_pt45 = tr.getVec<bool>("NonIsoMuonJets_pt45"+myVarSuffix_);
        const auto& Muons = tr.getVec<utility::LorentzVector>("Muons");
        const auto& MuonsCharge = tr.getVec<int>("Muons_charge");
        const auto& MuonsMiniIso = tr.getVec<float>("Muons_iso");
        const auto& GoodMuons = tr.getVec<bool>("GoodMuons"+myVarSuffix_);
        const auto& NonIsoMuons = tr.getVec<bool>("NonIsoMuons"+myVarSuffix_);
        const auto& NGoodMuons = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        const auto& GoodMuons_pt20 = tr.getVec<bool>("GoodMuons_pt20"+myVarSuffix_);
        const auto& Electrons = tr.getVec<utility::LorentzVector>("Electrons");
        const auto& ElectronsCharge = tr.getVec<int>("Electrons_charge");
        const auto& ElectronsMiniIso = tr.getVec<float>("Electrons_iso");
        const auto& GoodElectrons = tr.getVec<bool>("GoodElectrons"+myVarSuffix_);
        const auto& NGoodElectrons = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);
        const auto& GoodElectrons_pt20 = tr.getVec<bool>("GoodElectrons_pt20"+myVarSuffix_);
        const auto& etaCut = tr.getVar<double>("etaCut");
        const auto& NGoodJets_pt30       = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_);
        const auto& NGoodBJets_pt30 = tr.getVar<int>("NGoodBJets_pt30"+myVarSuffix_);
        const auto& NGoodBJets_pt45 = tr.getVar<int>("NGoodBJets_pt45"+myVarSuffix_);
        const auto& NNonIsoMuonJets_pt30 = tr.getVar<int>("NNonIsoMuonJets_pt30"+myVarSuffix_);

        // Define electron HEM15/16 veto
        bool vetoedHEMelectron = objectInHEM(Electrons, -3.00, -1.30, -1.57, -0.87, 20.0);
        tr.registerDerivedVar("vetoedHEMelectron"+myVarSuffix_, vetoedHEMelectron);

        // Define muon HEM15/16 veto
        bool vetoedHEMmuon = objectInHEM(Muons, -3.00, -1.30, -1.57, -0.87, 20.0);
        tr.registerDerivedVar("vetoedHEMmuon"+myVarSuffix_, vetoedHEMmuon);

        // Define jet HEM15/16 veto
        bool vetoedHEMjet = objectInHEM(Jets, -3.20, -1.10, -1.77, -0.67, 20.0);
        tr.registerDerivedVar("vetoedHEMjet"+myVarSuffix_, vetoedHEMjet);

        // Define b jet HEM15/16 veto
        bool vetoedHEMbjet = objectInHEM(Jets, GoodBJets, -3.20, -1.10, -1.77, -0.67, 20.0);
        tr.registerDerivedVar("vetoedHEMbjet"+myVarSuffix_, vetoedHEMbjet);

        // Define lepton HEM15/16 veto
        bool vetoedHEMlepton = vetoedHEMelectron || vetoedHEMmuon;
        tr.registerDerivedVar("vetoedHEMlepton"+myVarSuffix_, vetoedHEMlepton);

        // Create the eventCounter variable to keep track of processed events
        int w = 1;
        bool passMadHT = true;
        if(runType == "MC")
        {
            // For MC, a "Weight" branch is provided by TreeMaker in the ntuples
            const float Weight = tr.getVar<float>("Weight");
            w = (Weight >= 0.0) ? 1 : -1;

            // "weightAbsVal" calculated by samples.cc provides no sign information
            // So determine that here using the weight from TreeMaker, where
            // the value is irrelevant
            const auto& weightAbsVal = tr.getVar<double>("weightAbsVal");

            // Reregister "Weight" with the newly calculated weight coming
            // from samples.cc and save the original Weight in a new "WeightTM" field
            tr.registerDerivedVar<float>("Weight",   w*weightAbsVal);
            tr.registerDerivedVar<float>("WeightTM", Weight);

            // Register lumi * xsec as its own variable for users
            // For 2018, if an event would be vetoed, use 2018 pre-HEM lumi in lumi * xsec
            // otherwise, use the nominal 2018 or respective year lumi
            const auto& Lumi = tr.getVar<double>("Lumi");
            double FinalLumi = Lumi;
            if (runYear == "2018" && vetoedHEMelectron)
            {
                const auto& Lumi_preHEM = tr.getVar<double>("Lumi_preHEM");
                tr.registerDerivedVar<double>("LumiXsec", w*weightAbsVal*Lumi_preHEM);

                FinalLumi = Lumi_preHEM;

            } else
            {
                tr.registerDerivedVar<double>("LumiXsec", w*weightAbsVal*Lumi);
            }
            tr.registerDerivedVar<double>("FinalLumi", FinalLumi);

            // --------------------
            // Check passMadHT here
            // --------------------
            const auto& madHT  = tr.getVar<float>("madHT");

            // Exclude events with MadGraph HT > 70 from the WJets and DY inclusive samples
            // in order to avoid double counting with the HT-binned samples
            if(filetag.find("DYJetsToLL_M-50_Incl") != std::string::npos && madHT > 70.0) passMadHT = false;
            if(filetag.find("WJetsToLNu_Incl")      != std::string::npos && madHT > 70.0) passMadHT = false;

            // Stitch TTbar samples together
            // remove HT overlap
            if((filetag.find("TTJets_Incl")               != std::string::npos || 
                filetag.find("TTJets_SingleLeptFromT")    != std::string::npos || 
                filetag.find("TTJets_SingleLeptFromTbar") != std::string::npos || 
                filetag.find("TTJets_DiLept")             != std::string::npos) && madHT > 600) 
            {
                passMadHT = false;
            }

            // also remove lepton overlap from the inclusive sample
            const auto& GenElectrons = tr.getVec<utility::LorentzVector>("GenElectrons");
            const auto& GenMuons     = tr.getVec<utility::LorentzVector>("GenMuons");
            const auto& GenTaus      = tr.getVec<utility::LorentzVector>("GenTaus");
            int NGenLeptons          = GenElectrons.size() + GenMuons.size() + GenTaus.size();

            if (filetag.find("TTJets_Incl") != std::string::npos && NGenLeptons > 0) passMadHT = false;
        }
        else if(runType == "Data")
        {
            // For Data, no "Weight" branch is provided
            // Put in trivial 1.0 for it and "WeightTM"
            // for use later on in the code
            tr.registerDerivedVar<float>("Weight",   1.0);
            tr.registerDerivedVar<float>("WeightTM", 1.0);
        }

        tr.registerDerivedVar<bool>("passMadHT" + myVarSuffix_, passMadHT);

        tr.registerDerivedVar<int>("eventCounter",w);        

        // HT of jets
        double ht = 0.0, ht_pt30 = 0.0, ht_pt45 = 0.0;         
        double ht_NonIsoMuon_pt30 = 0.0, ht_NonIsoMuon_pt45 = 0.0;
        for(unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
        {            
            double pT = Jets.at(ijet).Pt();
            double absEta = abs(Jets.at(ijet).Eta());

            if(GoodJets[ijet] && pT > 30 && absEta < etaCut) ht_pt30 += pT;
            if(GoodJets[ijet] && pT > 40 && absEta < etaCut)      ht += pT;
            if(GoodJets[ijet] && pT > 45 && absEta < etaCut) ht_pt45 += pT;
            
            if(NonIsoMuonJets_pt30[ijet] && pT > 30 && absEta < etaCut) ht_NonIsoMuon_pt30 += pT;
            if(NonIsoMuonJets_pt45[ijet] && pT > 45 && absEta < etaCut) ht_NonIsoMuon_pt45 += pT;
        }
        tr.registerDerivedVar("HT_trigger"+myVarSuffix_, ht);
        tr.registerDerivedVar("HT_trigger_pt30"+myVarSuffix_, ht_pt30);
        tr.registerDerivedVar("HT_trigger_pt45"+myVarSuffix_, ht_pt45);
        tr.registerDerivedVar("HT_NonIsoMuon_pt30"+myVarSuffix_, ht_NonIsoMuon_pt30);
        tr.registerDerivedVar("HT_NonIsoMuon_pt45"+myVarSuffix_, ht_NonIsoMuon_pt45);

        // Put leptons together
        auto* GoodLeptons = new std::vector<std::pair<std::string, utility::LorentzVector>>();
        auto* GoodLeptons_pt20 = new std::vector<std::pair<std::string, utility::LorentzVector>>();
        auto* GoodLeptonsCharge = new std::vector<int>();
        auto* GoodLeptonsCharge_pt20 = new std::vector<int>();
        auto* GoodLeptonsMiniIso = new std::vector<double>();
        auto* GoodLeptonsMiniIso_pt20 = new std::vector<double>();

        auto* GoodNonIsoMuonsCharge = new std::vector<int>();
        auto* GoodNonIsoMuonsMiniIso = new std::vector<double>();

        int NGoodLeptons = 0;
        int NGoodLeptons_pt20 = 0;
        for(unsigned int imu = 0; imu < Muons.size(); ++imu)
        {
            if(GoodMuons[imu])
            {
                utility::LorentzVector muon = Muons.at(imu);
                GoodLeptons->push_back( std::make_pair("m", muon) );
                GoodLeptonsCharge->push_back( MuonsCharge.at(imu) );
                GoodLeptonsMiniIso->push_back( MuonsMiniIso.at(imu) );
                NGoodLeptons++;
            }
            if(GoodMuons_pt20[imu])
            {
                utility::LorentzVector muon = Muons.at(imu);
                GoodLeptons_pt20->push_back( std::make_pair("m", muon) );
                GoodLeptonsCharge_pt20->push_back( MuonsCharge.at(imu) );
                GoodLeptonsMiniIso_pt20->push_back( MuonsMiniIso.at(imu) );
                NGoodLeptons_pt20++;
            }

            if(NonIsoMuons[imu])
            {
                GoodNonIsoMuonsCharge->push_back( MuonsCharge.at(imu) );
                GoodNonIsoMuonsMiniIso->push_back( MuonsMiniIso.at(imu) );
            }
        }
        for(unsigned int iel = 0; iel < Electrons.size(); ++iel)
        {
            if(GoodElectrons[iel])
            {
                utility::LorentzVector electron = Electrons.at(iel);
                GoodLeptons->push_back( std::make_pair("e", electron) );
                GoodLeptonsCharge->push_back( ElectronsCharge.at(iel) );
                GoodLeptonsMiniIso->push_back( ElectronsMiniIso.at(iel) );
                NGoodLeptons++;
            }
            if(GoodElectrons_pt20[iel])
            {
                utility::LorentzVector electron = Electrons.at(iel);
                GoodLeptons_pt20->push_back( std::make_pair("e", electron) );
                GoodLeptonsCharge_pt20->push_back( ElectronsCharge.at(iel) );
                GoodLeptonsMiniIso_pt20->push_back( ElectronsMiniIso.at(iel) );
                NGoodLeptons_pt20++;
            }
        }

        tr.registerDerivedVec("GoodLeptons"+myVarSuffix_, GoodLeptons);
        tr.registerDerivedVar("NGoodLeptons"+myVarSuffix_, NGoodLeptons);
        tr.registerDerivedVec("GoodLeptonsCharge"+myVarSuffix_, GoodLeptonsCharge);
        tr.registerDerivedVec("GoodLeptonsMiniIso"+myVarSuffix_, GoodLeptonsMiniIso);
        tr.registerDerivedVec("GoodNonIsoMuonsCharge"+myVarSuffix_, GoodNonIsoMuonsCharge);
        tr.registerDerivedVec("GoodNonIsoMuonsMiniIso"+myVarSuffix_, GoodNonIsoMuonsMiniIso);
        tr.registerDerivedVec("GoodLeptons_pt20"+myVarSuffix_, GoodLeptons_pt20);
        tr.registerDerivedVar("NGoodLeptons_pt20"+myVarSuffix_, NGoodLeptons_pt20);
        tr.registerDerivedVec("GoodLeptonsCharge_pt20"+myVarSuffix_, GoodLeptonsCharge_pt20);
        tr.registerDerivedVec("GoodLeptonsMiniIso_pt20"+myVarSuffix_, GoodLeptonsMiniIso_pt20);
       
        // M(l,b); closest to 105 GeV if multiple combinations (halfway between 30 and 180 GeV)
        double Mbl = 0;
        double Mbldiff = 999.;
        for(const auto& pair : *GoodLeptons)
        {
            utility::LorentzVector lepton = pair.second;
            for(unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
            {            
                if(!GoodBJets_pt30[ijet]) continue;
                utility::LorentzVector bjet = Jets.at(ijet);                
                double mbl = (lepton+bjet).M();
                if( abs(mbl - 105) < Mbldiff)
                {
                    Mbl = mbl;
                    Mbldiff = abs(mbl - 105);
                }
            }
        }
        tr.registerDerivedVar("Mbl"+myVarSuffix_,Mbl);
        
        //Find single lepton for HistoContainer
        TLorentzVector singleLepton;
        for(unsigned int i = 0; i < Muons.size(); ++i)
        {
            if(!GoodMuons[i]) continue;
            if(Muons[i].Pt() > 20)
            {
                singleLepton = utility::convertLV<TLorentzVector, utility::LorentzVector>(Muons[i]);
                break;
            }
        }
        for(unsigned int i = 0; i < Electrons.size(); ++i)
        {
            if(!GoodElectrons[i]) continue;
            if(Electrons[i].Pt() > 20)
            {
                singleLepton = utility::convertLV<TLorentzVector, utility::LorentzVector>(Electrons[i]);
                break;
            }
        }
        tr.registerDerivedVar("singleLepton"+myVarSuffix_, singleLepton);            

        // 2 lepton onZ selection variables
        // Default to clear non-physical value in case of same sign or opposite flavor leptons
        bool onZ = false;
        double mll = -999.0;
        if( GoodLeptons->size() == 2 )
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

        // ---------------------------------------------------
        // -- Calculate DeltaR between 2 bjets for 0 lepton
        // ---------------------------------------------------
        // calculate dR_bjets for old baseline selection
        double dR_bjets_old = -1;
        if(NGoodBJets_pt45 >= 2)
        {
            std::vector<utility::LorentzVector> bjets;
            for(unsigned int ijet = 0; ijet < Jets.size(); ijet++)
            {
                if(!GoodBJets_pt45[ijet]) continue;
                bjets.push_back(Jets.at(ijet));
            }
            int n = bjets.size();
            std::vector<double> deltaRs( n*(n - 1)/2, 0.0 );
            for(int i = 0; i < n; i++) 
            {
                for(int j = i+1; j < n; j++) 
                {
                    deltaRs[i+j-1] = utility::DeltaR(bjets[i], bjets[j]);
                }
            }
            dR_bjets_old = *std::max_element(deltaRs.begin(), deltaRs.end());
        }
        tr.registerDerivedVar("dR_bjets_old"+myVarSuffix_, dR_bjets_old);

        // calculate dR_bjets for new baseline selection
        double dR_bjets = -1;
        if(NGoodBJets_pt30 >= 2)
        {
            std::vector<utility::LorentzVector> bjets;
            for(unsigned int ijet = 0; ijet < Jets.size(); ijet++)
            {
                if(!GoodBJets_pt30[ijet]) continue;
                bjets.push_back(Jets.at(ijet));
            }
            int n = bjets.size();
            std::vector<double> deltaRs( n*(n - 1)/2, 0.0 );
            for(int i = 0; i < n; i++) 
            {
                for(int j = i+1; j < n; j++) 
                {
                    deltaRs[i+j-1] = utility::DeltaR(bjets[i], bjets[j]);
                }
            }
            dR_bjets = *std::max_element(deltaRs.begin(), deltaRs.end());
        }
        tr.registerDerivedVar("dR_bjets"+myVarSuffix_, dR_bjets);

        // An event must "earn" its keep by passing any of the main 0L, 1L, 2L baselines.
        // If an event does not pass any of these main selections, time-intensive
        // modules can be skipped in the pipeline as long as the relevant analyzer is also informed and also skips to the next event
        // This boolean only informs modules about an event, it does not make the final decision to skip an event
        bool lostCauseEvent = true;
                
        // Now combine all relevant baselines into a single bool
        // Here we select the main baselines for all three channels
        // If all four booleans are false, then the event is considered a "lost cause"
        lostCauseEvent &= !(NGoodJets_pt30       >= 6 and ht_pt30            > 500.0 and NGoodBJets_pt30 >= 1) and
                          !(NNonIsoMuonJets_pt30 >= 7 and ht_NonIsoMuon_pt30 > 500.0);

        tr.registerDerivedVar<bool>("lostCauseEvent" + myVarSuffix_, lostCauseEvent);
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
