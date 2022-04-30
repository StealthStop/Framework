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
            const auto& GenParticles = tr.getVec<utility::LorentzVector>("GenParticles");
            const auto& Mbl_Index    = tr.getVar<int>("used_bjet_for_mbl"+myVarSuffix_);
            const auto& Jets         = tr.getVec<utility::LorentzVector>(("Jets"+myVarSuffix_));

            for(unsigned int gpi=0; gpi < GenParticles.size(); gpi++ ) 
            {
                if(Mbl_Index < 0) continue;
                double deltaR = utility::DeltaR(Jets.at(Mbl_Index), GenParticles.at(gpi) );
                if(deltaR < 0.4) used_gen_bjet_for_mbl = gpi;
            }            
        }
        tr.registerDerivedVar("used_gen_bjet_for_mbl"+myVarSuffix_, used_gen_bjet_for_mbl);
    }
    void getJetMatch(const std::vector<utility::LorentzVector>& JetVectors, const std::vector<bool>& GoodJets_pt30, const std::vector<bool>& GoodBJets_pt30, 
                     double& twoLep_mbldiff, int& used_jet, unsigned int other_used_jet, double& TwoLep_Mbl, const utility::LorentzVector& lep,  std::pair<int, int>& TwoLep_Mbl_Idx, const bool checkB)
    {
        for(unsigned int j=0; j < JetVectors.size(); j++)
        {
            bool pass_b = checkB ? GoodBJets_pt30.at(j):GoodJets_pt30.at(j);
            if(pass_b && j != other_used_jet)
            {
                double twoLep_mbl = (lep + JetVectors.at(j)).M();
                if (abs(twoLep_mbl - 105) < twoLep_mbldiff)
                {
                    TwoLep_Mbl = twoLep_mbl;
                    twoLep_mbldiff = abs(twoLep_mbl - 105);
                    used_jet = j;
                    TwoLep_Mbl_Idx.first = j;
                }                                           
            }
        }
    }
    
    bool objectInHEM(const std::vector<utility::LorentzVector>& objects, const double etalow, const double etahigh, const double philow, const double phihigh, const double ptcut, const std::string& runYear) const
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

    bool objectInHEM(const std::vector<utility::LorentzVector>& objects, const std::vector<bool>& filter, const double etalow, const double etahigh, const double philow, const double phihigh, const double ptcut, const std::string& runYear) const
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
        const auto& runType = tr.getVar<std::string>("runtype");
        const auto& runNumber = tr.getVar<unsigned int>("RunNum");
        const auto& Jets = tr.getVec<utility::LorentzVector>(("Jets"+myVarSuffix_));
        const auto& GoodJets = tr.getVec<bool>("GoodJets"+myVarSuffix_);
        const auto& GoodJets_pt30 = tr.getVec<bool>("GoodJets_pt30"+myVarSuffix_);
        const auto& GoodBJets = tr.getVec<bool>("GoodBJets"+myVarSuffix_);
        const auto& GoodBJets_pt30 = tr.getVec<bool>("GoodBJets_pt30"+myVarSuffix_);
        const auto& GoodBJets_pt45 = tr.getVec<bool>("GoodBJets_pt45"+myVarSuffix_);
        const auto& NonIsoMuonJets_pt30 = tr.getVec<bool>("NonIsoMuonJets_pt30"+myVarSuffix_);
        const auto& NonIsoMuonJets_pt45 = tr.getVec<bool>("NonIsoMuonJets_pt45"+myVarSuffix_);
        const auto& Muons = tr.getVec<utility::LorentzVector>("Muons");
        const auto& MuonsCharge = tr.getVec<int>("Muons_charge");
        const auto& MuonsMiniIso = tr.getVec<float>("Muons_iso");
        const auto& GoodMuons = tr.getVec<bool>("GoodMuons"+myVarSuffix_);
        const auto& NGoodMuons = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        const auto& GoodMuons_pt20 = tr.getVec<bool>("GoodMuons_pt20"+myVarSuffix_);
        const auto& Electrons = tr.getVec<utility::LorentzVector>("Electrons");
        const auto& ElectronsCharge = tr.getVec<int>("Electrons_charge");
        const auto& ElectronsMiniIso = tr.getVec<float>("Electrons_iso");
        const auto& GoodElectrons = tr.getVec<bool>("GoodElectrons"+myVarSuffix_);
        const auto& NGoodElectrons = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);
        const auto& GoodElectrons_pt20 = tr.getVec<bool>("GoodElectrons_pt20"+myVarSuffix_);
        const auto& etaCut = tr.getVar<double>("etaCut");
        const auto& NGoodBJets_pt30 = tr.getVar<int>("NGoodBJets_pt30"+myVarSuffix_);
        const auto& NGoodBJets_pt45 = tr.getVar<int>("NGoodBJets_pt45"+myVarSuffix_);
       
        // Define electron HEM15/16 veto
        bool vetoedHEMelectron = objectInHEM(Electrons, -3.00, -1.30, -1.57, -0.87, 20.0, runYear);
        tr.registerDerivedVar("vetoedHEMelectron"+myVarSuffix_, vetoedHEMelectron);

        // Define muon HEM15/16 veto
        bool vetoedHEMmuon = objectInHEM(Muons, -3.00, -1.30, -1.57, -0.87, 20.0, runYear);
        tr.registerDerivedVar("vetoedHEMmuon"+myVarSuffix_, vetoedHEMmuon);

        // Define jet HEM15/16 veto
        bool vetoedHEMjet = objectInHEM(Jets, -3.20, -1.10, -1.77, -0.67, 20.0, runYear);
        tr.registerDerivedVar("vetoedHEMjet"+myVarSuffix_, vetoedHEMjet);

        // Define b jet HEM15/16 veto
        bool vetoedHEMbjet = objectInHEM(Jets, GoodBJets, -3.20, -1.10, -1.77, -0.67, 20.0, runYear);
        tr.registerDerivedVar("vetoedHEMbjet"+myVarSuffix_, vetoedHEMbjet);

        // Define lepton HEM15/16 veto
        bool vetoedHEMlepton = vetoedHEMelectron || vetoedHEMmuon;
        tr.registerDerivedVar("vetoedHEMlepton"+myVarSuffix_, vetoedHEMlepton);

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
        tr.registerDerivedVec("GoodLeptons_pt20"+myVarSuffix_, GoodLeptons_pt20);
        tr.registerDerivedVar("NGoodLeptons_pt20"+myVarSuffix_, NGoodLeptons_pt20);
        tr.registerDerivedVec("GoodLeptonsCharge_pt20"+myVarSuffix_, GoodLeptonsCharge_pt20);
        tr.registerDerivedVec("GoodLeptonsMiniIso_pt20"+myVarSuffix_, GoodLeptonsMiniIso_pt20);
       
        // M(l,b); closest to 105 GeV if multiple combinations (halfway between 30 and 180 GeV)
        double Mbl = 0;
        auto*  MblVec = new std::vector<double>();
        double Mbldiff = 999.;
        int used_bjet_for_mbl = -1;
        for(const auto& pair : *GoodLeptons)
        {
            utility::LorentzVector lepton = pair.second;
            for(unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
            {            
                if(!GoodBJets_pt30[ijet]) continue;
                utility::LorentzVector bjet = Jets.at(ijet);                
                double mbl = (lepton+bjet).M();
                MblVec->push_back(mbl);
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
        tr.registerDerivedVec("MblVec"+myVarSuffix_,MblVec);
        genMatch(tr);
        
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
        bool onZ = false;
        double mll = 0;
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


        // Two Lepton Mbl definiton
        utility::LorentzVector lep1, lep2;
        int  used_jet1 = -1, used_jet2 =-1;
        double TwoLep_Mbl1=-1, TwoLep_Mbl2=-1;
        double twoLep_mbl1diff=999, twoLep_mbl2diff=999;
        std::pair<int, int> TwoLep_Mbl1_Idx(-1, -1), TwoLep_Mbl2_Idx(-1, -1);
       
        if(NGoodLeptons_pt20==2)
        {
            std::vector<std::pair<utility::LorentzVector,int>> GoodBJetsVec_pt30;

            if(GoodLeptons_pt20->at(0).second.Pt() >= GoodLeptons_pt20->at(1).second.Pt())
            {
                lep1 = GoodLeptons_pt20->at(0).second;
                lep2 = GoodLeptons_pt20->at(1).second;
                TwoLep_Mbl1_Idx.second = 0;
                TwoLep_Mbl2_Idx.second = 1;
            }
            else
            {
                lep1 = GoodLeptons_pt20->at(1).second;
                lep2 = GoodLeptons_pt20->at(0).second;
                TwoLep_Mbl1_Idx.second  = 1;
                TwoLep_Mbl2_Idx.second  = 0;
            }
            for(unsigned int j=0 ; j < Jets.size(); j++)
            {
                if(GoodBJets_pt30.at(j)) GoodBJetsVec_pt30.push_back(std::make_pair(Jets.at(j), j));
            }
            if(NGoodBJets_pt30 >= 2) //need three cases: N bjets >= 2, Nbjets == 1, Nbjets = 0. for Nbjets >=2 matches each lepton with best bjet
            {
                getJetMatch(Jets, GoodJets_pt30, GoodBJets_pt30, twoLep_mbl1diff, used_jet1, -1,  TwoLep_Mbl1, lep1, TwoLep_Mbl1_Idx, true);
                getJetMatch(Jets, GoodJets_pt30, GoodBJets_pt30, twoLep_mbl2diff, used_jet2, -1,  TwoLep_Mbl2, lep2, TwoLep_Mbl2_Idx, true);
                if(used_jet1 == used_jet2 && twoLep_mbl1diff > twoLep_mbl2diff) //if leptons are matched to same bjet, keep better one and rematch the other lepton
                {
                    twoLep_mbl1diff = 999;
                    getJetMatch(Jets, GoodJets_pt30, GoodBJets_pt30, twoLep_mbl1diff, used_jet1, used_jet2, TwoLep_Mbl1, lep1, TwoLep_Mbl1_Idx, true);
                }
                else if(used_jet1 == used_jet2 && twoLep_mbl1diff <= twoLep_mbl2diff) 
                {
                    twoLep_mbl2diff = 999;
                    getJetMatch(Jets, GoodJets_pt30, GoodBJets_pt30, twoLep_mbl2diff, used_jet2, used_jet1, TwoLep_Mbl2, lep2, TwoLep_Mbl2_Idx, true);
                }
            }
            else if(NGoodBJets_pt30 == 1) //in this case the better Mbl pair is made of the leptons, then the other lepton is paired with best goodjet
            {
                if( abs((lep1 + GoodBJetsVec_pt30.at(0).first).M() - 105) <= abs((lep2 + GoodBJetsVec_pt30.at(0).first).M() - 105) )
                {
                    TwoLep_Mbl1 = (lep1 + GoodBJetsVec_pt30.at(0).first).M();
                    TwoLep_Mbl1_Idx.first = GoodBJetsVec_pt30.at(0).second;
                    used_jet1 = GoodBJetsVec_pt30.at(0).second;
                    getJetMatch(Jets, GoodJets_pt30, GoodBJets_pt30, twoLep_mbl2diff, used_jet2, used_jet1, TwoLep_Mbl2, lep2, TwoLep_Mbl2_Idx, false);
                }
                else
                { 
                    TwoLep_Mbl2 = (lep2 + GoodBJetsVec_pt30.at(0).first).M();
                    TwoLep_Mbl2_Idx.first =GoodBJetsVec_pt30.at(0).second;
                    used_jet2 = GoodBJetsVec_pt30.at(0).second;
                    getJetMatch(Jets, GoodJets_pt30, GoodBJets_pt30, twoLep_mbl1diff, used_jet1, used_jet2, TwoLep_Mbl1, lep1, TwoLep_Mbl1_Idx, false);
                }
            }
            else           //with no bjets, leptons are paired with best goodjet
            {
                getJetMatch(Jets, GoodJets_pt30, GoodBJets_pt30, twoLep_mbl1diff, used_jet1, -1, TwoLep_Mbl1, lep1, TwoLep_Mbl1_Idx, false);
                getJetMatch(Jets, GoodJets_pt30, GoodBJets_pt30, twoLep_mbl2diff, used_jet2, -1, TwoLep_Mbl2, lep2, TwoLep_Mbl2_Idx, false);

                if(used_jet1 == used_jet2 && twoLep_mbl1diff >=  twoLep_mbl2diff) //if leptons are matched to same jet, keep better one and find a new match for other
                {
                    twoLep_mbl1diff = 999;
                    getJetMatch(Jets, GoodJets_pt30, GoodBJets_pt30, twoLep_mbl1diff, used_jet1, used_jet2, TwoLep_Mbl1, lep1, TwoLep_Mbl1_Idx, false);
                }
                else if(used_jet1 == used_jet2 && twoLep_mbl1diff <  twoLep_mbl2diff)
                {
                    twoLep_mbl2diff = 999;
                    getJetMatch(Jets, GoodJets_pt30, GoodBJets_pt30, twoLep_mbl2diff, used_jet2, used_jet1, TwoLep_Mbl2, lep2, TwoLep_Mbl2_Idx, false);
                }
            }
        }      

        tr.registerDerivedVar("TwoLep_Mbl1"+myVarSuffix_, TwoLep_Mbl1);
        tr.registerDerivedVar("TwoLep_Mbl2"+myVarSuffix_, TwoLep_Mbl2);
        tr.registerDerivedVar("TwoLep_Mbl1_Idx"+myVarSuffix_, TwoLep_Mbl1_Idx);
        tr.registerDerivedVar("TwoLep_Mbl2_Idx"+myVarSuffix_, TwoLep_Mbl2_Idx);

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
