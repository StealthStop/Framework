#ifndef BASELINE_H
#define BASELINE_H

class Baseline
{
private:
    std::string myVarSuffix_;

    void baseline(NTupleReader& tr)
    {
        const auto& runtype                = tr.getVar<std::string>("runtype");     
        const auto& filetag                = tr.getVar<std::string>("filetag");
        const auto& runYear                = tr.getVar<std::string>("runYear");
        const auto& RunNum                 = tr.getVar<unsigned int>("RunNum");
        const auto& blind                  = tr.getVar<bool>("blind");
        const auto& TriggerPass            = tr.getVec<int>("TriggerPass");
        const auto& TriggerNames           = utility::splitString(tr.getBranchTitle("TriggerPass"));
        // common variables for all channels 
        const auto& JetID                  = tr.getVar<bool>("JetID"                +myVarSuffix_);
        const auto& HT_trigger_pt30        = tr.getVar<double>("HT_trigger_pt30"    +myVarSuffix_);
        const auto& NNonIsoMuons           = tr.getVar<int>("NNonIsoMuons"          +myVarSuffix_);
        const auto& NGoodLeptons           = tr.getVar<int>("NGoodLeptons"          +myVarSuffix_);
        const auto& NGoodJets_pt30         = tr.getVar<int>("NGoodJets_pt30"        +myVarSuffix_);
        const auto& NGoodBJets_pt30        = tr.getVar<int>("NGoodBJets_pt30"       +myVarSuffix_);      
        // get variables for 0-Lepton
        const auto& HT_trigger_pt45        = tr.getVar<double>("HT_trigger_pt45"    +myVarSuffix_);
        const auto& NGoodJets_pt45         = tr.getVar<int>("NGoodJets_pt45"        +myVarSuffix_);
        const auto& NGoodBJets_pt45        = tr.getVar<int>("NGoodBJets_pt45"       +myVarSuffix_);
        const auto& dR_bjets_old           = tr.getVar<double>("dR_bjets_old"       +myVarSuffix_);
        const auto& dR_bjets               = tr.getVar<double>("dR_bjets"           +myVarSuffix_);
        const auto& ntops                  = tr.getVar<int>("ntops"                 +myVarSuffix_);
        // get variables for 1-Lepton
        const auto& Mbl                    = tr.getVar<double>("Mbl"                +myVarSuffix_);
        // get variables for 2-Lepton
        const auto& NGoodMuons             = tr.getVar<int>("NGoodMuons"            +myVarSuffix_);
        const auto& NGoodElectrons         = tr.getVar<int>("NGoodElectrons"        +myVarSuffix_);
        const auto& GoodLeptonsCharge      = tr.getVec<int>("GoodLeptonsCharge"     +myVarSuffix_);
        const auto& onZ                    = tr.getVar<bool>("onZ"                  +myVarSuffix_);
        // get variables for QCD CR        
        const auto& HT_NonIsoMuon_pt30     = tr.getVar<double>("HT_NonIsoMuon_pt30" +myVarSuffix_);
        const auto& NNonIsoMuonJets_pt30   = tr.getVar<int>("NNonIsoMuonJets_pt30"  +myVarSuffix_);
        // get variables for HEM veto 
        const auto& vetoedHEMelectron      = tr.getVar<bool>("vetoedHEMelectron"    +myVarSuffix_);
        bool passElectronHEMveto = !(vetoedHEMelectron && runtype == "Data" && runYear == "2018" && RunNum >= 319077);

        // --------------------
        // Data dependent stuff
        // --------------------
        bool passTriggerAllHad     = false, passTriggerMuonsRefAN = false; // 0-Lepton
        bool passTriggerMuon       = false, passTriggerElectron   = false; // 1-Lepton
        bool passTriggerNonIsoMuon = false;                                // QCD CR

        if (runYear.find("2016") != std::string::npos)
        {
            passTriggerAllHad     = PassTriggerAllHad2016(TriggerNames,     TriggerPass);
            passTriggerMuonsRefAN = PassTriggerMuonsRefAN(TriggerNames,     TriggerPass);
            passTriggerMuon       = PassTriggerMuon2016(TriggerNames,       TriggerPass);
            passTriggerElectron   = PassTriggerElectron2016(TriggerNames,   TriggerPass);
            passTriggerNonIsoMuon = PassTriggerNonIsoMuon2016(TriggerNames, TriggerPass);
            
        }
        else if (runYear == "2017")
        {
            passTriggerAllHad     = PassTriggerAllHad2017(TriggerNames,     TriggerPass);
            passTriggerMuonsRefAN = PassTriggerMuonsRefAN(TriggerNames,     TriggerPass);
            passTriggerMuon       = PassTriggerMuon2017(TriggerNames,       TriggerPass);
            passTriggerElectron   = PassTriggerElectron2017(TriggerNames,   TriggerPass);
            passTriggerNonIsoMuon = PassTriggerNonIsoMuon2017(TriggerNames, TriggerPass);           
        }
        else if (runYear == "2018pre" || runYear == "2018post" || runYear == "2018")
        {
            passTriggerAllHad     = PassTriggerAllHad2018(TriggerNames,     TriggerPass);
            passTriggerMuonsRefAN = PassTriggerMuonsRefAN(TriggerNames,     TriggerPass); 
            passTriggerMuon       = PassTriggerMuon2018(TriggerNames,       TriggerPass);
            passTriggerElectron   = PassTriggerElectron2018(TriggerNames,   TriggerPass);
            passTriggerNonIsoMuon = PassTriggerNonIsoMuon2018(TriggerNames, TriggerPass);
        }

        bool passTrigger         = true; // checking for jet, mu, el triggers
        bool passTriggerHadMC    = true; // checking emulators for jet triggers
        bool passTriggerMC       = true; // checking emulators for muon & electron triggers
        bool passNonIsoTrigger   = true; // checking for QCD CR triggers
        bool passNonIsoTriggerMC = true; // checking emulators forfor QCD CR
        bool passBlindHad_Good   = true;
        bool passBlindLep_Good   = true;        
        bool passBlind2Lep_Good  = true;
        if (runtype == "Data")
        {            
            // Pass the right trigger
            if (filetag.find("Data_JetHT")          != std::string::npos && !passTriggerAllHad)     passTrigger       = false;
            if (filetag.find("Data_SingleMuon")     != std::string::npos && !passTriggerMuon)       passTrigger       = false;
            if (filetag.find("Data_SingleElectron") != std::string::npos && !passTriggerElectron)   passTrigger       = false;
            if (filetag.find("Data_SingleMuon")     != std::string::npos && !passTriggerNonIsoMuon) passNonIsoTrigger = false;
        }
        
        // Blind data AND MC together, always
        if (NGoodJets_pt30 >= 9 && blind) passBlindHad_Good  = false;
        if (NGoodJets_pt30 >= 9 && blind) passBlindLep_Good  = false;
        if (NGoodJets_pt30 >= 8 && blind) passBlind2Lep_Good = false;

        // ------------------
        // MC dependent stuff
        // ------------------
        bool passMadHT = true;
        if(runtype == "MC")
        {
            const auto& madHT  = tr.getVar<float>("madHT");

            // Exclude events with MadGraph HT > 100 (70) from the WJets (DY) inclusive samples
            // in order to avoid double counting with the HT-binned samples
            // For DY, a special exception is made for 2016 and 2016APV where no HT-binned samples are present - 12 Oct 2022 JCH
            if(filetag.find("DYJetsToLL_M-50_Incl") != std::string::npos && madHT > 70 && runYear.find("2016") == std::string::npos) passMadHT = false;

            // For Wjets, 70to100 HT-binned samples are present for all years except 2016preVFP, which starts at 100to200
            if(filetag.find("WJetsToLNu_Incl") != std::string::npos && madHT > 70  && runYear.find("2016preVFP") == std::string::npos) passMadHT = false;
            if(filetag.find("WJetsToLNu_Incl") != std::string::npos && madHT > 100 && runYear.find("2016preVFP") != std::string::npos) passMadHT = false;

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
            
            // MC modeling of the trigger
            if( !passTriggerAllHad )                       passTriggerHadMC    = false;
            if( !passTriggerMuon && !passTriggerElectron ) passTriggerMC       = false;
            if( !passTriggerNonIsoMuon )                   passNonIsoTriggerMC = false;
        }

        // --------------------------
        // MET Filter dependent stuff
        // --------------------------
        const auto& globalSuperTightHalo2016Filter      = static_cast<bool>(tr.getVar<int>("globalSuperTightHalo2016Filter")    );
        const auto& PrimaryVertexFilter                 = static_cast<bool>(tr.getVar<int>("PrimaryVertexFilter")               );
        const auto& BadPFMuonFilter                     = static_cast<bool>(tr.getVar<int>("BadPFMuonFilter")                   );
        const auto& EcalDeadCellTriggerPrimitiveFilter  = static_cast<bool>(tr.getVar<int>("EcalDeadCellTriggerPrimitiveFilter"));
        const auto& HBHEIsoNoiseFilter                  = static_cast<bool>(tr.getVar<int>("HBHEIsoNoiseFilter")                );
        const auto& HBHENoiseFilter                     = static_cast<bool>(tr.getVar<int>("HBHENoiseFilter")                   );
        bool passMETFilters  = globalSuperTightHalo2016Filter && PrimaryVertexFilter && BadPFMuonFilter && EcalDeadCellTriggerPrimitiveFilter && HBHEIsoNoiseFilter && HBHENoiseFilter;
                
        // -----------------------------------
        // Define 0-Lepton Baseline Selections
        // -----------------------------------
        // full old baseline
        bool passBaseline0l_old = JetID                  &&
                                  passMETFilters         &&
                                  passMadHT              &&
                                  passTrigger            &&
                                  passTriggerHadMC       &&
                                  (runtype != "Data"     || filetag.find("Data_JetHT") != std::string::npos) &&
                                  NGoodLeptons == 0      &&
                                  HT_trigger_pt45 > 500  &&
                                  NGoodJets_pt45 >= 6    &&
                                  NGoodBJets_pt45 >= 2   &&
                                  ntops >= 2             &&
                                  dR_bjets_old >= 1.0    ;

        // proto baseline 
        bool passBaseline0l_pre = JetID                  &&
                                  passMETFilters         &&
                                  passMadHT              &&
                                  passTrigger            &&
                                  passTriggerHadMC       &&
                                  (runtype != "Data"     || filetag.find("Data_JetHT") != std::string::npos) &&
                                  NGoodLeptons == 0      &&
                                  HT_trigger_pt30 > 500  && 
                                  NGoodJets_pt45 >= 6    &&
                                  NGoodBJets_pt45 >= 1   ;

        // full baseline
        bool passBaseline0l_Good = passBaseline0l_pre   &&
                                   passElectronHEMveto  &&
                                   NNonIsoMuons == 0    &&
                                   NGoodJets_pt30 >= 7  &&
                                   ntops >= 2           &&
                                   NGoodBJets_pt30 >= 2 &&
                                   dR_bjets >= 1.0      ;
        
        // full baseline for bilind
        bool passBaseline0l_Good_blind = passBaseline0l_Good &&
                                         passBlindHad_Good;
       
        // baseline for HEM study 
        bool passBaseline0l_Good_noHEMveto = passBaseline0l_pre   &&
                                             NNonIsoMuons == 0    &&
                                             NGoodJets_pt30 >= 7  &&
                                             ntops >= 2           &&
                                             NGoodBJets_pt30 >= 2 &&
                                             dR_bjets >= 1.0      ;

        // baseline for trigger eff. & SF
        bool passBaseline0l_pt45 = JetID                 &&
                                   passMETFilters        &&
                                   passMadHT             &&
                                   passElectronHEMveto   &&
                                   NGoodMuons == 1       &&
                                   HT_trigger_pt45 > 500 &&
                                   NGoodJets_pt45 >= 6   ;

        // -----------------------------------
        // Define 1-Lepton Baseline Selections
        // -----------------------------------
        // proto baseline
        bool passBaselineGoodOffline1l = JetID                 &&
                                         passMETFilters        &&
                                         HT_trigger_pt30 > 500 &&
                                         passMadHT             &&
                                         NGoodBJets_pt30 >= 1  &&
                                         (50 < Mbl && Mbl < 250);
        // good muons for full baseline
        bool passBaseline1mu_Good = passBaselineGoodOffline1l &&
                                    passTrigger               &&
                                    passTriggerMC             &&
                                    (runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos) &&
                                    NGoodMuons == 1           && 
                                    NGoodElectrons == 0       &&
                                    NGoodJets_pt30 >= 7;
        
        // good lecetrons for full baseline
        bool passBaseline1el_Good = passBaselineGoodOffline1l &&
                                    passTrigger               &&
                                    passTriggerMC             &&
                                    (runtype != "Data" || filetag.find("Data_SingleElectron") != std::string::npos) &&
                                    NGoodElectrons == 1       &&
                                    NGoodMuons == 0           &&
                                    NGoodJets_pt30 >= 7;

        // full baseline 
        bool passBaseline1l_Good = (passBaseline1mu_Good || passBaseline1el_Good) &&
                                    passElectronHEMveto;

        // full baseline for blind
        bool passBaseline1l_Good_blind = (passBaseline1mu_Good || passBaseline1el_Good) &&
                                         passElectronHEMveto &&
                                         passBlindLep_Good;

        // baseline for HEM study  
        bool passBaseline1l_Good_noHEMveto = passBaseline1mu_Good || passBaseline1el_Good;

        // NOT sure what they are for !!!
        bool passBaseline1mu_Good_loose = passBaselineGoodOffline1l &&
                                          passTrigger               &&
                                          passTriggerMC             &&
                                          (runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos) &&
                                          NGoodMuons == 1           && 
                                          NGoodElectrons == 0       &&
                                          NGoodJets_pt30 >= 5;

        bool passBaseline1el_Good_loose = passBaselineGoodOffline1l &&
                                          passTrigger               &&
                                          passTriggerMC             &&
                                          (runtype != "Data" || filetag.find("Data_SingleElectron") != std::string::npos) &&
                                          NGoodElectrons == 1       &&
                                          NGoodMuons == 0           &&
                                          NGoodJets_pt30 >= 5;

        bool passBaseline1l_Good_loose = (passBaseline1mu_Good_loose || passBaseline1el_Good_loose) && 
                                         passElectronHEMveto &&
                                         passBlindLep_Good;

        bool passBaseline1l_HT500_Good = passBaseline1l_Good &&
                                         passElectronHEMveto &&
                                         HT_trigger_pt30 > 500;
                                         
        bool passBaseline1l_HT700_Good = passBaseline1l_Good &&
                                         passElectronHEMveto &&
                                         HT_trigger_pt30 > 700;

        // -----------------------------------
        // Define 2-Lepton Baseline Selections
        // ----------------------------------- 
        // proto baseline
        bool passBaseline2l_base = JetID                  &&
                                   passMETFilters         &&
                                   passMadHT              &&
                                   passTrigger            &&
                                   passTriggerMC          &&
                                   passElectronHEMveto    &&
                                   //passBlindLep_Good      &&                                  
                                   (runtype != "Data"  || (NGoodMuons >= 1 && filetag.find("Data_SingleMuon") != std::string::npos ) 
                                                       || (NGoodElectrons == 2 && filetag.find("Data_SingleElectron") != std::string::npos) ) &&
                                   NNonIsoMuons == 0      &&
                                   HT_trigger_pt30 > 500  &&
                                   NGoodBJets_pt30 >= 1   &&
                                   NGoodJets_pt30 >= 6    &&
                                   NGoodLeptons == 2 ? GoodLeptonsCharge[0]!=GoodLeptonsCharge[1] : false;
        
        // ful baseline
        bool passBaseline2l_Good       = passBaseline2l_base && !onZ;

        // full baseline for blind
        bool passBaseline2l_Good_blind = passBaseline2l_Good && passBlind2Lep_Good;

        // NOT sure what it is for !!! 
        bool passBaseline2lonZ_Good    = passBaseline2l_base &&  onZ;

        // -----------------------------------------
        // Define QCD CR Selections
        //  -- common for all 3 channels for now !!!
        // -----------------------------------------
        bool pass_qcdCR = JetID                    && 
                          passMETFilters           &&
                          passMadHT                &&
                          passNonIsoTrigger        &&
                          passNonIsoTriggerMC      &&
                          passElectronHEMveto      &&
                          (runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos) &&
                          HT_NonIsoMuon_pt30 > 500 &&
                          NNonIsoMuons == 1        &&
                          NGoodMuons == 0          &&
                          NGoodElectrons == 0      &&
                          NNonIsoMuonJets_pt30 >= 7;

        bool pass_qcdCR_1b    = pass_qcdCR && NGoodBJets_pt30 >= 1;
        bool pass_qcdCR_1t    = pass_qcdCR                          && ntops >= 1;
        bool pass_qcdCR_1b_1t = pass_qcdCR && NGoodBJets_pt30 >= 1  && ntops >= 1;
        bool pass_qcdCR_2b    = pass_qcdCR && NGoodBJets_pt30 >= 2;

        // -------------------
        // Register all things
        // -------------------
        // 0-lepton things 
        tr.registerDerivedVar<bool>("passBaseline0l_old"            +myVarSuffix_, passBaseline0l_old           ); // old baseline based on jet triggers      
        tr.registerDerivedVar<bool>("passBaseline0l_pre"            +myVarSuffix_, passBaseline0l_pre           ); 
        tr.registerDerivedVar<bool>("passBaseline0l_Good"           +myVarSuffix_, passBaseline0l_Good          );
        tr.registerDerivedVar<bool>("passBaseline0l_Good_blind"     +myVarSuffix_, passBaseline0l_Good_blind    ); // for blinding
        tr.registerDerivedVar<bool>("passBaseline0l_Good_noHEMveto" +myVarSuffix_, passBaseline0l_Good_noHEMveto); // for HEM study
        tr.registerDerivedVar<bool>("passBaseline0l_pt45"           +myVarSuffix_, passBaseline0l_pt45          ); // for trigger SF
        tr.registerDerivedVar<bool>("passBlindHad_Good"             +myVarSuffix_, passBlindHad_Good            );
        tr.registerDerivedVar<bool>("passTriggerAllHad"             +myVarSuffix_, passTriggerAllHad            );
        tr.registerDerivedVar<bool>("passTriggerHadMC"              +myVarSuffix_, passTriggerHadMC             );
        tr.registerDerivedVar<bool>("passTriggerMuonsRefAN"         +myVarSuffix_, passTriggerMuonsRefAN        ); // muon trigger for preselection of trigger SF 
        // 1-lepton things
        tr.registerDerivedVar<bool>("passBaselineGoodOffline1l"     +myVarSuffix_, passBaselineGoodOffline1l    );
        tr.registerDerivedVar<bool>("passBaseline1l_Good"           +myVarSuffix_, passBaseline1l_Good          );
        tr.registerDerivedVar<bool>("passBaseline1l_Good_blind"     +myVarSuffix_, passBaseline1l_Good_blind    );
        tr.registerDerivedVar<bool>("passBaseline1l_Good_noHEMveto" +myVarSuffix_, passBaseline1l_Good_noHEMveto);
        tr.registerDerivedVar<bool>("passBaseline1l_Good_loose"     +myVarSuffix_, passBaseline1l_Good_loose    );
        tr.registerDerivedVar<bool>("passBaseline1l_HT500_Good"     +myVarSuffix_, passBaseline1l_HT500_Good    );
        tr.registerDerivedVar<bool>("passBaseline1l_HT700_Good"     +myVarSuffix_, passBaseline1l_HT700_Good    );
        tr.registerDerivedVar<bool>("passBlindLep_Good"             +myVarSuffix_, passBlindLep_Good            );
        tr.registerDerivedVar<bool>("passTriggerMuon"               +myVarSuffix_, passTriggerMuon              );
        tr.registerDerivedVar<bool>("passTriggerElectron"           +myVarSuffix_, passTriggerElectron          );
        // 2-lepton things
        tr.registerDerivedVar<bool>("passBaseline2l_base "          +myVarSuffix_, passBaseline2l_base          ); 
        tr.registerDerivedVar<bool>("passBaseline2l_Good"           +myVarSuffix_, passBaseline2l_Good          );
        tr.registerDerivedVar<bool>("passBaseline2l_Good_blind"     +myVarSuffix_, passBaseline2l_Good_blind    );
        tr.registerDerivedVar<bool>("passBaseline2lonZ_Good"        +myVarSuffix_, passBaseline2lonZ_Good       );
        // QCD CR things        
        tr.registerDerivedVar<bool>("pass_qcdCR"                    +myVarSuffix_, pass_qcdCR                   );
        tr.registerDerivedVar<bool>("pass_qcdCR_1b"                 +myVarSuffix_, pass_qcdCR_1b                );
        tr.registerDerivedVar<bool>("pass_qcdCR_1t"                 +myVarSuffix_, pass_qcdCR_1t                );
        tr.registerDerivedVar<bool>("pass_qcdCR_1b_1t"              +myVarSuffix_, pass_qcdCR_1b_1t             );
        tr.registerDerivedVar<bool>("pass_qcdCR_2b"                 +myVarSuffix_, pass_qcdCR_2b                );
        tr.registerDerivedVar<bool>("passNonIsoTrigger"             +myVarSuffix_, passNonIsoTrigger            );
        tr.registerDerivedVar<bool>("passNonIsoTriggerMC"           +myVarSuffix_, passNonIsoTriggerMC          );
        // common things       
        tr.registerDerivedVar<bool>("passMadHT"                     +myVarSuffix_, passMadHT                    );
        tr.registerDerivedVar<bool>("passMETFilters"                +myVarSuffix_, passMETFilters               );
        tr.registerDerivedVar<bool>("passElectronHEMveto"           +myVarSuffix_, passElectronHEMveto          );
        tr.registerDerivedVar<bool>("passTrigger"                   +myVarSuffix_, passTrigger                  );
        tr.registerDerivedVar<bool>("passTriggerMC"                 +myVarSuffix_, passTriggerMC                );
        
    }

    bool PassTriggerGeneral(std::vector<std::string>& mytriggers, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        bool passTrigger = false;
        for(unsigned int i=0; i<TriggerNames.size(); ++i)
        {
            if(TriggerPass.at(i) != 1)
                continue;
            std::string trigname = TriggerNames.at(i);
            if( std::any_of(mytriggers.begin(), mytriggers.end(), [&] (std::string s) { return trigname.find(s)!=std::string::npos; }) )
            {
                passTrigger = true;
                break;
            }
        }
        return passTrigger;
    }

    // -------------------------
    // Jet Triggers for 0-Lepton
    // -------------------------
    bool PassTriggerAllHad2016(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {
            "HLT_PFJet450",
            "HLT_PFHT900",
            "HLT_PFHT450_SixJet40_BTagCSV_p056",
            "HLT_PFHT400_SixJet30_DoubleBTagCSV_p056",            
        };
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerAllHad2017(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {
            "HLT_PFJet500",
            "HLT_PFHT1050",
            "HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2",
            "HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2",
            "HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5",
            "HLT_PFHT380_SixJet32_DoubleBTagCSV_p075",
            "HLT_PFHT430_SixJet40_BTagCSV_p080",
        };
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerAllHad2018(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass) 
    {
        std::vector<std::string> mytriggers = {
            "HLT_PFJet500",
            "HLT_PFHT1050",
            "HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2",
            "HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2",
            "HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5",
            "HLT_PFHT430_SixPFJet40_PFBTagDeepCSV_1p5",
            "HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94",
            "HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59",
        };
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    // muon trigger for preselection of trigger SF (0-Lepton)
    bool PassTriggerMuonsRefAN(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_IsoMu24", "HLT_IsoMu22_eta2p1"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    // -------------------------------------
    // Muon & Electron Triggers for 1-Lepton
    // -------------------------------------
    bool PassTriggerMuon2016(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_IsoMu24","HLT_IsoTkMu24","HLT_Mu50","HLT_TkMu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerMuon2017(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        // IsoMu24 ended up prescaled for 3.6\fb, so also include the IsoMu27 here
        std::vector<std::string> mytriggers = {"HLT_IsoMu24","HLT_IsoMu27","HLT_Mu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerMuon2018(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_IsoMu24","HLT_IsoMu27","HLT_Mu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerElectron2016(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Ele27_WPTight_Gsf","HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon175"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerElectron2017(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        // Ele35 fine for whole year, gets complicated for the other triggers
        std::vector<std::string> mytriggers = {"HLT_Ele35_WPTight","HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon200"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerElectron2018(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Ele35_WPTight","HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon200"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerPhoton(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Photon175_v"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }
   
    // ----------------------------
    // NonIsoMu Triggers for QCD CR
    // ----------------------------
    bool PassTriggerNonIsoMuon2016(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Mu50", "HLT_TkMu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }
    
    bool PassTriggerNonIsoMuon2017(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Mu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerNonIsoMuon2018(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Mu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }



public:
    Baseline(std::string myVarSuffix = "")
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up Baseline"<<std::endl;
    }
    
    void operator()(NTupleReader& tr)
    {
        baseline(tr);
    }
};

#endif
