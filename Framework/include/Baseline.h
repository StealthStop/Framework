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
        const auto& TriggerNames           = utility::splitString(tr.getBranchTitle("TriggerPass"));
        const auto& TriggerPass            = tr.getVec<int>("TriggerPass");
        const auto& NNonIsoMuons           = tr.getVar<int>("NNonIsoMuons"+myVarSuffix_);
        const auto& NGoodLeptons           = tr.getVar<int>("NGoodLeptons"+myVarSuffix_);
        const auto& GoodLeptonsCharge      = tr.getVec<int>("GoodLeptonsCharge"+myVarSuffix_);
        const auto& NGoodMuons             = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        const auto& NGoodPlusMuons         = tr.getVar<int>("NGoodPlusMuons"+myVarSuffix_);
        const auto& NGoodMinusMuons        = tr.getVar<int>("NGoodMinusMuons"+myVarSuffix_);
        const auto& NGoodElectrons         = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);
        const auto& NGoodPlusElectrons     = tr.getVar<int>("NGoodPlusElectrons"+myVarSuffix_);
        const auto& NGoodMinusElectrons    = tr.getVar<int>("NGoodMinusElectrons"+myVarSuffix_);
        const auto& NGoodLeptons_pt20      = tr.getVar<int>("NGoodLeptons_pt20"+myVarSuffix_);
        const auto& GoodLeptonsCharge_pt20 = tr.getVec<int>("GoodLeptonsCharge_pt20"+myVarSuffix_);
        const auto& HT_trigger_pt30        = tr.getVar<double>("HT_trigger_pt30"+myVarSuffix_);
        const auto& HT_trigger_pt45        = tr.getVar<double>("HT_trigger_pt45"+myVarSuffix_);
        const auto& HT_NonIsoMuon_pt30     = tr.getVar<double>("HT_NonIsoMuon_pt30"+myVarSuffix_);
        const auto& onZ                    = tr.getVar<bool>("onZ"+myVarSuffix_); 
        const auto& JetID                  = tr.getVar<bool>("JetID"+myVarSuffix_);
        const auto& NGoodJets_pt40         = tr.getVar<int>("NGoodJets_pt40"+myVarSuffix_); 
        const auto& NGoodJets_pt45         = tr.getVar<int>("NGoodJets_pt45"+myVarSuffix_);
        const auto& NGoodBJets_pt45        = tr.getVar<int>("NGoodBJets_pt45"+myVarSuffix_);
        const auto& NGoodJets_pt30         = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_);
        const auto& NGoodBJets_pt30        = tr.getVar<int>("NGoodBJets_pt30"+myVarSuffix_); 
        const auto& NNonIsoMuonJets_pt30   = tr.getVar<int>("NNonIsoMuonJets_pt30"+myVarSuffix_);         
        const auto& NGoodPhotons           = tr.getVar<int>("NGoodPhotons"+myVarSuffix_);
        const auto& Mbl                    = tr.getVar<double>("Mbl"+myVarSuffix_);
        const auto& vetoedHEMelectron      = tr.getVar<bool>("vetoedHEMelectron"+myVarSuffix_);
        const auto& NGoodBJetsCSV_pt30     = tr.getVar<int>("NGoodBJetsCSV_pt30"+myVarSuffix_); 
        const auto& ntops                  = tr.getVar<int>("ntops"+myVarSuffix_);
        const auto& dR_bjets_old           = tr.getVar<double>("dR_bjets_old"+myVarSuffix_);
        const auto& dR_bjets               = tr.getVar<double>("dR_bjets"+myVarSuffix_); 
 
        bool passElectronHEMveto = !(vetoedHEMelectron && runtype == "Data" && runYear == "2018" && RunNum >= 319077);

        // ------------------------------
        // -- Data dependent stuff
        // ------------------------------
        bool passTriggerPhoton = PassTriggerPhoton(TriggerNames, TriggerPass);
        bool passTriggerAllHad = false, passTriggerMuon = false, passTriggerElectron = false, passTriggerNonIsoMuon = false, passTriggerIsoMu = false;
        bool passTriggerMuonsRefAN = false, passTriggerRefAN = false;  
        if (runYear.find("2016") != std::string::npos)
        {
            passTriggerAllHad = PassTriggerAllHad2016(TriggerNames, TriggerPass);
            passTriggerMuon = PassTriggerMuon2016(TriggerNames, TriggerPass);
            passTriggerElectron = PassTriggerElectron2016(TriggerNames, TriggerPass);
            passTriggerNonIsoMuon = PassTriggerNonIsoMuon2016(TriggerNames, TriggerPass);
            passTriggerIsoMu = PassTriggerIsoMu2016(TriggerNames, TriggerPass);
            passTriggerMuonsRefAN = PassTriggerMuonsRefAN(TriggerNames, TriggerPass);
            passTriggerRefAN      = PassTriggerRefAN(TriggerNames, TriggerPass);
            
        }
        else if (runYear == "2017")
        {
            passTriggerAllHad = PassTriggerAllHad2017(TriggerNames, TriggerPass);
            passTriggerMuon = PassTriggerMuon2017(TriggerNames, TriggerPass);
            passTriggerElectron = PassTriggerElectron2017(TriggerNames, TriggerPass);
            passTriggerNonIsoMuon = PassTriggerNonIsoMuon2017(TriggerNames, TriggerPass);           
            passTriggerIsoMu = PassTriggerIsoMu2017(TriggerNames, TriggerPass);
            passTriggerMuonsRefAN = PassTriggerMuonsRefAN(TriggerNames, TriggerPass); 
        }
        else if (runYear == "2018pre" || runYear == "2018post" || runYear == "2018")
        {
            passTriggerAllHad = PassTriggerAllHad2018(TriggerNames, TriggerPass); 
            passTriggerMuon = PassTriggerMuon2018(TriggerNames, TriggerPass);
            passTriggerElectron = PassTriggerElectron2018(TriggerNames, TriggerPass);
            passTriggerNonIsoMuon = PassTriggerNonIsoMuon2018(TriggerNames, TriggerPass);
            passTriggerIsoMu = PassTriggerIsoMu2018(TriggerNames, TriggerPass);
            passTriggerMuonsRefAN = PassTriggerMuonsRefAN(TriggerNames, TriggerPass);
        }

        bool passTrigger   = true;
        bool passTriggerMC = true;
        bool passTriggerHadMC = true; 
        bool passNonIsoTrigger = true;
        bool passIsoMuTrigger = true;
        bool passNonIsoTriggerMC = true;
        bool passIsoMuTriggerMC = true;
        bool passBlindHad_Good = true;
        bool passBlindLep_Good = true;        
        bool passBlind2Lep_Good = true;
        if (runtype == "Data")
        {            
            // Pass the right trigger
            if (filetag.find("Data_JetHT") != std::string::npos && !passTriggerAllHad) passTrigger = false;
            if (filetag.find("Data_SingleMuon") != std::string::npos && !passTriggerMuon) passTrigger = false;
            if (filetag.find("Data_SingleElectron") != std::string::npos && !passTriggerElectron) passTrigger = false;
            if (filetag.find("Data_SinglePhoton") != std::string::npos && !passTriggerPhoton) passTrigger = false;
            if (filetag.find("Data_SingleMuon") != std::string::npos && !passTriggerNonIsoMuon) passNonIsoTrigger = false;
            if (filetag.find("Data_SingleMuon") != std::string::npos && !passTriggerIsoMu) passIsoMuTrigger = false;
        }
        
        // Blind data AND MC together, always
        if (NGoodJets_pt30 >= 9 && blind) passBlindHad_Good = false;
        if (NGoodJets_pt30 >= 9 && blind) passBlindLep_Good = false;
        if (NGoodJets_pt30 >= 8 && blind) passBlind2Lep_Good = false;

        // ------------------------
        // -- MC dependent stuff - 
        // -----------------------
        bool passMadHT = true;
        if(runtype == "MC")
        {
            const auto& madHT  = tr.getVar<float>("madHT");

            // Exclude events with MadGraph HT > 100 (70) from the WJets (DY) inclusive samples
            // in order to avoid double counting with the HT-binned samples
            // For DY, a special exception is made for 2016 and 2016APV where no HT-binned samples are present - 12 Oct 2022 JCH
            if(filetag.find("DYJetsToLL_M-50_Incl") != std::string::npos && madHT > 70 && runYear.find("2016") == std::string::npos) passMadHT = false;
            if(filetag.find("WJetsToLNu_Incl") != std::string::npos      && madHT > 100) passMadHT = false;

            // Stitch TTbar samples together
            // remove HT overlap
            if((filetag.find("TTJets_Incl") != std::string::npos || 
                filetag.find("TTJets_SingleLeptFromT") != std::string::npos || 
                filetag.find("TTJets_SingleLeptFromTbar") != std::string::npos || 
                filetag.find("TTJets_DiLept") != std::string::npos) && madHT > 600) 
            {
                passMadHT = false;
            }
            // also remove lepton overlap from the inclusive sample
            const auto& GenElectrons        = tr.getVec<utility::LorentzVector>("GenElectrons");
            const auto& GenMuons            = tr.getVec<utility::LorentzVector>("GenMuons");
            const auto& GenTaus             = tr.getVec<utility::LorentzVector>("GenTaus");
            int NGenLeptons = GenElectrons.size() + GenMuons.size() + GenTaus.size();
            if (filetag.find("TTJets_Incl") != std::string::npos && NGenLeptons > 0) passMadHT = false;
            
            // MC modeling of the trigger
            if( !passTriggerAllHad ) passTriggerHadMC = false;
            if( !passTriggerMuon && !passTriggerElectron ) passTriggerMC = false;
            if( !passTriggerNonIsoMuon ) passNonIsoTriggerMC = false;
            if( !passTriggerIsoMu ) passIsoMuTriggerMC = false;
        }

        //------------------------
        // -- MET Filter dependent stuff
        // -----------------------
        const auto& globalSuperTightHalo2016Filter      = static_cast<bool>( tr.getVar<int>("globalSuperTightHalo2016Filter") );
        const auto& PrimaryVertexFilter                 = static_cast<bool>( tr.getVar<int>("PrimaryVertexFilter") );
        const auto& BadPFMuonFilter                     = static_cast<bool>( tr.getVar<int>("BadPFMuonFilter") );
        const auto& EcalDeadCellTriggerPrimitiveFilter  = static_cast<bool>( tr.getVar<int>("EcalDeadCellTriggerPrimitiveFilter") );
        const auto& HBHEIsoNoiseFilter                  = static_cast<bool>( tr.getVar<int>("HBHEIsoNoiseFilter") );
        const auto& HBHENoiseFilter                     = static_cast<bool>( tr.getVar<int>("HBHENoiseFilter") );
        bool passMETFilters  = globalSuperTightHalo2016Filter && PrimaryVertexFilter && BadPFMuonFilter && EcalDeadCellTriggerPrimitiveFilter && HBHEIsoNoiseFilter && HBHENoiseFilter;
                
        // --------------------------------
        // -- Define 0 Lepton Selections  
        // --------------------------------
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
        bool passBaseline0l_good = passBaseline0l_pre   &&
                                   passElectronHEMveto  &&
                                   NNonIsoMuons == 0    &&
                                   NGoodJets_pt30 >= 7  &&
                                   ntops >= 2           &&
                                   NGoodBJets_pt30 >= 2 &&
                                   dR_bjets >= 1.0      ;
       
        bool passBaseline0l_good_blind = passBaseline0l_good &&
                                         passBlindHad_Good;
 
        bool passBaseline0l_good_noHEMveto = passBaseline0l_pre   &&
                                             NNonIsoMuons == 0    &&
                                             NGoodJets_pt30 >= 7  &&
                                             ntops >= 2           &&
                                             NGoodBJets_pt30 >= 2 &&
                                             dR_bjets >= 1.0      ;

        bool passBaseline0l_pre_loose = JetID                  &&
                                        passMETFilters         &&
                                        passMadHT              &&
                                        passTrigger            &&
                                        passTriggerHadMC       &&
                                        (runtype != "Data"     || filetag.find("Data_JetHT") != std::string::npos) &&
                                        NGoodLeptons == 0      &&
                                        HT_trigger_pt30 > 500  && 
                                        NGoodJets_pt45 >= 5    &&
                                        NGoodBJets_pt45 >= 1   ;

        bool passBaseline0l_good_loose = passBaseline0l_pre_loose   &&
                                         passElectronHEMveto  &&
                                         NNonIsoMuons == 0    &&
                                         NGoodJets_pt30 >= 5  &&
                                         ntops >= 2           &&
                                         NGoodBJets_pt30 >= 2 &&
                                         dR_bjets >= 1.0      &&
                                         passBlindHad_Good;

        // QCD CR 
        bool pass_qcdCR = JetID                     &&
                          passMETFilters            &&
                          passMadHT                 &&
                          passNonIsoTrigger         &&
                          passNonIsoTriggerMC       &&
                          (runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos) &&
                          NGoodLeptons == 0         &&
                          NNonIsoMuons == 1         &&
                          HT_NonIsoMuon_pt30 > 300  &&
                          NNonIsoMuonJets_pt30 >= 7 ;

        // ----------------------------------------------------------
        // -- Define 0 Lepton Baseline for Trigger Efficiency & SF
        // ----------------------------------------------------------
        bool passBaseline0l_hadTrig= JetID                 &&
                                     passMETFilters        &&
                                     passMadHT             &&
                                     NGoodLeptons == 0     &&
                                     HT_trigger_pt45 > 300 &&
                                     NGoodJets_pt45 >= 5   &&
                                     NGoodBJets_pt45 >= 1  ; // >= 0

        bool passBaseline0l_hadMuTrig = JetID                 &&
                                        passMETFilters        &&
                                        passMadHT             &&
                                        passIsoMuTrigger      &&
                                        NGoodMuons == 1       &&
                                        HT_trigger_pt45 > 300 &&
                                        NGoodJets_pt45 >= 5   &&
                                        NGoodBJets_pt45 >= 1  ;

        // reference analysis pre-selections
        bool passBaseline0l_refAN = JetID                 &&
                                    passMETFilters        &&
                                    passMadHT             &&
                                    //passTriggerMuonsRefAN &&
                                    NGoodMuons == 1       &&
                                    HT_trigger_pt30 > 500 &&
                                    NGoodJets_pt40 >= 6   && 
                                    NGoodBJets_pt30 >= 2  ; 

        // reference analysis pre-selections with our pt > 45 cut
        bool passBaseline0l_refAN_pt45 = JetID                 &&
                                         passMETFilters        &&
                                         passMadHT             &&
                                         NGoodMuons == 1       &&
                                         HT_trigger_pt45 > 500 &&
                                         NGoodJets_pt45 >= 6   &&
                                         NGoodBJets_pt45 >= 2  ;

        // reference analysis pre-selections with CSV
        bool passBaseline0l_csv_refAN = JetID                 &&
                                        passMETFilters        &&
                                        passMadHT             &&
                                        NGoodMuons == 1       &&
                                        HT_trigger_pt30 > 500 &&
                                        NGoodJets_pt40 >= 6   &&
                                        NGoodBJetsCSV_pt30 >= 2; 

        // latest preselctions with HEMveto
        bool passBaseline0l_pt45 = JetID                 &&
                                   passMETFilters        &&
                                   passMadHT             &&
                                   NGoodMuons == 1       &&
                                   HT_trigger_pt45 > 500 &&
                                   NGoodJets_pt45 >= 6   &&
                                   NGoodBJets_pt45 >= 2  ;

        // ------------------------------------
        // -- Define 1-lepton proto-baseline
        // ------------------------------------
        bool passBaselineGoodOffline1l = JetID                 &&
                                         passMETFilters        &&
                                         HT_trigger_pt30 > 500 &&
                                         passMadHT             &&
                                         NGoodBJets_pt30 >= 1  &&
                                         (50 < Mbl && Mbl < 250);

        // -------------------------------
        // -- Define 1 Lepton Baseline
        // -------------------------------
        bool passBaseline1mu_Good = passBaselineGoodOffline1l &&
                                    passTrigger               &&
                                    passTriggerMC             &&
                                    (runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos) &&
                                    NGoodMuons == 1           && 
                                    NGoodElectrons == 0       &&
                                    NGoodJets_pt30 >= 7;

        bool passBaseline1el_Good = passBaselineGoodOffline1l &&
                                    passTrigger               &&
                                    passTriggerMC             &&
                                    (runtype != "Data" || filetag.find("Data_SingleElectron") != std::string::npos) &&
                                    NGoodElectrons == 1       &&
                                    NGoodMuons == 0           &&
                                    NGoodJets_pt30 >= 7;

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

        bool passBaseline1l_Good_noHEMveto = passBaseline1mu_Good || passBaseline1el_Good;

        bool passBaseline1l_Good = (passBaseline1mu_Good || passBaseline1el_Good) && 
                                    passElectronHEMveto;

        bool passBaseline1l_Good_blind = (passBaseline1mu_Good || passBaseline1el_Good) && 
                                         passElectronHEMveto &&
                                         passBlindLep_Good;

        bool passBaseline1l_Good_loose = (passBaseline1mu_Good_loose || passBaseline1el_Good_loose) && 
                                         passElectronHEMveto &&
                                         passBlindLep_Good;

        bool passBaseline1l_HT500_Good = passBaseline1l_Good &&
                                         passElectronHEMveto &&
                                         HT_trigger_pt30 > 500;
                                         
        bool passBaseline1l_HT700_Good = passBaseline1l_Good &&
                                         passElectronHEMveto &&
                                         HT_trigger_pt30 > 700;

        bool passBaseline1l_NonIsoMuon = HT_NonIsoMuon_pt30 > 500 &&
                                         passElectronHEMveto      &&
                                         passMETFilters           &&
                                         passMadHT                &&
                                         passNonIsoTrigger        &&
                                         passNonIsoTriggerMC      &&
                                         (runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos) &&
                                         NNonIsoMuons == 1        &&
                                         NGoodMuons == 0          &&
                                         NGoodElectrons == 0      &&
                                         JetID                    &&
                                         NNonIsoMuonJets_pt30 >= 7;

        // ----------------------------------
        // -- Define 2 Lepton Baseline
        // ----------------------------------  
        bool passBaseline2l_base = JetID                  &&
                                   passMETFilters         &&
                                   passMadHT              &&
                                   passTrigger            &&
                                   passTriggerMC          &&
                                   passElectronHEMveto    &&
                                   NNonIsoMuons == 0      &&
                                   HT_trigger_pt30 > 500  &&
                                   //passBlindLep_Good      &&                                  
                                   (runtype != "Data"  || (NGoodMuons >= 1 && filetag.find("Data_SingleMuon") != std::string::npos ) 
                                                       || (NGoodElectrons == 2 && filetag.find("Data_SingleElectron") != std::string::npos) ) &&
                                   NGoodBJets_pt30 >= 1   &&
                                   NGoodJets_pt30 >= 6    &&
                                   NGoodLeptons == 2 ? GoodLeptonsCharge[0]!=GoodLeptonsCharge[1] : false;
   
        bool passBaseline2lonZ_Good    = passBaseline2l_base &&  onZ;
        bool passBaseline2l_Good       = passBaseline2l_base && !onZ;
        bool passBaseline2l_Good_blind = passBaseline2l_Good && passBlind2Lep_Good;

        // -----------------------------------
        // -- Define 2 Lepton pt20 Baseline
        // -----------------------------------
        bool passBaseline2l_pt20 =    JetID              &&
                                      passMETFilters     &&
                                      passMadHT          &&
                                      //passTrigger        &&
                                      //passTriggerMC      &&
                                      (runtype != "Data" || (NGoodMuons == 2 && filetag.find("Data_SingleMuon") != std::string::npos ) 
                                                         || (NGoodElectrons == 2 && filetag.find("Data_SingleElectron") != std::string::npos) ) &&
                                      NGoodJets_pt30 >= 6 &&
                                      NGoodBJets_pt30 >= 1 && 
                                      NGoodLeptons_pt20 == 2 ? GoodLeptonsCharge_pt20[0]!=GoodLeptonsCharge_pt20[1] : false;

        // -----------------------------------
        // -- Define 2 Lepton pt30 Baseline
        // -----------------------------------
        bool passBaseline2l_pt30 =    JetID              &&
                                      passMETFilters     &&
                                      passMadHT          &&
                                      //passTrigger        &&
                                      //passTriggerMC      &&
                                      (runtype != "Data" || (NGoodMuons == 2 && filetag.find("Data_SingleMuon") != std::string::npos ) 
                                                         || (NGoodElectrons == 2 && filetag.find("Data_SingleElectron") != std::string::npos) ) &&
                                      NGoodJets_pt30  >= 4 &&
                                      NGoodBJets_pt30 >= 1 &&
                                      NGoodLeptons == 2 ? GoodLeptonsCharge[0]!=GoodLeptonsCharge[1] : false;

        // -----------------------------------
        // -- Define 1 e and 1 m Baseline
        // -----------------------------------        
        bool passBaseline1e1m_Good = JetID                 &&
                                     passMETFilters        &&
                                     HT_trigger_pt30 > 300 &&
                                     passMadHT             &&
                                     passTriggerMuon       &&
                                     NGoodMuons == 1       &&
                                     NGoodElectrons == 1   &&
                                     (runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos) &&
                                     (NGoodPlusMuons + NGoodPlusElectrons == 1)          &&
                                     (NGoodMinusMuons + NGoodMinusElectrons == 1)        &&
                                     NGoodJets_pt30 <= 6   && 
                                     NGoodBJets_pt30 >= 2;
        
        // -----------------------------------
        // -- Define 1 Photon Baseline
        // -----------------------------------
        bool passBaseline1photon_Good = passMadHT           &&
                                        passMETFilters      &&
                                        passTrigger         &&
                                        passTriggerMC       &&
                                        NGoodPhotons == 1   &&
                                        NGoodLeptons == 0   && 
                                        NGoodJets_pt30 >= 7; 

        tr.registerDerivedVar<bool>("passBaseline0l_old"+myVarSuffix_,            passBaseline0l_old);       
        tr.registerDerivedVar<bool>("passBaseline0l_pre"+myVarSuffix_,            passBaseline0l_pre); 
        tr.registerDerivedVar<bool>("passBaseline0l_good"+myVarSuffix_,           passBaseline0l_good);
        tr.registerDerivedVar<bool>("passBaseline0l_good_blind"+myVarSuffix_,     passBaseline0l_good_blind);
        tr.registerDerivedVar<bool>("passBaseline0l_good_loose"+myVarSuffix_,     passBaseline0l_good_loose);
        tr.registerDerivedVar<bool>("passBaseline0l_good_noHEMveto"+myVarSuffix_, passBaseline0l_good_noHEMveto);
        tr.registerDerivedVar<bool>("pass_qcdCR"+myVarSuffix_,                    pass_qcdCR);
        tr.registerDerivedVar<bool>("passBaseline0l_hadTrig"+myVarSuffix_,        passBaseline0l_hadTrig);    // 0l trigger study
        tr.registerDerivedVar<bool>("passBaseline0l_hadMuTrig"+myVarSuffix_,      passBaseline0l_hadMuTrig);  //
        tr.registerDerivedVar<bool>("passBaseline0l_refAN"+myVarSuffix_,          passBaseline0l_refAN);      //
        tr.registerDerivedVar<bool>("passBaseline0l_refAN_pt45"+myVarSuffix_,     passBaseline0l_refAN_pt45); //
        tr.registerDerivedVar<bool>("passBaseline0l_csv_refAN"+myVarSuffix_,      passBaseline0l_csv_refAN);  //
        tr.registerDerivedVar<bool>("passBaseline0l_pt45"+myVarSuffix_,           passBaseline0l_pt45);       // 0l trigger study
        tr.registerDerivedVar<bool>("passBaselineGoodOffline1l"+myVarSuffix_,     passBaselineGoodOffline1l);
        tr.registerDerivedVar<bool>("passBaseline1l_Good"+myVarSuffix_,           passBaseline1l_Good);
        tr.registerDerivedVar<bool>("passBaseline1l_Good_blind"+myVarSuffix_,     passBaseline1l_Good_blind);
        tr.registerDerivedVar<bool>("passBaseline1l_Good_loose"+myVarSuffix_,     passBaseline1l_Good_loose);
        tr.registerDerivedVar<bool>("passBaseline1l_Good_noHEMveto"+myVarSuffix_, passBaseline1l_Good_noHEMveto);
        tr.registerDerivedVar<bool>("passBaseline1l_HT500_Good"+myVarSuffix_,     passBaseline1l_HT500_Good);
        tr.registerDerivedVar<bool>("passBaseline1l_HT700_Good"+myVarSuffix_,     passBaseline1l_HT700_Good);
        tr.registerDerivedVar<bool>("passBaseline1l_NonIsoMuon"+myVarSuffix_,     passBaseline1l_NonIsoMuon);
        tr.registerDerivedVar<bool>("passBaseline2lonZ_Good"+myVarSuffix_,        passBaseline2lonZ_Good);
        tr.registerDerivedVar<bool>("passBaseline2l_Good"+myVarSuffix_,           passBaseline2l_Good);
        tr.registerDerivedVar<bool>("passBaseline2l_Good_blind"+myVarSuffix_,     passBaseline2l_Good_blind);
        tr.registerDerivedVar<bool>("passBaseline2l_pt20"+myVarSuffix_,           passBaseline2l_pt20);
        tr.registerDerivedVar<bool>("passBaseline2l_pt30"+myVarSuffix_,           passBaseline2l_pt30);
        tr.registerDerivedVar<bool>("passBaseline1photon_Good"+myVarSuffix_,      passBaseline1photon_Good);
        tr.registerDerivedVar<bool>("passBaseline1e1m_Good"+myVarSuffix_,         passBaseline1e1m_Good);
        tr.registerDerivedVar<bool>("passBlindHad_Good"+myVarSuffix_,             passBlindHad_Good);
        tr.registerDerivedVar<bool>("passBlindLep_Good"+myVarSuffix_,             passBlindLep_Good);
        tr.registerDerivedVar<bool>("passTriggerAllHad"+myVarSuffix_,             passTriggerAllHad);
        tr.registerDerivedVar<bool>("passTriggerMuon"+myVarSuffix_,               passTriggerMuon);
        tr.registerDerivedVar<bool>("passTriggerElectron"+myVarSuffix_,           passTriggerElectron);
        tr.registerDerivedVar<bool>("passTrigger"+myVarSuffix_,                   passTrigger);
        tr.registerDerivedVar<bool>("passNonIsoTrigger"+myVarSuffix_,             passNonIsoTrigger);
        tr.registerDerivedVar<bool>("passIsoMuTrigger"+myVarSuffix_,              passIsoMuTrigger);
        tr.registerDerivedVar<bool>("passTriggerMuonsRefAN"+myVarSuffix_,         passTriggerMuonsRefAN);
        tr.registerDerivedVar<bool>("passTriggerRefAN"+myVarSuffix_,              passTriggerRefAN); 
        tr.registerDerivedVar<bool>("passTriggerMC"+myVarSuffix_,                 passTriggerMC);
        tr.registerDerivedVar<bool>("passTriggerHadMC"+myVarSuffix_,              passTriggerHadMC);
        tr.registerDerivedVar<bool>("passNonIsoTriggerMC"+myVarSuffix_,           passNonIsoTriggerMC);
        tr.registerDerivedVar<bool>("passIsoMuTriggerMC"+myVarSuffix_,            passIsoMuTriggerMC);
        tr.registerDerivedVar<bool>("passMadHT"+myVarSuffix_,                     passMadHT);
        tr.registerDerivedVar<bool>("passMETFilters"+myVarSuffix_,                passMETFilters);
        tr.registerDerivedVar<bool>("passElectronHEMveto"+myVarSuffix_,           passElectronHEMveto);
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

    bool PassTriggerAllHad2016(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {
            "HLT_PFJet450",
            //"HLT PFJet500",
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
            //"HLT_PFJet550",
            "HLT_PFHT1050",
            "HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2",
            "HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2",
            "HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5",
            // -- newer ones based on Semra's study
            "HLT_PFHT380_SixJet32_DoubleBTagCSV_p075",
            "HLT_PFHT430_SixJet40_BTagCSV_p080",
        };
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerAllHad2018(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass) 
    {
        std::vector<std::string> mytriggers = {
            "HLT_PFJet500",
            //"HLT_PFJet550",
            "HLT_PFHT1050",
            "HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2",
            "HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2",
            "HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5",
            // -- newer ones based on Semra's study
            "HLT_PFHT430_SixPFJet40_PFBTagDeepCSV_1p5",
            "HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94",
            "HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59",
        };
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

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
        // Preliminary list for 2018 
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
        // Preliminary list for 2018 
        std::vector<std::string> mytriggers = {"HLT_Ele35_WPTight","HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon200"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerPhoton(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Photon175_v"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }
    
    bool PassTriggerNonIsoMuon2016(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Mu50", "HLT_TkMu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }
    
    bool PassTriggerNonIsoMuon2017(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        //No track based muon 50 trigger in 2017 in our tuples
        std::vector<std::string> mytriggers = {"HLT_Mu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerNonIsoMuon2018(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        //No track based muon 50 trigger in 2018 in our tuples
        std::vector<std::string> mytriggers = {"HLT_Mu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    // -------------------------------------------------------------
    // for 0-lepton trigger study, we will also include muon trigger
    // -------------------------------------------------------------
    bool PassTriggerIsoMu2016(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_IsoMu24"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerIsoMu2017(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_IsoMu24"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerIsoMu2018(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_IsoMu24"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    // --------------------------------------------------------------------------------------------------
    // for 0-lepton trigger study, review the AN-2016/411 which is our reference for choosing the tiggers
    // --------------------------------------------------------------------------------------------------
    bool PassTriggerMuonsRefAN(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {   
        std::vector<std::string> mytriggers = {"HLT_IsoMu24", "HLT_IsoMu22_eta2p1"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }    

    bool PassTriggerRefAN(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_PFHT450_SixJet40_BTagCSV_p056", "HLT_PFHT400_SixJet30_DoubleBTagCSV_p056", "HLT PFJet450"}; 
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
