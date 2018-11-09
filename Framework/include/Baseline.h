#ifndef BASELINE_H
#define BASELINE_H

class Baseline
{
private:
    std::string myVarSuffix_;

    void baseline(NTupleReader& tr)
    {
        const auto& runtype             = tr.getVar<std::string>("runtype");     
        const auto& filetag             = tr.getVar<std::string>("filetag");
        const auto& blind               = tr.getVar<bool>("blind");
        const auto& TriggerNames        = tr.getVec<std::string>("TriggerNames");
        const auto& TriggerPass         = tr.getVec<int>("TriggerPass");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons"+myVarSuffix_);
        const auto& NGoodMuons          = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons      = tr.getVar<int>("NGoodElectrons");
        const auto& NJets_pt45          = tr.getVar<int>("NJets_pt45"+myVarSuffix_);
        const auto& HT_trigger          = tr.getVar<double>("HT_trigger"+myVarSuffix_);
        const auto& HT_trigger_pt30     = tr.getVar<double>("HT_trigger_pt30"+myVarSuffix_);
        const auto& NBJets_pt45         = tr.getVar<int>("NBJets_pt45"+myVarSuffix_);
        const auto& NBJets_pt30         = tr.getVar<int>("NBJets_pt30"+myVarSuffix_); 
        const auto& NJets_pt30          = tr.getVar<int>("NJets_pt30"+myVarSuffix_); 
        const auto& onZ                 = tr.getVar<bool>("onZ"+myVarSuffix_); 
        const auto& JetID               = tr.getVar<bool>("JetID"+myVarSuffix_); 
        const auto& NGoodJets_pt45      = tr.getVar<int>("NGoodJets_pt45"+myVarSuffix_);
        const auto& NGoodBJets_pt45     = tr.getVar<int>("NGoodBJets_pt45"+myVarSuffix_);
        const auto& NGoodBJets_pt30     = tr.getVar<int>("NGoodBJets_pt30"+myVarSuffix_); 
        const auto& NGoodJets_pt30      = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_); 
        const auto& NGoodPhotons        = tr.getVar<int>("NGoodPhotons");
        
        // ------------------------------
        // -- Data dependent stuff
        // ------------------------------
        
        bool passTriggerAllHad   = PassTriggerAllHad2016(TriggerNames, TriggerPass) || PassTriggerAllHad2017(TriggerNames, TriggerPass);
        bool passTriggerMuon     = PassTriggerMuon(TriggerNames, TriggerPass);
        bool passTriggerElectron = PassTriggerElectron(TriggerNames, TriggerPass);
        bool passTriggerPhoton   = PassTriggerPhoton(TriggerNames, TriggerPass);

        bool passTrigger   = true;
        bool passTriggerMC = true;
        bool passBlindHad  = true;
        bool passBlindLep  = true;
        
        bool passBlindHad_Good = true;
        bool passBlindLep_Good = true;
        
        if (runtype == "Data")
        {
            // Pass the right trigger
            if (filetag == "Data_JetHT" && !passTriggerAllHad) passTrigger = false;
            if (filetag == "Data_SingleMuon" && !passTriggerMuon) passTrigger = false;
            if (filetag == "Data_SingleElectron" && !passTriggerElectron) passTrigger = false;
            if (filetag == "Data_SinglePhoton" && !passTriggerPhoton) passTrigger = false;

            // Blinding data 
            if (NJets_pt30 >= 9 && blind) passBlindHad = false;
            if (NJets_pt30 >= 8 && blind) passBlindLep = false;
            
            if (NGoodJets_pt30 >= 9 && blind) passBlindHad_Good = false;
            if (NGoodJets_pt30 >= 8 && blind) passBlindLep_Good = false;
        }
        
        // ------------------------
        // -- MC dependent stuff - 
        // -----------------------
        bool passMadHT = true;
        if(runtype == "MC")
        {
            const auto& madHT  = tr.getVar<double>("madHT");

            // Exclude events with MadGraph HT > 100 from the DY & WJets inclusive samples
            if(filetag == "DYJetsToLL_M-50_Incl" && madHT > 100) passMadHT = false;
            if(filetag == "WJetsToLNu_Incl" && madHT > 100) passMadHT = false;

            // Stitch TTbar samples together
            // remove HT overlap
            if( (filetag == "TTJets_Incl" || filetag == "TTJets_SingleLeptFromT" || filetag == "TTJets_SingleLeptFromTbar" || filetag == "TTJets_DiLept") 
                && madHT > 600) passMadHT = false;
            // also remove lepton overlap from the inclusive sample
            const auto& GenElectrons        = tr.getVec<TLorentzVector>("GenElectrons");
            const auto& GenMuons            = tr.getVec<TLorentzVector>("GenMuons");
            const auto& GenTaus             = tr.getVec<TLorentzVector>("GenTaus");
            int NGenLeptons = GenElectrons.size() + GenMuons.size() + GenTaus.size();
            if (filetag == "TTJets_Incl" && NGenLeptons > 0) passMadHT = false;
            
            // MC modeling of the trigger
            if( !passTriggerMuon && !passTriggerElectron ) passTriggerMC = false;
        }

        
        // -------------------------------
        // -- Define 0 Lepton Baseline
        // -------------------------------
        
        bool passBaseline0l = JetID              &&
                              passMadHT          &&
                              passTrigger        &&
                              passTriggerMC      &&
                              passBlindHad       &&
                              NGoodLeptons == 0  && 
                              NJets_pt45 >= 7    && 
                              HT_trigger > 500   && 
                              NBJets_pt45 >= 2;
        
        bool passBaseline0l_Good = JetID          &&
                              passMadHT           &&
                              passTrigger         &&
                              passTriggerMC       &&
                              passBlindHad_Good   &&
                              NGoodLeptons == 0   && 
                              NGoodJets_pt45 >= 7 && 
                              HT_trigger > 500    && 
                              NGoodBJets_pt45 >= 2;

        bool passBaseline0l_hadTrig = JetID      &&
                              passMadHT          &&
                              passTrigger        &&
                              passTriggerMC      &&
                              passBlindHad       &&
                              NJets_pt45 >= 7    && 
                              HT_trigger > 500   && 
                              NBJets_pt45 >= 2;
        
        bool passBaseline0l_hadTrig_Good = JetID  &&
                              passMadHT           &&
                              passTrigger         &&
                              passTriggerMC       &&
                              passBlindHad_Good   &&
                              NGoodJets_pt45 >= 7 && 
                              HT_trigger > 500    && 
                              NGoodBJets_pt45 >= 2;
        

        // -------------------------------
        // -- Define 1 Lepton Baseline
        // -------------------------------
        
        bool passBaseline1mu = JetID                 &&
                               HT_trigger_pt30 > 300 &&
                               passMadHT             &&
                               passTrigger           &&
                               passTriggerMC         &&            
                               (runtype != "Data" || filetag == "Data_SingleMuon") &&
                               passBlindLep          &&
                               NGoodMuons == 1       && 
                               NGoodElectrons == 0   &&
                               NJets_pt30 >= 7       && 
                               NBJets_pt30 >= 1;

        bool passBaseline1el = JetID                 &&
                               HT_trigger_pt30 > 300 &&
                               passMadHT             &&
                               passTrigger           &&
                               passTriggerMC         &&
                               (runtype != "Data" || filetag == "Data_SingleElectron") &&
                               passBlindLep          &&
                               NGoodElectrons == 1   &&
                               NGoodMuons == 0       &&
                               NJets_pt30 >= 7       && 
                               NBJets_pt30 >= 1;

        bool passBaseline1l = passBaseline1mu || passBaseline1el;
        
        bool passBaseline1mu_Good = JetID            &&
                               HT_trigger_pt30 > 300 &&
                               passMadHT             &&
                               passTrigger           &&
                               passTriggerMC         &&
                               (runtype != "Data" || filetag == "Data_SingleMuon") &&
                               passBlindLep_Good     &&
                               NGoodMuons == 1       && 
                               NGoodElectrons == 0   &&
                               NGoodJets_pt30 >= 7   && 
                               NGoodBJets_pt30 >= 1;

        bool passBaseline1el_Good = JetID            &&
                               HT_trigger_pt30 > 300 &&
                               passMadHT             &&
                               passTrigger           &&
                               passTriggerMC         &&
                               (runtype != "Data" || filetag == "Data_SingleElectron") &&
                               passBlindLep_Good     &&
                               NGoodElectrons == 1   &&
                               NGoodMuons == 0       &&
                               NGoodJets_pt30 >= 7   && 
                               NGoodBJets_pt30 >= 1;

        bool passBaseline1l_Good = passBaseline1mu_Good || passBaseline1el_Good;

        // ----------------------------------
        // -- Define 2 Lepton onZ Baseline
        // ----------------------------------
        
        bool passBaseline2lonZ = JetID              &&
                                 passMadHT          &&
                                 passTrigger        &&
                                 passTriggerMC      &&
                                 NGoodLeptons == 2  && 
                                 onZ                &&
                                 (runtype != "Data" || (NGoodMuons == 2 && filetag == "Data_SingleMuon" ) 
                                                    || (NGoodElectrons == 2 && filetag == "Data_SingleElectron") ) &&
                                 NJets_pt30 >= 7; 
        
        bool passBaseline2lonZ_Good = JetID         &&
                                 passMadHT          &&
                                 passTrigger        &&
                                 passTriggerMC      &&
                                 NGoodLeptons == 2  && 
                                 onZ                &&
                                 (runtype != "Data" || (NGoodMuons == 2 && filetag == "Data_SingleMuon" ) 
                                                    || (NGoodElectrons == 2 && filetag == "Data_SingleElectron") ) &&
                                 NGoodJets_pt30 >= 7; 
        
        // -----------------------------------
        // -- Define 2 Lepton offZ Baseline
        // -----------------------------------
        
        bool passBaseline2l = JetID              &&
                              passMadHT          &&
                              passTrigger        &&
                              passTriggerMC      &&
                              passBlindLep       &&
                              NGoodLeptons == 2  && 
                              !onZ               &&
                              (runtype != "Data" || (NGoodMuons >= 1 && filetag == "Data_SingleMuon" ) 
                                                 || (NGoodElectrons == 2 && filetag == "Data_SingleElectron") ) &&
                              NJets_pt30 >= 7    && 
                              NBJets_pt30 >= 1;
        
        bool passBaseline2l_Good = JetID          &&
                              passMadHT           &&
                              passTrigger         &&
                              passTriggerMC       &&
                              passBlindLep_Good   &&
                              NGoodLeptons == 2   && 
                              !onZ                &&
                              (runtype != "Data"  || (NGoodMuons >= 1 && filetag == "Data_SingleMuon" ) 
                                                  || (NGoodElectrons == 2 && filetag == "Data_SingleElectron") ) &&
                              NGoodJets_pt30 >= 7 && 
                              NGoodBJets_pt30 >= 1;
        
        // -----------------------------------
        // -- Define 1 Photon Baseline
        // -----------------------------------

        bool passBaseline1photon_Good = passMadHT           &&
                                        passTrigger         &&
                                        passTriggerMC       &&
                                        NGoodPhotons == 1   &&
                                        NGoodLeptons == 0   && 
                                        NGoodJets_pt30 >= 7; 
        
        tr.registerDerivedVar<bool>("passBaseline0l"+myVarSuffix_,              passBaseline0l);
        tr.registerDerivedVar<bool>("passBaseline0l_Good"+myVarSuffix_,         passBaseline0l_Good);
        tr.registerDerivedVar<bool>("passBaseline0l_hadTrig"+myVarSuffix_,      passBaseline0l_hadTrig);
        tr.registerDerivedVar<bool>("passBaseline0l_hadTrig_Good"+myVarSuffix_, passBaseline0l_hadTrig_Good);
        tr.registerDerivedVar<bool>("passBaseline1l"+myVarSuffix_,              passBaseline1l);
        tr.registerDerivedVar<bool>("passBaseline1l_Good"+myVarSuffix_,         passBaseline1l_Good);
        tr.registerDerivedVar<bool>("passBaseline1mu"+myVarSuffix_,             passBaseline1mu);
        tr.registerDerivedVar<bool>("passBaseline1el"+myVarSuffix_,             passBaseline1el);
        tr.registerDerivedVar<bool>("passBaseline2lonZ"+myVarSuffix_,           passBaseline2lonZ);
        tr.registerDerivedVar<bool>("passBaseline2lonZ_Good"+myVarSuffix_,      passBaseline2lonZ_Good);
        tr.registerDerivedVar<bool>("passBaseline2l"+myVarSuffix_,              passBaseline2l);
        tr.registerDerivedVar<bool>("passBaseline2l_Good"+myVarSuffix_,         passBaseline2l_Good);
        tr.registerDerivedVar<bool>("passBaseline1photon_Good"+myVarSuffix_,    passBaseline1photon_Good);
        tr.registerDerivedVar<bool>("passBlindHad"+myVarSuffix_,                passBlindHad);
        tr.registerDerivedVar<bool>("passBlindLep"+myVarSuffix_,                passBlindLep);
        tr.registerDerivedVar<bool>("passBlindHad_Good"+myVarSuffix_,           passBlindHad_Good);
        tr.registerDerivedVar<bool>("passBlindLep_Good"+myVarSuffix_,           passBlindLep_Good);
        tr.registerDerivedVar<bool>("passTriggerAllHad"+myVarSuffix_,           passTriggerAllHad);
        tr.registerDerivedVar<bool>("passTriggerMuon"+myVarSuffix_,             passTriggerMuon);
        tr.registerDerivedVar<bool>("passTriggerElectron"+myVarSuffix_,         passTriggerElectron);
        tr.registerDerivedVar<bool>("passTrigger"+myVarSuffix_,                 passTrigger);
        tr.registerDerivedVar<bool>("passTriggerMC"+myVarSuffix_,               passTriggerMC);
        tr.registerDerivedVar<bool>("passMadHT"+myVarSuffix_,                   passMadHT);
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
            //"HLT_PFHT900",
            "HLT_PFHT450_SixJet40_BTagCSV",
            "HLT_PFHT400_SixJet30_DoubleBTagCSV",            
        };
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerAllHad2017(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {
            //"HLT_PFHT1050",
            "HLT_PFHT380_SixPFJet32_DoublePFBTagCSV",
            "HLT_PFHT430_SixPFJet40_PFBTagCSV",
        };
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerMuon(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_IsoMu24","HLT_IsoTkMu24_v"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerElectron(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Ele27_WPTight_Gsf"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerPhoton(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Photon175_v"};
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
