#ifndef BASELINE_H
#define BASELINE_H

class Baseline
{
private:
    void baseline(NTupleReader& tr)
    {

        const auto& runtype      = tr.getVar<std::string>("runtype");     
        const auto& filetag      = tr.getVar<std::string>("filetag");
        const auto& blind        = tr.getVar<bool>("blind");
        const auto& TriggerNames = tr.getVec<std::string>("TriggerNames");
        const auto& TriggerPass  = tr.getVec<int>("TriggerPass");
        const auto& NGoodLeptons = tr.getVar<int>("NGoodLeptons");
        const auto& NJets_pt45   = tr.getVar<int>("NJets_pt45");
        const auto& HT_trigger   = tr.getVar<double>("HT_trigger");
        const auto& NBJets_pt45  = tr.getVar<int>("NBJets_pt45");
        const auto& NJets_pt30   = tr.getVar<int>("NJets_pt30"); 

        // ------------------------
        // -- MC dependent stuff
        // -----------------------
        bool passMadHT = true;
        if(runtype == "MC"){
            const auto& madHT  = tr.getVar<double>("madHT");
            // Exclude events with MadGraph HT > 100 from the DY inclusive sample
            if(filetag == "DYJetsToLL_M-50_Incl" && madHT > 100) passMadHT = false;
        }
        
        // ------------------------------
        // -- Data dependent stuff
        // ------------------------------
        
        bool passTriggerAllHad   = PassTriggerAllHad(TriggerNames, TriggerPass);
        bool passTriggerMuon     = PassTriggerMuon(TriggerNames, TriggerPass);
        bool passTriggerElectron = PassTriggerElectron(TriggerNames, TriggerPass);
        bool passTrigger = true;
        bool passBlind   = true;
        if (runtype == "Data")
        {
            // Pass the right trigger
            if (filetag == "Data_JetHT" && !passTriggerAllHad) passTrigger = false;
            if (filetag == "Data_SingleMuon" && !passTriggerMuon) passTrigger = false;
            if (filetag == "Data_SingleElectron" && !passTriggerElectron) passTrigger = false;

            // Blinding data 
            if (NJets_pt30 >= 7 && blind) passBlind = false;
        }        
        
        // -------------------------------
        // -- Define 0 Lepton Baseline
        // -------------------------------
        
        bool passBaseline0l = passMadHT          &&
                              passTrigger        &&
                              passBlind          &&
                              NGoodLeptons == 0  && 
                              NJets_pt45 >= 6    && 
                              HT_trigger > 500   && 
                              NBJets_pt45 >= 2;

        tr.registerDerivedVar<bool>("passBaseline0l",passBaseline0l);

    }

    bool PassTriggerGeneral(std::vector<std::string> &mytriggers, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
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

    bool PassTriggerAllHad(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers {
            //"HLT_PFHT1050", // 2017 trigger
            //"HLT_PFHT900"
            //"HLT_PFHT380_SixPFJet32_DoublePFBTagCSV", // 2017 trigger
            //"HLT_PFHT430_SixPFJet40_PFBTagCSV", // 2017 trigger
            "HLT_PFHT450_SixJet40_BTagCSV",
                "HLT_PFHT400_SixJet30_DoubleBTagCSV",            
                };
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerMuon(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers {"HLT_IsoMu24","HLT_IsoTkMu24_v"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerElectron(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers {"HLT_Ele27_WPTight_Gsf"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

public:
    Baseline() 
    {
    }
    
    void operator()(NTupleReader& tr)
    {
        baseline(tr);
    }
};

#endif
