#ifndef JET_H
#define JET_H

class Jet
{
private:
    std::string myVarSuffix_;
    const float puIDcut_medium[2][4][4]= {
        {{ 0.20,-0.56,-0.43,-0.38},
        { 0.62,-0.39,-0.32,-0.29},
        { 0.86,-0.10,-0.15,-0.08},
        { 0.93,0.19,0.04,0.12}},
        {{0.26,-0.33,-0.54,-0.37},
        {0.68,-0.04,-0.43,-0.30},
        {0.90,0.36,-0.16,-0.09},
        {0.96,0.61,0.14,0.12}}
    };
    const float puIDcut_tight[2][4][4] = {
        {{ 0.71,-0.32,-0.30,-0.22},
        { 0.87,-0.08,-0.16,-0.12},
        { 0.94,-0.24,0.05,0.10},
        { 0.97,0.48,0.26,0.29}},
        {{0.77,0.38,-0.31,-0.21},
        {0.90,0.60,-0.12,-0.13},
        {0.96,0.82,0.20,0.09},
        {0.98,0.92,0.47,0.29}}
    };

    bool passPuID(float pt, float eta, float disc_value, std::string runYear, const float cuts[2][4][4]){
        int yearidx = 0;
        if (runYear.find("2016") == std::string::npos) yearidx = 1;
        int ptidx = 0;
        int etaidx = 0;
        if( pt < 20) ptidx = 0;
        else if( pt < 30) ptidx = 1;
        else if( pt < 40) ptidx = 2;
        else if( pt < 50) ptidx = 3;
        if( std::abs(eta) < 2.5f) etaidx = 0;
        else if( std::abs(eta) < 2.75f) etaidx = 1;
        else if( std::abs(eta) < 3.0f) etaidx = 2;
        else if( std::abs(eta) < 5.0f) etaidx = 3;
        return cuts[yearidx][ptidx][etaidx] < disc_value || pt > 50.0f;
    }

    void setJetVar(bool pass, std::vector<bool>& boolVec, int& n)
    {
        if( pass )
        {
            boolVec.push_back(true);                
            n++;
        }
        else
        {
            boolVec.push_back(false);
        }            
    }

    void jet(NTupleReader& tr)
    {
        const auto& Jets          = tr.getVec<utility::LorentzVector>(("Jets"+myVarSuffix_));
        const auto& etaCut        = tr.getVar<double>("etaCut");
        const auto& Muons         = tr.getVec<utility::LorentzVector>("Muons");
        const auto& GoodMuons     = tr.getVec<bool>("GoodMuons"+myVarSuffix_);
        const auto& NMuons        = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        const auto& Electrons     = tr.getVec<utility::LorentzVector>("Electrons");
        const auto& GoodElectrons = tr.getVec<bool>("GoodElectrons"+myVarSuffix_);
        const auto& NElectrons    = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);
        const auto& NonIsoMuons   = tr.getVec<bool>("NonIsoMuons"+myVarSuffix_);
        const auto& puId          = tr.getVec<float>("Jets"+myVarSuffix_+"_pileupId");
        const auto& runYear       = tr.getVar<std::string>("runYear");

        //Adding code to create a vector of GoodJets -> defined as the jet collection that eliminates the closest jet to any good lepton (muon or electron) 
        //if that delta R is less than 0.4 and the pT of the jet and lepton is approximately the same
        auto tempGoodJets       = std::make_unique<std::vector<bool>>(Jets.size(), true);
        auto tempNonIsoMuonJets = std::make_unique<std::vector<bool>>(Jets.size(), true);

        if( NMuons > 0 ) 
        {
            for(unsigned int imu = 0; imu < Muons.size(); ++imu)
            {            
                utility::LorentzVector myMuon = Muons.at(imu);
                double         tempDeltaR = 10.0;
                int            tempJetIt  = -1;                
                for(unsigned int myJetIt = 0; myJetIt < Jets.size(); ++myJetIt ) 
                {                    
                    utility::LorentzVector myJet = Jets.at(myJetIt);
                    //First check pT matching between Jet and Muon
                    if( std::fabs( myJet.Pt() - myMuon.Pt() ) / myMuon.Pt() > 1.0 ) continue; 
                    double jetDeltaR = utility::DeltaR(myMuon, myJet);                    
                    if( jetDeltaR < tempDeltaR ) 
                    {
                        tempDeltaR = jetDeltaR;
                        tempJetIt  = myJetIt;
                    }                    
                }//END of looping through jets

                if( tempDeltaR < 0.4 && tempJetIt != -1) 
                {
                    if(GoodMuons[imu])                     tempGoodJets->at(tempJetIt)       = false;
                    if(GoodMuons[imu] || NonIsoMuons[imu]) tempNonIsoMuonJets->at(tempJetIt) = false;
                }
            }//END of looping through muons
        }//END of NMuons if statement

        if( NElectrons > 0 ) 
        {
            for(unsigned int iel = 0; iel < Electrons.size(); ++iel)
            {
                if(!GoodElectrons[iel]) continue;
                utility::LorentzVector myElectron = Electrons.at(iel);
                double         tempDeltaR = 10.0;
                int            tempJetIt  = -1;                
                for(unsigned int myJetIt = 0; myJetIt < Jets.size(); ++myJetIt ) 
                {               
                    utility::LorentzVector myJet = Jets.at(myJetIt);
                    //Check pT matching between Jet and Electron
                    if( std::fabs( myJet.Pt() - myElectron.Pt() ) / myElectron.Pt() > 1.0 ) continue;
                    double jetDeltaR = utility::DeltaR(myElectron, myJet);                    
                    if( jetDeltaR < tempDeltaR ) 
                    {
                        tempDeltaR = jetDeltaR;
                        tempJetIt  = myJetIt;
                    }
                }//END of looping through jets                
                if( tempDeltaR < 0.4 && tempJetIt != -1) 
                {
                    tempGoodJets->at(tempJetIt)       = false;
                    tempNonIsoMuonJets->at(tempJetIt) = false;
                }
            }//END of looping through electrons
        }//END of NElectrons if statement

        auto& jets_pt20_ = tr.createDerivedVec<bool>("Jets_pt20"+myVarSuffix_);
        auto& jets_pt30_ = tr.createDerivedVec<bool>("Jets_pt30"+myVarSuffix_);
        auto& jets_puID_ = tr.createDerivedVec<bool>("Jets_pt30_pileupId"+myVarSuffix_);
        auto& jets_pt40_ = tr.createDerivedVec<bool>("Jets_pt40"+myVarSuffix_);
        auto& jets_pt45_ = tr.createDerivedVec<bool>("Jets_pt45"+myVarSuffix_);
        int NJets_pt20 = 0, NJets_pt30 = 0, NJets_pt40 = 0, NJets_pt45 = 0, NJets_puID = 0;

        auto& goodjets_      = tr.createDerivedVec<bool>("GoodJets"+myVarSuffix_);
        auto& goodjets_puID_ = tr.createDerivedVec<bool>("GoodJets_pt30_PuIdMedium"+myVarSuffix_);
        auto& goodjets_pt20_ = tr.createDerivedVec<bool>("GoodJets_pt20"+myVarSuffix_);
        auto& goodjets_pt30_ = tr.createDerivedVec<bool>("GoodJets_pt30"+myVarSuffix_);
        auto& goodjets_pt40_ = tr.createDerivedVec<bool>("GoodJets_pt40"+myVarSuffix_);
        auto& goodjets_pt45_ = tr.createDerivedVec<bool>("GoodJets_pt45"+myVarSuffix_);
        int NGoodJets = 0, NGoodJetsPuIdMedium=0, NGoodJets_pt20 = 0, NGoodJets_pt30 = 0, NGoodJets_pt40 = 0, NGoodJets_pt45 = 0;

        auto& nonIsoMuonjets_      = tr.createDerivedVec<bool>("NonIsoMuonJets"+myVarSuffix_);
        auto& nonIsoMuonjets_pt20_ = tr.createDerivedVec<bool>("NonIsoMuonJets_pt20"+myVarSuffix_);
        auto& nonIsoMuonjets_pt30_ = tr.createDerivedVec<bool>("NonIsoMuonJets_pt30"+myVarSuffix_);
        auto& nonIsoMuonjets_pt40_ = tr.createDerivedVec<bool>("NonIsoMuonJets_pt40"+myVarSuffix_);
        auto& nonIsoMuonjets_pt45_ = tr.createDerivedVec<bool>("NonIsoMuonJets_pt45"+myVarSuffix_);

        int NNonIsoMuonJets = 0, NNonIsoMuonJets_pt20 = 0, NNonIsoMuonJets_pt30 = 0, NNonIsoMuonJets_pt40 = 0, NNonIsoMuonJets_pt45 = 0;

        for(unsigned int i = 0; i < Jets.size(); ++i ) 
        {               
            utility::LorentzVector lv = Jets.at(i);
         
            bool puID_medium = passPuID(lv.Pt(), lv.Eta(), puId.at(i), runYear, puIDcut_medium);
            bool puID_tight  = passPuID(lv.Pt(), lv.Eta(), puId.at(i), runYear, puIDcut_tight);
  
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 20, jets_pt20_, NJets_pt20); 
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 30, jets_pt30_, NJets_pt30);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 30 && puID_medium, jets_puID_, NJets_puID);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 40, jets_pt40_, NJets_pt40);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 45, jets_pt45_, NJets_pt45);

            setJetVar( abs(lv.Eta()) < etaCut &&             tempGoodJets->at(i), goodjets_,      NGoodJets     );
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 20 && goodjets_.at(i), goodjets_pt20_, NGoodJets_pt20);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 30 && goodjets_.at(i), goodjets_pt30_, NGoodJets_pt30);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 40 && goodjets_.at(i), goodjets_pt40_, NGoodJets_pt40);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 45 && goodjets_.at(i), goodjets_pt45_, NGoodJets_pt45);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 10 && puID_medium && goodjets_.at(i),
                    goodjets_puID_, NGoodJetsPuIdMedium );

            setJetVar( abs(lv.Eta()) < etaCut &&             tempNonIsoMuonJets->at(i), nonIsoMuonjets_,      NNonIsoMuonJets     );
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 20 && nonIsoMuonjets_.at(i), nonIsoMuonjets_pt20_, NNonIsoMuonJets_pt20);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 30 && nonIsoMuonjets_.at(i), nonIsoMuonjets_pt30_, NNonIsoMuonJets_pt30);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 40 && nonIsoMuonjets_.at(i), nonIsoMuonjets_pt40_, NNonIsoMuonJets_pt40);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 45 && nonIsoMuonjets_.at(i), nonIsoMuonjets_pt45_, NNonIsoMuonJets_pt45);
        }

        tr.registerDerivedVar("NJets_pt20"+myVarSuffix_,  NJets_pt20);
        tr.registerDerivedVar("NJets_pt30"+myVarSuffix_,  NJets_pt30);
        tr.registerDerivedVar("NJets_pt40"+myVarSuffix_,  NJets_pt40);
        tr.registerDerivedVar("NJets_pt45"+myVarSuffix_,  NJets_pt45);
        tr.registerDerivedVar("NJets_pt30_pileupId"+myVarSuffix_,  NJets_puID);

        tr.registerDerivedVar("NGoodJets"     +myVarSuffix_, NGoodJets);       
        tr.registerDerivedVar("NGoodJetsPuIdMedium"     +myVarSuffix_, NGoodJetsPuIdMedium);       
        tr.registerDerivedVar("NGoodJets_pt20"+myVarSuffix_, NGoodJets_pt20); 
        tr.registerDerivedVar("NGoodJets_pt20_double"+myVarSuffix_, static_cast<double>(NGoodJets_pt20)); 
        tr.registerDerivedVar("NGoodJets_pt30"+myVarSuffix_, NGoodJets_pt30);
        tr.registerDerivedVar("NGoodJets_pt30_double"+myVarSuffix_, static_cast<double>(NGoodJets_pt30));
        tr.registerDerivedVar("NGoodJets_pt40"+myVarSuffix_, NGoodJets_pt40);
        tr.registerDerivedVar("NGoodJets_pt45"+myVarSuffix_, NGoodJets_pt45);
        tr.registerDerivedVar("NGoodJets_pt45_double"+myVarSuffix_, static_cast<double>(NGoodJets_pt45));

        tr.registerDerivedVar("NNonIsoMuonJets"     +myVarSuffix_, NNonIsoMuonJets);        
        tr.registerDerivedVar("NNonIsoMuonJets_pt20"+myVarSuffix_, NNonIsoMuonJets_pt20);
        tr.registerDerivedVar("NNonIsoMuonJets_pt30"+myVarSuffix_, NNonIsoMuonJets_pt30);
        tr.registerDerivedVar("NNonIsoMuonJets_pt40"+myVarSuffix_, NNonIsoMuonJets_pt40);
        tr.registerDerivedVar("NNonIsoMuonJets_pt45"+myVarSuffix_, NNonIsoMuonJets_pt45);

        int NGoodJets_pt30_inclusive = NGoodJets_pt30 > 12 ? 12 : NGoodJets_pt30;
        tr.registerDerivedVar("NGoodJets_pt30_inclusive"+myVarSuffix_, NGoodJets_pt30_inclusive);
        tr.registerDerivedVar("NGoodJets_pt30_inclusive_shift"+myVarSuffix_, NGoodJets_pt30_inclusive - 7);

        int NNonIsoMuonJets_pt30_inclusive = NNonIsoMuonJets_pt30 > 12 ? 12 : NNonIsoMuonJets_pt30;
        tr.registerDerivedVar("NNonIsoMuonJets_pt30_inclusive"+myVarSuffix_, NNonIsoMuonJets_pt30_inclusive);
        tr.registerDerivedVar("NNonIsoMuonJets_pt30_inclusive_shift"+myVarSuffix_, NNonIsoMuonJets_pt30_inclusive - 7);
        tr.createDerivedVec<bool>("AllJets"+myVarSuffix_, Jets.size(), true);

    }

public:
    Jet(std::string myVarSuffix = "") 
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up Jet"<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        jet(tr);
    }
};

#endif
