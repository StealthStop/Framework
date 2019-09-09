#ifndef JET_H
#define JET_H

class Jet
{
private:
    std::string myVarSuffix_;

    void setJetVar(bool pass, std::vector<bool>* boolVec, int& n)
    {
        if( pass )
        {
            boolVec->push_back(true);                
            n++;
        }
        else
        {
            boolVec->push_back(false);
        }            
    }

    void jet(NTupleReader& tr)
    {
        const auto& Jets          = tr.getVec<TLorentzVector>(("Jets"+myVarSuffix_));
        const auto& etaCut        = tr.getVar<double>("etaCut");

        const auto& Muons         = tr.getVec<TLorentzVector>("Muons");
        const auto& GoodMuons     = tr.getVec<bool>("GoodMuons"+myVarSuffix_);
        const auto& NMuons        = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        const auto& Electrons     = tr.getVec<TLorentzVector>("Electrons");
        const auto& GoodElectrons = tr.getVec<bool>("GoodElectrons"+myVarSuffix_);
        const auto& NElectrons    = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);

        const auto& NonIsoMuons   = tr.getVec<bool>("NonIsoMuons"+myVarSuffix_);
        const auto& NNonIsoMuons  = tr.getVar<int>("NNonIsoMuons"+myVarSuffix_);

        const auto& doQCDCR       = tr.getVar<bool>("doQCDCR");

        //Adding code to create a vector of GoodJets -> defined as the jet collection that eliminates the closest jet to any good lepton (muon or electron) 
        //if that delta R is less than 0.4 and the pT of the jet and lepton is approximately the same

        auto tempGoodJets = std::make_unique<std::vector<bool>>(Jets.size(), true);

        if( NMuons > 0 ) 
        {
            for(unsigned int imu = 0; imu < Muons.size(); ++imu)
            {            
                if( doQCDCR ) {
                    if( !GoodMuons[imu] && !NonIsoMuons[imu] ) continue;
                    
                    TLorentzVector myMuon = Muons.at(imu);
                    double         tempDeltaR = 10.0;
                    int            tempJetIt  = -1;                
                    for( int myJetIt = 0; myJetIt < Jets.size(); ++myJetIt ) 
                    {                    
                        TLorentzVector myJet = Jets.at(myJetIt);
                        //First check pT matching between Jet and Muon
                        if( std::fabs( myJet.Pt() - myMuon.Pt() ) / myMuon.Pt() > 1.0 ) continue; 
                        double jetDeltaR = myMuon.DeltaR(myJet);                    
                        if( jetDeltaR < tempDeltaR ) 
                        {
                            tempDeltaR = jetDeltaR;
                            tempJetIt  = myJetIt;
                        }                    
                    }//END of looping through jets

                    if( tempDeltaR < 0.4 && tempJetIt != -1) 
                    {
                        tempGoodJets->at(tempJetIt) = false;
                    }
                }
                
                else {
                    if( !GoodMuons[imu] ) continue;

                    TLorentzVector myMuon = Muons.at(imu);
                    double         tempDeltaR = 10.0;
                    int            tempJetIt  = -1;                
                    for( int myJetIt = 0; myJetIt < Jets.size(); ++myJetIt ) 
                    {                    
                        TLorentzVector myJet = Jets.at(myJetIt);
                        //First check pT matching between Jet and Muon
                        if( std::fabs( myJet.Pt() - myMuon.Pt() ) / myMuon.Pt() > 1.0 ) continue; 
                        double jetDeltaR = myMuon.DeltaR(myJet);                    
                        if( jetDeltaR < tempDeltaR ) 
                        {
                            tempDeltaR = jetDeltaR;
                            tempJetIt  = myJetIt;
                        }                    
                    }//END of looping through jets
                    if( tempDeltaR < 0.4 && tempJetIt != -1) 
                    {
                        tempGoodJets->at(tempJetIt) = false;
                    }
                }
            }//END of looping through muons
        }//END of NMuons if statement

        if( NElectrons > 0 ) 
        {
            for(unsigned int iel = 0; iel < Electrons.size(); ++iel)
            {
                if(!GoodElectrons[iel]) continue;
                TLorentzVector myElectron = Electrons.at(iel);
                double         tempDeltaR = 10.0;
                int            tempJetIt  = -1;                
                for( int myJetIt = 0; myJetIt < Jets.size(); ++myJetIt ) 
                {               
                    TLorentzVector myJet = Jets.at(myJetIt);
                    //Check pT matching between Jet and Electron
                    if( std::fabs( myJet.Pt() - myElectron.Pt() ) / myElectron.Pt() > 1.0 ) continue;
                    double jetDeltaR = myElectron.DeltaR(myJet);                    
                    if( jetDeltaR < tempDeltaR ) 
                    {
                        tempDeltaR = jetDeltaR;
                        tempJetIt  = myJetIt;
                    }
                }//END of looping through jets                
                if( tempDeltaR < 0.4 && tempJetIt != -1) 
                {
                    tempGoodJets->at(tempJetIt) = false;
                }
            }//END of looping through electrons
        }//END of NElectrons if statement

        auto* jets_pt30_ = new std::vector<bool>();
        auto* jets_pt40_ = new std::vector<bool>();
        auto* jets_pt45_ = new std::vector<bool>();
        int NJets_pt30 = 0, NJets_pt40 = 0, NJets_pt45 = 0;

        auto* goodjets_       = new std::vector<bool>();
        auto* goodjets_pt30_  = new std::vector<bool>();
        auto* goodjets_pt40_  = new std::vector<bool>();
        auto* goodjets_pt45_  = new std::vector<bool>();
        int NGoodJets = 0, NGoodJets_pt30 = 0, NGoodJets_pt40 = 0, NGoodJets_pt45 = 0;

        for( int i = 0; i < Jets.size(); ++i ) 
        {               
            TLorentzVector lv = Jets.at(i);
            
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 30, jets_pt30_, NJets_pt30);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 40, jets_pt40_, NJets_pt40);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 45, jets_pt45_, NJets_pt45);

            setJetVar( abs(lv.Eta()) < etaCut &&                 tempGoodJets->at(i), goodjets_,      NGoodJets     );
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 30 && goodjets_->at(i),    goodjets_pt30_, NGoodJets_pt30);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 40 && goodjets_->at(i),    goodjets_pt40_, NGoodJets_pt40);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 45 && goodjets_->at(i),    goodjets_pt45_, NGoodJets_pt45);
        }

        tr.registerDerivedVec("Jets_pt30" +myVarSuffix_,  jets_pt30_);
        tr.registerDerivedVar("NJets_pt30"+myVarSuffix_,  NJets_pt30);
        tr.registerDerivedVec("Jets_pt40" +myVarSuffix_,  jets_pt40_);
        tr.registerDerivedVar("NJets_pt40"+myVarSuffix_,  NJets_pt40);
        tr.registerDerivedVec("Jets_pt45" +myVarSuffix_,  jets_pt45_);
        tr.registerDerivedVar("NJets_pt45"+myVarSuffix_,  NJets_pt45);
        
        tr.registerDerivedVec("GoodJets"      +myVarSuffix_, goodjets_);
        tr.registerDerivedVar("NGoodJets"     +myVarSuffix_, NGoodJets);
        tr.registerDerivedVec("GoodJets_pt30" +myVarSuffix_, goodjets_pt30_);
        tr.registerDerivedVar("NGoodJets_pt30"+myVarSuffix_, NGoodJets_pt30);
        tr.registerDerivedVec("GoodJets_pt40" +myVarSuffix_, goodjets_pt40_);
        tr.registerDerivedVar("NGoodJets_pt40"+myVarSuffix_, NGoodJets_pt40);
        tr.registerDerivedVec("GoodJets_pt45" +myVarSuffix_, goodjets_pt45_);
        tr.registerDerivedVar("NGoodJets_pt45"+myVarSuffix_, NGoodJets_pt45);

        int NGoodJets_pt30_inclusive = NGoodJets_pt30 > 12 ? 12 : NGoodJets_pt30;
        tr.registerDerivedVar("NGoodJets_pt30_inclusive"+myVarSuffix_, NGoodJets_pt30_inclusive);
        tr.registerDerivedVar("NGoodJets_pt30_inclusive_shift"+myVarSuffix_, NGoodJets_pt30_inclusive - 7);
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
