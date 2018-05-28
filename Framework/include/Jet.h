#ifndef JET_H
#define JET_H

class Jet
{
private:
    std::vector<TLorentzVector>* jets_pt30_;
    std::vector<TLorentzVector>* jets_pt40_;
    std::vector<TLorentzVector>* jets_pt45_;
    std::vector<double>*         jets_pt45_csv_;
    
    std::vector<TLorentzVector>* goodJets_;
    std::vector<TLorentzVector>* goodjets_pt30_;
    std::vector<TLorentzVector>* goodjets_pt40_;
    std::vector<TLorentzVector>* goodjets_pt45_;
    std::vector<double>*         goodjets_CSV_;
        

    void jet(NTupleReader& tr)
    {
        const auto& Jets        = tr.getVec<TLorentzVector>("Jets");
        const auto& etaCut      = tr.getVar<double>("etaCut");
        const auto& Jets_CSV    = tr.getVec<double>("Jets_bDiscriminatorCSV");

        const auto& Muons       = tr.getVec<TLorentzVector>("GoodMuons");
        const auto& NMuons      = tr.getVar<int>("NGoodMuons");
        const auto& Electrons   = tr.getVec<TLorentzVector>("GoodElectrons");
        const auto& NElectrons  = tr.getVar<int>("NGoodElectrons");

        jets_pt30_      = new std::vector<TLorentzVector>();
        jets_pt40_      = new std::vector<TLorentzVector>();
        jets_pt45_      = new std::vector<TLorentzVector>();
        jets_pt45_csv_  = new std::vector<double>();
        
        for (unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
        {
            TLorentzVector lv(Jets.at(ijet));
            if( abs(lv.Eta()) < etaCut && 
                lv.Pt() > 30 
                )
            {
                jets_pt30_->push_back(lv);
                
                if(lv.Pt() > 40) {
                    
                    jets_pt40_->push_back(lv);
                    
                    if(lv.Pt() > 45) {   
                        jets_pt45_->push_back(lv);
                        jets_pt45_csv_->push_back(Jets_CSV.at(ijet));
                    }
                }
            }
        }

        tr.registerDerivedVar("NJets",       Jets.size());
        tr.registerDerivedVec("Jets_pt30",   jets_pt30_);
        tr.registerDerivedVar("NJets_pt30", (jets_pt30_==nullptr)?0:jets_pt30_->size());
        tr.registerDerivedVec("Jets_pt40",   jets_pt40_);
        tr.registerDerivedVar("NJets_pt40", (jets_pt40_==nullptr)?0:jets_pt40_->size());
        tr.registerDerivedVec("Jets_pt45",   jets_pt45_);
        tr.registerDerivedVec("Jets_pt45_CSV",   jets_pt45_csv_);
        tr.registerDerivedVar("NJets_pt45", (jets_pt45_==nullptr)?0:jets_pt45_->size());

        //Adding code to create a vector of GoodJets -> defined as the jet collection that eliminates the closest jet to any good lepton (muon or electron) if that delta R is less than 0.4 and the pT of the jet and lepton is approximately the same
        
        goodJets_       = new std::vector<TLorentzVector>();
        goodjets_pt30_  = new std::vector<TLorentzVector>();
        goodjets_pt40_  = new std::vector<TLorentzVector>();
        goodjets_pt45_  = new std::vector<TLorentzVector>();
        goodjets_CSV_   = new std::vector<double>();

        std::vector<TLorentzVector>   goodjets_      = Jets;
        std::vector<double>           goodjets_csv_  = Jets_CSV;
        std::vector<int>              removedJetsIndices;

        if( NMuons > 0 ) {

            for( TLorentzVector myMuon : Muons ) {
                
                double          tempDeltaR = 10.0;
                int             tempJetIt  = -1;
                
                for( int myJetIt = 0; myJetIt < Jets.size(); ++myJetIt ) {
                    
                    TLorentzVector myJet = Jets.at(myJetIt);
                    //First check pT matching between Jet and Muon
                    if( std::fabs( myJet.Pt() - myMuon.Pt() ) / myMuon.Pt() > 1.0 ) continue; 
                    double jetDeltaR = myMuon.DeltaR(myJet);
                    
                    if( jetDeltaR < tempDeltaR ) {
                        tempDeltaR = jetDeltaR;
                        tempJetIt  = myJetIt;
                    }
                    
                }//END of looping through jets

                if( tempDeltaR < 0.4 ) {
                    if( std::find(removedJetsIndices.begin(), removedJetsIndices.end(), tempJetIt ) == removedJetsIndices.end() )  removedJetsIndices.push_back( tempJetIt );
                }
            }//END of looping through muons
        }//END of NMuons if statement
        
        if( NElectrons > 0 ) {

            for( TLorentzVector myElectron : Electrons ) {
                
                double      tempDeltaR = 10.0;
                int         tempJetIt  = -1;
                
                for( int myJetIt = 0; myJetIt < Jets.size(); ++myJetIt ) {
               
                    TLorentzVector myJet = Jets.at(myJetIt);
                    //Check pT matching between Jet and Electron
                    if( std::fabs( myJet.Pt() - myElectron.Pt() ) / myElectron.Pt() > 1.0 ) continue;
                    double jetDeltaR = myElectron.DeltaR(myJet);
                    
                    if( jetDeltaR < tempDeltaR ) {
                        tempDeltaR = jetDeltaR;
                        tempJetIt  = myJetIt;
                    }
                }//END of looping through jets
                
                if( tempDeltaR < 0.4 ) {
                    if( std::find(removedJetsIndices.begin(), removedJetsIndices.end(), tempJetIt ) == removedJetsIndices.end() )  removedJetsIndices.push_back( tempJetIt );
                }
            }//END of looping through electrons
        }//END of NElectrons if statement
        
        std::sort(removedJetsIndices.begin(), removedJetsIndices.end(), std::greater<>());
        
        for( int i = 0; i < removedJetsIndices.size(); i++ ) {
            goodjets_.erase( goodjets_.begin() + removedJetsIndices[i] );
            goodjets_csv_.erase( goodjets_csv_.begin() + removedJetsIndices[i] );
        }
        
        for( TLorentzVector myGoodJet : goodjets_ ){
            if( myGoodJet.Pt() > 30 && std::fabs( myGoodJet.Eta()) < etaCut ) {
                goodjets_pt30_->push_back(myGoodJet);
                if( myGoodJet.Pt() > 40 && std::fabs( myGoodJet.Eta() ) < etaCut ) {
                    goodjets_pt40_->push_back(myGoodJet);
                    if( myGoodJet.Pt() >  45 && std::fabs( myGoodJet.Eta() ) < etaCut )  {
                        goodjets_pt45_->push_back(myGoodJet);
                    }
                }
            }
        }

        for( double goodCsvValue : goodjets_csv_ ) {
            goodjets_CSV_->push_back( goodCsvValue );
        }
        
        for( TLorentzVector myGoodJet : goodjets_ ) {
            goodJets_->push_back( myGoodJet );
        }
        
        tr.registerDerivedVec("GoodJets",                   goodJets_);
        tr.registerDerivedVar("NGoodJets",                  goodjets_.size());
        tr.registerDerivedVec("GoodJets_bDiscriminatorCSV", goodjets_CSV_);
        tr.registerDerivedVec("GoodJets_pt30",              goodjets_pt30_);
        tr.registerDerivedVar("NGoodJets_pt30",            (goodjets_pt30_==nullptr)?0:goodjets_pt30_->size());
        tr.registerDerivedVec("GoodJets_pt40",              goodjets_pt40_);
        tr.registerDerivedVar("NGoodJets_pt40",            (goodjets_pt40_==nullptr)?0:goodjets_pt40_->size());
        tr.registerDerivedVec("GoodJets_pt45",              goodjets_pt45_);
        tr.registerDerivedVar("NGoodJets_pt45",            (goodjets_pt45_==nullptr)?0:goodjets_pt45_->size());
    }

public:
    Jet() : 
        jets_pt30_(nullptr)
      , jets_pt40_(nullptr)
      , jets_pt45_(nullptr)
      , jets_pt45_csv_(nullptr)
      , goodJets_(nullptr)
      , goodjets_pt30_(nullptr)
      , goodjets_pt40_(nullptr)
      , goodjets_pt45_(nullptr)
      , goodjets_CSV_(nullptr)
    {
    }

    void operator()(NTupleReader& tr)
    {
        jet(tr);
    }
};

#endif
