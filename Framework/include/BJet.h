#ifndef BJET_H
#define BJET_H

class BJet
{
private:
    std::vector<TLorentzVector>* bjets_;
    std::vector<TLorentzVector>* bjets_pt30_;
    std::vector<TLorentzVector>* bjets_pt40_;
    std::vector<TLorentzVector>* bjets_pt45_;

    std::vector<TLorentzVector>* bjets_tight_;
    std::vector<TLorentzVector>* bjets_pt30_tight_;
    std::vector<TLorentzVector>* bjets_pt45_tight_;
    
    std::vector<TLorentzVector>* goodbjets_;
    std::vector<TLorentzVector>* goodbjets_pt30_;
    std::vector<TLorentzVector>* goodbjets_pt40_;
    std::vector<TLorentzVector>* goodbjets_pt45_;
    
    std::vector<TLorentzVector>* goodbjets_tight_;
    std::vector<TLorentzVector>* goodbjets_pt30_tight_;
    std::vector<TLorentzVector>* goodbjets_pt45_tight_;
    
    void bjet(NTupleReader& tr)
    {
        const auto& Jets = tr.getVec<TLorentzVector>("Jets");
        const auto& Jets_bDiscriminatorCSV = tr.getVec<double>("Jets_bDiscriminatorCSV");
        const auto& etaCut = tr.getVar<double>("etaCut");
        const auto& JetsID = tr.getVec<bool>("Jets_ID");
        
        bjets_ = new std::vector<TLorentzVector>();
        bjets_pt30_ = new std::vector<TLorentzVector>();
        bjets_pt40_ = new std::vector<TLorentzVector>();
        bjets_pt45_ = new std::vector<TLorentzVector>();

        bjets_tight_ = new std::vector<TLorentzVector>();
        bjets_pt30_tight_ = new std::vector<TLorentzVector>();
        bjets_pt45_tight_ = new std::vector<TLorentzVector>();


        for (unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
        {
            if( !JetsID.at(ijet) ) continue;
            
            TLorentzVector lv(Jets.at(ijet));
            
            if( abs(lv.Eta()) < etaCut &&
                Jets_bDiscriminatorCSV.at(ijet) > 0.8484 
                )
            {
                bjets_->push_back(lv); 
                if( lv.Pt() > 30 )
                    bjets_pt30_->push_back(lv);
                if( lv.Pt() > 40 )
                    bjets_pt40_->push_back(lv);
                if( lv.Pt() > 45 )
                    bjets_pt45_->push_back(lv);
            }
            
            if( abs(lv.Eta()) < etaCut &&
                Jets_bDiscriminatorCSV.at(ijet) > 0.9535 
                )
            {
                bjets_tight_->push_back(lv); 
                if( lv.Pt() > 30 )
                    bjets_pt30_tight_->push_back(lv);
                if( lv.Pt() > 45 )
                    bjets_pt45_tight_->push_back(lv);
            }
        }


        tr.registerDerivedVec("BJets",        bjets_);
        tr.registerDerivedVar("NBJets",      static_cast<int>((bjets_==nullptr)?0:bjets_->size()));
        tr.registerDerivedVec("BJets_pt30",   bjets_pt30_);
        tr.registerDerivedVar("NBJets_pt30", static_cast<int>((bjets_pt30_==nullptr)?0:bjets_pt30_->size()));
        tr.registerDerivedVec("BJets_pt40",   bjets_pt40_);
        tr.registerDerivedVar("NBJets_pt40", static_cast<int>((bjets_pt40_==nullptr)?0:bjets_pt40_->size()));
        tr.registerDerivedVec("BJets_pt45",   bjets_pt45_);
        tr.registerDerivedVar("NBJets_pt45", static_cast<int>((bjets_pt45_==nullptr)?0:bjets_pt45_->size()));

        tr.registerDerivedVec("BJets_tight",        bjets_tight_);
        tr.registerDerivedVar("NBJets_tight",      static_cast<int>((bjets_tight_==nullptr)?0:bjets_tight_->size()));
        tr.registerDerivedVec("BJets_pt30_tight",   bjets_pt30_tight_);
        tr.registerDerivedVar("NBJets_pt30_tight", static_cast<int>((bjets_pt30_tight_==nullptr)?0:bjets_pt30_tight_->size()));
        tr.registerDerivedVec("BJets_pt45_tight",   bjets_pt45_tight_);
        tr.registerDerivedVar("NBJets_pt45_tight", static_cast<int>((bjets_pt45_tight_==nullptr)?0:bjets_pt45_tight_->size()));

        goodbjets_ = new std::vector<TLorentzVector>();
        goodbjets_pt30_ = new std::vector<TLorentzVector>();
        goodbjets_pt40_ = new std::vector<TLorentzVector>();
        goodbjets_pt45_ = new std::vector<TLorentzVector>();

        goodbjets_tight_ = new std::vector<TLorentzVector>();
        goodbjets_pt30_tight_ = new std::vector<TLorentzVector>();
        goodbjets_pt45_tight_ = new std::vector<TLorentzVector>();
        
        const auto& GoodJets        = tr.getVec<TLorentzVector>("GoodJets");
        const auto& GoodJets_CSV    = tr.getVec<double>("GoodJets_bDiscriminatorCSV");

        for( unsigned int igjet = 0; igjet < GoodJets.size(); ++igjet ) {

            TLorentzVector myGoodJet( GoodJets.at(igjet) );

            if( myGoodJet.Eta() < etaCut &&
                GoodJets_CSV.at(igjet) > 0.8484 ) {

                goodbjets_->push_back( myGoodJet );
                if( myGoodJet.Pt() > 30 )
                    goodbjets_pt30_->push_back( myGoodJet );
                if( myGoodJet.Pt() > 40 )
                    goodbjets_pt40_->push_back( myGoodJet );
                if( myGoodJet.Pt() > 45 )
                    goodbjets_pt45_->push_back( myGoodJet );
            }

            if( myGoodJet.Eta() < etaCut &&
                GoodJets_CSV.at(igjet) > 0.9535 ) {
            
                if( myGoodJet.Pt() > 30 )
                    goodbjets_pt30_tight_->push_back( myGoodJet );
                if( myGoodJet.Pt() > 45 )
                    goodbjets_pt45_tight_->push_back( myGoodJet );
            }

        }
        
        tr.registerDerivedVec("GoodBJets",        goodbjets_);
        tr.registerDerivedVar("NGoodBJets",      static_cast<int>((goodbjets_==nullptr)?0:goodbjets_->size()));
        tr.registerDerivedVec("GoodBJets_pt30",   goodbjets_pt30_);
        tr.registerDerivedVar("NGoodBJets_pt30", static_cast<int>((goodbjets_pt30_==nullptr)?0:goodbjets_pt30_->size()));
        tr.registerDerivedVec("GoodBJets_pt40",   goodbjets_pt40_);
        tr.registerDerivedVar("NGoodBJets_pt40", static_cast<int>((goodbjets_pt40_==nullptr)?0:goodbjets_pt40_->size()));
        tr.registerDerivedVec("GoodBJets_pt45",   goodbjets_pt45_);
        tr.registerDerivedVar("NGoodBJets_pt45", static_cast<int>((goodbjets_pt45_==nullptr)?0:goodbjets_pt45_->size()));
        
        tr.registerDerivedVec("GoodBJets_pt30_tight",   goodbjets_pt30_tight_);
        tr.registerDerivedVar("NGoodBJets_pt30_tight", static_cast<int>((goodbjets_pt30_tight_==nullptr)?0:goodbjets_pt30_tight_->size()));
        tr.registerDerivedVec("GoodBJets_pt45_tight",   goodbjets_pt45_tight_);
        tr.registerDerivedVar("NGoodBJets_pt45_tight", static_cast<int>((goodbjets_pt45_tight_==nullptr)?0:goodbjets_pt45_tight_->size()));

    }

public:
    BJet() 
        : bjets_(nullptr)
        , bjets_pt30_(nullptr)
        , bjets_pt40_(nullptr)
        , bjets_pt45_(nullptr)
        , bjets_tight_(nullptr)
        , bjets_pt30_tight_(nullptr)
        , bjets_pt45_tight_(nullptr)
        , goodbjets_(nullptr)
        , goodbjets_pt30_(nullptr)
        , goodbjets_pt40_(nullptr)
        , goodbjets_pt45_(nullptr)
        , goodbjets_pt30_tight_(nullptr)
        , goodbjets_pt45_tight_(nullptr)
    {}

    void operator()(NTupleReader& tr)
    {
        bjet(tr);
    }
};

#endif
