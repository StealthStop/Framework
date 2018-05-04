#ifndef BJET_H
#define BJET_H

class BJet
{
private:
    std::vector<TLorentzVector>* bjets_;
    std::vector<TLorentzVector>* bjets_pt30_;
    std::vector<TLorentzVector>* bjets_pt45_;

    std::vector<TLorentzVector>* bjets_tight_;
    std::vector<TLorentzVector>* bjets_pt30_tight_;
    std::vector<TLorentzVector>* bjets_pt45_tight_;
    
    std::vector<TLorentzVector>* bjets_clean_;
    std::vector<TLorentzVector>* bjets_pt30_clean_;
    std::vector<TLorentzVector>* bjets_pt45_clean_;

    std::vector<TLorentzVector>* bjets_tight_clean_;
    std::vector<TLorentzVector>* bjets_pt30_tight_clean_;
    std::vector<TLorentzVector>* bjets_pt45_tight_clean_;

    void bjet(NTupleReader& tr)
    {
        const auto& Jets = tr.getVec<TLorentzVector>("Jets");
        const auto& Jets_bDiscriminatorCSV = tr.getVec<double>("Jets_bDiscriminatorCSV");
        const auto& etaCut = tr.getVar<double>("etaCut");
        const auto& JetsID = tr.getVec<bool>("Jets_ID");
        const auto& cleanJets = tr.getVec<TLorentzVector>("Jetsclean");
        const auto& cleanJetsCSV = tr.getVec<double>("Jetsclean_bDiscriminatorCSV");
        const auto& cleanJetsID = tr.getVec<bool>("Jetsclean_ID");
        
        bjets_ = new std::vector<TLorentzVector>();
        bjets_pt30_ = new std::vector<TLorentzVector>();
        bjets_pt45_ = new std::vector<TLorentzVector>();

        bjets_tight_ = new std::vector<TLorentzVector>();
        bjets_pt30_tight_ = new std::vector<TLorentzVector>();
        bjets_pt45_tight_ = new std::vector<TLorentzVector>();

        bjets_clean_ = new std::vector<TLorentzVector>();
        bjets_pt30_clean_ = new std::vector<TLorentzVector>();
        bjets_pt45_clean_ = new std::vector<TLorentzVector>();

        bjets_tight_clean_ = new std::vector<TLorentzVector>();
        bjets_pt30_tight_clean_ = new std::vector<TLorentzVector>();
        bjets_pt45_tight_clean_ = new std::vector<TLorentzVector>();

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

        if( cleanJets.size() != 0 ) {
            for( unsigned int icjet = 0; icjet < cleanJets.size(); ++icjet ) 
            {
                if( !cleanJetsID.at(icjet) ) continue;
                TLorentzVector lv( cleanJets.at(icjet) );
            
                if( abs( lv.Eta() ) < etaCut &&
                    cleanJetsCSV.at(icjet) > 0.8484 
                    )
                {
                    bjets_clean_->push_back(lv);
                    if( lv.Pt() > 30 )
                        bjets_pt30_clean_->push_back(lv);
                    if( lv.Pt() > 45 )
                        bjets_pt45_clean_->push_back(lv);
                }

                if( abs( lv.Eta() ) < etaCut &&
                    cleanJetsCSV.at(icjet) > 0.9535
                    )
                {
                    bjets_tight_clean_->push_back(lv);
                    if( lv.Pt() > 30 )
                        bjets_pt30_tight_clean_->push_back(lv);
                    if( lv.Pt() > 45 )
                        bjets_pt45_tight_clean_->push_back(lv);
                }
            }
        }
        else if( cleanJets.size() == 0 ) {
            for( unsigned int icjet = 0; icjet < Jets.size(); ++icjet ) 
            {
                if( !jetsID.at(icjet) ) continue;
                TLorentzVector lv( Jets.at(icjet) );
            
                if( abs( lv.Eta() ) < etaCut &&
                    Jets_bDiscriminatorCSV.at(icjet) > 0.8484 
                    )
                {
                    bjets_clean_->push_back(lv);
                    if( lv.Pt() > 30 )
                        bjets_pt30_clean_->push_back(lv);
                    if( lv.Pt() > 45 )
                        bjets_pt45_clean_->push_back(lv);
                }

                if( abs( lv.Eta() ) < etaCut &&
                    Jets_bDiscriminatorCSV.at(icjet) > 0.9535
                    )
                {
                    bjets_tight_clean_->push_back(lv);
                    if( lv.Pt() > 30 )
                        bjets_pt30_tight_clean_->push_back(lv);
                    if( lv.Pt() > 45 )
                        bjets_pt45_tight_clean_->push_back(lv);
                }
            }   
        }

        tr.registerDerivedVec("BJets",        bjets_);
        tr.registerDerivedVar("NBJets",      (bjets_==nullptr)?0:bjets_->size());
        tr.registerDerivedVec("BJets_pt30",   bjets_pt30_);
        tr.registerDerivedVar("NBJets_pt30", (bjets_pt30_==nullptr)?0:bjets_pt30_->size());
        tr.registerDerivedVec("BJets_pt45",   bjets_pt45_);
        tr.registerDerivedVar("NBJets_pt45", (bjets_pt45_==nullptr)?0:bjets_pt45_->size());

        tr.registerDerivedVec("BJets_tight",        bjets_tight_);
        tr.registerDerivedVar("NBJets_tight",      (bjets_tight_==nullptr)?0:bjets_tight_->size());
        tr.registerDerivedVec("BJets_pt30_tight",   bjets_pt30_tight_);
        tr.registerDerivedVar("NBJets_pt30_tight", (bjets_pt30_tight_==nullptr)?0:bjets_pt30_tight_->size());
        tr.registerDerivedVec("BJets_pt45_tight",   bjets_pt45_tight_);
        tr.registerDerivedVar("NBJets_pt45_tight", (bjets_pt45_tight_==nullptr)?0:bjets_pt45_tight_->size());

        tr.registerDerivedVec("BJets_clean",        bjets_clean_);
        tr.registerDerivedVar("NBJets_clean",      (bjets_clean_==nullptr)?0:bjets_clean_->size());
        tr.registerDerivedVec("BJets_pt30_clean",   bjets_pt30_clean_);
        tr.registerDerivedVar("NBJets_pt30_clean", (bjets_pt30_clean_==nullptr)?0:bjets_pt30_clean_->size());
        tr.registerDerivedVec("BJets_pt45_clean",   bjets_pt45_clean_);
        tr.registerDerivedVar("NBJets_pt45_clean", (bjets_pt45_clean_==nullptr)?0:bjets_pt45_clean_->size());

        tr.registerDerivedVec("BJets_tight_clean",          bjets_tight_clean_);
        tr.registerDerivedVar("NBJets_tight_clean",        (bjets_tight_clean_==nullptr)?0:bjets_tight_clean_->size());
        tr.registerDerivedVec("BJets_pt30_tight_clean",     bjets_pt30_tight_clean_);
        tr.registerDerivedVar("NBJets_pt30_tight_clean",   (bjets_pt30_tight_clean_==nullptr)?0:bjets_pt30_tight_clean_->size());
        tr.registerDerivedVec("BJets_pt45_tight_clean",     bjets_pt45_tight_clean_);
        tr.registerDerivedVar("NBJets_pt45_tight_clean",   (bjets_pt45_tight_clean_==nullptr)?0:bjets_pt45_tight_clean_->size());
    }

public:
    BJet() 
        : bjets_(nullptr)
        , bjets_pt30_(nullptr)
        , bjets_pt45_(nullptr)
        , bjets_tight_(nullptr)
        , bjets_pt30_tight_(nullptr)
        , bjets_pt45_tight_(nullptr)
        , bjets_clean_(nullptr)
        , bjets_pt30_clean_(nullptr)
        , bjets_pt45_clean_(nullptr)
        , bjets_tight_clean_(nullptr)
        , bjets_pt30_tight_clean_(nullptr)
        , bjets_pt45_tight_clean_(nullptr)
    {}

    void operator()(NTupleReader& tr)
    {
        bjet(tr);
    }
};

#endif
