#ifndef JET_H
#define JET_H

class Jet
{
private:
    std::vector<TLorentzVector>* jets_pt30_;
    std::vector<TLorentzVector>* jets_pt45_;
    std::vector<TLorentzVector>* jets_pt30_clean_;
    std::vector<TLorentzVector>* jets_pt45_clean_;

    void jet(NTupleReader& tr)
    {
        const auto& Jets = tr.getVec<TLorentzVector>("Jets");
        const auto& etaCut = tr.getVar<double>("etaCut");
        const auto& jetsID = tr.getVec<bool>("Jets_ID");
        const auto& cleanJets = tr.getVec<TLorentzVector>("Jetsclean");
        const auto& cleanJetsID = tr.getVec<bool>("Jetsclean_ID");

        jets_pt30_ = new std::vector<TLorentzVector>();
        jets_pt45_ = new std::vector<TLorentzVector>();
        for (unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
        {
            if( !jetsID.at(ijet) ) continue;

            TLorentzVector lv(Jets.at(ijet));
            if( abs(lv.Eta()) < etaCut && 
                lv.Pt() > 30 
                )
            {
                jets_pt30_->push_back(lv);
                if(lv.Pt() > 45)
                    jets_pt45_->push_back(lv);
            }
        }

        jets_pt30_clean_ = new std::vector<TLorentzVector>();
        jets_pt45_clean_ = new std::vector<TLorentzVector>();

        if( cleanJets.size() != 0 ) {
            for (unsigned int icjet = 0; icjet < cleanJets.size(); ++icjet)
            {
                if( !cleanJetsID.at(icjet) ) continue;

                TLorentzVector lv(cleanJets.at(icjet));
                if( abs(lv.Eta()) < etaCut &&
                    lv.Pt() > 30
                    )
                {
                    jets_pt30_clean_->push_back(lv);
                    if( lv.Pt() > 45 )
                        jets_pt45_clean_->push_back(lv);
                }
            }
        }
        else if( cleanJets.size() == 0 ) {
            for (unsigned int icjet = 0; icjet < Jets.size(); ++icjet)
            {
                if( !jetsID.at(icjet) ) continue;

                TLorentzVector lv(Jets.at(icjet));
                if( abs(lv.Eta()) < etaCut &&
                    lv.Pt() > 30
                    )
                {
                    jets_pt30_clean_->push_back(lv);
                    if( lv.Pt() > 45 )
                        jets_pt45_clean_->push_back(lv);
                }
            }
        }

        tr.registerDerivedVar("NJets",       Jets.size());
        tr.registerDerivedVec("Jets_pt30",   jets_pt30_);
        tr.registerDerivedVar("NJets_pt30", (jets_pt30_==nullptr)?0:jets_pt30_->size());
        tr.registerDerivedVec("Jets_pt45",   jets_pt45_);
        tr.registerDerivedVar("NJets_pt45", (jets_pt45_==nullptr)?0:jets_pt45_->size());

        tr.registerDerivedVar("NCleanJets",         cleanJets.size());
        tr.registerDerivedVec("Jets_pt30_clean",    jets_pt30_clean_);
        tr.registerDerivedVar("NJets_pt30_clean",  (jets_pt30_clean_==nullptr)?0:jets_pt30_clean_->size());
        tr.registerDerivedVec("Jets_pt45_clean",    jets_pt45_clean_);
        tr.registerDerivedVar("NJets_pt45_clean",  (jets_pt45_clean_==nullptr)?0:jets_pt45_clean_->size());
    }

public:
    Jet() : 
        jets_pt30_(nullptr)
      , jets_pt45_(nullptr)
      , jets_pt30_clean_(nullptr)
      , jets_pt45_clean_(nullptr)
    {
    }

    void operator()(NTupleReader& tr)
    {
        jet(tr);
    }
};

#endif
