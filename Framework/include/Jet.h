#ifndef JET_H
#define JET_H

class Jet
{
private:
    std::vector<TLorentzVector>* jets_pt30_;
    std::vector<TLorentzVector>* jets_pt45_;
    std::vector<TLorentzVector>* clean_jets_pt30_;
    std::vector<TLorentzVector>* clean_jets_pt45_;

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

        clean_jets_pt30_ = new std::vector<TLorentzVector>();
        clean_jets_pt45_ = new std::vector<TLorentzVector>();

        for (unsigned int icjet = 0; icjet < cleanJets.size(); ++icjet)
        {
            if( !cleanJetsID.at(icljet) ) continue;

            TLorentzVector lv(cleanJets.at(icjet));
            if( abs(lv.Eta()) < etaCut &&
                lv.Pt() > 30
                )
            {
                clean_jets_pt30_->push_back(lv);
                if( lv.Pt() > 45 )
                    clean_jets_pt45_->push_back(lv);
            }
        }

        tr.registerDerivedVar("NJets",       Jets.size());
        tr.registerDerivedVec("Jets_pt30",   jets_pt30_);
        tr.registerDerivedVar("NJets_pt30", (jets_pt30_==nullptr)?0:jets_pt30_->size());
        tr.registerDerivedVec("Jets_pt45",   jets_pt45_);
        tr.registerDerivedVar("NJets_pt45", (jets_pt45_==nullptr)?0:jets_pt45_->size());

        tr.registerDerivedVar("NCleanJets",         cleanJets.size());
        tr.reigsterDerivedVec("Clean_Jets_pt30",    clean_jets_pt30_);
        tr.registerDerivedVar("NClean_Jets_pt30",  (clean_jets_pt30_==nullptr)?0:clean_jets_pt30->size());
        tr.registerDerivedVec("Clean_Jets_pt45",    clean_jets_pt45_);
        tr.registerDerivedVar("NClean_Jets_pt45",  (clean_jets_pt45_==nullptr)?0:clean_jets_pt45->size());
    }

public:
    Jet() : 
        jets_pt30_(nullptr),
        jets_pt45_(nullptr)
    {
    }

    void operator()(NTupleReader& tr)
    {
        jet(tr);
    }
};

#endif
