#ifndef BJET_H
#define BJET_H

class BJet
{
private:
    std::vector<TLorentzVector>* bjets_;
    std::vector<TLorentzVector>* bjets_pt30_;
    std::vector<TLorentzVector>* bjets_pt45_;

    void bjet(NTupleReader& tr)
    {
        const auto& Jets = tr.getVec<TLorentzVector>("Jets");
        const auto& Jets_bDiscriminatorCSV = tr.getVec<double>("Jets_bDiscriminatorCSV");
        const auto& etaCut = tr.getVar<double>("etaCut");
        
        bjets_ = new std::vector<TLorentzVector>();
        bjets_pt30_ = new std::vector<TLorentzVector>();
        bjets_pt45_ = new std::vector<TLorentzVector>();
        for (unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
        {
            TLorentzVector lv(Jets.at(ijet));
            if( abs(lv.Eta()) < etaCut &&
                Jets_bDiscriminatorCSV.at(ijet) > 0.8484 
                )
            {
                bjets_->push_back(lv); 
                if(lv.Pt() > 30)
                    bjets_pt30_->push_back(lv);
                if(lv.Pt() > 45)
                    bjets_pt45_->push_back(lv);
            }
        }

        tr.registerDerivedVec("BJets",        bjets_);
        tr.registerDerivedVar("NBJets",      (bjets_==nullptr)?0:bjets_->size());
        tr.registerDerivedVec("BJets_pt30",   bjets_pt30_);
        tr.registerDerivedVar("NBJets_pt30", (bjets_pt30_==nullptr)?0:bjets_pt30_->size());
        tr.registerDerivedVec("BJets_pt45",   bjets_pt45_);
        tr.registerDerivedVar("NBJets_pt45", (bjets_pt45_==nullptr)?0:bjets_pt45_->size());
    }

public:
    BJet() 
        : bjets_(nullptr)
        , bjets_pt30_(nullptr)
        , bjets_pt45_(nullptr)
    {}

    void operator()(NTupleReader& tr)
    {
        bjet(tr);
    }
};

#endif
