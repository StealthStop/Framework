#ifndef BJET_H
#define BJET_H

class BJet
{
private:
    std::vector<TLorentzVector>* bjets_;
    std::vector<TLorentzVector>* bjets_pt45_;
    void bjet(NTupleReader& tr)
    {
        const std::vector<TLorentzVector>& Jets = tr.getVec<TLorentzVector>("Jets");
        const std::vector<double>& Jets_bDiscriminatorCSV = tr.getVec<double>("Jets_bDiscriminatorCSV");
            
        for (unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
        {
            TLorentzVector lv(Jets.at(ijet));
            if( abs(lv.Eta()) < 2.4 && 
                lv.Pt() > 30 && 
                Jets_bDiscriminatorCSV.at(ijet) > 0.8484
                )
            {
                bjets_->push_back(lv);
                if(lv.Pt() > 45)
                    bjets_pt45_->push_back(lv);
            }
        }

        tr.registerDerivedVec("BJets", bjets_);
        tr.registerDerivedVar("NBJets", bjets_->size());
        tr.registerDerivedVec("BJets_pt45", bjets_pt45_);
        tr.registerDerivedVar("NBJets_pt45", bjets_pt45_->size());
    }

public:
    BJet() 
        : bjets_(nullptr)
    {}

    void operator()(NTupleReader& tr)
    {
        bjet(tr);
    }
};

#endif
