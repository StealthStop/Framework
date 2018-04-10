#ifndef JET_H
#define JET_H

class Jet
{
private:
    std::vector<TLorentzVector>* jets_pt30_;
    std::vector<TLorentzVector>* jets_pt45_;

    void jet(NTupleReader& tr)
    {
        const auto& Jets = tr.getVec<TLorentzVector>("Jets");
        const auto& etaCut = tr.getVar<double>("etaCut");

        jets_pt30_ = new std::vector<TLorentzVector>();
        jets_pt45_ = new std::vector<TLorentzVector>();
        for (unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
        {
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

        tr.registerDerivedVar("NJets",       Jets.size());
        tr.registerDerivedVec("Jets_pt30",   jets_pt30_);
        tr.registerDerivedVar("NJets_pt30", (jets_pt30_==nullptr)?0:jets_pt30_->size());
        tr.registerDerivedVec("Jets_pt45",   jets_pt45_);
        tr.registerDerivedVar("NJets_pt45", (jets_pt45_==nullptr)?0:jets_pt45_->size());
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
