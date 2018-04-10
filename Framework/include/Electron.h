#ifndef ELECTRON_H
#define ELECTRON_H

class Electron
{
private:
    std::vector<TLorentzVector>* good_electrons_;
    std::vector<int>* good_electrons_charge_;
    void electron(NTupleReader& tr)
    {
        const auto& allElectrons = tr.getVec<TLorentzVector>("Electrons");
        const auto& allElectrons_passIso = tr.getVec<bool>("Electrons_passIso");
        const auto& allElectrons_charge  = tr.getVec<int>("Electrons_charge");
        const auto& allElectrons_tightID = tr.getVec<bool>("Electrons_tightID");
        const auto& etaCut = tr.getVar<double>("etaCut");
            
        good_electrons_ = new std::vector<TLorentzVector>();
        good_electrons_charge_ = new std::vector<int>();
        for (unsigned int iel = 0; iel < allElectrons.size(); ++iel)
        {
            TLorentzVector lvel(allElectrons.at(iel));
            if( abs(lvel.Eta()) < etaCut && 
                lvel.Pt() > 30 && 
                allElectrons_passIso.at(iel) &&
                allElectrons_tightID.at(iel) 
                )
            {
                good_electrons_->push_back(lvel);
                good_electrons_charge_->push_back(allElectrons_charge.at(iel));
            }
        }


        tr.registerDerivedVec("GoodElectrons", good_electrons_);
        tr.registerDerivedVar("NGoodElectrons", (good_electrons_==nullptr)?0:good_electrons_->size());
        tr.registerDerivedVec("GoodElectronsCharge", good_electrons_charge_);
    }

public:
    Electron() 
        : good_electrons_(nullptr)
        , good_electrons_charge_(nullptr) 
    {}

    void operator()(NTupleReader& tr)
    {
        electron(tr);
    }
};

#endif
