#ifndef PHOTON_H
#define PHOTON_H

class Photon
{
private:
    void photon(NTupleReader& tr)
    {
        const auto& allPhotons = tr.getVec<TLorentzVector>("Photons");
        const auto& allPhotons_fullID = tr.getVec<bool>("Photons_fullID");
        const auto& etaCut = 1.4;//tr.getVar<double>("etaCut");

        auto* good_electrons_ = new std::vector<bool>();
        int NGoodPhotons = 0;
        for(unsigned int iel = 0; iel < allPhotons.size(); ++iel)
        {
            TLorentzVector lvel = allPhotons.at(iel);
            if( abs(lvel.Eta()) < etaCut && 
                lvel.Pt() > 200 && 
                allPhotons_fullID.at(iel) 
                )
            {
                good_electrons_->push_back(true);
                NGoodPhotons++;
            }
            else
            {
                good_electrons_->push_back(false);
            }
        }

        tr.registerDerivedVec("GoodPhotons", good_electrons_);
        tr.registerDerivedVar("NGoodPhotons", NGoodPhotons);
    }

public:
    Photon() 
    {}

    void operator()(NTupleReader& tr)
    {
        photon(tr);
    }
};

#endif
