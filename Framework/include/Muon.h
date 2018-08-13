#ifndef MUON_H
#define MUON_H

class Muon
{
private:
    void muon(NTupleReader& tr)
    {
        const auto& allMuons = tr.getVec<TLorentzVector>("Muons");
        const auto& allMuons_passIso = tr.getVec<bool>("Muons_passIso");
        const auto& allMuons_charge = tr.getVec<int>("Muons_charge");
        const auto& etaCut = tr.getVar<double>("etaCut");
        const auto& MET = tr.getVar<double>("MET"); 
        const auto& METPhi = tr.getVar<double>("METPhi");

        TLorentzVector lvMET;
        lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);

        auto* good_muons_ = new std::vector<bool>();
        auto* muons_mtw_ = new std::vector<double>;
        int NGoodMuons = 0;
        for(unsigned int imu = 0; imu < allMuons.size(); ++imu)
        {            
            TLorentzVector lvmu(allMuons.at(imu));
            double mtw = sqrt( 2*( lvMET.Pt()*lvmu.Pt() - (lvMET.Px()*lvmu.Px() + lvMET.Py()*lvmu.Py()) ) );
            muons_mtw_->push_back(mtw);
            if( abs(lvmu.Eta()) < etaCut && 
                lvmu.Pt() > 30 && 
                allMuons_passIso.at(imu)
                )
            {
                good_muons_->push_back(true);
                NGoodMuons++;
            }
            else
            {
                good_muons_->push_back(false);
            }
        }
        
        tr.registerDerivedVec("GoodMuons",  good_muons_);
        tr.registerDerivedVar("NGoodMuons", NGoodMuons );
        tr.registerDerivedVec("MuonsMTW",   muons_mtw_ );
    }

public:
    Muon() 
    {}

    void operator()(NTupleReader& tr)
    {
        muon(tr);
    }
};

#endif
