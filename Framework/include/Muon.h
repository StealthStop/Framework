#ifndef MUON_H
#define MUON_H

class Muon
{
private:
    std::string myVarSuffix_;

    void muon(NTupleReader& tr)
    {
        const auto& allMuons = tr.getVec<TLorentzVector>("Muons");
        const auto& allMuons_passIso = tr.getVec<bool>("Muons_passIso");
        const auto& allMuons_charge = tr.getVec<int>("Muons_charge");
        const auto& allMuons_medID = tr.getVec<bool>("Muons_mediumID");
        const auto& etaCut = tr.getVar<double>("etaCut");
        const auto& MET = tr.getVar<double>("MET"); 
        const auto& METPhi = tr.getVar<double>("METPhi");

        TLorentzVector lvMET;
        lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);

        auto* good_muons_ = new std::vector<bool>();
        auto* muons_mtw_ = new std::vector<double>;
        int NGoodMuons = 0;
        int NGoodPlusMuons = 0;
        int NGoodMinusMuons = 0;
        for(unsigned int imu = 0; imu < allMuons.size(); ++imu)
        {            
            TLorentzVector lvmu(allMuons.at(imu));
            double mtw = sqrt( 2*( lvMET.Pt()*lvmu.Pt() - (lvMET.Px()*lvmu.Px() + lvMET.Py()*lvmu.Py()) ) );
            muons_mtw_->push_back(mtw);
            if( abs(lvmu.Eta()) < etaCut && 
                lvmu.Pt() > 30 && 
                allMuons_passIso.at(imu) &&
                allMuons_medID.at(imu)
                )
            {
                good_muons_->push_back(true);
                NGoodMuons++;
                if( allMuons_charge.at(imu) ==  1 ) NGoodPlusMuons++;
                else if( allMuons_charge.at(imu) == -1 ) NGoodMinusMuons++;
                else std::cout<<"Charge values in nTuples are different"<<std::endl;
            }
            else
            {
                good_muons_->push_back(false);
            }
        }
        
        tr.registerDerivedVec("GoodMuons"+myVarSuffix_,  good_muons_);
        tr.registerDerivedVar("NGoodMuons"+myVarSuffix_, NGoodMuons );
        tr.registerDerivedVar("NGoodPlusMuons"+myVarSuffix_, NGoodPlusMuons );
        tr.registerDerivedVar("NGoodMinusMuons"+myVarSuffix_, NGoodMinusMuons );
        tr.registerDerivedVec("MuonsMTW"+myVarSuffix_,   muons_mtw_ );
    }

public:
    Muon(std::string myVarSuffix = "") 
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up Muon"<<std::endl;   
    }
    
    void operator()(NTupleReader& tr)
    {
        muon(tr);
    }
};

#endif
