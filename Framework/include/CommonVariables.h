#ifndef COMMONVARIABLES_H
#define COMMONVARIABLES_H

class CommonVariables
{
private:
    void commonVariables(NTupleReader& tr)
    {
        // Get needed branches
        const auto& Jets = tr.getVec<TLorentzVector>("Jets");
        const auto& BJets_pt30 = tr.getVec<TLorentzVector>("BJets_pt30");
        const auto& GoodMuons = tr.getVec<TLorentzVector>("GoodMuons");
        const auto& GoodMuonsMTW = tr.getVec<double>("GoodMuonsMTW");
        const auto& GoodElectrons = tr.getVec<TLorentzVector>("GoodElectrons");
        const auto& GoodElectronsMTW = tr.getVec<double>("GoodElectronsMTW");
        const auto& etaCut = tr.getVar<double>("etaCut");

        // HT with jets with pT>40
        double ht = 0;
        for (TLorentzVector jet : Jets)
        {
            if(jet.Pt() > 40 && abs(jet.Eta()) < etaCut)
                ht += jet.Pt();
        }
        tr.registerDerivedVar("HT_trigger", ht);

        // Put leptons together
        std::vector<TLorentzVector>* GoodLeptons = new std::vector<TLorentzVector>();
        for (TLorentzVector muon : GoodMuons)
            GoodLeptons->push_back(muon);
        for (TLorentzVector electron : GoodElectrons)
            GoodLeptons->push_back(electron);
        tr.registerDerivedVec("GoodLeptons", GoodLeptons);
        tr.registerDerivedVar("NGoodLeptons", GoodLeptons->size());
        
        // M(l,b); closest to 105 GeV if multiple combinations (halfway between 30 and 180 GeV)
        double Mbl = 0;
        double Mbldiff = 999.;
        for(TLorentzVector lepton : *GoodLeptons)
        {
            for(TLorentzVector bjet : BJets_pt30)
            {
                double mbl = (lepton+bjet).M();
                if( abs(mbl - 105) < Mbldiff)
                {
                    Mbl = mbl;
                    Mbldiff = abs(mbl - 105);
                }
            }
        }
        tr.registerDerivedVar("Mbl",Mbl);
        
        //Find single lepton for HistoContainer
        TLorentzVector singleLepton;
        double mTLep = 999.9;
        for(unsigned int i = 0; i < GoodMuons.size(); ++i)
        {
            if(GoodMuons[i].Pt() > 20)
            {
                singleLepton = GoodMuons[i];
                mTLep = GoodMuonsMTW[i];
                break;
            }
        }
        for(unsigned int i = 0; i < GoodElectrons.size(); ++i)
        {
            if(GoodElectrons[i].Pt() > 20)
            {
                singleLepton = GoodElectrons[i];
                mTLep = GoodElectronsMTW[i];
                break;
            }
        }
        tr.registerDerivedVar("singleLepton", singleLepton);            
    }

public:
    CommonVariables()
    {}

    void operator()(NTupleReader& tr)
    {
        commonVariables(tr);
    }
};

#endif