#ifndef MakeMT2Hemispheres_h
#define MakeMT2Hemispheres_h

#include "Framework/Framework/include/MT2Hemispheres.h"
#include "Framework/Framework/include/Davismt2.h"

#include <vector>
#include <iostream>
#include <cmath>

class MakeMT2Hemispheres
{
private:
    std::string jetMaskName_, nJetName_, myVarSuffix_;

    const double CalcMT2(const double testmass, const bool massive, const TLorentzVector& visible1, const TLorentzVector& visible2, const TLorentzVector& MET) const
    {
        double pa[3];
        double pb[3];
        double pmiss[3];

        pmiss[0] = 0;
        pmiss[1] = MET.Px();
        pmiss[2] = MET.Py();

        pa[0] = massive ? visible1.M() : 0;
        pa[1] = visible1.Px();
        pa[2] = visible1.Py();

        pb[0] = massive ? visible2.M() : 0;
        pb[1] = visible2.Px();
        pb[2] = visible2.Py();

        Davismt2 mt2;
        mt2.set_momenta(pa, pb, pmiss);
        mt2.set_mn(testmass);
        return mt2.get_mt2();
    }

    void getHemispheres(NTupleReader& tr) const
    {
        // these variables values the same as  MT2
        const double testmass      = 0.0;
        const bool massive         = false; 
        const int hemi_association = 3; // 3: 3th method, 'lund' used by MT2  

        const auto& met            = tr.getVar<double>("MET");
        const auto& metPhi         = tr.getVar<double>("METPhi");
        const auto& Jets           = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets       = tr.getVec<bool>(jetMaskName_);
        const auto& NGoodJets      = tr.getVar<int>(nJetName_);
        const auto& GoodLeptons    = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");

        double stopMass            = 0.0;
        double hemi1Mass           = -9999.9, hemi1Eta = -9999.9, hemi1Phi = -9999.9, hemi1Pt = -9999.9;
        double hemi2Mass           = -9999.9, hemi2Eta = -9999.9, hemi2Phi = -9999.9, hemi2Pt = -9999.9;
        double dR_hemi1hemi2       = -1;
        double dPhi_hemi1hemi2     = -1;

        if(NGoodJets >= 2)
        {
            TLorentzVector MET;
            MET.SetPtEtaPhiM(met, 0.0, metPhi, 0.0);

            vector<float> px, py, pz, E;
            for(int i=0; i < Jets.size(); ++i)
            {
                if(!GoodJets[i]) continue;
                px.push_back(Jets[i].Px());
                py.push_back(Jets[i].Py());
                pz.push_back(Jets[i].Pz());
                E .push_back(Jets[i].E ());
            }
            for(const auto& pair : GoodLeptons)
            {
                px.push_back(pair.second.Px());
                py.push_back(pair.second.Py());
                pz.push_back(pair.second.Pz());
                E .push_back(pair.second.E ());                
            }

            // Get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)  
            Hemisphere hemi(px, py, pz, E, 2, hemi_association);
            vector<int> grouping = hemi.getGrouping();
            TLorentzVector pseudojet1(0.0, 0.0, 0.0, 0.0);
            TLorentzVector pseudojet2(0.0, 0.0, 0.0, 0.0);

            for(int i=0; i < px.size(); ++i)
            {
                if(grouping[i] == 1)
                {
                    pseudojet1.SetPx(pseudojet1.Px() + px[i]);
                    pseudojet1.SetPy(pseudojet1.Py() + py[i]);
                    pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
                    pseudojet1.SetE( pseudojet1.E()  + E[i]);
                }
            
                else if(grouping[i] == 2)
                {
                    pseudojet2.SetPx(pseudojet2.Px() + px[i]);
                    pseudojet2.SetPy(pseudojet2.Py() + py[i]);
                    pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
                    pseudojet2.SetE( pseudojet2.E()  + E[i]);
                }
            }
        
            stopMass = CalcMT2(testmass, massive, pseudojet1, pseudojet2, MET);
        
            // --------------------------
            // -- stop MT2 hemispheres 
            // --------------------------
            hemi1Mass       = pseudojet1.M();
            hemi1Eta        = pseudojet1.Eta();
            hemi1Phi        = pseudojet1.Phi();
            hemi1Pt         = pseudojet1.Pt();
            hemi2Mass       = pseudojet2.M();
            hemi2Eta        = pseudojet2.Eta();
            hemi2Phi        = pseudojet2.Phi();
            hemi2Pt         = pseudojet2.Pt(); 
            dR_hemi1hemi2   = pseudojet1.DeltaR(pseudojet2);    
            dPhi_hemi1hemi2 = pseudojet1.DeltaPhi(pseudojet2);

        }
        tr.registerDerivedVar("stopMass"+myVarSuffix_,stopMass);
        tr.registerDerivedVar("hemi1Mass"+myVarSuffix_,hemi1Mass);
        tr.registerDerivedVar("hemi1Eta"+myVarSuffix_,hemi1Eta);
        tr.registerDerivedVar("hemi1Phi"+myVarSuffix_,hemi1Phi);
        tr.registerDerivedVar("hemi1Pt"+myVarSuffix_,hemi1Pt);
        tr.registerDerivedVar("hemi2Mass"+myVarSuffix_,hemi2Mass);
        tr.registerDerivedVar("hemi2Eta"+myVarSuffix_,hemi2Eta);
        tr.registerDerivedVar("hemi2Phi"+myVarSuffix_,hemi2Phi);
        tr.registerDerivedVar("hemi2Pt"+myVarSuffix_,hemi2Pt);
        tr.registerDerivedVar("dR_hemi1hemi2"+myVarSuffix_,dR_hemi1hemi2);
        tr.registerDerivedVar("dPhi_hemi1hemi2"+myVarSuffix_,dPhi_hemi1hemi2);
    }

public:    
    MakeMT2Hemispheres(const std::string& jetMaskName = "GoodJets_pt45", const std::string& nJetName = "NGoodJets_pt45", const std::string& myVarSuffix = "")
        : jetMaskName_(jetMaskName)
        , nJetName_(nJetName)
        , myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up MT2Hemispheres with jet collection: \""<<jetMaskName<<"\""<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        getHemispheres(tr);
    }
};

#endif
