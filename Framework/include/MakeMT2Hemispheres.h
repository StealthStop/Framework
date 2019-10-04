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
        //-----------------------------------
        //Calculate/find the folowing variables
        //-----------------------------------
        const double testmass = 0.0;
        const bool massive = true; 
        const int hemi_association = 1;

        const auto& met = tr.getVar<double>("MET");
        const auto& metPhi = tr.getVar<double>("METPhi");
        const auto& Jets = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets = tr.getVec<bool>(jetMaskName_);
        const auto& NGoodJets = tr.getVar<int>(nJetName_);
        const auto& GoodLeptons = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");

        double stopMass = 0.0;
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
            for(int i=0; i<px.size(); ++i)
            {
                if(grouping[i]==1)
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
        }
        tr.registerDerivedVar("stopMass"+myVarSuffix_,stopMass);
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
