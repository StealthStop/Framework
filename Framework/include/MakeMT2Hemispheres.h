#ifndef MakeMT2Hemispheres_h
#define MakeMT2Hemispheres_h

#include "Framework/Framework/include/MT2Hemispheres.h"
#include "Framework/Framework/include/MT2HemispheresUtilities.h"
#include "Framework/Framework/include/Davismt2.h"
#include "Framework/Framework/include/TMctLib.h"
#include "Framework/Framework/include/mctlib.h"

#include <vector>
#include <iostream>
#include <cmath>


class MakeMT2Hemispheres
{
private:
    std::string myVarSuffix_;

    double CalcMT2(double testmass, bool massive, TLorentzVector visible1, TLorentzVector visible2, TLorentzVector MET)
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

        Davismt2 *mt2 = new Davismt2();
        mt2->set_momenta(pa, pb, pmiss);
        mt2->set_mn(testmass);
        double MT2 = mt2->get_mt2();
        delete mt2;
        return MT2;
    }

    void getHemispheres(NTupleReader& tr) 
    {
        //-----------------------------------
        //Calculate/find the folowing variables
        //-----------------------------------
        double testmass = 0.0;
        bool massive = true; 
        int hemi_association = 1;

        const auto& met = tr.getVar<double>("MET");
        const auto& metPhi = tr.getVar<double>("METPhi");
        const auto& jet = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt45 = tr.getVec<bool>("GoodJets_pt45");
        const auto& NGoodJets_pt45 = tr.getVar<int>("NGoodJets_pt45");

        double stopMass = 0.0;
        if(NGoodJets_pt45 >= 2)
        {
            TLorentzVector MET;
            MET.SetPtEtaPhiM(met, 0.0, metPhi, 0.0);

            vector<float> px, py, pz, E;
            for(int i=0; i<jet.size(); ++i)
            {
                if(!GoodJets_pt45[i]) continue;
                px.push_back(jet[i].Px());
                py.push_back(jet[i].Py());
                pz.push_back(jet[i].Pz());
                E .push_back(jet[i].E ());
            }

            // Get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)  
            Hemisphere hemi(px, py, pz, E, 2, hemi_association);
            vector<int> grouping = hemi.getGrouping();

            TLorentzVector pseudojet1(0.,0.,0.,0.);
            TLorentzVector pseudojet2(0.,0.,0.,0.);
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
    MakeMT2Hemispheres(std::string myVarSuffix = "")
        : myVarSuffix_       (myVarSuffix)
    {
        std::cout<<"Setting up MT2Hemispheres"<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        getHemispheres(tr);
    }
};

#endif
