#ifndef MakeStopHemispheres_h
#define MakeStopHemispheres_h

#include "Framework/Framework/include/Hemispheres.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include <vector>
#include <iostream>
#include <cmath>

class MakeStopHemispheres
{
private:
    std::string jetName_, jetMaskName_, nJetName_, myVarSuffix_;
    Hemisphere::SeedMethod seedMethod_;

    template<typename T> void orderVars(T& stop1, T& stop2, const T pseudo1, const T pseudo2, const bool pseudo1Tostop1) const
    {
        if(pseudo1Tostop1) 
        {            
            stop1 = pseudo1;
            stop2 = pseudo2;      
        } 
        else
        {
            stop1 = pseudo2;
            stop2 = pseudo1;      
        }
    }
    
    void getHemispheres(NTupleReader& tr) const
    {
        const auto& met                   = tr.getVar<double>("MET");
        const auto& metPhi                = tr.getVar<double>("METPhi");
        const auto& Jets                  = tr.getVec<TLorentzVector>(jetName_);
        const auto& GoodJets              = tr.getVec<bool>(jetMaskName_);
        const auto& NGoodJets             = tr.getVar<int>(nJetName_);
        const auto& GoodLeptons           = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        const auto& event_beta_z          = tr.getVar<double>("event_beta_z"); // for NN
        const auto& event_phi_rotate      = tr.getVar<double>("event_phi_rotate"); // for NN      

        static const int hemi_association = 3; // 3: 3th method, 'lund' used by MT2  
        TLorentzVector stop1_PtRank,       stop2_PtRank;
        TLorentzVector stop1_MassRank,     stop2_MassRank;
        TLorentzVector stop1_ScalarPtRank, stop2_ScalarPtRank;

        TLorentzVector stop1_PtRank_cm,       stop2_PtRank_cm;
        TLorentzVector stop1_MassRank_cm,     stop2_MassRank_cm;
        TLorentzVector stop1_ScalarPtRank_cm, stop2_ScalarPtRank_cm;

        double stop1ScalarPt_ScalarPtRank = -9999.9, stop2ScalarPt_ScalarPtRank = -9999.9;
        double MT2                        = 0.0;
        double dR_stop1stop2              = -1;
        double dPhi_stop1stop2            = -1;
        double difference_stopMasses      = -9999.9;
        double average_stopMasses         = -9999.9;
        double relativeDiff_stopMasses    = -9999.9;
 
        if(NGoodJets >= 2)
        {
            TLorentzVector MET;
            MET.SetPtEtaPhiM(met, 0.0, metPhi, 0.0);

            std::vector<float> px, py, pz, E;
            for(unsigned int i=0; i < Jets.size(); ++i)
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
            asymm_mt2_lester_bisect::disableCopyrightMessage();
            Hemisphere hemi(px, py, pz, E, seedMethod_, hemi_association); // to get MT2 hemisphere jets
            std::vector<int> grouping = hemi.getGrouping();
            TLorentzVector pseudojet1, pseudojet2;
            TLorentzVector pseudojet1_cm, pseudojet2_cm;
            double pseudojet1ScalarPt = 0.0, pseudojet2ScalarPt = 0.0;

            // Perform vector and scalar sum of mega jets
            for(unsigned int i=0; i < px.size(); ++i)
            {
                TLorentzVector obj;
                obj.SetPxPyPzE(px[i], py[i], pz[i], E[i]);

                if(grouping[i] == 1)
                {
                    //vector sum
                    pseudojet1 += obj;

                    //scalar sum
                    pseudojet1ScalarPt += obj.Pt();
                }
                else if(grouping[i] == 2)
                {
                    //vector sum
                    pseudojet2 += obj;

                    //scalar sum
                    pseudojet2ScalarPt += obj.Pt();
                }
            }

            // ---------------------------------------
            // -- Rank mega jets (Mass, Pt, Scalar Pt) 
            // ---------------------------------------
            orderVars(stop1_PtRank,               stop2_PtRank,               pseudojet1,         pseudojet2,         pseudojet1.Pt()    > pseudojet2.Pt()   );
            orderVars(stop1_MassRank,             stop2_MassRank,             pseudojet1,         pseudojet2,         pseudojet1.M()     > pseudojet2.M()    );
            orderVars(stop1_ScalarPtRank,         stop2_ScalarPtRank,         pseudojet1,         pseudojet2,         pseudojet1ScalarPt > pseudojet2ScalarPt);
            orderVars(stop1ScalarPt_ScalarPtRank, stop2ScalarPt_ScalarPtRank, pseudojet1ScalarPt, pseudojet2ScalarPt, pseudojet1ScalarPt > pseudojet2ScalarPt);                

            MT2                     = ttUtility::coreMT2calc(pseudojet1, pseudojet2, MET);
            dR_stop1stop2           = pseudojet1.DeltaR(pseudojet2);    
            dPhi_stop1stop2         = pseudojet1.DeltaPhi(pseudojet2);
            difference_stopMasses   = abs( pseudojet1.M() - pseudojet2.M() );
            average_stopMasses      = ( 0.5 * ( pseudojet1.M() + pseudojet2.M() ) );
            relativeDiff_stopMasses = difference_stopMasses / average_stopMasses;

            // ----------------------------------------------------
            // -- Boost to CM frame & Rotate the hemispheres for NN 
            // ----------------------------------------------------
            pseudojet1_cm = pseudojet1;
            pseudojet2_cm = pseudojet2;
        
            TVector3 boostVec(0.0, 0.0, -event_beta_z);  
    
            pseudojet1_cm.Boost(boostVec);
            pseudojet1_cm.RotateZ(-event_phi_rotate);

            pseudojet2_cm.Boost(boostVec);
            pseudojet2_cm.RotateZ(-event_phi_rotate); 

            // Rank the hemispheres (Mass, Pt, Scalar Pt) 
            orderVars(stop1_PtRank_cm,       stop2_PtRank_cm,       pseudojet1_cm, pseudojet2_cm, pseudojet1_cm.Pt() > pseudojet2_cm.Pt());
            orderVars(stop1_MassRank_cm,     stop2_MassRank_cm,     pseudojet1_cm, pseudojet2_cm, pseudojet1_cm.M()  > pseudojet2_cm.M() );
            orderVars(stop1_ScalarPtRank_cm, stop2_ScalarPtRank_cm, pseudojet1_cm, pseudojet2_cm, pseudojet1ScalarPt > pseudojet2ScalarPt);

        }
        tr.registerDerivedVar("stop1_PtRank"+myVarSuffix_,stop1_PtRank);
        tr.registerDerivedVar("stop2_PtRank"+myVarSuffix_,stop2_PtRank);
        tr.registerDerivedVar("stop1_PtRank"+myVarSuffix_+"_mass",stop1_PtRank.M()); // for NN
        tr.registerDerivedVar("stop2_PtRank"+myVarSuffix_+"_mass",stop2_PtRank.M());
        tr.registerDerivedVar("stop1_PtRank"+myVarSuffix_+"_pt",stop1_PtRank.Pt());
        tr.registerDerivedVar("stop2_PtRank"+myVarSuffix_+"_pt",stop2_PtRank.Pt());
        tr.registerDerivedVar("stop1_PtRank"+myVarSuffix_+"_phi",stop1_PtRank.Phi());
        tr.registerDerivedVar("stop2_PtRank"+myVarSuffix_+"_phi",stop2_PtRank.Phi());
        tr.registerDerivedVar("stop1_PtRank"+myVarSuffix_+"_eta",stop1_PtRank.Eta());
        tr.registerDerivedVar("stop2_PtRank"+myVarSuffix_+"_eta",stop2_PtRank.Eta());
        tr.registerDerivedVar("stop1_PtRank_cm"+myVarSuffix_+"_mass",stop1_PtRank_cm.M()); // for NN
        tr.registerDerivedVar("stop2_PtRank_cm"+myVarSuffix_+"_mass",stop2_PtRank_cm.M());
        tr.registerDerivedVar("stop1_PtRank_cm"+myVarSuffix_+"_pt",stop1_PtRank_cm.Pt());
        tr.registerDerivedVar("stop2_PtRank_cm"+myVarSuffix_+"_pt",stop2_PtRank_cm.Pt());
        tr.registerDerivedVar("stop1_PtRank_cm"+myVarSuffix_+"_phi",stop1_PtRank_cm.Phi());
        tr.registerDerivedVar("stop2_PtRank_cm"+myVarSuffix_+"_phi",stop2_PtRank_cm.Phi());
        tr.registerDerivedVar("stop1_PtRank_cm"+myVarSuffix_+"_eta",stop1_PtRank_cm.Eta());
        tr.registerDerivedVar("stop2_PtRank_cm"+myVarSuffix_+"_eta",stop2_PtRank_cm.Eta());
        tr.registerDerivedVar("stop1_MassRank"+myVarSuffix_,stop1_MassRank);
        tr.registerDerivedVar("stop2_MassRank"+myVarSuffix_,stop2_MassRank);
        tr.registerDerivedVar("stop1_MassRank"+myVarSuffix_+"_mass",stop1_MassRank.M()); // for NN
        tr.registerDerivedVar("stop2_MassRank"+myVarSuffix_+"_mass",stop2_MassRank.M());
        tr.registerDerivedVar("stop1_MassRank"+myVarSuffix_+"_pt",stop1_MassRank.Pt());
        tr.registerDerivedVar("stop2_MassRank"+myVarSuffix_+"_pt",stop2_MassRank.Pt());
        tr.registerDerivedVar("stop1_MassRank"+myVarSuffix_+"_phi",stop1_MassRank.Phi());
        tr.registerDerivedVar("stop2_MassRank"+myVarSuffix_+"_phi",stop2_MassRank.Phi());
        tr.registerDerivedVar("stop1_MassRank"+myVarSuffix_+"_eta",stop1_MassRank.Eta());
        tr.registerDerivedVar("stop2_MassRank"+myVarSuffix_+"_eta",stop2_MassRank.Eta());
        tr.registerDerivedVar("stop1_MassRank_cm"+myVarSuffix_+"_mass",stop1_MassRank_cm.M()); // for NN
        tr.registerDerivedVar("stop2_MassRank_cm"+myVarSuffix_+"_mass",stop2_MassRank_cm.M());
        tr.registerDerivedVar("stop1_MassRank_cm"+myVarSuffix_+"_pt",stop1_MassRank_cm.Pt());
        tr.registerDerivedVar("stop2_MassRank_cm"+myVarSuffix_+"_pt",stop2_MassRank_cm.Pt());
        tr.registerDerivedVar("stop1_MassRank_cm"+myVarSuffix_+"_phi",stop1_MassRank_cm.Phi());
        tr.registerDerivedVar("stop2_MassRank_cm"+myVarSuffix_+"_phi",stop2_MassRank_cm.Phi());
        tr.registerDerivedVar("stop1_MassRank_cm"+myVarSuffix_+"_eta",stop1_MassRank_cm.Eta());
        tr.registerDerivedVar("stop2_MassRank_cm"+myVarSuffix_+"_eta",stop2_MassRank_cm.Eta());
        tr.registerDerivedVar("stop1_ScalarPtRank"+myVarSuffix_,stop1_ScalarPtRank);
        tr.registerDerivedVar("stop2_ScalarPtRank"+myVarSuffix_,stop2_ScalarPtRank);
        tr.registerDerivedVar("stop1_ScalarPtRank"+myVarSuffix_+"_mass",stop1_ScalarPtRank.M()); // for NN
        tr.registerDerivedVar("stop2_ScalarPtRank"+myVarSuffix_+"_mass",stop2_ScalarPtRank.M());
        tr.registerDerivedVar("stop1_ScalarPtRank"+myVarSuffix_+"_pt",stop1_ScalarPtRank.Pt());
        tr.registerDerivedVar("stop2_ScalarPtRank"+myVarSuffix_+"_pt",stop2_ScalarPtRank.Pt());
        tr.registerDerivedVar("stop1_ScalarPtRank"+myVarSuffix_+"_phi",stop1_ScalarPtRank.Phi());
        tr.registerDerivedVar("stop2_ScalarPtRank"+myVarSuffix_+"_phi",stop2_ScalarPtRank.Phi());
        tr.registerDerivedVar("stop1_ScalarPtRank"+myVarSuffix_+"_eta",stop1_ScalarPtRank.Eta());
        tr.registerDerivedVar("stop2_ScalarPtRank"+myVarSuffix_+"_eta",stop2_ScalarPtRank.Eta());
        tr.registerDerivedVar("stop1_ScalarPtRank_cm"+myVarSuffix_+"_mass",stop1_ScalarPtRank_cm.M()); // for NN
        tr.registerDerivedVar("stop2_ScalarPtRank_cm"+myVarSuffix_+"_mass",stop2_ScalarPtRank_cm.M());
        tr.registerDerivedVar("stop1_ScalarPtRank_cm"+myVarSuffix_+"_pt",stop1_ScalarPtRank_cm.Pt());
        tr.registerDerivedVar("stop2_ScalarPtRank_cm"+myVarSuffix_+"_pt",stop2_ScalarPtRank_cm.Pt());
        tr.registerDerivedVar("stop1_ScalarPtRank_cm"+myVarSuffix_+"_phi",stop1_ScalarPtRank_cm.Phi());
        tr.registerDerivedVar("stop2_ScalarPtRank_cm"+myVarSuffix_+"_phi",stop2_ScalarPtRank_cm.Phi());
        tr.registerDerivedVar("stop1_ScalarPtRank_cm"+myVarSuffix_+"_eta",stop1_ScalarPtRank_cm.Eta());
        tr.registerDerivedVar("stop2_ScalarPtRank_cm"+myVarSuffix_+"_eta",stop2_ScalarPtRank_cm.Eta());
        tr.registerDerivedVar("stop1ScalarPt_ScalarPtRank"+myVarSuffix_,stop1ScalarPt_ScalarPtRank);
        tr.registerDerivedVar("stop2ScalarPt_ScalarPtRank"+myVarSuffix_,stop1ScalarPt_ScalarPtRank); 
        tr.registerDerivedVar("MT2"+myVarSuffix_,MT2);
        tr.registerDerivedVar("dR_stop1stop2"+myVarSuffix_,dR_stop1stop2);
        tr.registerDerivedVar("dPhi_stop1stop2"+myVarSuffix_,dPhi_stop1stop2);
        tr.registerDerivedVar("difference_stopMasses"+myVarSuffix_,difference_stopMasses);
        tr.registerDerivedVar("average_stopMasses"+myVarSuffix_,average_stopMasses);
        tr.registerDerivedVar("relativeDiff_stopMasses"+myVarSuffix_,relativeDiff_stopMasses);
    }

public:    
    MakeStopHemispheres(const std::string& jetName = "Jets", const std::string& jetMaskName = "GoodJets_pt45", const std::string& nJetName = "NGoodJets_pt45", const std::string& myVarSuffix = "", const Hemisphere::SeedMethod& seedMethod = Hemisphere::InvMassSeed)
        : jetName_(jetName) 
        , jetMaskName_(jetMaskName)
        , nJetName_(nJetName)
        , myVarSuffix_(myVarSuffix)
        , seedMethod_(seedMethod)
    {
        std::cout<<"Setting up StopHemispheres with jet collection: \""<<jetMaskName<<"\""<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        getHemispheres(tr);
    }
};

#endif
