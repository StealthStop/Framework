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

    template<typename T> void orderVars(T& stop1, T& stop2, const T hemi1, const T hemi2, const bool hemi1ToStop1) const
    {
        if(hemi1ToStop1) 
        {            
            stop1 = hemi1;
            stop2 = hemi2;      
        } 
        else
        {
            stop1 = hemi2;
            stop2 = hemi1;      
        }
    }
    
    void getHemispheres(NTupleReader& tr) const
    {
        const auto& met               = tr.getVar<double>("MET");
        const auto& metPhi            = tr.getVar<double>("METPhi");
        const auto& Jets              = tr.getVec<TLorentzVector>(jetName_);
        const auto& GoodJets          = tr.getVec<bool>(jetMaskName_);
        const auto& NGoodJets         = tr.getVar<int>(nJetName_);
        const auto& GoodLeptons       = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        const auto& event_beta_z      = tr.getVar<double>("event_beta_z"); // for NN
        const auto& event_phi_rotate  = tr.getVar<double>("event_phi_rotate"); // for NN      

        static const int hemi_association = 3; // 3: 3th method, 'lund' used by MT2  

        TLorentzVector Stop1,                 Stop2;
        TLorentzVector Stop1_PtRank,          Stop2_PtRank;
        TLorentzVector Stop1_MassRank,        Stop2_MassRank;
        TLorentzVector Stop1_ScalarPtRank,    Stop2_ScalarPtRank;
        TLorentzVector Stop1_cm,              Stop2_cm;
        TLorentzVector Stop1_PtRank_cm,       Stop2_PtRank_cm;
        TLorentzVector Stop1_MassRank_cm,     Stop2_MassRank_cm;
        TLorentzVector Stop1_ScalarPtRank_cm, Stop2_ScalarPtRank_cm;

        double Stop1ScalarPt              = -9999.9, Stop2ScalarPt = -9999.9;
        double Stop1ScalarPt_ScalarPtRank = -9999.9, Stop2ScalarPt_ScalarPtRank = -9999.9;
        double MT2                        = 0.0;
        double dR_Stop1Stop2              = -1;
        double dPhi_Stop1Stop2            = -1;
        double difference_stopMasses      = -9999.9;
        double average_stopMasses         = -9999.9;
        double relativeDiff_stopMasses    = -9999.9;
        double MT2_cm                     = 0.0;
        double dR_Stop1Stop2_cm           = -1;
        double dPhi_Stop1Stop2_cm         = -1;
 
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
            Hemisphere hemi(px, py, pz, E, seedMethod_, hemi_association); // to get hemisphere jets
            std::vector<int> grouping = hemi.getGrouping();

            // Perform vector and scalar sum of mega jets
            for(unsigned int i=0; i < px.size(); ++i)
            {
                TLorentzVector obj;
                obj.SetPxPyPzE(px[i], py[i], pz[i], E[i]);

                if(grouping[i] == 1)
                {
                    //vector sum
                    Stop1 += obj;

                    //scalar sum
                    Stop1ScalarPt += obj.Pt();
                }
                else if(grouping[i] == 2)
                {
                    //vector sum
                    Stop2 += obj;

                    //scalar sum
                    Stop2ScalarPt += obj.Pt();
                }
            }

            // ---------------------------------------
            // -- Rank mega jets (Mass, Pt, Scalar Pt) 
            // ---------------------------------------
            orderVars(Stop1_PtRank,               Stop2_PtRank,               Stop1,         Stop2,         Stop1.Pt()    > Stop2.Pt()   );
            orderVars(Stop1_MassRank,             Stop2_MassRank,             Stop1,         Stop2,         Stop1.M()     > Stop2.M()    );
            orderVars(Stop1_ScalarPtRank,         Stop2_ScalarPtRank,         Stop1,         Stop2,         Stop1ScalarPt > Stop2ScalarPt);
            orderVars(Stop1ScalarPt_ScalarPtRank, Stop2ScalarPt_ScalarPtRank, Stop1ScalarPt, Stop2ScalarPt, Stop1ScalarPt > Stop2ScalarPt);                

            MT2                     = ttUtility::coreMT2calc(Stop1, Stop2, MET);
            dR_Stop1Stop2           = Stop1.DeltaR(Stop2);    
            dPhi_Stop1Stop2         = Stop1.DeltaPhi(Stop2);
            difference_stopMasses   = abs( Stop1.M() - Stop2.M() );
            average_stopMasses      = ( 0.5 * ( Stop1.M() + Stop2.M() ) );
            relativeDiff_stopMasses = difference_stopMasses / average_stopMasses;

            // ----------------------------------------------------
            // -- Boost to CM frame & Rotate the hemispheres for NN 
            // ----------------------------------------------------
            Stop1_cm = Stop1;
            Stop2_cm = Stop2;
        
            TVector3 boostVec(0.0, 0.0, -event_beta_z);  
    
            Stop1_cm.Boost(boostVec);
            Stop1_cm.RotateZ(-event_phi_rotate);
            
            Stop2_cm.Boost(boostVec);
            Stop2_cm.RotateZ(-event_phi_rotate); 

            // Rank the hemispheres (Mass, Pt, Scalar Pt) 
            orderVars(Stop1_PtRank_cm,       Stop2_PtRank_cm,       Stop1_cm, Stop2_cm, Stop1_cm.Pt() > Stop2_cm.Pt());
            orderVars(Stop1_MassRank_cm,     Stop2_MassRank_cm,     Stop1_cm, Stop2_cm, Stop1_cm.M()  > Stop2_cm.M() );
            orderVars(Stop1_ScalarPtRank_cm, Stop2_ScalarPtRank_cm, Stop1_cm, Stop2_cm, Stop1ScalarPt > Stop2ScalarPt);

            MT2_cm             = ttUtility::coreMT2calc(Stop1_cm, Stop2_cm, MET);
            dR_Stop1Stop2_cm   = Stop1_cm.DeltaR(Stop2_cm);
            dPhi_Stop1Stop2_cm = Stop1_cm.DeltaPhi(Stop2_cm);

        }
        // without any rank
        tr.registerDerivedVar("MT2_cm"+myVarSuffix_,MT2_cm);
        tr.registerDerivedVar("dR_Stop1Stop2_cm"+myVarSuffix_,dR_Stop1Stop2_cm);
        tr.registerDerivedVar("dPhi_Stop1Stop2_cm"+myVarSuffix_,dPhi_Stop1Stop2_cm);
        tr.registerDerivedVar("Stop1_mass_cm"+myVarSuffix_,Stop1_cm.M());
        tr.registerDerivedVar("Stop2_mass_cm"+myVarSuffix_,Stop2_cm.M());
        tr.registerDerivedVar("Stop1_pt_cm"+myVarSuffix_,Stop1_cm.Pt());
        tr.registerDerivedVar("Stop2_pt_cm"+myVarSuffix_,Stop2_cm.Pt());
        tr.registerDerivedVar("Stop1_phi_cm"+myVarSuffix_,Stop1_cm.Phi());
        tr.registerDerivedVar("Stop2_phi_cm"+myVarSuffix_,Stop2_cm.Phi());
        tr.registerDerivedVar("Stop1_eta_cm"+myVarSuffix_,Stop1_cm.Eta());
        tr.registerDerivedVar("Stop2_eta_cm"+myVarSuffix_,Stop2_cm.Eta());
        tr.registerDerivedVar("Stop1_scalarPt_cm"+myVarSuffix_,Stop1ScalarPt);
        tr.registerDerivedVar("Stop2_scalarPt_cm"+myVarSuffix_,Stop2ScalarPt);
        // Pt Rank
        tr.registerDerivedVar("Stop1_PtRank"+myVarSuffix_,Stop1_PtRank);
        tr.registerDerivedVar("Stop2_PtRank"+myVarSuffix_,Stop2_PtRank);
        tr.registerDerivedVar("Stop1_mass_PtRank_cm"+myVarSuffix_,Stop1_PtRank_cm.M()); 
        tr.registerDerivedVar("Stop2_mass_PtRank_cm"+myVarSuffix_,Stop2_PtRank_cm.M());
        tr.registerDerivedVar("Stop1_pt_PtRank_cm"+myVarSuffix_,Stop1_PtRank_cm.Pt());
        tr.registerDerivedVar("Stop2_pt_PtRank_cm"+myVarSuffix_,Stop2_PtRank_cm.Pt());
        tr.registerDerivedVar("Stop1_phi_PtRank_cm"+myVarSuffix_,Stop1_PtRank_cm.Phi());
        tr.registerDerivedVar("Stop2_phi_PtRank_cm"+myVarSuffix_,Stop2_PtRank_cm.Phi());
        tr.registerDerivedVar("Stop1_eta_PtRank_cm"+myVarSuffix_,Stop1_PtRank_cm.Eta());
        tr.registerDerivedVar("Stop2_eta_PtRank_cm"+myVarSuffix_,Stop2_PtRank_cm.Eta());
        // Mask Rank
        tr.registerDerivedVar("Stop1_MassRank"+myVarSuffix_,Stop1_MassRank);
        tr.registerDerivedVar("Stop2_MassRank"+myVarSuffix_,Stop2_MassRank);
        tr.registerDerivedVar("Stop1_mass_MassRank_cm"+myVarSuffix_,Stop1_MassRank_cm.M()); 
        tr.registerDerivedVar("Stop2_mass_MassRank_cm"+myVarSuffix_,Stop2_MassRank_cm.M());
        tr.registerDerivedVar("Stop1_pt_MassRank_cm"+myVarSuffix_,Stop1_MassRank_cm.Pt());
        tr.registerDerivedVar("Stop2_pt_MassRank_cm"+myVarSuffix_,Stop2_MassRank_cm.Pt());
        tr.registerDerivedVar("Stop1_phi_MassRank_cm"+myVarSuffix_,Stop1_MassRank_cm.Phi());
        tr.registerDerivedVar("Stop2_phi_MassRank_cm"+myVarSuffix_,Stop2_MassRank_cm.Phi());
        tr.registerDerivedVar("Stop1_eta_MassRank_cm"+myVarSuffix_,Stop1_MassRank_cm.Eta());
        tr.registerDerivedVar("Stop2_eta_MassRank_cm"+myVarSuffix_,Stop2_MassRank_cm.Eta());
        // ScalarPt Rank
        tr.registerDerivedVar("Stop1_ScalarPtRank"+myVarSuffix_,Stop1_ScalarPtRank);
        tr.registerDerivedVar("Stop2_ScalarPtRank"+myVarSuffix_,Stop2_ScalarPtRank);
        tr.registerDerivedVar("Stop1_mass_ScalarPtRank_cm"+myVarSuffix_,Stop1_ScalarPtRank_cm.M());
        tr.registerDerivedVar("Stop2_mass_ScalarPtRank_cm"+myVarSuffix_,Stop2_ScalarPtRank_cm.M());
        tr.registerDerivedVar("Stop1_pt_ScalarPtRank_cm"+myVarSuffix_,Stop1_ScalarPtRank_cm.Pt());
        tr.registerDerivedVar("Stop2_pt_ScalarPtRank_cm"+myVarSuffix_,Stop2_ScalarPtRank_cm.Pt());
        tr.registerDerivedVar("Stop1_phi_ScalarPtRank_cm"+myVarSuffix_,Stop1_ScalarPtRank_cm.Phi());
        tr.registerDerivedVar("Stop2_phi_ScalarPtRank_cm"+myVarSuffix_,Stop2_ScalarPtRank_cm.Phi());
        tr.registerDerivedVar("Stop1_eta_ScalarPtRank_cm"+myVarSuffix_,Stop1_ScalarPtRank_cm.Eta());
        tr.registerDerivedVar("Stop2_eta_ScalarPtRank_cm"+myVarSuffix_,Stop2_ScalarPtRank_cm.Eta());
        // others
        tr.registerDerivedVar("MT2"+myVarSuffix_,MT2);
        tr.registerDerivedVar("dR_Stop1Stop2"+myVarSuffix_,dR_Stop1Stop2);
        tr.registerDerivedVar("dPhi_Stop1Stop2"+myVarSuffix_,dPhi_Stop1Stop2);
        tr.registerDerivedVar("difference_stopMasses"+myVarSuffix_,difference_stopMasses);
        tr.registerDerivedVar("average_stopMasses"+myVarSuffix_,average_stopMasses);
        tr.registerDerivedVar("relativeDiff_stopMasses"+myVarSuffix_,relativeDiff_stopMasses);
        tr.registerDerivedVar("Stop1ScalarPt_ScalarPtRank"+myVarSuffix_,Stop1ScalarPt_ScalarPtRank);
        tr.registerDerivedVar("Stop2ScalarPt_ScalarPtRank"+myVarSuffix_,Stop2ScalarPt_ScalarPtRank);
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
