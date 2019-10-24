#ifndef MakeMT2Hemispheres_h
#define MakeMT2Hemispheres_h

#include "Framework/Framework/include/MT2Hemispheres.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h" 
#include <vector>
#include <iostream>
#include <cmath>

class MakeMT2Hemispheres
{
private:
    std::string jetMaskName_, nJetName_, myVarSuffix_;

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
        const auto& met                = tr.getVar<double>("MET");
        const auto& metPhi             = tr.getVar<double>("METPhi");
        const auto& Jets               = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets           = tr.getVec<bool>(jetMaskName_);
        const auto& NGoodJets          = tr.getVar<int>(nJetName_);
        const auto& GoodLeptons        = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");

        static const int hemi_association = 3; // 3: 3th method, 'lund' used by MT2  
        TLorentzVector stop1_PtRank,       stop2_PtRank;
        TLorentzVector stop1_MassRank,     stop2_MassRank;
        TLorentzVector stop1_ScalarPtRank, stop2_ScalarPtRank;
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
            Hemisphere hemi(px, py, pz, E, 2, hemi_association);
            std::vector<int> grouping = hemi.getGrouping();
            TLorentzVector pseudojet1, pseudojet2;
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
        }
        tr.registerDerivedVar("stop1_PtRank"+myVarSuffix_,stop1_PtRank);
        tr.registerDerivedVar("stop2_PtRank"+myVarSuffix_,stop2_PtRank);
        tr.registerDerivedVar("stop1_MassRank"+myVarSuffix_,stop1_MassRank);
        tr.registerDerivedVar("stop2_MassRank"+myVarSuffix_,stop2_MassRank);
        tr.registerDerivedVar("stop1_ScalarPtRank"+myVarSuffix_,stop1_ScalarPtRank);
        tr.registerDerivedVar("stop2_ScalarPtRank"+myVarSuffix_,stop2_ScalarPtRank);
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
