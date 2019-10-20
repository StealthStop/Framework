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

    void getHemispheres(NTupleReader& tr) const
    {
        // these variables values the same as  MT2
        const int hemi_association = 3; // 3: 3th method, 'lund' used by MT2  

        const auto& met                         = tr.getVar<double>("MET");
        const auto& metPhi                      = tr.getVar<double>("METPhi");
        const auto& Jets                        = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets                    = tr.getVec<bool>(jetMaskName_);
        const auto& NGoodJets                   = tr.getVar<int>(nJetName_);
        const auto& GoodLeptons                 = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        double MT2                              = 0.0;
        double stop1Mass                        = -9999.9, stop1Eta = -9999.9, stop1Phi = -9999.9, stop1Pt = -9999.9;
        double stop2Mass                        = -9999.9, stop2Eta = -9999.9, stop2Phi = -9999.9, stop2Pt = -9999.9;
        double dR_stop1stop2                    = -1;
        double dPhi_stop1stop2                  = -1;
        double stop1Mass_PtRank                 = -9999.9, stop2Mass_PtRank   = -9999.9; 
        double stop1Mass_MassRank               = -9999.9, stop2Mass_MassRank = -9999.9;
        double difference_stopMasses_PtRank     = -9999.9;
        double relativeDiff_stopMasses_PtRank   = -9999.9;
        double difference_stopMasses_MassRank   = -9999.9;
        double relativeDiff_stopMasses_MassRank = -9999.9;        
        double difference_stopMasses            = -9999.9;
        double average_stopMasses               = -9999.9;
        double relativeDiff_stopMasses          = -9999.9;
 
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
                
            // --------------------------
            // -- stop MT2 stopspheres 
            // --------------------------
            MT2             = ttUtility::coreMT2calc(pseudojet1, pseudojet2, MET);
            stop1Mass       = pseudojet1.M();
            stop1Eta        = pseudojet1.Eta();
            stop1Phi        = pseudojet1.Phi();
            stop1Pt         = pseudojet1.Pt();
            stop2Mass       = pseudojet2.M();
            stop2Eta        = pseudojet2.Eta();
            stop2Phi        = pseudojet2.Phi();
            stop2Pt         = pseudojet2.Pt(); 
            dR_stop1stop2   = pseudojet1.DeltaR(pseudojet2);    
            dPhi_stop1stop2 = pseudojet1.DeltaPhi(pseudojet2);
            
            // stopMass Pt rank
            if (pseudojet1.Pt() > pseudojet2.Pt()) 
            {
                stop1Mass_PtRank               = pseudojet1.M();
                stop2Mass_PtRank               = pseudojet2.M();      
                difference_stopMasses_PtRank   = ( pseudojet1.M() - pseudojet2.M() );
                relativeDiff_stopMasses_PtRank = ( pseudojet1.M() - pseudojet2.M() ) / ( 0.5 * ( pseudojet1.M() + pseudojet2.M() ) );     
            } 
            else
            {
                stop1Mass_PtRank               = pseudojet2.M();
                stop2Mass_PtRank               = pseudojet1.M();
                difference_stopMasses_PtRank   = ( pseudojet2.M() - pseudojet1.M() );
                relativeDiff_stopMasses_PtRank = ( pseudojet2.M() - pseudojet1.M() ) / ( 0.5 * ( pseudojet1.M() + pseudojet2.M() ) );

            }
    
            // stopMass Mass rank
            if (pseudojet1.M() > pseudojet2.M())
            {
                stop1Mass_MassRank               = pseudojet1.M();
                stop2Mass_MassRank               = pseudojet2.M();
                difference_stopMasses_MassRank   = ( pseudojet1.M() - pseudojet2.M() );
                relativeDiff_stopMasses_MassRank = ( pseudojet1.M() - pseudojet2.M() ) / ( 0.5 * ( pseudojet1.M() + pseudojet2.M() ) );
            } 
            else
            {
                stop1Mass_MassRank               = pseudojet2.M();
                stop2Mass_MassRank               = pseudojet1.M();
                difference_stopMasses_MassRank   = ( pseudojet2.M() - pseudojet1.M() );
                relativeDiff_stopMasses_MassRank = ( pseudojet2.M() - pseudojet1.M() ) / ( 0.5 * ( pseudojet1.M() + pseudojet2.M() ) );
            }
        
            // stop Masses difference & average & relative difference
            difference_stopMasses   = ( pseudojet1.M() - pseudojet2.M() );
            average_stopMasses      = ( 0.5 * ( pseudojet1.M() + pseudojet2.M() ) );
            relativeDiff_stopMasses = ( pseudojet1.M() - pseudojet2.M() ) / ( 0.5 * ( pseudojet1.M() + pseudojet2.M() ) );

        }
        tr.registerDerivedVar("MT2"+myVarSuffix_,MT2);
        tr.registerDerivedVar("stop1Mass"+myVarSuffix_,stop1Mass);
        tr.registerDerivedVar("stop1Eta"+myVarSuffix_,stop1Eta);
        tr.registerDerivedVar("stop1Phi"+myVarSuffix_,stop1Phi);
        tr.registerDerivedVar("stop1Pt"+myVarSuffix_,stop1Pt);
        tr.registerDerivedVar("stop2Mass"+myVarSuffix_,stop2Mass);
        tr.registerDerivedVar("stop2Eta"+myVarSuffix_,stop2Eta);
        tr.registerDerivedVar("stop2Phi"+myVarSuffix_,stop2Phi);
        tr.registerDerivedVar("stop2Pt"+myVarSuffix_,stop2Pt);
        tr.registerDerivedVar("dR_stop1stop2"+myVarSuffix_,dR_stop1stop2);
        tr.registerDerivedVar("dPhi_stop1stop2"+myVarSuffix_,dPhi_stop1stop2);
        tr.registerDerivedVar("stop1Mass_PtRank"+myVarSuffix_,stop1Mass_PtRank);
        tr.registerDerivedVar("stop2Mass_PtRank"+myVarSuffix_,stop2Mass_PtRank);
        tr.registerDerivedVar("stop1Mass_MassRank"+myVarSuffix_,stop1Mass_MassRank);
        tr.registerDerivedVar("stop2Mass_MassRank"+myVarSuffix_,stop2Mass_MassRank);
        tr.registerDerivedVar("difference_stopMasses_PtRank"+myVarSuffix_,difference_stopMasses_PtRank);
        tr.registerDerivedVar("relativeDiff_stopMasses_PtRank"+myVarSuffix_,relativeDiff_stopMasses_PtRank);
        tr.registerDerivedVar("difference_stopMasses_MassRank"+myVarSuffix_,difference_stopMasses_MassRank);
        tr.registerDerivedVar("relativeDiff_stopMasses_MassRank"+myVarSuffix_,relativeDiff_stopMasses_MassRank);
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
