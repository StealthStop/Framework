#ifndef FATJETCOMBINE_H
#define FATJETCOMBINE_H

#include "Framework/Framework/include/Utility.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/interface/lester_mt2_bisect.h"

class FatJetCombine
{
private:
    std::string myVarSuffix_;

    void fatjetcombine(NTupleReader& tr)
    {
        const auto& JetsAK8               = tr.getVec<TLorentzVector>("JetsAK8");
        const auto& Tau1                  = tr.getVec<double>("JetsAK8_NsubjettinessTau1");
        const auto& Tau2                  = tr.getVec<double>("JetsAK8_NsubjettinessTau2");
        const auto& Tau3                  = tr.getVec<double>("JetsAK8_NsubjettinessTau3");
        const auto& softDropMass          = tr.getVec<double>("JetsAK8_softDropMass");
        const auto& Muons                 = tr.getVec<TLorentzVector>("Muons");
        const auto& GoodMuons             = tr.getVec<bool>("GoodMuons"+myVarSuffix_);
        const auto& Electrons             = tr.getVec<TLorentzVector>("Electrons");
        const auto& GoodElectrons         = tr.getVec<bool>("GoodElectrons"+myVarSuffix_);
        const auto& NMuons                = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        const auto& NElectrons            = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);

        const auto& TwoLep_Mbl1_Idx      = tr.getVar<std::pair<unsigned int, unsigned int>>("TwoLep_Mbl1_Idx"+myVarSuffix_);
        const auto& TwoLep_Mbl2_Idx      = tr.getVar<std::pair<unsigned int, unsigned int>>("TwoLep_Mbl2_Idx"+myVarSuffix_);
        
        const auto& Jets                 = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
        const auto& GoodLeptons          = tr.getVec<std::pair<std::string,TLorentzVector>>("GoodLeptons"+myVarSuffix_);
        const auto& NGoodLeptons         = tr.getVar<int>("NGoodLeptons"+myVarSuffix_);

        const auto& subjets              = tr.getVec<std::vector<TLorentzVector>>("JetsAK8_subjets"+myVarSuffix_);
        const auto& subjets_CSV          = tr.getVec<std::vector<double>>("JetsAK8_subjets_bDiscriminatorCSV"+myVarSuffix_);
//        const auto& medium_CSV           = tr.getVar<double>("deepCSV_WP_medium"+myVarSuffix_);
//        const auto& loose_CSV           = tr.getVar<double>("deepCSV_WP_loose"+myVarSuffix_);

        const auto& MET                  = tr.getVar<double>("MET"+myVarSuffix_);
        const auto& METPhi               = tr.getVar<double>("METPhi"+myVarSuffix_);
        
        //First clean out leptons from JetsAK8 collection

        std::vector<bool> initGoodJets(JetsAK8.size(), true);
        auto& GoodJetsAK8         = tr.createDerivedVec<bool>("GoodJetsAK8"+myVarSuffix_);
        GoodJetsAK8 = initGoodJets;
        
        if(NMuons > 0)
        {
            for(unsigned int mu=0; mu < Muons.size(); mu++)
            {
                double minDeltaR = 0.8;
                TLorentzVector myMuon = Muons.at(mu);
                int muonCand = -1;
                for(unsigned int j=0; j < JetsAK8.size(); j++)
                {
                    TLorentzVector myJet = JetsAK8.at(j);
                    if( std::fabs(myMuon.Pt() - myJet.Pt()) / myMuon.Pt() < 1 && myMuon.DeltaR(myJet) < minDeltaR)
                    {
                        minDeltaR = myMuon.DeltaR(myJet);
                        muonCand = j;
                    }
                }
                if(GoodMuons.at(mu) && muonCand != -1) GoodJetsAK8.at(muonCand) = false;
            }
        }

        if(NElectrons > 0)
        {
            for(unsigned int el=0; el < Electrons.size(); el++)
            {
                double minDeltaR = 0.8;
                TLorentzVector myElec = Electrons.at(el);
                int elecCand = -1;
                for(unsigned int j=0; j < JetsAK8.size(); j++)
                {
                    TLorentzVector myJet = JetsAK8.at(j);
                    if( std::fabs(myElec.Pt() - myJet.Pt()) / myElec.Pt() < 1 && myElec.DeltaR(myJet) < minDeltaR)
                    {
                        minDeltaR = myElec.DeltaR(myJet);
                        elecCand = j;                  
                    }
                }
                if(GoodElectrons.at(el) && elecCand != -1) GoodJetsAK8.at(elecCand) = false;
            }
        }
        int NGoodJetsAK8 = 0;

        //Now apply eta and pt cuts
        for(unsigned int j=0; j < JetsAK8.size(); j++)
        {
            if(abs(JetsAK8.at(j).Eta()) > 2.4) GoodJetsAK8.at(j) = false;
            if(JetsAK8.at(j).Pt() < 30) GoodJetsAK8.at(j) = false;
            if(GoodJetsAK8.at(j)) NGoodJetsAK8 += 1;
        }
        
        //Begin fat jet reconstruction for 2L
        TLorentzVector CombinedJet1,CombinedJet2;
        double CombinedJet1_T3T1=0,CombinedJet1_T2T1=0,CombinedJet1_SDM=0;
        double CombinedJet2_T3T1=0,CombinedJet2_T2T1=0,CombinedJet2_SDM=0;

        TLorentzVector NlinoJet1,NlinoJet2;

        bool B1Found = false, B2Found = false;
        bool B1_Ovl = false, B2_Ovl = false;
        int nBottom1Ovl=0,nBottom2Ovl=0;
        int nNlinoCand = 0;
        TLorentzVector lvMET;
        lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);


        auto& BOvlDeltaR        = tr.createDerivedVec<double>("BOvlDeltaR"+myVarSuffix_);
        auto& BOvlPt            = tr.createDerivedVec<double>("BOvlPt"+myVarSuffix_);
        auto& BOvlDCSV          = tr.createDerivedVec<double>("BOvlDCSV"+myVarSuffix_);
         
        if (NGoodLeptons == 2 && NGoodJetsAK8 <= 4 )
        {
            TLorentzVector Bottom1 = Jets[TwoLep_Mbl1_Idx.first];
            TLorentzVector Bottom2 = Jets[TwoLep_Mbl2_Idx.first];
            TLorentzVector BLVec_1 = Bottom1 + GoodLeptons[TwoLep_Mbl1_Idx.second].second;
            TLorentzVector BLVec_2 = Bottom2 + GoodLeptons[TwoLep_Mbl2_Idx.second].second;
            std::vector<TLorentzVector> GoodJetsAK8Vec;
            std::vector<double> GoodJetsAK8SDM,GoodJetsAK8T3T1,GoodJetsAK8T2T1;
            std::vector<bool> NlinoCandidate = GoodJetsAK8;

            int Bottom1AK8Cand,Bottom2AK8Cand;
            
            if (NGoodJetsAK8 > 2)
            {
                double minDR1 = 0.2,minDR2 = 0.2;
                for (unsigned int j=0; j < JetsAK8.size(); j++)
                {
                    if (GoodJetsAK8.at(j) && JetsAK8.at(j).DeltaR(Bottom1) < minDR1 && abs(1 - JetsAK8.at(j).Pt()/Bottom1.Pt()) < 0.2)
                    {
                        Bottom1AK8Cand = j;
                        minDR1 = JetsAK8.at(j).DeltaR(Bottom1);
                        B1Found = true;
                    }
                    else if (GoodJetsAK8.at(j) && JetsAK8.at(j).DeltaR(Bottom2) < minDR2 && abs(1 - JetsAK8.at(j).Pt()/Bottom2.Pt()) < 0.2)
                    {
                        Bottom2AK8Cand = j;
                        minDR2 = JetsAK8.at(j).DeltaR(Bottom2);
                        B2Found = true;
                    }
                }
                if (B1Found) NlinoCandidate.at(Bottom1AK8Cand) = false;
                if (B2Found) NlinoCandidate.at(Bottom2AK8Cand) = false;
            }
            
            for (unsigned int n = 0; n < NlinoCandidate.size(); n++)
            {
                if (NlinoCandidate.at(n)) nNlinoCand += 1;
            }

            if (nNlinoCand == 2)
            {
                double minDR1 = 0.2, minDR2 = 0.2;
                std::pair<int, int> likelyB1 (-1, -1) , likelyB2 (-1, -1);            
             
                for(unsigned int j=0; j < JetsAK8.size(); j++)
                {
                    if(NlinoCandidate.at(j))
                    {
                        GoodJetsAK8Vec.push_back(JetsAK8.at(j));
                        GoodJetsAK8SDM.push_back(softDropMass.at(j));
                        double Tau3Tau1 = Tau3.at(j) / Tau1.at(j);
                        double Tau2Tau1 = Tau2.at(j) / Tau1.at(j);
                        GoodJetsAK8T3T1.push_back(Tau3Tau1);
                        GoodJetsAK8T2T1.push_back(Tau2Tau1);
                        
                        for(unsigned int s=0; s < subjets.at(j).size(); s++)
                        {
                            TLorentzVector sj = subjets.at(j).at(s);
                            if (!B1Found && sj.DeltaR(Bottom1) < minDR1 && abs( 1 - sj.Pt()/Bottom1.Pt()) < 0.2 && subjets_CSV.at(j).at(s) > 0) 
                            {
                                nBottom1Ovl += 1;
                                likelyB1.first = j;
                                likelyB1.second = s;
                                minDR1 = sj.DeltaR(Bottom1);
                            }
                            if (!B2Found && sj.DeltaR(Bottom2) < 0.2 && abs( 1 - sj.Pt()/Bottom2.Pt()) < 0.2 && subjets_CSV.at(j).at(s) > 0)
                            {
                                likelyB2.first = j;
                                likelyB2.second = s;
                                nBottom2Ovl += 1;
                                minDR2 = sj.DeltaR(Bottom2);
                            }                                                 
                        }
                        
                        if (likelyB1.first != -1)
                        {
                            B1Found = true;
                            BOvlDeltaR.push_back(minDR1);                      
                            BOvlPt.push_back( subjets.at(likelyB1.first).at(likelyB1.second).Pt() / Bottom1.Pt());
                            BOvlDCSV.push_back( subjets_CSV.at(likelyB1.first).at(likelyB1.second));
                        }
                        if (likelyB2.first != -1)
                        {
                            B2Found = true;
                            BOvlDeltaR.push_back(minDR2);                      
                            BOvlPt.push_back( subjets.at(likelyB2.first).at(likelyB2.second).Pt() / Bottom2.Pt());
                            BOvlDCSV.push_back( subjets_CSV.at(likelyB2.first).at(likelyB2.second));
                        }
                    }
                }
                if (likelyB1.first != -1)
                {
                    B1_Ovl = true;                
                    BOvlDeltaR.push_back(minDR1);                      
                    BOvlPt.push_back( subjets.at(likelyB1.first).at(likelyB1.second).Pt() / Bottom1.Pt());
                    BOvlDCSV.push_back( subjets_CSV.at(likelyB1.first).at(likelyB1.second));
                }
                if (likelyB2.first != -1)
                {
                    B2_Ovl = true;
                    BOvlDeltaR.push_back(minDR2);                      
                    BOvlPt.push_back( subjets.at(likelyB2.first).at(likelyB2.second).Pt() / Bottom2.Pt());
                    BOvlDCSV.push_back( subjets_CSV.at(likelyB2.first).at(likelyB2.second));
                }
            
                TLorentzVector poss1_CombinedJet1 = BLVec_1 + GoodJetsAK8Vec.at(0);
                TLorentzVector poss1_CombinedJet2 = BLVec_2 + GoodJetsAK8Vec.at(1);
                TLorentzVector poss2_CombinedJet1 = BLVec_1 + GoodJetsAK8Vec.at(1);
                TLorentzVector poss2_CombinedJet2 = BLVec_2 + GoodJetsAK8Vec.at(0);
                
                TLorentzVector poss1_NlinoJet1 = GoodJetsAK8Vec.at(0);
                TLorentzVector poss1_NlinoJet2 = GoodJetsAK8Vec.at(1);
                TLorentzVector poss2_NlinoJet1 = GoodJetsAK8Vec.at(1);
                TLorentzVector poss2_NlinoJet2 = GoodJetsAK8Vec.at(0);
                
                if (B1_Ovl)
                {
                    poss1_CombinedJet1 = GoodLeptons[TwoLep_Mbl1_Idx.second].second + GoodJetsAK8Vec.at(0);
                    poss2_CombinedJet1 =  GoodLeptons[TwoLep_Mbl1_Idx.second].second + GoodJetsAK8Vec.at(1);
                    poss1_NlinoJet1 = GoodJetsAK8Vec.at(0) - Bottom1;
                    poss2_NlinoJet1 = GoodJetsAK8Vec.at(1) - Bottom1;
                }
                if (B2_Ovl)
                {
                    poss1_CombinedJet2 = GoodLeptons[TwoLep_Mbl2_Idx.second].second + GoodJetsAK8Vec.at(1);
                    poss2_CombinedJet2 = GoodLeptons[TwoLep_Mbl2_Idx.second].second + GoodJetsAK8Vec.at(0);
                    poss1_NlinoJet2 = GoodJetsAK8Vec.at(1) - Bottom2;
                    poss2_NlinoJet2 = GoodJetsAK8Vec.at(0) - Bottom2;                    
                }

                if(abs(poss1_CombinedJet1.M() - poss1_CombinedJet2.M()) < abs(poss2_CombinedJet1.M() - poss2_CombinedJet2.M()))
                {
                    CombinedJet1 = poss1_CombinedJet1;
                    CombinedJet1_T3T1 = GoodJetsAK8T3T1.at(0);
                    CombinedJet1_T2T1 = GoodJetsAK8T2T1.at(0);
                    CombinedJet1_SDM = GoodJetsAK8SDM.at(0);
                    
                    NlinoJet1 = poss1_NlinoJet1;

                    CombinedJet2 = poss1_CombinedJet2;
                    CombinedJet2_T3T1 = GoodJetsAK8T3T1.at(1);
                    CombinedJet2_T2T1 = GoodJetsAK8T2T1.at(1);
                    CombinedJet2_SDM = GoodJetsAK8SDM.at(1);
                    
                    NlinoJet2 = poss1_NlinoJet2;
            
                }
                else
                {
                    CombinedJet1 = poss2_CombinedJet1;
                    CombinedJet1_T3T1 = GoodJetsAK8T3T1.at(1);
                    CombinedJet1_T2T1 = GoodJetsAK8T2T1.at(1);
                    CombinedJet1_SDM = GoodJetsAK8SDM.at(1);
                    
                    NlinoJet1 = poss2_NlinoJet1;
                    
                    CombinedJet2 = poss2_CombinedJet2;
                    CombinedJet2_T3T1 = GoodJetsAK8T3T1.at(0);
                    CombinedJet2_T2T1 = GoodJetsAK8T2T1.at(0);
                    CombinedJet2_SDM = GoodJetsAK8SDM.at(0);
                    
                    NlinoJet2 = poss2_NlinoJet2;

                }
            }
        }
        
        asymm_mt2_lester_bisect::disableCopyrightMessage();
        
        tr.registerDerivedVar("FatJetCombined1"+myVarSuffix_, CombinedJet1);
        tr.registerDerivedVar("FatJetCombined2"+myVarSuffix_, CombinedJet2);
        tr.registerDerivedVar("FatJetMT2"+myVarSuffix_, ttUtility::coreMT2calc(CombinedJet1,CombinedJet2,lvMET));
        tr.registerDerivedVar("FatJetCombined1_SDM"+myVarSuffix_, CombinedJet1_SDM);
        tr.registerDerivedVar("FatJetCombined2_SDM"+myVarSuffix_, CombinedJet2_SDM);
        tr.registerDerivedVar("FatJetCombined1_T3T1"+myVarSuffix_, CombinedJet1_T3T1);
        tr.registerDerivedVar("FatJetCombined2_T3T1"+myVarSuffix_, CombinedJet2_T3T1);
        tr.registerDerivedVar("FatJetCombined1_T2T1"+myVarSuffix_, CombinedJet1_T2T1);
        tr.registerDerivedVar("FatJetCombined2_T2T1"+myVarSuffix_, CombinedJet2_T2T1);

        tr.registerDerivedVar("FatJetNlino1"+myVarSuffix_, NlinoJet1);
        tr.registerDerivedVar("FatJetNlino2"+myVarSuffix_, NlinoJet2);

        tr.registerDerivedVar("nNlinoCand"+myVarSuffix_, nNlinoCand);
        
        tr.registerDerivedVar("BOvl_B1Filter"+myVarSuffix_, B1Found);
        tr.registerDerivedVar("BOvl_B2Filter"+myVarSuffix_, B2Found);
        
        tr.registerDerivedVar("NGoodJetsAK8"+myVarSuffix_, NGoodJetsAK8);
        tr.registerDerivedVar("nBottom1Ovl"+myVarSuffix_, nBottom1Ovl);
        tr.registerDerivedVar("nBottom2Ovl"+myVarSuffix_, nBottom2Ovl);
        
    }
    

public:
    FatJetCombine(std::string myVarSuffix = "") 
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up FatJetCombine"<<std::endl;   
    }
    
    void operator()(NTupleReader& tr)
    {
        fatjetcombine(tr);
    }
};

#endif
