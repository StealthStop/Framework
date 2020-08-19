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
        const auto& prunedMass            = tr.getVec<double>("JetsAK8_prunedMass");
        const auto& Muons                 = tr.getVec<TLorentzVector>("Muons");
//        const auto& GoodMuons             = tr.getVec<bool>("GoodMuons"+myVarSuffix_);
        const auto& Electrons             = tr.getVec<TLorentzVector>("Electrons");
//        const auto& GoodElectrons         = tr.getVec<bool>("GoodElectrons"+myVarSuffix_);
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
        for (unsigned int j = 0; j < JetsAK8.size(); j++)
        {
            GoodJetsAK8.push_back(true);
        }                
        
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
//                    if( std::fabs(myMuon.Pt() - myJet.Pt()) / myMuon.Pt() < 1 && myMuon.DeltaR(myJet) < minDeltaR)
                    if(myMuon.DeltaR(myJet) < minDeltaR)
                    {
                        minDeltaR = myMuon.DeltaR(myJet);
                        muonCand = j;
                    }
                }
                if(muonCand != -1) GoodJetsAK8.at(muonCand) = false;
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
                    // if( std::fabs(myElec.Pt() - myJet.Pt()) / myElec.Pt() < 1 && myElec.DeltaR(myJet) < minDeltaR)
                    if(myElec.DeltaR(myJet) < minDeltaR)
                    {
                        minDeltaR = myElec.DeltaR(myJet);
                        elecCand = j;                  
                    }
                }
                if(elecCand != -1) GoodJetsAK8.at(elecCand) = false;
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
        double CombinedJet1_T3T1=0,CombinedJet1_T2T1=0,CombinedJet1_T3T2=0,CombinedJet1_SDM=0,CombinedJet1_Pruned = 0;
        double CombinedJet2_T3T1=0,CombinedJet2_T2T1=0,CombinedJet2_T3T2=0,CombinedJet2_SDM=0,CombinedJet2_Pruned = 0;
        TLorentzVector CombinedJet1_MaxSubjet,CombinedJet2_MaxSubjet,CombinedJet1_MinSubjet,CombinedJet2_MinSubjet;

        TLorentzVector NlinoJet1,NlinoJet2;

        bool B1Found = false, B2Found = false;
        bool B1_Ovl = false, B2_Ovl = false;
        int nBottom1Ovl=0,nBottom2Ovl=0;
        int nNlinoCand = 0;
        TLorentzVector lvMET;
        lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);
        int FatJet1_Idx = -1, FatJet2_Idx = -1;
        //std::vector<int> FatJet_Idx;
        auto& BOvlDeltaR        = tr.createDerivedVec<double>("BOvlDeltaR"+myVarSuffix_);
        auto& BOvlPt            = tr.createDerivedVec<double>("BOvlPt"+myVarSuffix_);
        auto& BOvlDCSV          = tr.createDerivedVec<double>("BOvlDCSV"+myVarSuffix_);
        auto& NlinoNSubjets     = tr.createDerivedVec<int>("NlinoNSubjets"+myVarSuffix_);        
         
        if (NGoodLeptons == 2)// && NGoodJetsAK8 <= 4 )
        {
            TLorentzVector Bottom1 = Jets[TwoLep_Mbl1_Idx.first];
            TLorentzVector Bottom2 = Jets[TwoLep_Mbl2_Idx.first];
            TLorentzVector BLVec_1 = Bottom1 + GoodLeptons[TwoLep_Mbl1_Idx.second].second;
            TLorentzVector BLVec_2 = Bottom2 + GoodLeptons[TwoLep_Mbl2_Idx.second].second;
            std::vector<TLorentzVector> GoodJetsAK8Vec;
            std::vector<double> GoodJetsAK8SDM,GoodJetsAK8Pruned,GoodJetsAK8T3T1,GoodJetsAK8T2T1,GoodJetsAK8T3T2;
            std::vector<TLorentzVector> GoodJetsAK8MaxSubjet,GoodJetsAK8MinSubjet;
            std::vector<bool> NlinoCandidate = GoodJetsAK8;

            int Bottom1AK8Cand = -1,Bottom2AK8Cand = -1;
            
            if (NGoodJetsAK8 >= 2)
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
                    if (softDropMass.at(j) < 20) NlinoCandidate.at(j) = false;
                    if (subjets.at(j).size() == 0) NlinoCandidate.at(j) = false;
                }
                if (B1Found) NlinoCandidate.at(Bottom1AK8Cand) = false;
                if (B2Found) NlinoCandidate.at(Bottom2AK8Cand) = false;

            }
            
            for (unsigned int n = 0; n < NlinoCandidate.size(); n++)
            {
                if (NlinoCandidate.at(n)) nNlinoCand += 1;
                if (NlinoCandidate.at(n)) NlinoNSubjets.push_back(subjets.at(n).size());
            }
//            std::cout << nNlinoCand << std::endl;

            if (nNlinoCand == 2)
            {
                double minDR1 = 0.2, minDR2 = 0.2;
                std::pair<int, int> likelyB1 (-1, -1) , likelyB2 (-1, -1);
             
                for(unsigned int j=0; j < JetsAK8.size(); j++)
                {
                    if(NlinoCandidate.at(j))
                    {
                        
                        if (FatJet1_Idx != -1 && FatJet2_Idx == -1) 
                        {
                            FatJet2_Idx = j;
                        }
                        if (FatJet1_Idx == -1) 
                        {
                            FatJet1_Idx = j;
                        }
                        GoodJetsAK8Vec.push_back(JetsAK8.at(j));
                        GoodJetsAK8SDM.push_back(softDropMass.at(j));
                        GoodJetsAK8T3T1.push_back(Tau3.at(j) / Tau1.at(j));
                        GoodJetsAK8T2T1.push_back(Tau2.at(j) / Tau1.at(j));
                        GoodJetsAK8T3T2.push_back(Tau3.at(j) / Tau2.at(j));
                        GoodJetsAK8Pruned.push_back(prunedMass.at(j));

                        double minMass = 999;
                        double maxMass = 0;
                        int maxsj=-1,minsj=-1;
                        for(unsigned int s=0; s < subjets.at(j).size(); s++)
                        {
                            if (subjets.at(j).at(s).M() > maxMass)
                            {
                                maxMass = subjets.at(j).at(s).M();
                                maxsj = s;
                            }
                            if (subjets.at(j).at(s).M() < minMass)
                            {
                                minMass = subjets.at(j).at(s).M();
                                minsj = s;
                            }
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
                        if (maxsj == -1) GoodJetsAK8MaxSubjet.push_back(CombinedJet1_MaxSubjet);
                        else GoodJetsAK8MaxSubjet.push_back(subjets.at(j).at(maxsj));
                        if (minsj == -1) GoodJetsAK8MinSubjet.push_back(CombinedJet1_MinSubjet);
                        else GoodJetsAK8MinSubjet.push_back(subjets.at(j).at(minsj));
                        

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
                    CombinedJet1_T3T2 = GoodJetsAK8T3T2.at(0);
                    CombinedJet1_SDM = GoodJetsAK8SDM.at(0);
                    CombinedJet1_Pruned = GoodJetsAK8Pruned.at(0);
                    CombinedJet1_MaxSubjet = GoodJetsAK8MaxSubjet.at(0);
                    CombinedJet1_MinSubjet = GoodJetsAK8MinSubjet.at(0);
                    
                    NlinoJet1 = poss1_NlinoJet1;

                    CombinedJet2 = poss1_CombinedJet2;
                    CombinedJet2_T3T1 = GoodJetsAK8T3T1.at(1);
                    CombinedJet2_T2T1 = GoodJetsAK8T2T1.at(1);
                    CombinedJet2_T3T2 = GoodJetsAK8T3T2.at(1);
                    CombinedJet2_SDM = GoodJetsAK8SDM.at(1);
                    CombinedJet2_Pruned = GoodJetsAK8Pruned.at(1);
                    CombinedJet2_MaxSubjet = GoodJetsAK8MaxSubjet.at(1);
                    CombinedJet2_MinSubjet = GoodJetsAK8MinSubjet.at(1);

                    
                    NlinoJet2 = poss1_NlinoJet2;
            
                }
                else
                {
                    CombinedJet1 = poss2_CombinedJet1;
                    CombinedJet1_T3T1 = GoodJetsAK8T3T1.at(1);
                    CombinedJet1_T2T1 = GoodJetsAK8T2T1.at(1);
                    CombinedJet1_T3T2 = GoodJetsAK8T3T2.at(1);
                    CombinedJet1_SDM = GoodJetsAK8SDM.at(1);
                    CombinedJet1_Pruned = GoodJetsAK8Pruned.at(1);
                    CombinedJet1_MaxSubjet = GoodJetsAK8MaxSubjet.at(1);
                    CombinedJet1_MinSubjet = GoodJetsAK8MinSubjet.at(1);
                    
                    NlinoJet1 = poss2_NlinoJet1;
                    
                    CombinedJet2 = poss2_CombinedJet2;
                    CombinedJet2_T3T1 = GoodJetsAK8T3T1.at(0);
                    CombinedJet2_T2T1 = GoodJetsAK8T2T1.at(0);
                    CombinedJet2_T3T2 = GoodJetsAK8T3T2.at(0);
                    CombinedJet2_SDM = GoodJetsAK8SDM.at(0);
                    CombinedJet2_Pruned = GoodJetsAK8Pruned.at(0);
                    CombinedJet2_MaxSubjet = GoodJetsAK8MaxSubjet.at(0);
                    CombinedJet2_MinSubjet = GoodJetsAK8MinSubjet.at(0);

                    
                    NlinoJet2 = poss2_NlinoJet2;                                        
                }
//                std::cout << "--------------------------" << std::endl;
//                std::cout << CombinedJet1_SDM << " " << CombinedJet2_SDM << std::endl;
            }
            
        }

        asymm_mt2_lester_bisect::disableCopyrightMessage();
        
        tr.registerDerivedVar("FatJetCombined1"+myVarSuffix_, CombinedJet1);
        tr.registerDerivedVar("FatJetCombined2"+myVarSuffix_, CombinedJet2);
        tr.registerDerivedVar("FatJetMT2"+myVarSuffix_, ttUtility::coreMT2calc(CombinedJet1,CombinedJet2, lvMET));
        tr.registerDerivedVar("FatJetCombined1_SDM"+myVarSuffix_, CombinedJet1_SDM);
        tr.registerDerivedVar("FatJetCombined2_SDM"+myVarSuffix_, CombinedJet2_SDM);
        tr.registerDerivedVar("FatJetCombined1_Pruned"+myVarSuffix_, CombinedJet1_Pruned);
        tr.registerDerivedVar("FatJetCombined2_Pruned"+myVarSuffix_, CombinedJet2_Pruned);
        tr.registerDerivedVar("FatJetCombined1_T3T1"+myVarSuffix_, CombinedJet1_T3T1);
        tr.registerDerivedVar("FatJetCombined2_T3T1"+myVarSuffix_, CombinedJet2_T3T1);
        tr.registerDerivedVar("FatJetCombined1_T2T1"+myVarSuffix_, CombinedJet1_T2T1);
        tr.registerDerivedVar("FatJetCombined2_T2T1"+myVarSuffix_, CombinedJet2_T2T1);
        tr.registerDerivedVar("FatJetCombined1_T3T2"+myVarSuffix_, CombinedJet1_T3T2);
        tr.registerDerivedVar("FatJetCombined2_T3T2"+myVarSuffix_, CombinedJet2_T3T2);
        tr.registerDerivedVar("FatJetCombined1_MaxSubjet"+myVarSuffix_, CombinedJet1_MaxSubjet);
        tr.registerDerivedVar("FatJetCombined2_MaxSubjet"+myVarSuffix_, CombinedJet2_MaxSubjet);
        tr.registerDerivedVar("FatJetCombined1_MinSubjet"+myVarSuffix_, CombinedJet1_MinSubjet);
        tr.registerDerivedVar("FatJetCombined2_MinSubjet"+myVarSuffix_, CombinedJet2_MinSubjet);

        tr.registerDerivedVar("FatJet1_Idx"+myVarSuffix_, FatJet1_Idx);
        tr.registerDerivedVar("FatJet2_Idx"+myVarSuffix_, FatJet2_Idx);

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
