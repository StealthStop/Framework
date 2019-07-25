#ifndef MEGAJETCOMBINE_H
#define MEGAJEETCOMBINE_H

#include <iostream>
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/interface/lester_mt2_bisect.h"


class MegaJetCombine
{
private:
    typedef std::pair<std::vector<int>, std::vector<int>> pairListInt;
    std::string myVarSuffix_;
    std::map<int, const std::vector<pairListInt>> comboMap;
    
    const std::vector<std::pair< std::vector<int>, std::vector<int>>> getCombineList(const int N = 3, const int MinNJetsTotal = 2, const int MinNJetsPerCombo = 1)
    {
        auto start = std::chrono::high_resolution_clock::now();

        std::vector<pairListInt>  MegaJets;
        std::vector<int> JetIndex(N, -1);
        for (int i = 0; i<N; i++) JetIndex[i] = i;
        for (int i=0; i < (1<<N); i++)
        {
            int len = 0;
            std::vector<int> pairedCombo(N, -1);
            for (int j=0; j<N; j++) if ( (i&(1<<j)) ) pairedCombo[len++] = JetIndex[j];
            if (len >= MinNJetsTotal)
            {
                for (int k=1; k<(1<<(len-1)); k++)
                {
                    std::vector<int> firstOfPair, secOfPair;
                    // std::cout << "{ {";
                    bool  sepHere = false;
                    for (int d = 0; d < len; d++)
                    {
                        if ( (k&(1<<d)) ) 
                        { 
                            sepHere = true;
                            firstOfPair.emplace_back(pairedCombo[d]);
                            // std::cout << pairedCombo[d] << " ";
                        }
                    }
                    // std::cout << ",";
                    sepHere = false;
                    for (int d = 0; d< len; d++)
                    {  
                        if ( (k&(1<<d)) == 0)
                        {
                            sepHere = true;
                            secOfPair.emplace_back(pairedCombo[d]);
                            // std::cout << pairedCombo[d] <<" ";
                        }                
                    }
                
                    // std::cout << "} }" << std::endl;
                    if (firstOfPair.size() >= MinNJetsPerCombo && secOfPair.size() >= MinNJetsPerCombo)
                    {
                        MegaJets.emplace_back( std::move(firstOfPair), std::move(secOfPair) );
                        //MegaJets.push_back(std::move({std::move(firstOfPair), std::move(secOfPair)}));
                    }
                }
            }                       
        }
        return MegaJets;
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
//        std::cout << duration.count() << std::endl;
        //       std::cout <<"NJets: " << N << " Number of combinations: " <<  MegaJets.size() << std::endl;
    }
    
    void megaJetCombine(NTupleReader& tr)
    {
        const auto& NGoodJets_pt30       = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_);
        const auto& Jets                 = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
        const auto& GoodJets_pt30        = tr.getVec<bool>("GoodJets_pt30"+myVarSuffix_);
        const auto& GoodLeptons          = tr.getVec<std::pair<std::string,TLorentzVector>>("GoodLeptons"+myVarSuffix_);
        const auto& NGoodLeptons         = tr.getVar<int>("NGoodLeptons"+myVarSuffix_);
        const auto& GoodBJets_pt30       = tr.getVec<bool>("GoodBJets_pt30"+myVarSuffix_);
        const auto& NGoodBJets_pt30      = tr.getVar<int>("NGoodBJets_pt30"+myVarSuffix_);
        
        const auto& TwoLep_Mbl1_Idx      = tr.getVar<std::pair<int,int>>("TwoLep_Mbl1_Idx"+myVarSuffix_);
        const auto& TwoLep_Mbl2_Idx      = tr.getVar<std::pair<int,int>>("TwoLep_Mbl2_Idx"+myVarSuffix_);

        const auto& MET                  = tr.getVar<double>("MET"+myVarSuffix_);
        const auto& METPhi               = tr.getVar<double>("METPhi"+myVarSuffix_);
        TLorentzVector lvMET;
        lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);
      
        TLorentzVector RecoStop1, RecoStop2;
        if (NGoodLeptons == 2 && NGoodJets_pt30 >= 2)
        {
            TLorentzVector BLVec_1 = Jets[TwoLep_Mbl1_Idx.first] + GoodLeptons[TwoLep_Mbl1_Idx.second].second;
            TLorentzVector BLVec_2 = Jets[TwoLep_Mbl2_Idx.first] + GoodLeptons[TwoLep_Mbl2_Idx.second].second;
            double StopDiff = 9999;
            
            std::vector<TLorentzVector> JetsToCombine;
            
            for (int j =0; j < Jets.size(); j++)
            {
                if (GoodJets_pt30[j] && j != TwoLep_Mbl1_Idx.first && j != TwoLep_Mbl2_Idx.first) JetsToCombine.emplace_back(Jets[j]);
            }
            if( comboMap.find(JetsToCombine.size()) == comboMap.end() )
            {
//                std::cout<<"Adding this combo vector in map for size: "<<NGoodJets_pt30<<std::endl;
                comboMap.insert( std::move(std::pair<int, const std::vector<pairListInt>>( JetsToCombine.size(), std::move(getCombineList(JetsToCombine.size(), 6, 3)))) );
            }
            //          std::cout<<comboMap[NGoodJets_pt30].size()<<std::endl;
            for (int c=0; c < comboMap[JetsToCombine.size()].size(); c++)
            {
                TLorentzVector FirstOfPairSum, SecOfPairSum;
//                std::cout <<"First pair: ";
                for (int j=0; j < comboMap[JetsToCombine.size()][c].first.size(); j++)
                {
                    FirstOfPairSum +=  JetsToCombine[comboMap[JetsToCombine.size()][c].first[j]];
                    //  std::cout << comboMap[JetsToCombine.size()][c].first[j] << " ";
                }
                // std::cout  << std::endl;
                //std::cout << "Second pair: ";
                for (int j=0; j < comboMap[JetsToCombine.size()][c].second.size(); j++)
                {
                    SecOfPairSum += JetsToCombine[comboMap[JetsToCombine.size()][c].second[j]];
                    //     std::cout << comboMap[JetsToCombine.size()][c].second[j] << " ";
                }
                //std::cout << std::endl;
                if ( abs( (FirstOfPairSum + BLVec_1).M() - (SecOfPairSum + BLVec_2).M()) < abs( (FirstOfPairSum + BLVec_2).M() - (SecOfPairSum + BLVec_1).M() ))
                {
                    FirstOfPairSum += BLVec_1;
                    SecOfPairSum += BLVec_2;
                }
                else
                {
                    FirstOfPairSum += BLVec_2;
                    SecOfPairSum += BLVec_1;
                }
                
                double massDiff = abs(FirstOfPairSum.M() - SecOfPairSum.M());
                if ( massDiff < StopDiff) 
                {
                    RecoStop1 = FirstOfPairSum;
                    RecoStop2 = SecOfPairSum;
                    StopDiff = massDiff;
                }
            }
//            std::cout << RecoStop1.M() << " "<< RecoStop2.M() << "NJets: "<< JetsToCombine.size() << std::endl;
            //          std::cout << "-----------------------------" << std::endl;           
        }

        asymm_mt2_lester_bisect::disableCopyrightMessage();

        tr.registerDerivedVar("RecoStop1"+myVarSuffix_, RecoStop1);
        tr.registerDerivedVar("RecoStop2"+myVarSuffix_, RecoStop2);
        tr.registerDerivedVar("RecoStopMT2"+myVarSuffix_, ttUtility::coreMT2calc(RecoStop1,RecoStop2,lvMET));
    }    

public:
    MegaJetCombine(std::string myVarSuffix = "")
        : myVarSuffix_      (myVarSuffix)
    {
        std::cout<<"Setting up MegaJetCombine"<<std::endl;        
    }
    void operator()(NTupleReader& tr)
    {
        megaJetCombine(tr);
    }
};

#endif
