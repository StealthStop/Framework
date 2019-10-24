#ifndef MEGAJETCOMBINE_H
#define MEGAJEETCOMBINE_H

#include <iostream>
#include <algorithm>
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/interface/lester_mt2_bisect.h"


class MegaJetCombine
{
private:
    typedef std::pair<std::vector<int>, std::vector<int>> pairListInt;
    std::string myVarSuffix_;
    std::map<int, const std::vector<pairListInt>> comboMap;
    
    const std::vector<std::pair< std::vector<int>, std::vector<int>>> getCombineList(const int N = 3, const int MinNJetsTotal = 2, const unsigned int MinNJetsPerCombo = 1) //outputs possible combos of any number up to N jets into two megajets 
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
                    // std::cout << "{ {"; //these comments out std::cout options will print out the combinations in a more readable way
                    //bool  sepHere = false;
                    for (int d = 0; d < len; d++)
                    {
                        if ( (k&(1<<d)) ) 
                        { 
                            //sepHere = true;
                            firstOfPair.emplace_back(pairedCombo[d]);
                            // std::cout << pairedCombo[d] << " ";
                        }
                    }
                    // std::cout << ",";
                    //sepHere = false;
                    for (int d = 0; d< len; d++)
                    {  
                        if ( (k&(1<<d)) == 0)
                        {
                            //sepHere = true;
                            secOfPair.emplace_back(pairedCombo[d]);
                            // std::cout << pairedCombo[d] <<" ";
                        }                
                    }
                
                    // std::cout << "} }" << std::endl;
                    if (firstOfPair.size() >= MinNJetsPerCombo && secOfPair.size() >= MinNJetsPerCombo)
                    {
                        MegaJets.emplace_back( std::move(firstOfPair), std::move(secOfPair) );
                    }
                }
            }                       
        }
        return MegaJets;
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); //can be used to measure time it takes to compute the list
//        std::cout << duration.count() << std::endl;
        //       std::cout <<"NJets: " << N << " Number of combinations: " <<  MegaJets.size() << std::endl;
    }
    
    void megaJetCombine(NTupleReader& tr) //Main function for using the list of combos to compute megajets
    {
        //const auto& NGoodJets_pt30       = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_);
        const auto& Jets                 = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
        const auto& GoodJets_pt30        = tr.getVec<bool>("GoodJets_pt30"+myVarSuffix_);
        const auto& GoodLeptons          = tr.getVec<std::pair<std::string,TLorentzVector>>("GoodLeptons"+myVarSuffix_);
        const auto& NGoodLeptons         = tr.getVar<int>("NGoodLeptons"+myVarSuffix_);
        //const auto& GoodBJets_pt30       = tr.getVec<bool>("GoodBJets_pt30"+myVarSuffix_);
        //const auto& NGoodBJets_pt30      = tr.getVar<int>("NGoodBJets_pt30"+myVarSuffix_);
        const auto& TwoLep_Mbl1_Idx      = tr.getVar<std::pair<unsigned int, unsigned int>>("TwoLep_Mbl1_Idx"+myVarSuffix_);
        const auto& TwoLep_Mbl2_Idx      = tr.getVar<std::pair<unsigned int, unsigned int>>("TwoLep_Mbl2_Idx"+myVarSuffix_);
        const auto& MET                  = tr.getVar<double>("MET"+myVarSuffix_);
        const auto& METPhi               = tr.getVar<double>("METPhi"+myVarSuffix_);

        double massDiffThresh = 100; //specify threshold of mass difference between megajets to generate shortlist of candidates
        TLorentzVector lvMET;
        lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);
      
        TLorentzVector RecoStop1, RecoStop2;
        std::pair<std::vector<TLorentzVector>,std::vector<TLorentzVector>> RecoStopCands;
        
        auto& FirstComboCandidates         = tr.createDerivedVec<TLorentzVector>("FirstComboCandidates"+myVarSuffix_);
        auto& SecComboCandidates           = tr.createDerivedVec<TLorentzVector>("SecComboCandidates"+myVarSuffix_);
        auto& FirstStopMassSums            = tr.createDerivedVec<double>("FirstStopMassSums"+myVarSuffix_);
        auto& SecStopMassSums              = tr.createDerivedVec<double>("SecStopMassSums"+myVarSuffix_);
        auto& StopMassDiffs                = tr.createDerivedVec<double>("StopMassDiffs"+myVarSuffix_);
        auto& StopMT2s                     = tr.createDerivedVec<double>("StopMT2s"+myVarSuffix_);


        if (NGoodLeptons == 2)
        {
            TLorentzVector BLVec_1 = Jets[TwoLep_Mbl1_Idx.first] + GoodLeptons[TwoLep_Mbl1_Idx.second].second;
            TLorentzVector BLVec_2 = Jets[TwoLep_Mbl2_Idx.first] + GoodLeptons[TwoLep_Mbl2_Idx.second].second;
            //double StopDiff = 9999;
            
            std::vector<TLorentzVector> JetsToCombine;
            
            for (unsigned int j =0; j < Jets.size(); j++)
            {
                if (GoodJets_pt30[j] && j != TwoLep_Mbl1_Idx.first && j != TwoLep_Mbl2_Idx.first) JetsToCombine.emplace_back(Jets[j]);
            }
                                    
            if( comboMap.find(JetsToCombine.size()) == comboMap.end() ) //fills a map with the combo possibliites and checks to see if the combos corresponding to NJets already exist before computing
            {
                comboMap.insert( std::move(std::pair<int, const std::vector<pairListInt>>( JetsToCombine.size(), std::move(getCombineList(JetsToCombine.size(), 4, 2)))) );
            }
                                   
            for (unsigned int c=0; c < comboMap[JetsToCombine.size()].size(); c++)
            {
                TLorentzVector FirstOfPairSum = BLVec_1, SecOfPairSum = BLVec_2;
                for (unsigned int j=0; j < comboMap[JetsToCombine.size()][c].first.size(); j++)
                {
                    FirstOfPairSum +=  JetsToCombine[comboMap[JetsToCombine.size()][c].first[j]];
                }
                for (unsigned int j=0; j < comboMap[JetsToCombine.size()][c].second.size(); j++)
                {
                    SecOfPairSum += JetsToCombine[comboMap[JetsToCombine.size()][c].second[j]];
                }
                  
                if ( abs(FirstOfPairSum.M() - SecOfPairSum.M()) < massDiffThresh) 
                {
                    RecoStopCands.first.emplace_back(FirstOfPairSum);
                    RecoStopCands.second.emplace_back(SecOfPairSum);
                }
            }

            std::vector<double> TotalStopCandMassSums; //some extra variables are defined for neural network input
        
           for (unsigned int p=0; p < RecoStopCands.first.size(); p++)
           {
               FirstStopMassSums.emplace_back(RecoStopCands.first[p].M());
               SecStopMassSums.emplace_back(RecoStopCands.second[p].M());
               StopMT2s.emplace_back(ttUtility::coreMT2calc(RecoStopCands.first[p], RecoStopCands.second[p], lvMET));
               TotalStopCandMassSums.emplace_back((RecoStopCands.first[p] + RecoStopCands.second[p]).M());
           }
           std::vector<double>::iterator Max_mass = std::max_element(TotalStopCandMassSums.begin(), TotalStopCandMassSums.end());
           if (TotalStopCandMassSums.size() > 0)
           {
               RecoStop1 = RecoStopCands.first[std::distance(TotalStopCandMassSums.begin(), Max_mass)]; //selects megajet combos with max mass from shortlist as best option
               RecoStop2 = RecoStopCands.second[std::distance(TotalStopCandMassSums.begin(), Max_mass)];
           }
           
           for (unsigned int s=0; s < FirstStopMassSums.size(); s++)
           {
               StopMassDiffs.emplace_back(abs(FirstStopMassSums[s] - SecStopMassSums[s]));
            
           }
        } //End 2L Restriction

        FirstComboCandidates = RecoStopCands.first;
        SecComboCandidates = RecoStopCands.second;
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
