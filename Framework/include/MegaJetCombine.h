#ifndef MEGAJETCOMBINE_H
#define MEGAJEETCOMBINE_H

#include <iostream>
#include <algorithm>

class MegaJetCombine
{
private:
    std::string myVarSuffix_;

    inline std::vector<std::pair< std::vector<int>,std::vector<int>> > Combine(NTupleReader& tr) const
    {
        auto start = std::chrono::high_resolution_clock::now();
        const auto& NGoodJets_pt30 = tr.getVar<int>("NGoodJets_pt30");
        
        int N = 3;
        N = NGoodJets_pt30;
        int JetIndex[N];
        int MinNJetsTotal = 0;
        int MinNJetsPerCombo = 0;
        auto& MegaJets  = tr.createDerivedVec<std::pair< std::vector<int>,std::vector<int>> >("MegaJets");
        for (int i = 0; i<N; i++) JetIndex[i] = i+1;
        for (int i=0; i < (1<<N); i++)
        {
            int len = 0;
            int pairedCombo[N];
            for (int j=0; j<N; j++) if ( (i&(1<<j)) ) pairedCombo[len++] = JetIndex[j];
            if (len >= MinNJetsTotal)
            {
                std::pair< std::vector<int>, std::vector<int> > dummyPair;
               
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
                            firstOfPair.push_back(pairedCombo[d]);
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
                            secOfPair.push_back(pairedCombo[d]);
                            // std::cout << pairedCombo[d] <<" ";
                        }
                
                    }
                
                    // std::cout << "} }" << std::endl;
                    if (firstOfPair.size() >= MinNJetsPerCombo && secOfPair.size() >= MinNJetsPerCombo)
                    {
                        dummyPair.first = firstOfPair;
                        dummyPair.second = secOfPair;
                        MegaJets.push_back(dummyPair);
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

    void filterJets_2l(NTupleReader& tr)
    {
        const auto& Jets                         = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt30                = tr.getVec<bool>("GoodJets_pt30");
        const auto& GoodBJets_pt30               = tr.getVec<bool>("GoodBJets_pt30");
        const auto& GoodLeptons                  = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        const auto& NGoodJets_pt30               = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30              = tr.getVar<int>("NGoodBJets_pt30");
        const auto& NGoodLeptons                 = tr.getVar<int>("NGoodLeptons");

        const auto& TwoLep_Mbl1                  = tr.getVar<double>("TwoLep_Mbl1");
        const auto& TwoLep_Mbl2                  = tr.getVar<double>("TwoLep_Mbl2");
        const auto& TwoLep_Mbl1_Idx              = tr.getVar<std::pair<int,int>>("TwoLep_Mbl1_Idx");
        const auto& TwoLep_Mbl2_Idx              = tr.getVar<std::pair<int,int>>("TwoLep_Mbl2_Idx");




        /* std::vector<bool> filteredJets_2l = GoodJets_pt30;
        std::vector<TLorentzVector> GoodBJetsVec_pt30;
        TLorentzVector lep1, lep2;
        if (NGoodLeptons==2)
        {
            if (GoodLeptons.at(0).second.Pt() >= GoodLeptons.at(1).second.Pt())
            {
                lep1 = GoodLeptons.at(0).second;
                lep2 = GoodLeptons.at(1).second;
            }
            else
            {
                lep1 = GoodLeptons.at(1).second;
                lep2 = GoodLeptons.at(0).second;
            }
            if (NGoodBJets_pt30 >= 2)
            {
                for (int b=0; b < Jets.size(); b++)
                {
                    if (GoodBJets_pt30.at(b))
                    {
                        double mbl1 = (lep1 + Jets.at(b)).M(), mbl2 = (lep2 + Jets.at(b)).M();
                        if ( (mbl1 <25 || mbl1 > 250) && (mbl2 < 25 || mbl2 >250)) filteredJets_2l.at(b) = false;
                    }


                }
            }
        }
        int  NfilteredJets_2l=0;
        for (int g = 0; g < filteredJets_2l.size(); g++) if (filteredJets_2l.at(g)) NfilteredJets_2l += 1;
       
        if(NfilteredJets_2l != NGoodJets_pt30) 
        {
            std::cout<< "Something was actually filtered!";
            std::cout << "Number of good jets: " << NGoodJets_pt30 << " Number of good bjets: " << NGoodBJets_pt30 << " Number of filtered jets: " << NfilteredJets_2l << std::endl;
        }
        */
    }
    public:
    MegaJetCombine(std::string myVarSuffix = "")
        : myVarSuffix_      (myVarSuffix)
    {
        std::cout<<"Setting up MegaJetCombine"<<std::endl;
    }
    void operator()(NTupleReader& tr)
    {
        Combine(tr);
        filterJets_2l(tr);
    }
};

#endif
