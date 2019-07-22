#ifndef MEGAJETCOMBINE_H
#define MEGAJEETCOMBINE_H

#include <iostream>
#include <algorithm>

class MegaJetCombine
{
private:
    std::string myVarSuffix_;

    void Combine(NTupleReader& tr) const
    {
        auto start = std::chrono::high_resolution_clock::now();
        const auto& NGoodJets_pt30 = tr.getVar<int>("NGoodJets_pt30");
        
        int N = 13;
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
                        if ( (k&(1<<d)) ) { 
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
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        std::cout << duration.count() << std::endl;
        std::cout <<"NJets: " << N << " Number of combinations: " <<  MegaJets.size() << std::endl;
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
    }
};

#endif
