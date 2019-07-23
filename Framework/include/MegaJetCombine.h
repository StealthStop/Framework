#ifndef MEGAJETCOMBINE_H
#define MEGAJEETCOMBINE_H

#include <iostream>
#include <algorithm>



class MegaJetCombine
{
private:
    typedef std::pair<std::vector<int>, std::vector<int>> pairListInt;
    std::string myVarSuffix_;
    std::map<int, const std::vector<pairListInt>&&> comboMap;
    
    const std::vector<std::pair< std::vector<int>, std::vector<int>>> getCombineList(const int N = 3, const int MinNJetsTotal = 0, const int MinNJetsPerCombo = 0)
    {
        auto start = std::chrono::high_resolution_clock::now();

        std::vector<std::pair< std::vector<int>,std::vector<int>>>  MegaJets;
        std::vector<int> JetIndex(N, -1);
        for (int i = 0; i<N; i++) JetIndex[i] = i+1;
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
        const auto& NGoodJets_pt30 = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_);
        if( comboMap.find(NGoodJets_pt30) == comboMap.end() )
        {
            std::cout<<"Adding this combo vector in map for size: "<<NGoodJets_pt30<<std::endl;
            comboMap.insert( std::move(std::pair<int, const std::vector<pairListInt>>( NGoodJets_pt30, getCombineList(NGoodJets_pt30, 0, 0))) );
        }
        std::cout<<comboMap[NGoodJets_pt30].size()<<std::endl;
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
