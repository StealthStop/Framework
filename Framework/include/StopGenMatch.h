#ifndef STOPGENMATCH_H
#define STOPGENMATCH_H
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/interface/lester_mt2_bisect.h"

class StopGenMatch
{
private:
    std::string myVarSuffix_;

    inline int findParent(const int p, const int idx, const std::vector<int>& GenParticles_ParentId, const std::vector<int>& GenParticles_ParentIdx) const
    {
        if (idx == -1)
        {
            return -1;
        }
        else if(abs(GenParticles_ParentId[idx]) == p)
        {
            return GenParticles_ParentId[idx];
        }
        else
        {
            return findParent(p, GenParticles_ParentIdx[idx], GenParticles_ParentId, GenParticles_ParentIdx);
        }
    }

    //function to match all reco particles with all appropriate gen particles
    inline std::vector<std::tuple< int , int , double>> findAllDR(const std::vector<TLorentzVector>& GenParticles, const std::vector<TLorentzVector>& RecoParticles, 
                                                                  const std::vector<bool>& GoodGenParticles, int resPartID, const std::vector<int>& GenParticles_ParentId, 
                                                                  const std::vector<int>& GenParticles_ParentIdx, int nlino1_Idx, int nlino2_Idx) const
    {
        bool check_neutralinos = true;
        std::vector<std::tuple< int , int, double>> AllDR;
        std::tuple< int , int, double> DRtup;
        int check_resPartID = resPartID;
        for (unsigned int g=0; g < GenParticles.size(); g++)
        {
            if (resPartID == 1000022)
            {
                check_neutralinos = GenParticles_ParentIdx.at(g) == nlino1_Idx;
            }
            if (resPartID == -1000022)
            {
                check_resPartID = 1000022;
                check_neutralinos = GenParticles_ParentIdx.at(g) == nlino2_Idx;
            }
            for (unsigned int r=0; r < RecoParticles.size(); r++)
            {
                if (findParent(abs(check_resPartID), g, GenParticles_ParentId, GenParticles_ParentIdx) == check_resPartID && GoodGenParticles.at(g) && check_neutralinos)
                {
                    std::get<0>(DRtup) = g;
                    std::get<1>(DRtup) = r;
                    std::get<2>(DRtup) = GenParticles.at(g).DeltaR(RecoParticles.at(r));
                    AllDR.push_back(DRtup);
                }
            }
        }
        return AllDR;
    }       

    //function to sort for best matches and construct stop mass
    inline std::pair<std::vector<TLorentzVector>,int> getMatchedSum(const std::vector<std::tuple< int , int , double>>& AllDR, const std::vector<TLorentzVector>& RecoParticles, TLorentzVector& MatchedSum, std::vector<bool> availableDR, TLorentzVector& GenMatchedSum,const std::vector<TLorentzVector>& GenParticles, int NGenMatched = 0) const
    {
        double minDR = 999;
        std::tuple< int, int, double> bestDR;
        for (unsigned int d=0; d < AllDR.size(); d++)
        {
            if ( std::get<2>(AllDR.at(d)) < minDR && availableDR.at(d))
            {
                bestDR = AllDR.at(d);
                minDR = std::get<2>(AllDR.at(d));
            }
        }
        bool allgone = true;
        for (const auto& u : availableDR)
        {
            if (u) allgone = false;
        }
        if (!allgone)
        {
            if (std::get<2>(bestDR) < 1)
            {
                MatchedSum +=  RecoParticles.at(std::get<1>(bestDR));
                GenMatchedSum += GenParticles.at(std::get<0>(bestDR));
                NGenMatched += 1;
            }
            for (unsigned int d=0;  d < AllDR.size(); d++)
            {
                if (std::get<0>(AllDR.at(d)) == std::get<0>(bestDR) || std::get<1>(AllDR.at(d)) == std::get<1>(bestDR))
                {
                    availableDR.at(d) = false;
                }
            }
            return getMatchedSum( AllDR, RecoParticles, MatchedSum, availableDR, GenMatchedSum, GenParticles, NGenMatched);
        }
        else
        {
            std::pair<std::vector<TLorentzVector>,int> Sums;
            Sums.first.push_back(MatchedSum);
            Sums.first.push_back(GenMatchedSum);
            Sums.second = NGenMatched;
            return Sums;
        }
    }

    void genMatch(NTupleReader& tr)
    {
        const auto& runtype = tr.getVar<std::string>("runtype");

        if(runtype != "Data")
        {
            const auto& GenParticles            = tr.getVec<TLorentzVector>("GenParticles"+myVarSuffix_);
            const auto& GenParticles_PdgId      = tr.getVec<int>("GenParticles_PdgId"+myVarSuffix_);
            const auto& GenParticles_ParentId   = tr.getVec<int>("GenParticles_ParentId"+myVarSuffix_);
            const auto& GenParticles_ParentIdx  = tr.getVec<int>("GenParticles_ParentIdx"+myVarSuffix_);
            const auto& GenParticles_Status     = tr.getVec<int>("GenParticles_Status"+myVarSuffix_);            
            const auto& GoodJets_pt30           = tr.getVec<bool>("GoodJets_pt30"+myVarSuffix_);
            const auto& GoodLeptons             = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons"+myVarSuffix_);
            const auto& Jets                    = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
            const auto& MET                     = tr.getVar<double>("MET"+myVarSuffix_);
            const auto& METPhi                  = tr.getVar<double>("METPhi"+myVarSuffix_);
            const auto& GenMET                  = tr.getVar<double>("GenMET"+myVarSuffix_);
            const auto& GenMETPhi               = tr.getVar<double>("GenMETPhi"+myVarSuffix_);

            TLorentzVector lvMET;
            lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);
            
            TLorentzVector lvGenMET;
            lvGenMET.SetPtEtaPhiM(GenMET, 0.0, GenMETPhi, 0.0);

            std::vector<TLorentzVector> RecoParticles;
            for(unsigned int j=0; j < Jets.size(); j++)
            {
                if(GoodJets_pt30.at(j)) RecoParticles.push_back(Jets.at(j));
            }
            for(const auto& l : GoodLeptons)
            {
                RecoParticles.push_back(l.second);
            }
            std::vector<int> neutralinos_Idx;
            int nlino1_Idx = -1, nlino2_Idx = -1;
            //int nUsedGen = 0;
            //int totalGens = 0;

            //Define OkayParticles, which allows leptons/jets by status code and parent
            std::vector<bool> OkayGenParticles(GenParticles.size(), false);
            //int WJetCounter = 0, WPlusLepCounter = 0, WMinusLepCounter = 0;
            int TauLepCounter = 0;
            std::vector<int> WPlusLeps, WMinusLeps;
            for (unsigned int p=0; p < GenParticles.size(); p++)
            {
                int pdgid = GenParticles_PdgId.at(p);
                int momid = GenParticles_ParentId.at(p);
                int status = GenParticles_Status.at(p);
                bool is_lepton = ( abs(pdgid) == 11 || abs(pdgid) == 13 || abs(pdgid) == 15);
                bool is_jet = ( abs(pdgid) <= 5 || abs(pdgid) == 21);
                int WId = findParent(24, p, GenParticles_ParentId, GenParticles_ParentIdx);
                bool pass_lepton = is_lepton ? (status == 1) && (abs(momid) == 24 || abs(momid) == 15): false; //leptons must be status 1 and come from either a W or a tau
                bool pass_jet = is_jet ? status == 23 : false; //jets must have status 23
                int stopId = findParent(1000006, p, GenParticles_ParentId, GenParticles_ParentIdx);
                bool pass_stop = stopId != -1; //all gen particles must come from a stop
                bool filter = (pass_lepton || pass_jet) && pass_stop;

                if (filter)
                {
                    if (pass_lepton && WId == 24 ) WPlusLeps.push_back(p);
                    if (pass_lepton && WId == -24) WMinusLeps.push_back(p);
                    if ((abs(pdgid) == 11 || abs(pdgid) ==13) && pass_stop &&  abs(momid) == 15) TauLepCounter += 1;
                    OkayGenParticles.at(p) = true;
                }
            }
            //Define GoodGenParticles, which has no W radiation (for most part) and allows undecayed taus
            std::vector<bool> GoodGenParticles = OkayGenParticles;
            bool w_eplus = false, w_eminus = false, w_muplus = false, w_muminus =false;
            int em = 0, ep = 0, mm = 0, mp = 0;
            
            for ( const auto& w : WPlusLeps) 
            {
                if (GenParticles_PdgId.at(w) == 13)
                {
                    w_eminus = true;
                    em = w;
                }
                if (GenParticles_PdgId.at(w) == -13)
                {
                    w_eplus = true;
                    ep = w;
                }
                if (GenParticles_PdgId.at(w) == 11)
                {
                    w_muminus = true;
                    mm = w;
                }
                if (GenParticles_PdgId.at(w) == -11)
                {
                    w_muplus = true;
                    mp = w;
                }
                    
            } //removes pair produced leptons from radation. does not remove multiple pairs, need to fix for future studies
            if (w_eplus && w_eminus) 
            {
                GoodGenParticles.at(ep) = false;
                GoodGenParticles.at(em) = false;
            }
            if (w_muplus && w_muminus)
            {
                GoodGenParticles.at(mp) = false;
                GoodGenParticles.at(mm) = false;
            }
                
            bool keepTau = TauLepCounter == 0;
            for (unsigned int g=0; g < OkayGenParticles.size(); g++)
            {
                if (abs(GenParticles_PdgId.at(g)) == 15 && GenParticles_Status.at(g) == 2 && keepTau && findParent(1000006, g, GenParticles_ParentId, GenParticles_ParentIdx) != -1) GoodGenParticles.at(g) = true;
            }
            //int counter = 0;
            for (unsigned int g=0; g < GenParticles.size();g++)
            {
                if (abs(GenParticles_PdgId.at(g)) == 1000022 && findParent(1000006, g, GenParticles_ParentId, GenParticles_ParentIdx) != -1)
                {
                    neutralinos_Idx.push_back(g);
                }

            }
            
            if (neutralinos_Idx.size() == 2) //need to sort neutralinos by Pt
            {
                if (GenParticles.at(neutralinos_Idx.at(0)).Pt() > GenParticles.at(neutralinos_Idx.at(1)).Pt())
                {
                    nlino1_Idx = neutralinos_Idx.at(0);
                    nlino2_Idx = neutralinos_Idx.at(1);
                }
                else 
                {
                    nlino1_Idx = neutralinos_Idx.at(1);
                    nlino2_Idx = neutralinos_Idx.at(0);
                }
            }
            
            std::vector<int>  resParticleList{1000006, -1000006, 1000022, -1000022,  5000002, -5000002};
            std::vector<TLorentzVector> RecoSumList;
            std::vector<TLorentzVector> GenSumList;
            std::vector<int> NMatched(2,0);
            float NGenTotal = 0;
            double fracGenMatched = -1.0;
            for(unsigned int g=0; g < GenParticles.size(); g++)
            {
                if (GoodGenParticles[g]) NGenTotal += 1;
            }
            for(unsigned int p=0; p < resParticleList.size(); p++)
            {
                TLorentzVector initMatchedSum;
                TLorentzVector GenMatchedSum;
                
                std::vector<std::tuple< int , int , double>> AllDR = findAllDR(GenParticles, RecoParticles, GoodGenParticles, resParticleList[p], GenParticles_ParentId, GenParticles_ParentIdx, nlino1_Idx, nlino2_Idx);
                std::vector<bool> availableDR(AllDR.size(), true);
                std::pair<std::vector<TLorentzVector>,int> Sums = getMatchedSum( AllDR, RecoParticles, initMatchedSum, availableDR, GenMatchedSum, GenParticles);
                RecoSumList.push_back(Sums.first.at(0));
                GenSumList.push_back(Sums.first.at(1));
                
                if (resParticleList[p] == 1000006) NMatched[0] = Sums.second;
                else if (resParticleList[p] == -1000006) NMatched[1] = Sums.second;
            }
            int NTotMatched = NMatched[0] + NMatched[1];
            if (NGenTotal != 0) fracGenMatched = std::round((NTotMatched / NGenTotal)*1000)/1000;
//            std::cout << NTotMatched << " " << NGenTotal << " " << fracGenMatched << std::endl;
       

            
            asymm_mt2_lester_bisect::disableCopyrightMessage();
            
            tr.registerDerivedVar("GM_StopMT2"+myVarSuffix_,        ttUtility::coreMT2calc(RecoSumList.at(0),RecoSumList.at(1),lvMET));
            tr.registerDerivedVar("GM_StopGenMT2"+myVarSuffix_,     ttUtility::coreMT2calc(GenSumList.at(0),GenSumList.at(1),lvGenMET));
            tr.registerDerivedVar("GM_Stop1Mass"+myVarSuffix_,      RecoSumList.at(0).M());
            tr.registerDerivedVar("GM_Stop2Mass"+myVarSuffix_,      RecoSumList.at(1).M());
            tr.registerDerivedVar("GM_Stop1GenMass"+myVarSuffix_,   GenSumList.at(0).M());
            tr.registerDerivedVar("GM_Stop2GenMass"+myVarSuffix_,   GenSumList.at(1).M());            
            tr.registerDerivedVar("GM_Nlino1Mass"+myVarSuffix_,     RecoSumList.at(2).M());
            tr.registerDerivedVar("GM_Nlino2Mass"+myVarSuffix_,     RecoSumList.at(3).M());
            tr.registerDerivedVar("GM_Nlino1GenMass"+myVarSuffix_,  GenSumList.at(2).M());
            tr.registerDerivedVar("GM_Nlino2GenMass"+myVarSuffix_,  GenSumList.at(3).M());
            tr.registerDerivedVar("GM_Single1Mass"+myVarSuffix_,    RecoSumList.at(4).M());
            tr.registerDerivedVar("GM_Single2Mass"+myVarSuffix_,    RecoSumList.at(5).M());
            tr.registerDerivedVar("GM_Single1GenMass"+myVarSuffix_, GenSumList.at(4).M());
            tr.registerDerivedVar("GM_Single2GenMass"+myVarSuffix_, GenSumList.at(5).M());

            tr.registerDerivedVar("fracGenMatched"+myVarSuffix_,  fracGenMatched);
        }
    }

public:
    StopGenMatch(std::string myVarSuffix = "") 
        : myVarSuffix_       (myVarSuffix)
    {                
        std::cout<<"Setting up StopGenMatch"<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        genMatch(tr);
    }
};

#endif
