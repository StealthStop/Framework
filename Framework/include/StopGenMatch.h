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

    //function to generate all possible matches between gen and reco particles if they pass DR and pT cut
    inline std::vector<std::tuple< int , int , double>> findAllDR(const std::vector<TLorentzVector>& GenParticles, const std::vector<TLorentzVector>& RecoParticles, 
                                                                  const std::vector<bool>& GoodGenParticles, const int resPartID, const std::vector<int>& GenParticles_ParentId, 
                                                                  const std::vector<int>& GenParticles_ParentIdx, const int& nlino1_Idx, const int& nlino2_Idx, const double maxDR,const double maxPTratio) const
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
                bool passDR = GenParticles.at(g).DeltaR(RecoParticles.at(r)) < maxDR;
                bool passPT = abs(1 - RecoParticles.at(r).Pt()/GenParticles.at(g).Pt()) < maxPTratio;
                if (findParent(abs(check_resPartID), g, GenParticles_ParentId, GenParticles_ParentIdx) == check_resPartID && GoodGenParticles.at(g) && check_neutralinos && passDR && passPT)
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

    //function to sort for best matches
    void getMatches(const std::vector<std::tuple< int , int , double>>& AllDR, std::vector<std::pair<int,int>>& Matches, std::vector<bool> availableDR) const
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
            Matches.push_back(std::make_pair(std::get<0>(bestDR), std::get<1>(bestDR)));
            for (unsigned int d=0;  d < AllDR.size(); d++)
            {
                if (std::get<0>(AllDR.at(d)) == std::get<0>(bestDR) || std::get<1>(AllDR.at(d)) == std::get<1>(bestDR))
                {
                    availableDR.at(d) = false;
                }
            }
            getMatches( AllDR, Matches, availableDR);
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
            for(unsigned int j=0; j < GoodJets_pt30.size(); j++)
            {
                if(GoodJets_pt30.at(j)) RecoParticles.push_back(Jets.at(j));  //can replace this with any jet collection
            }

            for (const auto& l : GoodLeptons)
            {
                RecoParticles.push_back(l.second);
            }
            std::vector<int> neutralinos_Idx;
            int nlino1_Idx = -1, nlino2_Idx = -1;

            //Define OkayParticles, which allows leptons/jets by status code and parent
            std::vector<bool> OkayGenParticles(GenParticles.size(), false);
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
                bool pass_jet = is_jet ? status == 23 : false; //pre-radiation jets must have status 23, post-radiation have status 71
                int stopId = findParent(1000006, p, GenParticles_ParentId, GenParticles_ParentIdx);
                bool pass_stop = stopId != -1; //all gen particles must come from a stop
                bool in_acceptance = GenParticles.at(p).Pt() > 30 && abs(GenParticles.at(p).Eta()) < 2.4; //only gen match to jets that can be reconstructed
                bool filter = (pass_lepton || pass_jet) && pass_stop && in_acceptance;

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
            
            std::vector<int>  resParticleList{1000006, -1000006, 1000022, -1000022,  5000002, -5000002}; //gen match for stop, neutralino, and singlet
            std::vector<TLorentzVector> RecoSumList;
            std::vector<TLorentzVector> GenSumList;
            std::vector<std::pair<std::vector<double>, std::vector<double>>> DRandPTSumList;

            TLorentzVector AllGenSum,AllGenSumNlino;
            for (unsigned int g = 0; g < GenParticles.size(); g++)
            {
                if (GoodGenParticles.at(g) && findParent(1000006, g, GenParticles_ParentId, GenParticles_ParentIdx) == 1000006) AllGenSum += GenParticles.at(g);
                if (GoodGenParticles.at(g) && GenParticles_ParentIdx.at(g) == nlino1_Idx) AllGenSumNlino += GenParticles.at(g);
            }

//     For tracking the DR and pT of matched particles to check matching

            auto& GM_Stop1_DR     = tr.createDerivedVec<double>("GM_Stop1_DR"+myVarSuffix_);
            auto& GM_Stop1_PT     = tr.createDerivedVec<double>("GM_Stop1_PT"+myVarSuffix_);
            auto& GM_Stop2_DR     = tr.createDerivedVec<double>("GM_Stop2_DR"+myVarSuffix_);
            auto& GM_Stop2_PT     = tr.createDerivedVec<double>("GM_Stop2_PT"+myVarSuffix_);
            auto& GM_Nlino1_DR     = tr.createDerivedVec<double>("GM_Nlino1_DR"+myVarSuffix_);
            auto& GM_Nlino1_PT     = tr.createDerivedVec<double>("GM_Nlino1_PT"+myVarSuffix_);
            auto& GM_Nlino2_DR     = tr.createDerivedVec<double>("GM_Nlino2_DR"+myVarSuffix_);
            auto& GM_Nlino2_PT     = tr.createDerivedVec<double>("GM_Nlino2_PT"+myVarSuffix_);
//begin matching
            for(unsigned int p=0; p < resParticleList.size(); p++)
            {
                std::vector<std::pair<int, int>> Matches;
                double maxDR = 0.1; //set max DR allowed for matching
                double maxPTratio = 0.5; // set max pT allowed for matching
                
                std::vector<std::tuple< int , int , double>> AllDR = findAllDR(GenParticles, RecoParticles, GoodGenParticles, resParticleList[p], GenParticles_ParentId, GenParticles_ParentIdx, nlino1_Idx, nlino2_Idx, maxDR, maxPTratio);

                std::vector<bool> availableDR(AllDR.size(), true);

                getMatches( AllDR, Matches, availableDR);

                TLorentzVector RecoMatchedSum;
                TLorentzVector GenMatchedSum;
                std::vector<double> DRvec;
                std::vector<double> PTvec;
                for (unsigned int match = 0; match < Matches.size(); match++) //track whatever info about matches is needed
                {
                    GenMatchedSum  += GenParticles.at(Matches.at(match).first);
                    RecoMatchedSum += RecoParticles.at(Matches.at(match).second);
                    DRvec.push_back(GenParticles.at(Matches.at(match).first).DeltaR(RecoParticles.at(Matches.at(match).second)));
                    PTvec.push_back(abs(1 - GenParticles.at(Matches.at(match).first).Pt()/RecoParticles.at(Matches.at(match).second).Pt()));
                }
                //save info for all resonance particles in vector         
                RecoSumList.push_back(RecoMatchedSum);
                GenSumList.push_back(GenMatchedSum);
                DRandPTSumList.push_back(std::make_pair(DRvec, PTvec));
            }
            
            GM_Stop1_DR = DRandPTSumList.at(0).first;
            GM_Stop1_PT = DRandPTSumList.at(0).second;
            GM_Stop2_DR = DRandPTSumList.at(1).first;
            GM_Stop2_PT = DRandPTSumList.at(1).second;
            GM_Nlino1_DR = DRandPTSumList.at(2).first;
            GM_Nlino1_PT = DRandPTSumList.at(2).second;
            GM_Nlino2_DR = DRandPTSumList.at(3).first;
            GM_Nlino2_PT = DRandPTSumList.at(3).second;

            asymm_mt2_lester_bisect::disableCopyrightMessage();
            
            tr.registerDerivedVar("GM_StopMT2"+myVarSuffix_,        ttUtility::coreMT2calc(RecoSumList.at(0),RecoSumList.at(1),lvMET));
            tr.registerDerivedVar("GM_StopGenMT2"+myVarSuffix_,     ttUtility::coreMT2calc(GenSumList.at(0),GenSumList.at(1),lvGenMET));
            tr.registerDerivedVar("GM_Stop1"+myVarSuffix_,      RecoSumList.at(0));
            tr.registerDerivedVar("GM_Stop2"+myVarSuffix_,      RecoSumList.at(1));
            tr.registerDerivedVar("GM_Stop1Gen"+myVarSuffix_,   GenSumList.at(0));
            tr.registerDerivedVar("GM_Stop2Gen"+myVarSuffix_,   GenSumList.at(1));            
            tr.registerDerivedVar("GM_Nlino1"+myVarSuffix_,     RecoSumList.at(2));
            tr.registerDerivedVar("GM_Nlino2"+myVarSuffix_,     RecoSumList.at(3));
            tr.registerDerivedVar("GM_Nlino1Gen"+myVarSuffix_,  GenSumList.at(2));
            tr.registerDerivedVar("GM_Nlino2Gen"+myVarSuffix_,  GenSumList.at(3));
            tr.registerDerivedVar("GM_Single1"+myVarSuffix_,    RecoSumList.at(4));
            tr.registerDerivedVar("GM_Single2"+myVarSuffix_,    RecoSumList.at(5));
            tr.registerDerivedVar("GM_Single1Gen"+myVarSuffix_, GenSumList.at(4));
            tr.registerDerivedVar("GM_Single2Gen"+myVarSuffix_, GenSumList.at(5));
            tr.registerDerivedVar("GM_AllGen"+myVarSuffix_,     AllGenSum);
            tr.registerDerivedVar("GM_AllGenNlino"+myVarSuffix_,AllGenSumNlino);

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
