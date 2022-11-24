#ifndef STOPGENMATCH_H
#define STOPGENMATCH_H

// This class does gen matching to reconstruct the stops in an event both on gen level and reco level
// In the event, that there are no stops e.g. in just ttbar+jets scenario, the "stop" should be close to the mass of the top
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
    inline std::vector<std::tuple< int , int , double>> findAllDR(const std::vector<utility::LorentzVector>& GenParticles, const std::vector<utility::LorentzVector>& RecoParticles, 
                                                                  const std::vector<bool>& GoodGenParticles, const int resPartID, const std::vector<int>& GenParticles_ParentId, 
                                                                  const std::vector<int>& GenParticles_ParentIdx, const double maxDR,const double maxPTratio) const
    {
        std::vector<std::tuple< int , int, double>> AllDR;
        std::tuple< int , int, double> DRtup;
        int check_resPartID = resPartID;
        for (unsigned int r=0; r < RecoParticles.size(); r++)
        {
            for (unsigned int g=0; g < GenParticles.size(); g++)
            {
                bool passDR = utility::DeltaR(GenParticles.at(g), RecoParticles.at(r)) < maxDR;
                bool passPT = (GenParticles.at(g).Pt()/RecoParticles.at(r).Pt() > 1.0-maxPTratio and GenParticles.at(g).Pt()/RecoParticles.at(r).Pt() < 1.0+maxPTratio);
                if (findParent(abs(check_resPartID), g, GenParticles_ParentId, GenParticles_ParentIdx) == check_resPartID && GoodGenParticles.at(g) && passDR && passPT)
                {
                    std::get<0>(DRtup) = g;
                    std::get<1>(DRtup) = r;
                    std::get<2>(DRtup) = utility::DeltaR(GenParticles.at(g), RecoParticles.at(r));
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
                
                // We only take care when a gen is used, thus it cannot be matched to multiple reco
                // only matched to best reco. On the other hand, multiple gen are free to be matched
                // to a single reco.
                if (std::get<0>(AllDR.at(d)) == std::get<0>(bestDR) or (std::get<1>(AllDR.at(d)) == std::get<1>(bestDR)))
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
            const auto& GenParticles        = tr.getVec<utility::LorentzVector>("GenParticles");
            const auto& GenParticles_PdgId      = tr.getVec<int>("GenParticles_PdgId");
            const auto& GenParticles_ParentId   = tr.getVec<int>("GenParticles_ParentId");
            const auto& GenParticles_ParentIdx  = tr.getVec<int>("GenParticles_ParentIdx");
            const auto& GenParticles_Status     = tr.getVec<int>("GenParticles_Status");            
            const auto& Jets                    = tr.getVec<utility::LorentzVector>("Jets"+myVarSuffix_);
            const auto& Electrons               = tr.getVec<utility::LorentzVector>("Electrons");
            const auto& Muons                   = tr.getVec<utility::LorentzVector>("Muons");
            const auto& MET                     = tr.getVar<float>("MET");
            const auto& METPhi                  = tr.getVar<float>("METPhi");

            const auto& filetag                 = tr.getVar<std::string>("filetag");

            int commonAncestor = 6;
            if (filetag.find("mStop") != std::string::npos)
                commonAncestor = 1000006;

            utility::LorentzVector lvMET;
            lvMET.SetPt(MET); lvMET.SetEta(0.0); lvMET.SetPhi(METPhi); lvMET.SetE(MET);
            
            std::vector<utility::LorentzVector> RecoParticles;
            for(unsigned int j=0; j < Jets.size(); j++)
            {
                RecoParticles.push_back(Jets.at(j));  //can replace this with any jet collection
            }

            for (const auto& e : Electrons)
            {
                RecoParticles.push_back(e);
            }

            for (const auto& m : Muons)
            {
                RecoParticles.push_back(m);
            }

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
                int topId = findParent(6, p, GenParticles_ParentId, GenParticles_ParentIdx);
                bool pass_stop = stopId != -1; //all gen particles must come from a stop
                bool pass_top = topId != -1; //for ttbar, gen particles can come from top
                bool filter = (pass_lepton || pass_jet) && (pass_stop || pass_top) && true;//in_acceptance;

                if (filter)
                {
                    if (pass_lepton && WId == 24 ) WPlusLeps.push_back(p);
                    if (pass_lepton && WId == -24) WMinusLeps.push_back(p);
                    if ((abs(pdgid) == 11 || abs(pdgid) ==13) && (pass_stop || pass_top) &&  abs(momid) == 15) TauLepCounter += 1;
                    OkayGenParticles.at(p) = true;
                }
            }
            //Define GoodGenParticles, which has no W radiation (for most part) and allows undecayed taus
            std::vector<bool> GoodGenParticles = OkayGenParticles;
            bool wplus_eplus = false, wplus_eminus = false, wplus_muplus = false, wplus_muminus =false;
            bool wminus_eplus = false, wminus_eminus = false, wminus_muplus = false, wminus_muminus =false;
            int wpem = 0, wpep = 0, wpmm = 0, wpmp = 0;
            int wmem = 0, wmep = 0, wmmm = 0, wmmp = 0;
           
            for ( const auto& w : WPlusLeps) 
            {
                if (GenParticles_PdgId.at(w) == 13)
                {
                    wplus_eminus = true;
                    wpem = w;
                }
                if (GenParticles_PdgId.at(w) == -13)
                {
                    wplus_eplus = true;
                    wpep = w;
                }
                if (GenParticles_PdgId.at(w) == 11)
                {
                    wplus_muminus = true;
                    wpmm = w;
                }
                if (GenParticles_PdgId.at(w) == -11)
                {
                    wplus_muplus = true;
                    wpmp = w;
                }
                    
            } //removes pair produced leptons from radation. does not remove multiple pairs, need to fix for future studies
            for (const auto& w : WMinusLeps)
            {
                if (GenParticles_PdgId.at(w) == 13)
                {
                    wminus_eminus = true;
                    wmem = w;
                }
                if (GenParticles_PdgId.at(w) == -13)
                {
                    wminus_eplus = true;
                    wmep = w;
                }
                if (GenParticles_PdgId.at(w) == 11)
                {
                    wminus_muminus = true;
                    wmmm = w;
                }
                if (GenParticles_PdgId.at(w) == -11)
                {
                    wminus_muplus = true;
                    wmmp = w;
                }

            }
            if (wplus_eplus && wplus_eminus) 
            {
                GoodGenParticles.at(wpep) = false;
                GoodGenParticles.at(wpem) = false;
            }
            if (wplus_muplus && wplus_muminus)
            {
                GoodGenParticles.at(wpmp) = false;
                GoodGenParticles.at(wpmm) = false;
            }
            if (wminus_eminus && wminus_eplus) 
            {
                GoodGenParticles.at(wmep) = false;
                GoodGenParticles.at(wmem) = false;
            }
            if (wminus_muminus && wminus_muplus)
            {
                GoodGenParticles.at(wmmp) = false;
                GoodGenParticles.at(wmmm) = false;
            }
               
            bool keepTau = TauLepCounter == 0;
            for (unsigned int g=0; g < OkayGenParticles.size(); g++)
            {
                if (abs(GenParticles_PdgId.at(g)) == 15 && GenParticles_Status.at(g) == 2 && keepTau && findParent(commonAncestor, g, GenParticles_ParentId, GenParticles_ParentIdx) != -1) GoodGenParticles.at(g) = true;
            }

            std::vector<int> resParticleList{commonAncestor, -commonAncestor}; //gen match for stop, neutralino, and singlet

            std::vector<utility::LorentzVector> RecoSumList;

            for(unsigned int p=0; p < resParticleList.size(); p++)
            {
                std::vector<std::pair<int, int>> Matches;
                double maxDR = 0.4; //set max DR allowed for matching
                double maxPTratio = 999.0; // set max pT allowed for matching
                
                std::vector<std::tuple< int , int , double>> AllDR = findAllDR(GenParticles, RecoParticles, GoodGenParticles, resParticleList[p], GenParticles_ParentId, GenParticles_ParentIdx, maxDR, maxPTratio);

                std::vector<bool> availableDR(AllDR.size(), true);

                getMatches(AllDR, Matches, availableDR);

                utility::LorentzVector RecoMatchedSum;
            
                // Matches may have multiple GEN matched to a single RECO
                // So let's do some sneaky processing
                unsigned int skip = 0;
                while (skip < Matches.size())
                {
                    auto theRec = RecoParticles.at(Matches.at(skip).second);

                    RecoMatchedSum += theRec;

                    skip++;
                    if (skip == Matches.size()) {
                        break;
                    }
                }
                //save info for all resonance particles in vector         
                RecoSumList.push_back(RecoMatchedSum);
            }
            
            // For NN ntuples
            if (RecoSumList.at(0).Pt() > RecoSumList.at(1).Pt())
            {
                tr.registerDerivedVar("stop1_ptrank_mass"+myVarSuffix_, RecoSumList.at(0).M());
                tr.registerDerivedVar("stop2_ptrank_mass"+myVarSuffix_, RecoSumList.at(1).M());
            } else {
                tr.registerDerivedVar("stop1_ptrank_mass"+myVarSuffix_, RecoSumList.at(1).M());
                tr.registerDerivedVar("stop2_ptrank_mass"+myVarSuffix_, RecoSumList.at(0).M());
            }
            if (RecoSumList.at(0).M() > RecoSumList.at(1).M())
            {
                tr.registerDerivedVar("stop1_mrank_mass"+myVarSuffix_, RecoSumList.at(0).M());
                tr.registerDerivedVar("stop2_mrank_mass"+myVarSuffix_, RecoSumList.at(1).M());
            } else {
                tr.registerDerivedVar("stop1_mrank_mass"+myVarSuffix_, RecoSumList.at(1).M());
                tr.registerDerivedVar("stop2_mrank_mass"+myVarSuffix_, RecoSumList.at(0).M());
            }

            tr.registerDerivedVar("stop_avemass"+myVarSuffix_, (RecoSumList.at(0).M()+RecoSumList.at(1).M())/2.0);
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
        const auto& lostCauseEvent = tr.getVar<bool>("lostCauseEvent" + myVarSuffix_);
        const auto& fastMode       = tr.getVar<bool>("fastMode");

        if (!lostCauseEvent or !fastMode or tr.isFirstEvent())
            genMatch(tr);
    }
};

#endif
