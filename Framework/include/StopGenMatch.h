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
                //if (std::get<0>(AllDR.at(d)) == std::get<0>(bestDR))
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
            const auto& GenMET                  = tr.getVar<double>("GenMET");
            const auto& GenMETPhi               = tr.getVar<double>("GenMETPhi");

            const auto& filetag                 = tr.getVar<std::string>("filetag");

            int commonAncestor = 6;
            if (filetag.find("mStop") != std::string::npos)
                commonAncestor = 1000006;

            utility::LorentzVector lvMET;
            lvMET.SetPt(MET); lvMET.SetEta(0.0); lvMET.SetPhi(METPhi); lvMET.SetE(MET);
            
            utility::LorentzVector lvGenMET;
            lvGenMET.SetPt(GenMET); lvGenMET.SetEta(0.0); lvGenMET.SetPhi(GenMETPhi); lvGenMET.SetE(GenMET);
            
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
            std::vector<utility::LorentzVector> GenSumList;
            std::vector<std::pair<std::vector<double>, std::vector<double>>> DRandPTSumList;
            std::vector<std::vector<int> > pdgsList;
            std::vector<std::vector<int> > momList;
            std::vector<std::vector<double> > genEtaList;
            std::vector<std::vector<double> > recEtaList;
            std::vector<std::vector<double> > genPhiList;
            std::vector<std::vector<double> > recPhiList;
            std::vector<std::vector<double> > genPtList;
            std::vector<std::vector<double> > recPtList;

            utility::LorentzVector AllGenSum,AllGenSumNlino;
            for (unsigned int g = 0; g < GenParticles.size(); g++)
            {
                if (GoodGenParticles.at(g) && findParent(commonAncestor, g, GenParticles_ParentId, GenParticles_ParentIdx) == commonAncestor) AllGenSum += GenParticles.at(g);
            }

//     For tracking the DR and pT of matched particles to check matching

            auto& GM_Stop1_DR     = tr.createDerivedVec<double>("GM_Stop1_DR"+myVarSuffix_);
            auto& GM_Stop1_PT     = tr.createDerivedVec<double>("GM_Stop1_PT"+myVarSuffix_);
            auto& GM_Stop2_DR     = tr.createDerivedVec<double>("GM_Stop2_DR"+myVarSuffix_);
            auto& GM_Stop2_PT     = tr.createDerivedVec<double>("GM_Stop2_PT"+myVarSuffix_);
            auto& GM_Stop1_pdgs    = tr.createDerivedVec<int>("GM_Stop1_pdgs"+myVarSuffix_);
            auto& GM_Stop2_pdgs    = tr.createDerivedVec<int>("GM_Stop2_pdgs"+myVarSuffix_);
            auto& GM_Stop1_mom    = tr.createDerivedVec<int>("GM_Stop1_mom"+myVarSuffix_);
            auto& GM_Stop2_mom    = tr.createDerivedVec<int>("GM_Stop2_mom"+myVarSuffix_);
            auto& GM_Stop1_genpts    = tr.createDerivedVec<double>("GM_Stop1_genpts"+myVarSuffix_);
            auto& GM_Stop2_genpts    = tr.createDerivedVec<double>("GM_Stop2_genpts"+myVarSuffix_);
            auto& GM_Stop1_genphis    = tr.createDerivedVec<double>("GM_Stop1_genphis"+myVarSuffix_);
            auto& GM_Stop2_genphis    = tr.createDerivedVec<double>("GM_Stop2_genphis"+myVarSuffix_);
            auto& GM_Stop1_genetas    = tr.createDerivedVec<double>("GM_Stop1_genetas"+myVarSuffix_);
            auto& GM_Stop2_genetas    = tr.createDerivedVec<double>("GM_Stop2_genetas"+myVarSuffix_);
            auto& GM_Stop1_recpts    = tr.createDerivedVec<double>("GM_Stop1_recpts"+myVarSuffix_);
            auto& GM_Stop2_recpts    = tr.createDerivedVec<double>("GM_Stop2_recpts"+myVarSuffix_);
            auto& GM_Stop1_recphis    = tr.createDerivedVec<double>("GM_Stop1_recphis"+myVarSuffix_);
            auto& GM_Stop2_recphis    = tr.createDerivedVec<double>("GM_Stop2_recphis"+myVarSuffix_);
            auto& GM_Stop1_recetas    = tr.createDerivedVec<double>("GM_Stop1_recetas"+myVarSuffix_);
            auto& GM_Stop2_recetas    = tr.createDerivedVec<double>("GM_Stop2_recetas"+myVarSuffix_);

            for(unsigned int p=0; p < resParticleList.size(); p++)
            {
                std::vector<std::pair<int, int>> Matches;
                double maxDR = 0.4; //set max DR allowed for matching
                double maxPTratio = 999.0; // set max pT allowed for matching
                
                std::vector<std::tuple< int , int , double>> AllDR = findAllDR(GenParticles, RecoParticles, GoodGenParticles, resParticleList[p], GenParticles_ParentId, GenParticles_ParentIdx, maxDR, maxPTratio);

                std::vector<bool> availableDR(AllDR.size(), true);

                getMatches( AllDR, Matches, availableDR);

                utility::LorentzVector RecoMatchedSum;
                utility::LorentzVector GenMatchedSum;
                std::vector<double> DRvec;
                std::vector<double> PTvec;
                std::vector<int> pdgs;
                std::vector<int> mom;
                std::vector<double> genpts;
                std::vector<double> genetas;
                std::vector<double> genphis;
                std::vector<double> recpts;
                std::vector<double> recetas;
                std::vector<double> recphis;
            
                // Matches may have multiple GEN matched to a single RECO
                // So let's do some sneaky processing
                unsigned int skip = 0;
                while (skip < Matches.size())
                {
                    auto theGen = GenParticles.at(Matches.at(skip).first);
                    auto theRec = RecoParticles.at(Matches.at(skip).second);
                    //int recIdx = Matches.at(skip).second;
                    int genIdx = Matches.at(skip).first;

                    // Starting from this element onward we may have the same reco appearing in a row matched
                    // to multiple gen, so try and add all those gen together.
                    //for (unsigned int j = skip+1; j < Matches.size(); j++)
                    //{

                    //    int currentRecIdx = Matches.at(j).second;
                    //    skip = j;
   
                    //    // Here find another gen with the same reco match
                    //    if (recIdx == currentRecIdx)
                    //    {
                    //        theGen += GenParticles.at(Matches.at(j).first);
                    //    } else {
                    //        break;
                    //    }
                    //}

                    //std::cout << "idx: " << genIdx << " sta: " << GenParticles_Status.at(genIdx) << " mom: " << GenParticles_ParentId.at(genIdx) << " dau: " << GenParticles_PdgId.at(genIdx) << " dR: " << theGen.DeltaR(theRec) << " PtR: " << theGen.Pt()/theRec.Pt() << " eta: " << theRec.Eta() << " Phi: " << theRec.Phi() << std::endl; 

                    GenMatchedSum  += theGen;
                    RecoMatchedSum += theRec;
                    DRvec.push_back(utility::DeltaR(theGen, theRec));
                    PTvec.push_back(theGen.Pt()/theRec.Pt());
                    pdgs.push_back(GenParticles_PdgId.at(genIdx));
                    mom.push_back(GenParticles_ParentId.at(genIdx));
                    genetas.push_back(theGen.Eta());
                    genphis.push_back(theGen.Phi());
                    genpts.push_back(theGen.Pt());
                    recetas.push_back(theRec.Eta());
                    recphis.push_back(theRec.Phi());
                    recpts.push_back(theRec.Pt());

                    skip++;
                    if (skip == Matches.size()) {
                        break;
                    }
                }
                //std::cout << std::endl;
                //save info for all resonance particles in vector         
                RecoSumList.push_back(RecoMatchedSum);
                GenSumList.push_back(GenMatchedSum);
                DRandPTSumList.push_back(std::make_pair(DRvec, PTvec));
                pdgsList.push_back(pdgs);
                momList.push_back(mom);
                genEtaList.push_back(genetas);
                genPhiList.push_back(genphis);
                genPtList.push_back(genpts);
                recEtaList.push_back(recetas);
                recPhiList.push_back(recphis);
                recPtList.push_back(recpts);
            }
            
            // Let's do some weird stuff
            // For either signal or ttbar, we will have stops or "stops"
            GM_Stop1_DR = DRandPTSumList.at(0).first;
            GM_Stop1_PT = DRandPTSumList.at(0).second;

            GM_Stop2_DR = DRandPTSumList.at(1).first;
            GM_Stop2_PT = DRandPTSumList.at(1).second;

            GM_Stop1_pdgs = pdgsList.at(0);
            GM_Stop2_pdgs = pdgsList.at(1);

            GM_Stop1_mom = momList.at(0);
            GM_Stop2_mom = momList.at(1);

            GM_Stop1_genetas = genEtaList.at(0);
            GM_Stop2_genetas = genEtaList.at(1);

            GM_Stop1_recetas = recEtaList.at(0);
            GM_Stop2_recetas = recEtaList.at(1);

            GM_Stop1_genpts = genPtList.at(0);
            GM_Stop2_genpts = genPtList.at(1);

            GM_Stop1_recpts = recPtList.at(0);
            GM_Stop2_recpts = recPtList.at(1);

            GM_Stop1_genphis = genPhiList.at(0);
            GM_Stop2_genphis = genPhiList.at(1);

            GM_Stop1_recphis = recPhiList.at(0);
            GM_Stop2_recphis = recPhiList.at(1);

            tr.registerDerivedVar("GM_StopMT2"+myVarSuffix_,        ttUtility::coreMT2calc(utility::convertLV<TLorentzVector, utility::LorentzVector>(RecoSumList.at(0)),utility::convertLV<TLorentzVector, utility::LorentzVector>(RecoSumList.at(1)),utility::convertLV<TLorentzVector, utility::LorentzVector>(lvMET)));
            tr.registerDerivedVar("GM_StopGenMT2"+myVarSuffix_,     ttUtility::coreMT2calc(utility::convertLV<TLorentzVector, utility::LorentzVector>(GenSumList.at(0)), utility::convertLV<TLorentzVector, utility::LorentzVector>(GenSumList.at(1)),utility::convertLV<TLorentzVector, utility::LorentzVector>(lvGenMET)));
            tr.registerDerivedVar("GM_Stop1"+myVarSuffix_,      RecoSumList.at(0));
            tr.registerDerivedVar("GM_Stop2"+myVarSuffix_,      RecoSumList.at(1));
            tr.registerDerivedVar("GM_Stop1Gen"+myVarSuffix_,   GenSumList.at(0));
            tr.registerDerivedVar("GM_Stop2Gen"+myVarSuffix_,   GenSumList.at(1));            
            tr.registerDerivedVar("GM_AllGen"+myVarSuffix_,     AllGenSum);
            tr.registerDerivedVar("GM_AllGenNlino"+myVarSuffix_,AllGenSumNlino);

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
        genMatch(tr);
    }
};

#endif
