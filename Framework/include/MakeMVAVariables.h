#ifndef MAKEMVAVARIABLES_H
#define MAKEMVAVARIABLES_H 

#include "Framework/Framework/include/EventShapeVariables.h"
#include "Framework/Framework/include/Utility.h"

class MakeMVAVariables
{
private:
    class ComboLV
    {
    public:
        TLorentzVector v1;
        TLorentzVector v2;
        std::vector<int> jetCombo;
    };

    bool verb_;
    std::string myVarSuffix_;
    bool doGenMatch_;
    unsigned int nTopJets_, nLeptons_;
    std::string GoodJetsName_, NGoodJetsName_, GoodLeptonsName_, NGoodLeptonsName_, MVAJetName_, MVALeptonName_, ESVarName_;

    std::vector<int> decToBinary(int n, int max) const
    {
        std::vector<int> binaryNum(max, 0);
        int i = 0;
        while(n > 0) 
        {
            binaryNum[i] = n % 2;
            n = n / 2;
            i++;
        }
        return binaryNum;
    }

    bool genMatch(const NTupleReader& tr, const std::vector<TLorentzVector>& lv_all , const std::vector<int>& jetCombo) const
    {
        const auto& runtype = tr.getVar<std::string>("runtype");

        if(runtype != "Data" && jetCombo.size() > 0)
        {
            const auto& hadtopdaughters = tr.getVec<std::vector<const TLorentzVector*>>("hadtopdaughters"+myVarSuffix_);

            bool genMatched = true;
            int topID = 1, megaJetID = -1;
            for(const auto& daughters : hadtopdaughters)
            {
                int numMatchedJets = 1;
                for(const auto* d : daughters)
                {
                    for(unsigned int ijet = 0; ijet < lv_all.size(); ijet++)
                    {
                        double deltaR = d->DeltaR( lv_all.at(ijet) );
                        if(deltaR < 0.4)
                        {
                            if(numMatchedJets != 1 && genMatched) genMatched = megaJetID == jetCombo[ijet]; 
                            //std::cout<<"TopID: "<<topID<<" NDaughter: "<<daughters.size()<<" Mega Jet: "<<jetCombo[ijet]<<" GenMatch "<<genMatched<<std::endl;
                            if(topID != 1 && numMatchedJets == 1 && megaJetID == jetCombo[ijet]) genMatched = false;  
                            megaJetID = jetCombo[ijet];
                            numMatchedJets++;
                        }
                    }
                }
                topID++;
            }
            return genMatched;
        }
        else 
            return false;
    }

    void setJetCollection(const std::string& jetColl)
    {
        if(jetColl == "NonIsoMuonJets_pt30")
        {
            GoodJetsName_ = "NonIsoMuonJets_pt30";
            NGoodJetsName_ = "NNonIsoMuonJets_pt30";
            GoodLeptonsName_ = "GoodNonIsoMuons";
            NGoodLeptonsName_ = "NNonIsoMuons";
            MVAJetName_ = "JetNonIsoMuons";
            MVALeptonName_ ="GoodNonIsoMuons";
            ESVarName_ = "NonIsoMuons_";
        }
        else
        {
            GoodJetsName_ = "GoodJets_pt30";
            NGoodJetsName_ = "NGoodJets_pt30";
            GoodLeptonsName_ = "GoodLeptons";
            NGoodLeptonsName_ = "NGoodLeptons";
            MVAJetName_ = "Jet";
            MVALeptonName_ ="GoodLeptons";
            ESVarName_ = "";
        }
    }

    void makeMVAVariables(NTupleReader& tr)
    {
        const auto& Jets = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
        const auto& GoodJets = tr.getVec<bool>(GoodJetsName_+myVarSuffix_);
        //const auto& NGoodJets = tr.getVar<int>(NGoodJetsName_+myVarSuffix_);
        const auto& GoodLeptons = tr.getVec<std::pair<std::string, TLorentzVector>>(GoodLeptonsName_+myVarSuffix_);
        const auto& NGoodLeptons = tr.getVar<int>(NGoodLeptonsName_+myVarSuffix_);
        const auto& MET = tr.getVar<double>("MET"); 
        const auto& METPhi = tr.getVar<double>("METPhi");

        //--- Get the 4-vec for the MET
        TLorentzVector lvMET;
        lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);

        //--- Sum all jets
        TLorentzVector rlv_all;
        TLorentzVector rlv_pt20;
        for(auto jlv : Jets)
        {
            rlv_all += jlv;
            if (jlv.Pt() > 20.) rlv_pt20 += jlv;
        }

        //--- Fill vector of jet momenta in CM frame.
        //    Boost to the CM frame
        //TLorentzVector jetSum = rlv_all;
        //TVector3 cm_boost_beta_vec( -jetSum.BoostVector() );

        //    Boost to the CM frame (only in z)
        double event_beta_z = -9.0;
        double event_beta_z_pt20 = -9.0;
        double reco_jets_beta = rlv_all.Pz() / rlv_all.E();
        event_beta_z = reco_jets_beta;
        event_beta_z_pt20 = rlv_pt20.Pz() / rlv_pt20.E();
        TVector3 rec_boost_beta_vec( 0.0, 0.0, -reco_jets_beta );
        auto& cm_jets = tr.createDerivedVec<math::RThetaPhiVector>(ESVarName_+"cm_jets"+myVarSuffix_);
        auto& Jets_cm = tr.createDerivedVec<TLorentzVector>(ESVarName_+"Jets_cm"+myVarSuffix_);
        std::vector<TLorentzVector> Jets_;

        //--- Boost the GoodJets, Goodleptons, and the MET in the event 
        for(unsigned int j = 0; j < Jets.size(); j++)
        {
            if(!GoodJets[j]) continue;
            TLorentzVector jlvcm = Jets.at(j);            
            Jets_.push_back( jlvcm );
            
            jlvcm.Boost( rec_boost_beta_vec );
            Jets_cm.push_back( jlvcm );

            math::RThetaPhiVector cmvec( jlvcm.P(), jlvcm.Theta(), jlvcm.Phi() );
            cm_jets.push_back( cmvec );
        }
        auto GoodLeptons_cm = std::make_unique<std::vector<TLorentzVector>>();
        for(auto pair : GoodLeptons)
        {
            pair.second.Boost( rec_boost_beta_vec );
            GoodLeptons_cm->push_back( pair.second );
        }
        TLorentzVector lvMET_cm = lvMET;
        lvMET_cm.Boost( rec_boost_beta_vec );           

        //--- Try using only the 7 highest-P jets in the CM frame in the event shape vars.
        //    First, need to make a new input vector of jets containing only those jets.
        auto cm_jets_psort = cm_jets ;
        std::sort( cm_jets_psort.begin(), cm_jets_psort.end(), utility::compare_p ) ;
        std::vector<math::RThetaPhiVector> cm_jets_top6 ;

        auto Jets_cm_psort = Jets_cm;
        auto Jets_psort = Jets_;
        std::sort( Jets_cm_psort.begin(), Jets_cm_psort.end(), [](TLorentzVector v1, TLorentzVector v2){return v1.P() > v2.P();} );
        std::sort( Jets_psort.begin(), Jets_psort.end(), [](TLorentzVector v1, TLorentzVector v2){return v1.Pt() > v2.Pt();} );
        auto& Jets_cm_top6 = tr.createDerivedVec<TLorentzVector>(ESVarName_+"Jets_cm_top6"+myVarSuffix_);
        auto& Jets_top6 = tr.createDerivedVec<TLorentzVector>(ESVarName_+"Jets_top6"+myVarSuffix_);

        for(unsigned int ji=0; ji<cm_jets.size(); ji++ ) 
        {
            if ( ji < nTopJets_ ) 
            {
                cm_jets_top6.push_back( cm_jets_psort.at(ji) ) ;
                Jets_cm_top6.push_back( Jets_cm_psort.at(ji) ) ;
                Jets_top6.push_back( Jets_psort.at(ji) ) ;
            }
        } // ji

        if( verb_ ) 
        {
            printf("\n\n Unsorted and sorted CM jet lists.\n") ;
            for ( unsigned int ji=0; ji<cm_jets.size(); ji++ ) 
            {
                printf("  %2d :  (%7.1f, %7.3f, %7.3f) | (%7.1f, %7.3f, %7.3f)\n", ji,
                       cm_jets.at(ji).R(), cm_jets.at(ji).Theta(), cm_jets.at(ji).Phi(),
                       cm_jets_psort.at(ji).R(), cm_jets_psort.at(ji).Theta(), cm_jets_psort.at(ji).Phi() ) ;
                
                printf("  %2d :  (%7.1f, %7.3f, %7.3f) | (%7.1f, %7.3f, %7.3f)\n", ji,
                       Jets_cm.at(ji).P(), Jets_cm.at(ji).Theta(), Jets_cm.at(ji).Phi(),
                       Jets_cm_top6.at(ji).P(), Jets_cm_top6.at(ji).Theta(), Jets_cm_top6.at(ji).Phi() ) ;
            } // ji
            printf("\n\n") ;
        }

        //--- Make and get the event shape variables for the 6 highest-P jets in the CM frame
        EventShapeVariables esv_top6( cm_jets_top6 ) ;
        TVectorD eigen_vals_norm_top6 = esv_top6.getEigenValues() ;

        double fwm2_top6 = esv_top6.getFWmoment( 2 ) ;
        double fwm3_top6 = esv_top6.getFWmoment( 3 ) ;
        double fwm4_top6 = esv_top6.getFWmoment( 4 ) ;
        double fwm5_top6 = esv_top6.getFWmoment( 5 ) ;
        double fwm6_top6 = esv_top6.getFWmoment( 6 ) ;
        double fwm7_top6 = esv_top6.getFWmoment( 7 ) ;
        double fwm8_top6 = esv_top6.getFWmoment( 8 ) ;
        double fwm9_top6 = esv_top6.getFWmoment( 9 ) ;
        double fwm10_top6 = esv_top6.getFWmoment( 10 ) ;
        double jmt_ev0_top6 = eigen_vals_norm_top6[0] ;
        double jmt_ev1_top6 = eigen_vals_norm_top6[1] ;
        double jmt_ev2_top6 = eigen_vals_norm_top6[2] ;

        //--- register Variables
        //
        for(unsigned int i = 0; i < nTopJets_; i++)
        {
            tr.registerDerivedVar(MVAJetName_+"_pt_"+std::to_string(i+1)+myVarSuffix_,  static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).Pt()  : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_eta_"+std::to_string(i+1)+myVarSuffix_, static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).Eta() : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_phi_"+std::to_string(i+1)+myVarSuffix_, static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).Phi() : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_m_"+std::to_string(i+1)+myVarSuffix_,   static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).M()   : 0.0));
        }
        for(unsigned int i = 0; i < nLeptons_; i++)
        {
            tr.registerDerivedVar(MVALeptonName_+"_pt_"+std::to_string(i+1)+myVarSuffix_,  static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).Pt()  : 0.0));
            tr.registerDerivedVar(MVALeptonName_+"_eta_"+std::to_string(i+1)+myVarSuffix_, static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).Eta() : 0.0));
            tr.registerDerivedVar(MVALeptonName_+"_phi_"+std::to_string(i+1)+myVarSuffix_, static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).Phi() : 0.0));
            tr.registerDerivedVar(MVALeptonName_+"_m_"+std::to_string(i+1)+myVarSuffix_,   static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).M()   : 0.0));
        }

        tr.registerDerivedVar(ESVarName_+"lvMET_cm_pt"+myVarSuffix_,  static_cast<double>( lvMET_cm.Pt() ));
        tr.registerDerivedVar(ESVarName_+"lvMET_cm_eta"+myVarSuffix_, static_cast<double>( lvMET_cm.Eta()));
        tr.registerDerivedVar(ESVarName_+"lvMET_cm_phi"+myVarSuffix_, static_cast<double>( lvMET_cm.Phi()));
        tr.registerDerivedVar(ESVarName_+"lvMET_cm_m"+myVarSuffix_,   static_cast<double>( lvMET_cm.M()  ));
        tr.registerDerivedVar(ESVarName_+"fwm2_top6"+myVarSuffix_, fwm2_top6);
        tr.registerDerivedVar(ESVarName_+"fwm3_top6"+myVarSuffix_, fwm3_top6);
        tr.registerDerivedVar(ESVarName_+"fwm4_top6"+myVarSuffix_, fwm4_top6);
        tr.registerDerivedVar(ESVarName_+"fwm5_top6"+myVarSuffix_, fwm5_top6);
        tr.registerDerivedVar(ESVarName_+"fwm6_top6"+myVarSuffix_, fwm6_top6);
        tr.registerDerivedVar(ESVarName_+"fwm7_top6"+myVarSuffix_, fwm7_top6);
        tr.registerDerivedVar(ESVarName_+"fwm8_top6"+myVarSuffix_, fwm8_top6);
        tr.registerDerivedVar(ESVarName_+"fwm9_top6"+myVarSuffix_, fwm9_top6);
        tr.registerDerivedVar(ESVarName_+"fwm10_top6"+myVarSuffix_, fwm10_top6);
        tr.registerDerivedVar(ESVarName_+"jmt_ev0_top6"+myVarSuffix_, jmt_ev0_top6);
        tr.registerDerivedVar(ESVarName_+"jmt_ev1_top6"+myVarSuffix_, jmt_ev1_top6);
        tr.registerDerivedVar(ESVarName_+"jmt_ev2_top6"+myVarSuffix_, jmt_ev2_top6);
        tr.registerDerivedVar(ESVarName_+"event_beta_z"+myVarSuffix_, event_beta_z);
        tr.registerDerivedVar(ESVarName_+"event_beta_z_pt20"+myVarSuffix_, event_beta_z_pt20);

        // Sum jets, leptons, and MET in the CM frame to reco the SUSY particles
        std::pair<TLorentzVector, TLorentzVector> BestCombo, genBestCombo;
        bool genMatched = false;
        if(NGoodLeptons == 1 && doGenMatch_)
        {
            //--- Making a vector of all Jets, leptons, and MET
            std::vector<TLorentzVector> lv_all;
            for(unsigned int j = 0; j < Jets.size(); j++)
            {
                if(!GoodJets[j]) continue;
                lv_all.push_back( Jets.at(j) );
            }
            //adding the lepton and MET for 1 lepton selection
            lv_all.push_back( GoodLeptons[0].second + lvMET );
            // This wont work now that GoodLeptons is a std::vector<std::pair>
            //lv_all.insert( lv_all.end(), GoodLeptons.begin(), GoodLeptons.end() );
            //lv_all.push_back(lvMET);

            //--- Form all possible combos of summing lv into two lv
            //std::cout<<"----------------------"<<"Jets size: "<<NGoodJets<<" Total lv in event:"<<lv_all.size()<<"-----------------------"<<std::endl;
            int NAll = lv_all.size();
            int maxDec = 0;
            std::vector<ComboLV> combinedLV;
            for(int i = 1; i <= NAll - 1; i++) maxDec += pow(2,NAll - i);
            for(int i = 1; i <= maxDec; i++) 
            {
                std::vector<int> jetCombo = decToBinary( i, NAll);
                TLorentzVector v1, v2;
                int tempCount1 = 0, tempCount2 = 0;
                for(unsigned int ijet = 0; ijet < lv_all.size(); ijet++)
                {
                    if(jetCombo[ijet] == 0) 
                    {
                        v1 += lv_all.at(ijet);
                        tempCount1++;
                    }
                    else if(jetCombo[ijet] == 1)
                    {
                        v2 += lv_all.at(ijet);
                        tempCount2++;
                    }
                }
                //for (int j = lv_all.size() - 1; j >= 0; j--) std::cout << jetCombo[j];                    
                //std::cout<<" "<<tempCount1<<" "<<tempCount2<<" ("<<v1.M()<<", "<<v1.Pt()<<", "<<v1.Eta()<<", "<<v1.Phi()<<") ("
                //         <<v2.M()<<", "<<v2.Pt()<<", "<<v2.Eta()<<", "<<v2.Phi()<<")"<<std::endl;;
                combinedLV.push_back( {v1, v2, jetCombo} );
            }
            //std::cout<<combinedLV.size()<<std::endl;

            //--- Find the best combo of lv pair, looks more like the pair production
            double massDiff = 999999999999;
            std::vector<int> BestJetCombo;
            bool matched = false;
            for(const auto& cLV : combinedLV)
            {
                double mD = abs( cLV.v1.M() - cLV.v2.M() );
                if(mD < massDiff)
                {
                    massDiff = mD;
                    BestCombo = std::make_pair(cLV.v1, cLV.v2);
                    BestJetCombo = cLV.jetCombo;
                }

                bool m = genMatch(tr, lv_all, cLV.jetCombo);
                if(m && !matched)
                {
                    genBestCombo = std::make_pair(cLV.v1, cLV.v2);
                    matched = true;
                }
            }
            genMatched = genMatch(tr, lv_all, BestJetCombo);
            if(!matched) genBestCombo = BestCombo;
            //std::cout<<"Best mass diff Jets: ("<<BestCombo.first.M()<<", "<<BestCombo.first.Pt()<<", "<<BestCombo.first.Eta()<<", "<<BestCombo.first.Phi()<<") ("
            //         <<BestCombo.second.M()<<", "<<BestCombo.second.Pt()<<", "<<BestCombo.second.Eta()<<", "<<BestCombo.second.Phi()<<"):"
            //         <<" GenMatched: "<<genMatched<<std::endl;
            //std::cout<<"Best mass diff Jets: ("<<genBestCombo.first.M()<<", "<<genBestCombo.first.Pt()<<", "<<genBestCombo.first.Eta()<<", "<<genBestCombo.first.Phi()<<") ("
            //         <<genBestCombo.second.M()<<", "<<genBestCombo.second.Pt()<<", "<<genBestCombo.second.Eta()<<", "<<genBestCombo.second.Phi()<<"):"
            //         <<" GenMatched: "<<genMatched<<" Number of Jets: "<<NGoodJets<<std::endl;
        }
        tr.registerDerivedVar(ESVarName_+"BestCombo"+myVarSuffix_, BestCombo);        
        tr.registerDerivedVar(ESVarName_+"genBestCombo"+myVarSuffix_, genBestCombo);
        tr.registerDerivedVar(ESVarName_+"MegaJetsTopsGenMatched"+myVarSuffix_, genMatched);
        tr.registerDerivedVar(ESVarName_+"BestComboAvgMass"+myVarSuffix_, static_cast<double>(( BestCombo.first.M() + BestCombo.second.M() )/2));
    }
    
public:
    MakeMVAVariables(const bool verb = false, const std::string& myVarSuffix = "", const std::string& jetColl = "GoodJets_pt30", bool doGenMatch = false, bool printStatus = true, int nTopJets = 7, int nLeptons = 1)
        : verb_(verb)
        , myVarSuffix_(myVarSuffix)
        , doGenMatch_(doGenMatch)
        , nTopJets_(nTopJets)
        , nLeptons_(nLeptons)
    {
        if(printStatus) std::cout<<"Setting up MakeMVAVariables: using the \""<<jetColl<<"\" jet collection"<<std::endl;
        setJetCollection(jetColl);
    }

    void operator()(NTupleReader& tr)
    {
        makeMVAVariables(tr);
    }
};

#endif
