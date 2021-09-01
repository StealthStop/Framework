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

    class TLV
    {
    public:
        TLorentzVector tlv;
        double nEF;
        double cEF;
        double nHF;
        double cHF;
        double flavb;
        double flavg;
        double flavc;
        double flavuds;
        double flavq;
        double ptD;
        double axismajor;
        double axisminor;
        double multiplicity;
    };

    bool verb_;
    std::string myVarSuffix_;
    bool doGenMatch_;
    unsigned int nTopJets_, nLeptons_;
    std::string channel_;
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
        if (jetColl == "NonIsoMuonJets_pt30")
        {
            GoodJetsName_     = "NonIsoMuonJets_pt30";
            NGoodJetsName_    = "NNonIsoMuonJets_pt30";
            GoodLeptonsName_  = "GoodNonIsoMuons";
            NGoodLeptonsName_ = "NNonIsoMuons";
            MVAJetName_       = "JetNonIsoMuons";
            MVALeptonName_    = "GoodNonIsoMuons";
            ESVarName_        = "NonIsoMuons_";
        }
        else if (jetColl == "NonIsoMuonJets_pt45")
        {
            GoodJetsName_     = "NonIsoMuonJets_pt45";
            NGoodJetsName_    = "NNonIsoMuonJets_pt45";
            GoodLeptonsName_  = "GoodNonIsoMuons";
            NGoodLeptonsName_ = "NNonIsoMuons";
            MVAJetName_       = "JetNonIsoMuons";
            MVALeptonName_    = "GoodNonIsoMuons";
            ESVarName_        = "NonIsoMuons_";
        }    
        else if (jetColl == "GoodJets_pt30")
        {
            GoodJetsName_     = "GoodJets_pt30";
            NGoodJetsName_    = "NGoodJets_pt30";
            GoodLeptonsName_  = "GoodLeptons";
            NGoodLeptonsName_ = "NGoodLeptons";
            MVAJetName_       = "Jet";
            MVALeptonName_    = "GoodLeptons";
            ESVarName_        = "";
        }    
        else if (jetColl == "GoodJets_pt30_GoodLeptons_pt20")
        {
            GoodJetsName_     = "GoodJets_pt30";
            NGoodJetsName_    = "NGoodJets_pt30";
            GoodLeptonsName_  = "GoodLeptons_pt20";
            NGoodLeptonsName_ = "NGoodLeptons_pt20";
            MVAJetName_       = "Jet";
            MVALeptonName_    = "GoodLeptons";
            ESVarName_        = "";
        }        
        else if (jetColl == "GoodJets_pt45")
        {   
            GoodJetsName_     = "GoodJets_pt45" ;
            NGoodJetsName_    = "NGoodJets_pt45";
            GoodLeptonsName_  = "GoodLeptons"   ;
            NGoodLeptonsName_ = "NGoodLeptons"  ;
            MVAJetName_       = "Jet"           ;
            MVALeptonName_    = "GoodLeptons"   ;
            ESVarName_        = ""              ;
        }
    }

    void makeMVAVariables(NTupleReader& tr)
    {
        const auto& Jets                             = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_                       );
        const auto& Jets_bJetTagDeepFlavourtotb      = tr.getVec<double>("Jets"+myVarSuffix_+"_bJetTagDeepFlavourtotb"     );
        const auto& Jets_bJetTagDeepFlavourprobg     = tr.getVec<double>("Jets"+myVarSuffix_+"_bJetTagDeepFlavourprobg"    );
        const auto& Jets_bJetTagDeepFlavourtotq      = tr.getVec<double>("Jets"+myVarSuffix_+"_bJetTagDeepFlavourtotq"     );
        const auto& Jets_bJetTagDeepFlavourprobc     = tr.getVec<double>("Jets"+myVarSuffix_+"_bJetTagDeepFlavourprobc"    );
        const auto& Jets_bJetTagDeepFlavourprobuds   = tr.getVec<double>("Jets"+myVarSuffix_+"_bJetTagDeepFlavourprobuds"  );
        const auto& Jets_ptD                         = tr.getVec<double>("Jets"+myVarSuffix_+"_ptD"                        );
        const auto& Jets_axismajor                   = tr.getVec<double>("Jets"+myVarSuffix_+"_axismajor"                  );
        const auto& Jets_axisminor                   = tr.getVec<double>("Jets"+myVarSuffix_+"_axisminor"                  );
        const auto& Jets_multiplicity                = tr.getVec<int>(   "Jets"+myVarSuffix_+"_multiplicity"               );
        const auto& Jets_neutralEmEnergyFraction     = tr.getVec<double>("Jets"+myVarSuffix_+"_neutralEmEnergyFraction"    );
        const auto& Jets_chargedEmEnergyFraction     = tr.getVec<double>("Jets"+myVarSuffix_+"_chargedEmEnergyFraction"    );
        const auto& Jets_neutralHadronEnergyFraction = tr.getVec<double>("Jets"+myVarSuffix_+"_neutralHadronEnergyFraction");
        const auto& Jets_chargedHadronEnergyFraction = tr.getVec<double>("Jets"+myVarSuffix_+"_chargedHadronEnergyFraction");

        const auto& GoodJets           = tr.getVec<bool>(GoodJetsName_+myVarSuffix_                                     );
        const auto& NGoodJets          = tr.getVar<int>(NGoodJetsName_+myVarSuffix_                                     );
        const auto& GoodLeptons        = tr.getVec<std::pair<std::string, TLorentzVector>>(GoodLeptonsName_+myVarSuffix_);
        const auto& NGoodLeptons       = tr.getVar<int>(NGoodLeptonsName_+myVarSuffix_                                  );
        const auto& MET                = tr.getVar<double>("MET"                                                        ); 
        const auto& METPhi             = tr.getVar<double>("METPhi"                                                     );

        // Get the 4-vec for the MET
        TLorentzVector lvMET;
        lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);

        // Sum all jets
        TLorentzVector rlv_all, rlv_pt20;
        for(auto jlv : Jets)
        {
            rlv_all += jlv;
            if (jlv.Pt() > 20.) rlv_pt20 += jlv;
        }

        // Boost to the CM frame (only in z)
        double event_beta_z      = rlv_all.Pz() / rlv_all.E();
        double event_beta_z_pt20 = rlv_pt20.Pz() / rlv_pt20.E();
        TVector3 rec_boost_beta_vec( 0.0, 0.0, -event_beta_z );

        auto& cm_jets = tr.createDerivedVec<math::RThetaPhiVector>(ESVarName_+"cm_jets"+myVarSuffix_);
        std::vector<TLV> Jets_cm;
        std::vector<TLorentzVector> Jets_;

        // Boost the GoodJets, Goodleptons, and the MET in the event 
        for(unsigned int j = 0; j < Jets.size(); j++)
        {
            if(!GoodJets[j]) continue;
            TLorentzVector jlvcm = Jets.at(j);            
            Jets_.push_back( jlvcm );
            
            jlvcm.Boost( rec_boost_beta_vec );
            Jets_cm.push_back( {jlvcm, Jets_neutralEmEnergyFraction.at(j), Jets_chargedEmEnergyFraction.at(j), Jets_neutralHadronEnergyFraction.at(j),
                                       Jets_chargedHadronEnergyFraction.at(j), Jets_bJetTagDeepFlavourtotb.at(j), Jets_bJetTagDeepFlavourprobg.at(j),
                                       Jets_bJetTagDeepFlavourprobc.at(j), Jets_bJetTagDeepFlavourprobuds.at(j), Jets_bJetTagDeepFlavourtotq.at(j),
                                       Jets_ptD.at(j), Jets_axismajor.at(j), Jets_axisminor.at(j), static_cast<double>(Jets_multiplicity.at(j))} );

            math::RThetaPhiVector cmvec( jlvcm.P(), jlvcm.Theta(), jlvcm.Phi() );
            cm_jets.push_back( cmvec );
        }

        // Try using only the 7 highest-P jets in the CM frame in the event shape vars
        // First, need to make a new input vector of jets containing only those jets
        auto cm_jets_psort = cm_jets ;
        std::sort( cm_jets_psort.begin(), cm_jets_psort.end(), utility::compare_p ) ;
        std::vector<math::RThetaPhiVector> cm_jets_top6 ;

        auto Jets_cm_psort = Jets_cm;
        auto Jets_psort    = Jets_;
        std::sort( Jets_cm_psort.begin(), Jets_cm_psort.end(), [](TLV v1, TLV v2){return v1.tlv.P() > v2.tlv.P();} );
        std::sort( Jets_psort.begin(), Jets_psort.end(), utility::compare_pt_TLV );

        auto& Jets_cm_top6 = tr.createDerivedVec<TLorentzVector>(ESVarName_+"Jets_cm_top6"+channel_+myVarSuffix_);
        std::vector<double> Jets_cm_top6_flavb, Jets_cm_top6_flavg, Jets_cm_top6_flavc, Jets_cm_top6_flavuds, Jets_cm_top6_flavq;
        std::vector<double> Jets_cm_top6_ptD, Jets_cm_top6_axismajor, Jets_cm_top6_axisminor, Jets_cm_top6_multiplicity;
        std::vector<double> Jets_cm_top6_nEF, Jets_cm_top6_cEF, Jets_cm_top6_nHF, Jets_cm_top6_cHF;

        auto& Jets_top6 = tr.createDerivedVec<TLorentzVector>(ESVarName_+"Jets_top6"+myVarSuffix_);

        double phiMax = (NGoodJets > 0) ? Jets_cm_psort[0].tlv.Phi() : 0.0;

        for(unsigned int ji=0; ji<cm_jets.size(); ji++ ) 
        {
            if ( ji < nTopJets_ ) 
            {
                cm_jets_top6.push_back( cm_jets_psort.at(ji) ) ;

                TLorentzVector Jet_cm_psort = Jets_cm_psort.at(ji).tlv;                
                if(ji == 0)
                    Jet_cm_psort.SetPhi(0.0);
                else
                    Jet_cm_psort.RotateZ(-phiMax);                
                
                Jets_cm_top6.push_back             ( Jet_cm_psort                      );
                Jets_cm_top6_flavb.push_back       ( Jets_cm_psort.at(ji).flavb        );
                Jets_cm_top6_flavg.push_back       ( Jets_cm_psort.at(ji).flavg        );
                Jets_cm_top6_flavc.push_back       ( Jets_cm_psort.at(ji).flavc        );
                Jets_cm_top6_flavuds.push_back     ( Jets_cm_psort.at(ji).flavuds      );
                Jets_cm_top6_flavq.push_back       ( Jets_cm_psort.at(ji).flavq        );
                Jets_cm_top6_nEF.push_back         ( Jets_cm_psort.at(ji).nEF          );
                Jets_cm_top6_cEF.push_back         ( Jets_cm_psort.at(ji).cEF          );
                Jets_cm_top6_nHF.push_back         ( Jets_cm_psort.at(ji).nHF          );
                Jets_cm_top6_cHF.push_back         ( Jets_cm_psort.at(ji).cHF          );
                Jets_cm_top6_ptD.push_back         ( Jets_cm_psort.at(ji).ptD          );
                Jets_cm_top6_axismajor.push_back   ( Jets_cm_psort.at(ji).axismajor    );
                Jets_cm_top6_axisminor.push_back   ( Jets_cm_psort.at(ji).axisminor    );
                Jets_cm_top6_multiplicity.push_back( Jets_cm_psort.at(ji).multiplicity );
                Jets_top6.push_back                ( Jets_psort.at(ji)                 );
            }
        } // ji

        auto GoodLeptons_cm = std::make_unique<std::vector<TLorentzVector>>();
        for(unsigned int ilep = 0; ilep < GoodLeptons.size(); ilep++)
        {
            // Boost and rotate the good leptons into the same frame as the AK4 jets
            GoodLeptons_cm->push_back        ( GoodLeptons[ilep].second    );
            GoodLeptons_cm->at(ilep).Boost   ( rec_boost_beta_vec          );
            GoodLeptons_cm->at(ilep).RotateZ ( -phiMax                     );

        }

        TLorentzVector lvMET_cm = lvMET;
       
        // Boost and rotate the MET 4-vector into the same frame as the AK4 jets and good leptons
        lvMET_cm.Boost( rec_boost_beta_vec );           
        lvMET_cm.RotateZ( -phiMax );           

        if( verb_ ) 
        {
            printf("\n\n Unsorted and sorted CM jet lists.\n") ;
            for ( unsigned int ji=0; ji<cm_jets.size(); ji++ ) 
            {
                printf("  %2d :  (%7.1f, %7.3f, %7.3f) | (%7.1f, %7.3f, %7.3f)\n", ji,
                       Jets_cm.at(ji).tlv.P(), Jets_cm.at(ji).tlv.Theta(), Jets_cm.at(ji).tlv.Phi(),
                       Jets_cm_top6.at(ji).P(), Jets_cm_top6.at(ji).Theta(), Jets_cm_top6.at(ji).Phi() ) ;
            } // ji
            printf("\n\n") ;
        }
        
        // Make and get the event shape variables for the 6 highest-P jets in the CM frame
        EventShapeVariables esv_top6( cm_jets_top6 ) ;
        TVectorD eigen_vals_norm_top6 = esv_top6.getEigenValues() ;

        double fwm2_top6    = esv_top6.getFWmoment( 2 )  ;
        double fwm3_top6    = esv_top6.getFWmoment( 3 )  ;
        double fwm4_top6    = esv_top6.getFWmoment( 4 )  ;
        double fwm5_top6    = esv_top6.getFWmoment( 5 )  ;
        double fwm6_top6    = esv_top6.getFWmoment( 6 )  ;
        double fwm7_top6    = esv_top6.getFWmoment( 7 )  ;
        double fwm8_top6    = esv_top6.getFWmoment( 8 )  ;
        double fwm9_top6    = esv_top6.getFWmoment( 9 )  ;
        double fwm10_top6   = esv_top6.getFWmoment( 10 ) ;
        double jmt_ev0_top6 = eigen_vals_norm_top6[0]    ;
        double jmt_ev1_top6 = eigen_vals_norm_top6[1]    ;
        double jmt_ev2_top6 = eigen_vals_norm_top6[2]    ;

        // AK8 jet variables
        const auto& GoodJetsAK8        = tr.getVec<bool>("GoodJetsAK8"+myVarSuffix_                              );
        const auto& JetsAK8            = tr.getVec<TLorentzVector>("JetsAK8"+myVarSuffix_                        );
        const auto& Tau1               = tr.getVec<double>("JetsAK8"+myVarSuffix_+"_NsubjettinessTau1"           );
        const auto& Tau2               = tr.getVec<double>("JetsAK8"+myVarSuffix_+"_NsubjettinessTau2"           );
        const auto& Tau3               = tr.getVec<double>("JetsAK8"+myVarSuffix_+"_NsubjettinessTau3"           );
        const auto& softDropMass       = tr.getVec<double>("JetsAK8"+myVarSuffix_+"_softDropMass"                );
        const auto& prunedMass         = tr.getVec<double>("JetsAK8"+myVarSuffix_+"_prunedMass"                  );
        const auto& axismajor_AK8      = tr.getVec<double>("JetsAK8"+myVarSuffix_+"_axismajor"                   );
        const auto& axisminor_AK8      = tr.getVec<double>("JetsAK8"+myVarSuffix_+"_axisminor"                   );
        const auto& subjets            = tr.getVec<std::vector<TLorentzVector>>("JetsAK8"+myVarSuffix_+"_subjets");
        const auto& tDiscriminator_AK8 = tr.getVec<double>("JetsAK8"+myVarSuffix_+"_tDiscriminatorDeep"          );
        const auto& wDiscriminator_AK8 = tr.getVec<double>("JetsAK8"+myVarSuffix_+"_wDiscriminatorDeep"          );
        const auto& hDiscriminator_AK8 = tr.getVec<double>("JetsAK8"+myVarSuffix_+"_hDiscriminatorDeep"          );
        const auto& multiplicity_AK8   = tr.getVec<int>("JetsAK8"+myVarSuffix_+"_multiplicity"                   );
        
        std::vector<TLorentzVector> JetsAK8_TLV_cm;
        std::vector<double> JetsAK8_SDM, JetsAK8_Pruned, JetsAK8_Tau1, JetsAK8_Tau2, JetsAK8_Tau3, JetsAK8_axismajor, JetsAK8_axisminor, JetsAK8_tDiscriminator, JetsAK8_wDiscriminator, JetsAK8_hDiscriminator; 
        std::vector<int> JetsAK8_nsubjets, JetsAK8_multiplicity;
        for (unsigned int j = 0; j < JetsAK8.size(); j++)
        {
            TLorentzVector ijet = JetsAK8.at(j);

            // Boost and rotate the AK8 jets to same frame as the AK4 jets, MET, and leptons
            ijet.Boost( rec_boost_beta_vec );           
            ijet.RotateZ( -phiMax );           
            if (GoodJetsAK8.at(j))
            {
                JetsAK8_TLV_cm.push_back        ( ijet                     );
                JetsAK8_SDM.push_back           ( softDropMass.at(j)       );
                JetsAK8_Pruned.push_back        ( prunedMass.at(j)         );
                JetsAK8_Tau1.push_back          ( Tau1.at(j)               );
                JetsAK8_Tau2.push_back          ( Tau2.at(j)               );
                JetsAK8_Tau3.push_back          ( Tau3.at(j)               );
                JetsAK8_axismajor.push_back     ( axismajor_AK8.at(j)      );
                JetsAK8_axisminor.push_back     ( axisminor_AK8.at(j)      );
                JetsAK8_nsubjets.push_back      ( subjets.at(j).size()     );
                JetsAK8_tDiscriminator.push_back( tDiscriminator_AK8.at(j) );
                JetsAK8_wDiscriminator.push_back( wDiscriminator_AK8.at(j) );
                JetsAK8_hDiscriminator.push_back( hDiscriminator_AK8.at(j) );
                JetsAK8_multiplicity.push_back  ( multiplicity_AK8.at(j)   );
            }
        }

        // Need to zip up vectors so they are sorted simultaneously
        std::vector<std::tuple<TLorentzVector, double, double, double, double, double, double, double, int, double, double, double, int>> zipped_JetsAK8;
        for (unsigned int j = 0; j < JetsAK8_TLV_cm.size(); j++)
        {
            zipped_JetsAK8.push_back(std::make_tuple(JetsAK8_TLV_cm.at(j), JetsAK8_SDM.at(j), JetsAK8_Pruned.at(j), JetsAK8_Tau1.at(j), JetsAK8_Tau2.at(j), JetsAK8_Tau3.at(j), 
                                                     JetsAK8_axismajor.at(j), JetsAK8_axisminor.at(j), JetsAK8_nsubjets.at(j), 
                                                     JetsAK8_tDiscriminator.at(j), JetsAK8_wDiscriminator.at(j), JetsAK8_hDiscriminator.at(j), JetsAK8_multiplicity.at(j)));
        }

        // Now unzip
        std::sort(std::begin(zipped_JetsAK8), std::end(zipped_JetsAK8), [&](const auto& a, const auto& b){return std::get<0>(a).M() > std::get<0>(b).M();});

        std::vector<TLorentzVector> JetsAK8_sorted_TLV_cm;
        std::vector<double> JetsAK8_sorted_SDM, JetsAK8_sorted_Pruned, JetsAK8_sorted_Tau1, JetsAK8_sorted_Tau2, JetsAK8_sorted_Tau3;
        std::vector<double> JetsAK8_sorted_axismajor, JetsAK8_sorted_axisminor, JetsAK8_sorted_tDiscriminator, JetsAK8_sorted_wDiscriminator, JetsAK8_sorted_hDiscriminator;
        std::vector<int> JetsAK8_sorted_nsubjets, JetsAK8_sorted_multiplicity;
        for (unsigned int j = 0; j < zipped_JetsAK8.size(); j++)
        {
            JetsAK8_sorted_TLV_cm.push_back        ( std::get<0>(zipped_JetsAK8.at(j))  );
            JetsAK8_sorted_SDM.push_back           ( std::get<1>(zipped_JetsAK8.at(j))  );
            JetsAK8_sorted_Pruned.push_back        ( std::get<2>(zipped_JetsAK8.at(j))  );
            JetsAK8_sorted_Tau1.push_back          ( std::get<3>(zipped_JetsAK8.at(j))  );
            JetsAK8_sorted_Tau2.push_back          ( std::get<4>(zipped_JetsAK8.at(j))  );
            JetsAK8_sorted_Tau3.push_back          ( std::get<5>(zipped_JetsAK8.at(j))  );
            JetsAK8_sorted_axisminor.push_back     ( std::get<6>(zipped_JetsAK8.at(j))  );
            JetsAK8_sorted_axismajor.push_back     ( std::get<7>(zipped_JetsAK8.at(j))  );
            JetsAK8_sorted_nsubjets.push_back      ( std::get<8>(zipped_JetsAK8.at(j))  );
            JetsAK8_sorted_tDiscriminator.push_back( std::get<9>(zipped_JetsAK8.at(j))  );
            JetsAK8_sorted_wDiscriminator.push_back( std::get<10>(zipped_JetsAK8.at(j)) );
            JetsAK8_sorted_hDiscriminator.push_back( std::get<11>(zipped_JetsAK8.at(j)) );
            JetsAK8_sorted_multiplicity.push_back  ( std::get<12>(zipped_JetsAK8.at(j)) );
        }

        double phiMaxAK8 = (JetsAK8_sorted_TLV_cm.size() > 0) ? JetsAK8_sorted_TLV_cm.at(0).Phi() : 0.0;
        for(unsigned int j = 0; j < JetsAK8_sorted_TLV_cm.size(); j++)
        {
            if(j == 0)
                JetsAK8_sorted_TLV_cm.at(j).SetPhi(0.0);
            else
                JetsAK8_sorted_TLV_cm.at(j).RotateZ(-phiMaxAK8);
        }

        // Register Variables
        for(unsigned int i = 0; i < nTopJets_; i++)
        {
            tr.registerDerivedVar(MVAJetName_+"_pt_"+std::to_string(i+1)+channel_+myVarSuffix_,          static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).Pt()         : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_eta_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).Eta()        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_phi_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).Phi()        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_m_"+std::to_string(i+1)+channel_+myVarSuffix_,           static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).M()          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_flavb_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_flavb.at(i)        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_flavg_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_flavg.at(i)        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_flavc_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_flavc.at(i)        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_flavuds_"+std::to_string(i+1)+channel_+myVarSuffix_,     static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_flavuds.at(i)      : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_flavq_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_flavq.at(i)        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_nEF_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_nEF.at(i)          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_cEF_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_cEF.at(i)          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_nHF_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_nHF.at(i)          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_cHF_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_cHF.at(i)          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_ptD_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_ptD.at(i)          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_axismajor_"+std::to_string(i+1)+channel_+myVarSuffix_,   static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_axismajor.at(i)    : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_axisminor_"+std::to_string(i+1)+channel_+myVarSuffix_,   static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_axisminor.at(i)    : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_multiplicity_"+std::to_string(i+1)+channel_+myVarSuffix_,static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_multiplicity.at(i) : 0.0));
        }

        for(unsigned int i = 0; i < 5; i++)
        {
            tr.registerDerivedVar(MVAJetName_+"sAK8_pt_"+std::to_string(i+1)+channel_+myVarSuffix_,             static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_TLV_cm.at(i).Pt()    : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_eta_"+std::to_string(i+1)+channel_+myVarSuffix_,            static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_TLV_cm.at(i).Eta()   : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_phi_"+std::to_string(i+1)+channel_+myVarSuffix_,            static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_TLV_cm.at(i).Phi()   : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_m_"+std::to_string(i+1)+channel_+myVarSuffix_,              static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_TLV_cm.at(i).M()     : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_SDM_"+std::to_string(i+1)+channel_+myVarSuffix_,            static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_SDM.at(i)            : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_Pruned_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_Pruned.at(i)         : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_Tau1_"+std::to_string(i+1)+channel_+myVarSuffix_,           static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_Tau1.at(i)           : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_Tau2_"+std::to_string(i+1)+channel_+myVarSuffix_,           static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_Tau2.at(i)           : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_Tau3_"+std::to_string(i+1)+channel_+myVarSuffix_,           static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_Tau3.at(i)           : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_axismajor_"+std::to_string(i+1)+channel_+myVarSuffix_,      static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_axismajor.at(i)      : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_axisminor_"+std::to_string(i+1)+channel_+myVarSuffix_,      static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_axisminor.at(i)      : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_nsubjets_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_nsubjets.at(i)       : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_tDiscriminator_"+std::to_string(i+1)+channel_+myVarSuffix_, static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_tDiscriminator.at(i) : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_wDiscriminator_"+std::to_string(i+1)+channel_+myVarSuffix_, static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_wDiscriminator.at(i) : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_hDiscriminator_"+std::to_string(i+1)+channel_+myVarSuffix_, static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_hDiscriminator.at(i) : 0.0));
            tr.registerDerivedVar(MVAJetName_+"sAK8_multiplicity_"+std::to_string(i+1)+channel_+myVarSuffix_,   static_cast<double>( (JetsAK8_sorted_TLV_cm.size() >= i+1) ? JetsAK8_sorted_multiplicity.at(i)   : 0.0));
        }

        for(unsigned int i = 0; i < nLeptons_; i++)
        {
            tr.registerDerivedVar(MVALeptonName_+"_pt_"+std::to_string(i+1)+channel_+myVarSuffix_,      static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).Pt()            : 0.0));
            tr.registerDerivedVar(MVALeptonName_+"_eta_"+std::to_string(i+1)+channel_+myVarSuffix_,     static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).Eta()           : 0.0));
            tr.registerDerivedVar(MVALeptonName_+"_phi_"+std::to_string(i+1)+channel_+myVarSuffix_,     static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).Phi()           : 0.0));
            tr.registerDerivedVar(MVALeptonName_+"_m_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).M()             : 0.0));
        }

        tr.registerDerivedVar(ESVarName_+"lvMET_cm_pt"+channel_+myVarSuffix_,              static_cast<double>( lvMET_cm.Pt() ));
        tr.registerDerivedVar(ESVarName_+"lvMET_cm_eta"+channel_+myVarSuffix_,             static_cast<double>( lvMET_cm.Eta()));
        tr.registerDerivedVar(ESVarName_+"lvMET_cm_phi"+channel_+myVarSuffix_,             static_cast<double>( lvMET_cm.Phi()));
        tr.registerDerivedVar(ESVarName_+"lvMET_cm_m"+channel_+myVarSuffix_,               static_cast<double>( lvMET_cm.M()  ));
        tr.registerDerivedVar(ESVarName_+"fwm2_top6"+channel_+myVarSuffix_,       fwm2_top6                           );
        tr.registerDerivedVar(ESVarName_+"fwm3_top6"+channel_+myVarSuffix_,       fwm3_top6                           );
        tr.registerDerivedVar(ESVarName_+"fwm4_top6"+channel_+myVarSuffix_,       fwm4_top6                           );
        tr.registerDerivedVar(ESVarName_+"fwm5_top6"+channel_+myVarSuffix_,       fwm5_top6                           );
        tr.registerDerivedVar(ESVarName_+"fwm6_top6"+channel_+myVarSuffix_,       fwm6_top6                           );
        tr.registerDerivedVar(ESVarName_+"fwm7_top6"+channel_+myVarSuffix_,       fwm7_top6                           );
        tr.registerDerivedVar(ESVarName_+"fwm8_top6"+channel_+myVarSuffix_,       fwm8_top6                           );
        tr.registerDerivedVar(ESVarName_+"fwm9_top6"+channel_+myVarSuffix_,       fwm9_top6                           );
        tr.registerDerivedVar(ESVarName_+"fwm10_top6"+channel_+myVarSuffix_,      fwm10_top6                          );
        tr.registerDerivedVar(ESVarName_+"jmt_ev0_top6"+channel_+myVarSuffix_,    jmt_ev0_top6                        );
        tr.registerDerivedVar(ESVarName_+"jmt_ev1_top6"+channel_+myVarSuffix_,    jmt_ev1_top6                        );
        tr.registerDerivedVar(ESVarName_+"jmt_ev2_top6"+channel_+myVarSuffix_,    jmt_ev2_top6                        );
        tr.registerDerivedVar(ESVarName_+"event_beta_z"+myVarSuffix_,             event_beta_z                        );
        tr.registerDerivedVar(ESVarName_+"event_beta_z_pt20"+myVarSuffix_,        event_beta_z_pt20                   );
        tr.registerDerivedVar(ESVarName_+"event_phi_rotate"+myVarSuffix_,         phiMax                              );
        tr.registerDerivedVar(ESVarName_+"nMVAJets"+channel_+myVarSuffix_,        nTopJets_                           );

        // Sum jets, leptons, and MET in the CM frame to reco the SUSY particles
        std::pair<TLorentzVector, TLorentzVector> BestCombo, genBestCombo;
        bool genMatched = false;
        if(NGoodLeptons == 1 && doGenMatch_)
        {
            // Making a vector of all Jets, leptons, and MET
            std::vector<TLorentzVector> lv_all;
            for(unsigned int j = 0; j < Jets.size(); j++)
            {
                if(!GoodJets[j]) continue;
                lv_all.push_back( Jets.at(j) );
            }
            // Adding the lepton and MET for 1 lepton selection
            lv_all.push_back( GoodLeptons[0].second + lvMET );

            // Form all possible combos of summing lv into two lv
            int NAll   = lv_all.size();
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
                combinedLV.push_back( {v1, v2, jetCombo} );
            }

            // Find the best combo of lv pair, looks more like the pair production
            double massDiff = 999999999999;
            std::vector<int> BestJetCombo;
            bool matched = false;
            for(const auto& cLV : combinedLV)
            {
                double mD = abs( cLV.v1.M() - cLV.v2.M() );
                if(mD < massDiff)
                {
                    massDiff     = mD;
                    BestCombo    = std::make_pair(cLV.v1, cLV.v2);
                    BestJetCombo = cLV.jetCombo;
                }

                bool m = genMatch(tr, lv_all, cLV.jetCombo);
                if(m && !matched)
                {
                    genBestCombo = std::make_pair(cLV.v1, cLV.v2);
                    matched      = true;
                }
            }
            genMatched = genMatch(tr, lv_all, BestJetCombo);
            if(!matched) genBestCombo = BestCombo;
        }
        tr.registerDerivedVar(ESVarName_+"BestCombo"+myVarSuffix_,              BestCombo                                                            );        
        tr.registerDerivedVar(ESVarName_+"genBestCombo"+myVarSuffix_,           genBestCombo                                                         );
        tr.registerDerivedVar(ESVarName_+"MegaJetsTopsGenMatched"+myVarSuffix_, genMatched                                                           );
        tr.registerDerivedVar(ESVarName_+"BestComboAvgMass"+myVarSuffix_,       static_cast<double>(( BestCombo.first.M() + BestCombo.second.M() )/2));
    }
    
public:
    MakeMVAVariables(const bool verb = false, const std::string& myVarSuffix = "", const std::string& jetColl = "GoodJets_pt30", 
                     bool doGenMatch = false, bool printStatus = true, int nTopJets = 12, int nLeptons = 2, const std::string& channel = "")
        : verb_(verb)
        , myVarSuffix_(myVarSuffix)
        , doGenMatch_(doGenMatch)
        , nTopJets_(nTopJets)
        , nLeptons_(nLeptons)
        , channel_(channel)
    {
        if(printStatus) std::cout<<"Setting up MakeMVAVariables: using the \""<<jetColl+myVarSuffix<<"\" jet collection"<<std::endl;
        setJetCollection(jetColl);
    }

    void operator()(NTupleReader& tr)
    {
        makeMVAVariables(tr);
    }
};

#endif
