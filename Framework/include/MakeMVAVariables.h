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
        utility::LorentzVector v1;
        utility::LorentzVector v2;
        std::vector<int> jetCombo;
    };

    class TLV
    {
    public:
        utility::LorentzVector tlv;
        double nEF;
        double cEF;
        double nHF;
        double cHF;
        double flavb;
        double flavg;
        double flavc;
        double flavuds;
        double flavq;
        double CSVb;
        double CSVc;
        double CSVudsg;
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

    bool genMatch(const NTupleReader& tr, const std::vector<utility::LorentzVector>& lv_all , const std::vector<int>& jetCombo) const
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
                        double deltaR = utility::DeltaR(*d, lv_all.at(ijet));
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
        else if (jetColl == "GoodJetsPuIdMedium")
        {
            GoodJetsName_ = "GoodJetsPuIdMedium";
            NGoodJetsName_ = "NGoodJetsPuIdMedium";
            GoodLeptonsName_  = "GoodLeptons"   ;
            NGoodLeptonsName_ = "NGoodLeptons"  ;
            MVAJetName_       = "JetPuId"           ;
            MVALeptonName_    = "GoodLeptons"   ;
            ESVarName_        = "PuId_"              ;
        }
    }

    void makeMVAVariables(NTupleReader& tr)
    {
        const auto& Jets                             = tr.getVec<utility::LorentzVector>("Jets"+myVarSuffix_              );
        const auto& Jets_bJetTagDeepFlavourtotb      = tr.getVec<float>("Jets"+myVarSuffix_+"_bJetTagDeepFlavourtotb"     );
        const auto& Jets_bJetTagDeepFlavourprobg     = tr.getVec<float>("Jets"+myVarSuffix_+"_bJetTagDeepFlavourprobg"    );
        const auto& Jets_bJetTagDeepFlavourtotq      = tr.getVec<float>("Jets"+myVarSuffix_+"_bJetTagDeepFlavourtotq"     );
        const auto& Jets_bJetTagDeepFlavourprobc     = tr.getVec<float>("Jets"+myVarSuffix_+"_bJetTagDeepFlavourprobc"    );
        const auto& Jets_bJetTagDeepFlavourprobuds   = tr.getVec<float>("Jets"+myVarSuffix_+"_bJetTagDeepFlavourprobuds"  );
        const auto& Jets_bJetTagDeepCSVtotb          = tr.getVec<float>("Jets"+myVarSuffix_+"_bJetTagDeepCSVtotb"         );
        const auto& Jets_bJetTagDeepCSVprobc         = tr.getVec<float>("Jets"+myVarSuffix_+"_bJetTagDeepCSVprobc"        );
        const auto& Jets_bJetTagDeepCSVprobudsg      = tr.getVec<float>("Jets"+myVarSuffix_+"_bJetTagDeepCSVprobudsg"     );
        const auto& Jets_ptD                         = tr.getVec<float>("Jets"+myVarSuffix_+"_ptD"                        );
        const auto& Jets_axismajor                   = tr.getVec<float>("Jets"+myVarSuffix_+"_axismajor"                  );
        const auto& Jets_axisminor                   = tr.getVec<float>("Jets"+myVarSuffix_+"_axisminor"                  );
        const auto& Jets_multiplicity                = tr.getVec<int>(   "Jets"+myVarSuffix_+"_multiplicity"              );
        const auto& Jets_neutralEmEnergyFraction     = tr.getVec<float>("Jets"+myVarSuffix_+"_neutralEmEnergyFraction"    );
        const auto& Jets_chargedEmEnergyFraction     = tr.getVec<float>("Jets"+myVarSuffix_+"_chargedEmEnergyFraction"    );
        const auto& Jets_neutralHadronEnergyFraction = tr.getVec<float>("Jets"+myVarSuffix_+"_neutralHadronEnergyFraction");
        const auto& Jets_chargedHadronEnergyFraction = tr.getVec<float>("Jets"+myVarSuffix_+"_chargedHadronEnergyFraction");

        const auto& GoodJets           = tr.getVec<bool>(GoodJetsName_+myVarSuffix_                                     );
        const auto& NGoodJets          = tr.getVar<int>(NGoodJetsName_+myVarSuffix_                                     );
        const auto& GoodLeptons        = tr.getVec<std::pair<std::string, utility::LorentzVector>>(GoodLeptonsName_+myVarSuffix_);
        const auto& NGoodLeptons       = tr.getVar<int>(NGoodLeptonsName_+myVarSuffix_                                 );
        const auto& MET                = tr.getVar<float>("MET"                                                        ); 
        const auto& METPhi             = tr.getVar<float>("METPhi"                                                     );
        const auto& HT_Trigger_pt30    = tr.getVar<double>("HT_trigger_pt30"+myVarSuffix_                               );

        // Get the 4-vec for the MET
        utility::LorentzVector lvMET;
        lvMET.SetPt(MET); lvMET.SetEta(0.0); lvMET.SetPhi(METPhi); lvMET.SetE(MET);

        // Sum all jets
        utility::LorentzVector rlv_all, rlv_pt20;
        for(auto jlv : Jets)
        {
            rlv_all += jlv;
        }

        // Boost to the CM frame (only in z)
        double event_beta_z      = rlv_all.Pz() / rlv_all.E();
        utility::BoostVector rec_boost_beta_vec( 0.0, 0.0, -event_beta_z );

        auto& cm_jets = tr.createDerivedVec<math::RThetaPhiVector>(ESVarName_+"cm_jets"+myVarSuffix_);
        std::vector<TLV> Jets_cm;
        std::vector<utility::LorentzVector> Jets_;

        // Boost the GoodJets, Goodleptons, and the MET in the event 
        for(unsigned int j = 0; j < Jets.size(); j++)
        {
            if(!GoodJets[j]) continue;
            utility::LorentzVector jlv = Jets.at(j);            
            Jets_.push_back( jlv );
            utility::LorentzVector jlvcm = utility::Boost(jlv, rec_boost_beta_vec );

            Jets_cm.push_back( {jlvcm, Jets_neutralEmEnergyFraction.at(j), Jets_chargedEmEnergyFraction.at(j), Jets_neutralHadronEnergyFraction.at(j),
                                       Jets_chargedHadronEnergyFraction.at(j), Jets_bJetTagDeepFlavourtotb.at(j), Jets_bJetTagDeepFlavourprobg.at(j),
                                       Jets_bJetTagDeepFlavourprobc.at(j), Jets_bJetTagDeepFlavourprobuds.at(j), Jets_bJetTagDeepFlavourtotq.at(j),
                                       Jets_bJetTagDeepCSVtotb.at(j), Jets_bJetTagDeepCSVprobc.at(j), Jets_bJetTagDeepCSVprobudsg.at(j),
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
        std::sort( Jets_psort.begin(), Jets_psort.end(), utility::compare_pt_TLV<utility::LorentzVector, utility::LorentzVector> );

        auto& Jets_cm_top6 = tr.createDerivedVec<utility::LorentzVector>(ESVarName_+"Jets_cm_top6"+channel_+myVarSuffix_);
        std::vector<double> Jets_cm_top6_flavb, Jets_cm_top6_flavg, Jets_cm_top6_flavc, Jets_cm_top6_flavuds, Jets_cm_top6_flavq;
        std::vector<double> Jets_cm_top6_CSVb, Jets_cm_top6_CSVc, Jets_cm_top6_CSVudsg;
        std::vector<double> Jets_cm_top6_ptD, Jets_cm_top6_axismajor, Jets_cm_top6_axisminor, Jets_cm_top6_multiplicity;
        std::vector<double> Jets_cm_top6_nEF, Jets_cm_top6_cEF, Jets_cm_top6_nHF, Jets_cm_top6_cHF;

        auto& Jets_top6 = tr.createDerivedVec<utility::LorentzVector>(ESVarName_+"Jets_top6"+myVarSuffix_);

        double phiMax = (NGoodJets > 0) ? Jets_cm_psort[0].tlv.Phi() : 0.0;

        utility::LorentzVector combinedNp1thJetTLV;
        utility::LorentzVector combinedNthJetTLV;
        utility::LorentzVector combinedN1thJetTLV;
        for(unsigned int ji=0; ji<cm_jets.size(); ji++ ) 
        {
            utility::LorentzVector Jet_cm_psort = Jets_cm_psort.at(ji).tlv;                
            if(ji == 0)
                Jet_cm_psort.SetPhi(0.0);
            else
                Jet_cm_psort = utility::RotateZ(Jet_cm_psort, -phiMax);                

            if ( ji < nTopJets_ ) 
            {
                cm_jets_top6.push_back( cm_jets_psort.at(ji) ) ;

                Jets_cm_top6.push_back             ( Jet_cm_psort                      );
                Jets_cm_top6_flavb.push_back       ( Jets_cm_psort.at(ji).flavb        );
                Jets_cm_top6_flavg.push_back       ( Jets_cm_psort.at(ji).flavg        );
                Jets_cm_top6_flavc.push_back       ( Jets_cm_psort.at(ji).flavc        );
                Jets_cm_top6_flavuds.push_back     ( Jets_cm_psort.at(ji).flavuds      );
                Jets_cm_top6_flavq.push_back       ( Jets_cm_psort.at(ji).flavq        );
                Jets_cm_top6_CSVb.push_back        ( Jets_cm_psort.at(ji).CSVb         );
                Jets_cm_top6_CSVc.push_back        ( Jets_cm_psort.at(ji).CSVc         );
                Jets_cm_top6_CSVudsg.push_back     ( Jets_cm_psort.at(ji).CSVudsg      );
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

            // Add all jets after and including the nth jet (ji >= n-1)
            if ( ji >= nTopJets_-2 )
            {
                combinedN1thJetTLV += Jet_cm_psort;
            }
            if ( ji >= nTopJets_-1 )
            {
                combinedNthJetTLV += Jet_cm_psort;
            }
            if ( ji >= nTopJets_ )
            {
                combinedNp1thJetTLV += Jet_cm_psort;
            }
        } // ji

        tr.registerDerivedVar("combined" + std::to_string(nTopJets_-1) + "thToLast"+MVAJetName_+"_pt_cm"+channel_+myVarSuffix_,    static_cast<double>(combinedN1thJetTLV.Pt())    );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_-1) + "thToLast"+MVAJetName_+"_ptrHT_cm"+channel_+myVarSuffix_, static_cast<double>(combinedN1thJetTLV.Pt() / HT_Trigger_pt30)    );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_-1) + "thToLast"+MVAJetName_+"_eta_cm"+channel_+myVarSuffix_,   static_cast<double>(combinedN1thJetTLV.Eta())   );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_-1) + "thToLast"+MVAJetName_+"_phi_cm"+channel_+myVarSuffix_,   static_cast<double>(combinedN1thJetTLV.Phi())   );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_-1) + "thToLast"+MVAJetName_+"_m_cm"+channel_+myVarSuffix_,     static_cast<double>(combinedN1thJetTLV.M())     );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_-1) + "thToLast"+MVAJetName_+"_E_cm"+channel_+myVarSuffix_,     static_cast<double>(combinedN1thJetTLV.E())     );

        tr.registerDerivedVar("combined" + std::to_string(nTopJets_) + "thToLast"+MVAJetName_+"_pt_cm"+channel_+myVarSuffix_,    static_cast<double>(combinedNthJetTLV.Pt())    );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_) + "thToLast"+MVAJetName_+"_ptrHT_cm"+channel_+myVarSuffix_, static_cast<double>(combinedNthJetTLV.Pt() / HT_Trigger_pt30)    );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_) + "thToLast"+MVAJetName_+"_eta_cm"+channel_+myVarSuffix_,   static_cast<double>(combinedNthJetTLV.Eta())   );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_) + "thToLast"+MVAJetName_+"_phi_cm"+channel_+myVarSuffix_,   static_cast<double>(combinedNthJetTLV.Phi())   );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_) + "thToLast"+MVAJetName_+"_m_cm"+channel_+myVarSuffix_,     static_cast<double>(combinedNthJetTLV.M())     );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_) + "thToLast"+MVAJetName_+"_E_cm"+channel_+myVarSuffix_,     static_cast<double>(combinedNthJetTLV.E())     );

        tr.registerDerivedVar("combined" + std::to_string(nTopJets_ + 1) + "thToLast"+MVAJetName_+"_pt_cm"+channel_+myVarSuffix_,    static_cast<double>(combinedNp1thJetTLV.Pt())    );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_ + 1) + "thToLast"+MVAJetName_+"_ptrHT_cm"+channel_+myVarSuffix_, static_cast<double>(combinedNp1thJetTLV.Pt() / HT_Trigger_pt30)    );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_ + 1) + "thToLast"+MVAJetName_+"_eta_cm"+channel_+myVarSuffix_,   static_cast<double>(combinedNp1thJetTLV.Eta())   );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_ + 1) + "thToLast"+MVAJetName_+"_phi_cm"+channel_+myVarSuffix_,   static_cast<double>(combinedNp1thJetTLV.Phi())   );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_ + 1) + "thToLast"+MVAJetName_+"_m_cm"+channel_+myVarSuffix_,     static_cast<double>(combinedNp1thJetTLV.M())     );
        tr.registerDerivedVar("combined" + std::to_string(nTopJets_ + 1) + "thToLast"+MVAJetName_+"_E_cm"+channel_+myVarSuffix_,     static_cast<double>(combinedNp1thJetTLV.E())     );

        auto GoodLeptons_cm = std::make_unique<std::vector<utility::LorentzVector>>();
        for(unsigned int ilep = 0; ilep < GoodLeptons.size(); ilep++)
        {
            // Boost and rotate the good leptons into the same frame as the AK4 jets
            GoodLeptons_cm->push_back        ( GoodLeptons[ilep].second    );
            GoodLeptons_cm->at(ilep) = utility::Boost(GoodLeptons_cm->at(ilep), rec_boost_beta_vec    );
            GoodLeptons_cm->at(ilep) = utility::RotateZ(GoodLeptons_cm->at(ilep), -phiMax                     );

        }

        utility::LorentzVector lvMET_cm = lvMET;
       
        // Boost and rotate the MET 4-vector into the same frame as the AK4 jets and good leptons
        lvMET_cm = utility::Boost(lvMET_cm, rec_boost_beta_vec );           
        lvMET_cm = utility::RotateZ(lvMET_cm, -phiMax );           

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
        
        // Get the tops and boost and rotate them !
        auto& topsLV = tr.getVec<utility::LorentzVector>("topsLV"+myVarSuffix_);

        double top1mass = 0.0; double top2mass = 0.0;
        double top1pt = 0.0;   double top2pt = 0.0;
        double top1phi = 0.0;  double top2phi = 0.0;
        double top1eta = 0.0;  double top2eta = 0.0;

        if (topsLV.size() >= 1)
        {
             auto Top1 = topsLV.at(0);
             Top1 = utility::Boost(Top1, rec_boost_beta_vec);
             Top1 = utility::RotateZ(Top1, -phiMax);

             top1pt = Top1.Pt();
             top1eta = Top1.Eta();
             top1phi = Top1.Phi();
             top1mass = Top1.M();

             if (topsLV.size() >= 2)
             {
                 auto Top2 = topsLV.at(1);
                 Top2 = utility::Boost(Top2, rec_boost_beta_vec);
                 Top2 = utility::RotateZ(Top2, -phiMax);

                 top2pt = Top2.Pt();
                 top2eta = Top2.Eta();
                 top2phi = Top2.Phi();
                 top2mass = Top2.M();
             }
        }

        tr.registerDerivedVar("top1_pt_cm"+myVarSuffix_,   static_cast<double>( top1pt ));
        tr.registerDerivedVar("top1_eta_cm"+myVarSuffix_,  static_cast<double>( top1eta ));
        tr.registerDerivedVar("top1_phi_cm"+myVarSuffix_,  static_cast<double>( top1phi ));
        tr.registerDerivedVar("top1_mass_cm"+myVarSuffix_, static_cast<double>( top1mass ));

        tr.registerDerivedVar("top2_pt_cm"+myVarSuffix_,   static_cast<double>( top2pt ));
        tr.registerDerivedVar("top2_eta_cm"+myVarSuffix_,  static_cast<double>( top2eta ));
        tr.registerDerivedVar("top2_phi_cm"+myVarSuffix_,  static_cast<double>( top2phi ));
        tr.registerDerivedVar("top2_mass_cm"+myVarSuffix_, static_cast<double>( top2mass ));

        // Make and get the event shape variables for the 6 highest-P jets in the CM frame
        EventShapeVariables esv_top6( cm_jets_top6 ) ;
        TVectorD eigen_vals_norm_top6 = esv_top6.getEigenValues() ;

        double fwm2_top6    = esv_top6.getFWmoment( 2 ) ;
        double fwm3_top6    = esv_top6.getFWmoment( 3 ) ;
        double fwm4_top6    = esv_top6.getFWmoment( 4 ) ;
        double fwm5_top6    = esv_top6.getFWmoment( 5 ) ;
        double fwm6_top6    = esv_top6.getFWmoment( 6 ) ;
        double fwm7_top6    = esv_top6.getFWmoment( 7 ) ;
        double fwm8_top6    = esv_top6.getFWmoment( 8 ) ;
        double fwm9_top6    = esv_top6.getFWmoment( 9 ) ;
        double fwm10_top6   = esv_top6.getFWmoment( 10) ;
        double jmt_ev0_top6 = eigen_vals_norm_top6[0]   ;
        double jmt_ev1_top6 = eigen_vals_norm_top6[1]   ;
        double jmt_ev2_top6 = eigen_vals_norm_top6[2]   ;

        // Register Variables
        for(unsigned int i = 0; i < nTopJets_; i++)
        {
            tr.registerDerivedVar(MVAJetName_+"_pt_"+std::to_string(i+1)+channel_+myVarSuffix_,          static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).Pt()         : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_ptrHT_"+std::to_string(i+1)+channel_+myVarSuffix_,          static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).Pt() / HT_Trigger_pt30        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_eta_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).Eta()        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_phi_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).Phi()        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_m_"+std::to_string(i+1)+channel_+myVarSuffix_,           static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).M()          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_E_"+std::to_string(i+1)+channel_+myVarSuffix_,           static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6.at(i).E()          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_flavb_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_flavb.at(i)        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_flavg_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_flavg.at(i)        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_flavc_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_flavc.at(i)        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_flavuds_"+std::to_string(i+1)+channel_+myVarSuffix_,     static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_flavuds.at(i)      : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_flavq_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_flavq.at(i)        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_CSVb_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_CSVb.at(i)        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_CSVc_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_CSVc.at(i)        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_CSVudsg_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_CSVudsg.at(i)        : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_nEF_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_nEF.at(i)          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_cEF_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_cEF.at(i)          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_nHF_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_nHF.at(i)          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_cHF_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_cHF.at(i)          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_ptD_"+std::to_string(i+1)+channel_+myVarSuffix_,         static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_ptD.at(i)          : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_axismajor_"+std::to_string(i+1)+channel_+myVarSuffix_,   static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_axismajor.at(i)    : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_axisminor_"+std::to_string(i+1)+channel_+myVarSuffix_,   static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_axisminor.at(i)    : 0.0));
            tr.registerDerivedVar(MVAJetName_+"_multiplicity_"+std::to_string(i+1)+channel_+myVarSuffix_,static_cast<double>( (Jets_cm_top6.size() >= i+1) ? Jets_cm_top6_multiplicity.at(i) : 0.0));
        }

        for(unsigned int i = 0; i < nLeptons_; i++)
        {
            tr.registerDerivedVar(MVALeptonName_+"_pt_"+std::to_string(i+1)+channel_+myVarSuffix_,      static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).Pt()            : 0.0));
            tr.registerDerivedVar(MVALeptonName_+"_eta_"+std::to_string(i+1)+channel_+myVarSuffix_,     static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).Eta()           : 0.0));
            tr.registerDerivedVar(MVALeptonName_+"_phi_"+std::to_string(i+1)+channel_+myVarSuffix_,     static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).Phi()           : 0.0));
            tr.registerDerivedVar(MVALeptonName_+"_m_"+std::to_string(i+1)+channel_+myVarSuffix_,       static_cast<double>( (GoodLeptons_cm->size() >= i+1) ? GoodLeptons_cm->at(i).M()             : 0.0));
        }

        tr.registerDerivedVar(ESVarName_+"lvMET_cm_pt"+channel_+myVarSuffix_,  static_cast<double>( lvMET_cm.Pt() ));
        tr.registerDerivedVar(ESVarName_+"lvMET_cm_eta"+channel_+myVarSuffix_, static_cast<double>( lvMET_cm.Eta()));
        tr.registerDerivedVar(ESVarName_+"lvMET_cm_phi"+channel_+myVarSuffix_, static_cast<double>( lvMET_cm.Phi()));
        tr.registerDerivedVar(ESVarName_+"lvMET_cm_m"+channel_+myVarSuffix_,   static_cast<double>( lvMET_cm.M()  ));
        tr.registerDerivedVar(ESVarName_+"fwm2_top6"+channel_+myVarSuffix_,    static_cast<double>( fwm2_top6     ));
        tr.registerDerivedVar(ESVarName_+"fwm3_top6"+channel_+myVarSuffix_,    static_cast<double>( fwm3_top6     ));
        tr.registerDerivedVar(ESVarName_+"fwm4_top6"+channel_+myVarSuffix_,    static_cast<double>( fwm4_top6     ));
        tr.registerDerivedVar(ESVarName_+"fwm5_top6"+channel_+myVarSuffix_,    static_cast<double>( fwm5_top6     ));
        tr.registerDerivedVar(ESVarName_+"fwm6_top6"+channel_+myVarSuffix_,    static_cast<double>( fwm6_top6     ));
        tr.registerDerivedVar(ESVarName_+"fwm7_top6"+channel_+myVarSuffix_,    static_cast<double>( fwm7_top6     ));
        tr.registerDerivedVar(ESVarName_+"fwm8_top6"+channel_+myVarSuffix_,    static_cast<double>( fwm8_top6     ));
        tr.registerDerivedVar(ESVarName_+"fwm9_top6"+channel_+myVarSuffix_,    static_cast<double>( fwm9_top6     )); 
        tr.registerDerivedVar(ESVarName_+"fwm10_top6"+channel_+myVarSuffix_,   static_cast<double>( fwm10_top6    ));
        tr.registerDerivedVar(ESVarName_+"jmt_ev0_top6"+channel_+myVarSuffix_, static_cast<double>( jmt_ev0_top6  ));
        tr.registerDerivedVar(ESVarName_+"jmt_ev1_top6"+channel_+myVarSuffix_, static_cast<double>( jmt_ev1_top6  ));
        tr.registerDerivedVar(ESVarName_+"jmt_ev2_top6"+channel_+myVarSuffix_, static_cast<double>( jmt_ev2_top6  ));
        tr.registerDerivedVar(ESVarName_+"event_beta_z"+myVarSuffix_,          static_cast<double>( event_beta_z  ));
        tr.registerDerivedVar(ESVarName_+"event_phi_rotate"+myVarSuffix_,      static_cast<double>( phiMax        ));
        tr.registerDerivedVar(ESVarName_+"nMVAJets"+channel_+myVarSuffix_,                          nTopJets_      );
        tr.registerDerivedVar(ESVarName_+"nMVALeptons"+channel_+myVarSuffix_,                       nLeptons_      );

        // Sum jets, leptons, and MET in the CM frame to reco the SUSY particles
        std::pair<utility::LorentzVector, utility::LorentzVector> BestCombo, genBestCombo;
        bool genMatched = false;
        if(NGoodLeptons == 1 && doGenMatch_)
        {
            // Making a vector of all Jets, leptons, and MET
            std::vector<utility::LorentzVector> lv_all;
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
                utility::LorentzVector v1, v2;
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
        const auto& lostCauseEvent = tr.getVar<bool>("lostCauseEvent" + myVarSuffix_);
        const auto& fastMode       = tr.getVar<bool>("fastMode");

        if (!lostCauseEvent or !fastMode)
            makeMVAVariables(tr);
    }
};

#endif
