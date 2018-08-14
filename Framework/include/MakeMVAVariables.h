#ifndef MAKEMVAVARIABLES_H
#define MAKEMVAVARIABLES_H 

#include "Framework/Framework/src/EventShapeVariables.cc"
#include "Framework/Framework/src/get_cmframe_jets.c"

class MakeMVAVariables
{
private:
    bool verb_;
    std::string myVarSuffix_;

    std::vector<int> decToBinary(int n, int max)
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

    void makeMVAVariables(NTupleReader& tr)
    {
        const auto& Jets = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
        const auto& GoodLeptons = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        const auto& NGoodLeptons = tr.getVar<int>("NGoodLeptons");
        const auto& MET = tr.getVar<double>("MET"); 
        const auto& METPhi = tr.getVar<double>("METPhi");

        //--- Get the 4-vec for the MET
        TLorentzVector lvMET;
        lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);

        //--- Sum all jets, leptons, and MET together
        TLorentzVector rlv_all;
        for(auto jlv : Jets) rlv_all += jlv;
        for(auto pair : GoodLeptons) rlv_all += pair.second;
        rlv_all += lvMET;

        //--- Fill vector of jet momenta in CM frame.
        //    Boost to the CM frame
        //TLorentzVector jetSum = rlv_all;
        //TVector3 cm_boost_beta_vec( -jetSum.BoostVector() );

        //    Boost to the CM frame (only in z)
        double event_beta_z = -9.0;
        double reco_jets_beta = rlv_all.Pz() / rlv_all.E();
        event_beta_z = reco_jets_beta;
        TVector3 rec_boost_beta_vec( 0.0, 0.0, -reco_jets_beta );
        auto* cm_jets = new std::vector<math::RThetaPhiVector>();
        auto* Jets_cm = new std::vector<TLorentzVector>();
 
        for(auto jlvcm : Jets)
        {
            jlvcm.Boost( rec_boost_beta_vec );
            Jets_cm->push_back( jlvcm );

            math::RThetaPhiVector cmvec( jlvcm.P(), jlvcm.Theta(), jlvcm.Phi() );
            cm_jets->push_back( cmvec );
        }

        //--- Try using only the 6 highest-P jets in the CM frame in the event shape vars.
        //    First, need to make a new input vector of jets containing only those jets.
        auto cm_jets_psort = *cm_jets ;
        std::sort( cm_jets_psort.begin(), cm_jets_psort.end(), compare_p ) ;
        std::vector<math::RThetaPhiVector> cm_jets_top6 ;

        auto Jets_cm_psort = *Jets_cm;
        std::sort( Jets_cm_psort.begin(), Jets_cm_psort.end(), [](TLorentzVector v1, TLorentzVector v2){return v1.P() > v2.P();} );
        auto* Jets_cm_top6 = new std::vector<TLorentzVector>();

        for( unsigned int ji=0; ji<cm_jets->size(); ji++ ) 
        {
            if ( ji < 6 ) 
            {
                cm_jets_top6.push_back( cm_jets_psort.at(ji) ) ;
                Jets_cm_top6->push_back( Jets_cm_psort.at(ji) ) ;
            }
        } // ji

        if ( verb_ ) 
        {
            printf("\n\n Unsorted and sorted CM jet lists.\n") ;
            for ( unsigned int ji=0; ji<cm_jets->size(); ji++ ) 
            {
                printf("  %2d :  (%7.1f, %7.3f, %7.3f) | (%7.1f, %7.3f, %7.3f)\n", ji,
                       cm_jets->at(ji).R(), cm_jets->at(ji).Theta(), cm_jets->at(ji).Phi(),
                       cm_jets_psort.at(ji).R(), cm_jets_psort.at(ji).Theta(), cm_jets_psort.at(ji).Phi() ) ;
                
                printf("  %2d :  (%7.1f, %7.3f, %7.3f) | (%7.1f, %7.3f, %7.3f)\n", ji,
                       Jets_cm->at(ji).P(), Jets_cm->at(ji).Theta(), Jets_cm->at(ji).Phi(),
                       Jets_cm_top6->at(ji).P(), Jets_cm_top6->at(ji).Theta(), Jets_cm_top6->at(ji).Phi() ) ;
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
        tr.registerDerivedVec("cm_jets", cm_jets);
        tr.registerDerivedVec("Jets_cm", Jets_cm);
        tr.registerDerivedVec("Jets_cm_top6", Jets_cm_top6);
        tr.registerDerivedVar("fwm2_top6", fwm2_top6);
        tr.registerDerivedVar("fwm3_top6", fwm3_top6);
        tr.registerDerivedVar("fwm4_top6", fwm4_top6);
        tr.registerDerivedVar("fwm5_top6", fwm5_top6);
        tr.registerDerivedVar("fwm6_top6", fwm6_top6);
        tr.registerDerivedVar("fwm7_top6", fwm7_top6);
        tr.registerDerivedVar("fwm8_top6", fwm8_top6);
        tr.registerDerivedVar("fwm9_top6", fwm9_top6);
        tr.registerDerivedVar("fwm10_top6", fwm10_top6);
        tr.registerDerivedVar("jmt_ev0_top6", jmt_ev0_top6);
        tr.registerDerivedVar("jmt_ev1_top6", jmt_ev1_top6);
        tr.registerDerivedVar("jmt_ev2_top6", jmt_ev2_top6);
        tr.registerDerivedVar("event_beta_z", event_beta_z);

        // Sum jets, leptons, and MET in the CM frame to reco the SUSY particles
        if(NGoodLeptons == 1)
        {
            //--- Boost the leptons and MET in the event 
            std::vector<TLorentzVector> GoodLeptons_cm;
            for(auto pair : GoodLeptons)
            {
                pair.second.Boost( rec_boost_beta_vec );
                GoodLeptons_cm.push_back( pair.second );
            }
            TLorentzVector lvMET_cm = lvMET;
            lvMET_cm.Boost( rec_boost_beta_vec );
            
            //--- Making a vector of all Jets, leptons, and MET
            std::vector<TLorentzVector> lv_all_cm = *Jets_cm;
            //adding the lepton and MET for 1 lepton selection
            lv_all_cm.push_back( GoodLeptons_cm[0] + lvMET_cm );
            //lv_all_cm.insert( lv_all_cm.end(), GoodLeptons_cm.begin(), GoodLeptons_cm.end() );
            //lv_all_cm.push_back(lvMET_cm);

            //--- Form all possible combos of summing lv into two lv
            //std::cout<<"Jets size: "<<Jets_cm->size()<<" Total lv in event:"<<lv_all_cm.size()<<std::endl;
            int NJets_cm = lv_all_cm.size();
            int maxDec = 0;
            std::vector<std::pair<TLorentzVector, TLorentzVector>> combinedJets;
            for(int i = 1; i <= NJets_cm - 1; i++) maxDec += pow(2,NJets_cm - i);
            for(int i = 1; i <= maxDec; i++) 
            {                
                std::vector<int> jetCombo = decToBinary( i, NJets_cm);
                TLorentzVector v1, v2;
                int tempCount1 = 0, tempCount2 = 0;
                for(int ijet = 0; ijet < lv_all_cm.size(); ijet++)
                {
                    jetCombo[ijet];
                    if(jetCombo[ijet] == 0) 
                    {
                        v1 += lv_all_cm.at(ijet);
                        tempCount1++;
                    }
                    else if(jetCombo[ijet] == 1)
                    {
                        v2 += lv_all_cm.at(ijet);
                        tempCount2++;
                    }
                }
                //for (int j = lv_all_cm.size() - 1; j >= 0; j--)
                //    std::cout << jetCombo[j];                    
                //std::cout<<" "<<tempCount1<<" "<<tempCount2<<" ("<<v1.M()<<", "<<v1.Pt()<<", "<<v1.Eta()<<", "<<v1.Phi()<<") ("
                //         <<v2.M()<<", "<<v2.Pt()<<", "<<v2.Eta()<<", "<<v2.Phi()<<")"<<std::endl;;
                combinedJets.push_back( std::make_pair(v1, v2) );
            }
            //std::cout<<combinedJets.size()<<std::endl;

            //--- Find the best combo of lv pair, looks more like the pair production
            std::pair<TLorentzVector, TLorentzVector> BestComboJets;
            double massDiff = 999999999999;
            for(const auto& pair : combinedJets)
            {
                double mD = abs( pair.first.M() - pair.second.M() );
                if(mD < massDiff)
                {
                    massDiff = mD;
                    BestComboJets = pair;
                }
            }
            //std::cout<<"Best mass diff Jets_cm: ("<<BestComboJets.first.M()<<", "<<BestComboJets.first.Pt()<<", "<<BestComboJets.first.Eta()<<", "<<BestComboJets.first.Phi()<<") ("
            //         <<BestComboJets.second.M()<<", "<<BestComboJets.second.Pt()<<", "<<BestComboJets.second.Eta()<<", "<<BestComboJets.second.Phi()<<")"<<std::endl;
        }
    }

public:
    MakeMVAVariables(const bool verb = false, std::string myVarSuffix = "")
        : verb_(verb)
        , myVarSuffix_(myVarSuffix)
    {
    }

    void operator()(NTupleReader& tr)
    {
        makeMVAVariables(tr);
    }
};

#endif
