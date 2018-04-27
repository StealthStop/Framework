#ifndef MAKEMVAVARIABLES_H
#define MAKEMVAVARIABLES_H 

#include "Framework/Framework/src/EventShapeVariables.cc"
#include "Framework/Framework/src/get_cmframe_jets.c"

class MakeMVAVariables
{
private:
    bool verb_;
    Long64_t ei_;

    void makeMVAVariables(NTupleReader& tr)
    {
        const auto& Jets                         = tr.getVec<TLorentzVector>("Jets");
        const auto& Jets_bDiscriminatorCSV       = tr.getVec<double>        ("Jets_bDiscriminatorCSV");
        const auto& Jets_muonEnergyFraction      = tr.getVec<double>        ("Jets_muonEnergyFraction");
        const auto& Jets_neutralEmEnergyFraction = tr.getVec<double>        ("Jets_neutralEmEnergyFraction");
        const auto& JetsAK8                      = tr.getVec<TLorentzVector>("JetsAK8");
        const auto& JetsAK8_NsubjettinessTau1    = tr.getVec<double>        ("JetsAK8_NsubjettinessTau1");
        const auto& JetsAK8_NsubjettinessTau2    = tr.getVec<double>        ("JetsAK8_NsubjettinessTau2");
        const auto& JetsAK8_NsubjettinessTau3    = tr.getVec<double>        ("JetsAK8_NsubjettinessTau3");
        const auto& JetsAK8_softDropMass         = tr.getVec<double>        ("JetsAK8_softDropMass");
        const auto& Muons                        = tr.getVec<TLorentzVector>("Muons");
        const auto& Muons_charge                 = tr.getVec<int>           ("Muons_charge");
        const auto& Electrons                    = tr.getVec<TLorentzVector>("Electrons");
        const auto& Electrons_charge             = tr.getVec<int>           ("Electrons_charge");
        const auto& RunNum                       = tr.getVar<UInt_t>        ("RunNum");
        const auto& LumiBlockNum                 = tr.getVar<UInt_t>        ("LumiBlockNum");
        const auto& EvtNum                       = tr.getVar<ULong64_t>     ("EvtNum");
        const auto& nevts_ttree                  = tr.getVar<Long64_t>      ("nevts_ttree");
        const auto& modnum                       = tr.getVar<int>           ("modnum");

        ei_++;

        //--- Initialize all derived ntuple variables.
        double fwm2_top6 = 0.0;
        double fwm3_top6 = 0.0;
        double fwm4_top6 = 0.0;
        double fwm5_top6 = 0.0;
        double fwm6_top6 = 0.0;
        double jmt_ev0_top6 = 0.0;
        double jmt_ev1_top6 = 0.0;
        double jmt_ev2_top6 = 0.0;
        double fwm2_top6_tr_v3pt30 = 0.0;
        double fwm3_top6_tr_v3pt30 = 0.0;
        double fwm4_top6_tr_v3pt30 = 0.0;
        double fwm5_top6_tr_v3pt30 = 0.0;
        double fwm6_top6_tr_v3pt30 = 0.0;
        double jmt_ev0_top6_tr_v3pt30 = 0.0;
        double jmt_ev1_top6_tr_v3pt30 = 0.0;
        double jmt_ev2_top6_tr_v3pt30 = 0.0;
        double event_beta_z = -9.0;
        int    njets_pt45_eta24 = 0;
        int    njets_pt30_eta24 = 0;
        int    njets_pt20_eta24 = 0;
        int    njets_pt45_eta50 = 0;
        int    njets_pt30_eta50 = 0;
        int    njets_pt20_eta50 = 0;
        double pfht_pt40_eta24 = 0.0;
        double pfht_pt45_eta24 = 0.0;
        int    nleptons = 0;
        int    nbtag_csv85_pt30_eta24 = 0;
        double leppt1 = 0.0;
        double m_lep1_b = 0.0;
        double leppt2 = 0.0;
        double m_lep2_b = 0.0;
        int    evt_count = ei_;
        int    run = RunNum;
        int    lumi = LumiBlockNum;
        ULong64_t event = EvtNum;

        if ( verb_ ) 
        {
            printf("\n") ;
            printf(" =============== Count : %lld, run=%10d, lumi=%10d, event=%10llu\n", ei_, run, lumi, event ) ;
            printf("\n") ;
            printf("   AK4 jets:\n" ) ;
            for ( unsigned int rji=0; rji<Jets.size() ; rji ++ ) 
            {
                printf("  %3d : Pt = %7.1f, Eta = %7.3f, Phi = %7.3f | CSV = %8.3f , mu fr = %7.3f , neut EM fr = %7.3f\n",
                       rji, Jets.at(rji).Pt(), Jets.at(rji).Eta(), Jets.at(rji).Phi(),
                       Jets_bDiscriminatorCSV.at(rji),
                       Jets_muonEnergyFraction.at(rji),
                       Jets_neutralEmEnergyFraction.at(rji)
                    ) ;
            } // rji
            printf("\n") ;
        }

        if ( verb_ ) 
        {
            printf("\n") ;
            printf("   AK8 jets:\n" ) ;
            for ( unsigned int fji=0; fji<JetsAK8.size() ; fji++ ) 
            {
                double t31(0) ;
                double t21(0) ;
                if ( JetsAK8_NsubjettinessTau1.at(fji) > 0 ) 
                {
                    t31 = (JetsAK8_NsubjettinessTau3.at(fji)) / (JetsAK8_NsubjettinessTau1.at(fji)) ;
                    t21 = (JetsAK8_NsubjettinessTau2.at(fji)) / (JetsAK8_NsubjettinessTau1.at(fji)) ;
                }
                printf("  %3d : Pt = %7.1f, Eta = %7.3f, Phi = %7.3f | SD mass = %7.1f , tau3/tau1 = %7.3f , tau2/tau1 = %7.3f\n",
                       fji, JetsAK8.at(fji).Pt(), JetsAK8.at(fji).Eta(), JetsAK8.at(fji).Phi(),
                       JetsAK8_softDropMass.at(fji), t31, t21
                    ) ;
            } // ji
            printf("\n") ;
        }

        TLorentzVector rlv_all_jets ;

        int ngood_pt30_eta24(0) ;
        for ( unsigned int rji=0; rji < Jets.size() ; rji++ ) 
        {

            TLorentzVector jlv( Jets.at(rji) ) ;
            TLorentzVector rj_tlv( Jets.at(rji) ) ;

            rlv_all_jets += jlv ;

            double abseta = fabs( jlv.Eta() ) ;
            double pt = jlv.Pt() ;

            if ( abseta < 2.4 ) 
            {
                if ( pt > 45 ) njets_pt45_eta24 ++ ;
                if ( pt > 30 ) njets_pt30_eta24 ++ ;
                if ( pt > 20 ) njets_pt20_eta24 ++ ;
                if ( pt > 30 ) 
                {
                    ngood_pt30_eta24++ ;
                    if ( ngood_pt30_eta24 >= 7 ) 
                    {
                    }
                    if ( Jets_bDiscriminatorCSV.at(rji) > 0.85 ) 
                    {
                        nbtag_csv85_pt30_eta24++ ;
                    }
                }
                if ( pt > 40 ) 
                {
                    pfht_pt40_eta24 += pt ;
                }
                if ( pt > 45 ) 
                {
                    pfht_pt45_eta24 += pt ;
                }
            }
            if ( abseta < 5.0 ) 
            {
                if ( pt > 45 ) njets_pt45_eta50++ ;
                if ( pt > 30 ) njets_pt30_eta50++ ;
                if ( pt > 20 ) njets_pt20_eta50++ ;
            }

        } // rji

        double reco_jets_beta = rlv_all_jets.Pz() / rlv_all_jets.E() ;

        event_beta_z = reco_jets_beta ;

        TVector3 rec_boost_beta_vec( 0., 0., -1.*reco_jets_beta ) ;

        //--- Fill vector of jet momenta in CM frame.
        std::vector<math::RThetaPhiVector>* cm_jets = new std::vector<math::RThetaPhiVector>();

        for ( unsigned int rji=0; rji < Jets.size() ; rji++ ) 
        {
            TLorentzVector jlvcm( Jets.at(rji) ) ;
            jlvcm.Boost( rec_boost_beta_vec ) ;

            math::RThetaPhiVector cmvec( jlvcm.P(), jlvcm.Theta(), jlvcm.Phi() ) ;
            cm_jets->push_back( cmvec ) ;

        } // rji

        //--- Try using only the 6 highest-P jets in the CM frame in the event shape vars.
        //    First, need to make a new input vector of jets containing only those jets.

        std::vector<math::RThetaPhiVector> cm_jets_psort = *cm_jets ;
        std::sort( cm_jets_psort.begin(), cm_jets_psort.end(), compare_p ) ;
        std::vector<math::RThetaPhiVector> cm_jets_top6 ;
        for ( unsigned int ji=0; ji<cm_jets->size(); ji++ ) 
        {
            if ( ji < 6 ) cm_jets_top6.push_back( cm_jets_psort.at(ji) ) ;
        } // ji

        if ( verb_ ) 
        {
            printf("\n\n Unsorted and sorted CM jet lists.\n") ;
            for ( unsigned int ji=0; ji<cm_jets->size(); ji++ ) 
            {
                printf("  %2d :  (%7.1f, %7.3f, %7.3f) | (%7.1f, %7.3f, %7.3f)\n", ji,
                       cm_jets->at(ji).R(), cm_jets->at(ji).Theta(), cm_jets->at(ji).Phi(),
                       cm_jets_psort.at(ji).R(), cm_jets_psort.at(ji).Theta(), cm_jets_psort.at(ji).Phi() ) ;
            } // ji
            printf("\n\n") ;
        }

        EventShapeVariables esv_top6( cm_jets_top6 ) ;

        TVectorD eigen_vals_norm_top6 = esv_top6.getEigenValues() ;

        fwm2_top6 = esv_top6.getFWmoment( 2 ) ;
        fwm3_top6 = esv_top6.getFWmoment( 3 ) ;
        fwm4_top6 = esv_top6.getFWmoment( 4 ) ;
        fwm5_top6 = esv_top6.getFWmoment( 5 ) ;
        fwm6_top6 = esv_top6.getFWmoment( 6 ) ;

        jmt_ev0_top6 = eigen_vals_norm_top6[0] ;
        jmt_ev1_top6 = eigen_vals_norm_top6[1] ;
        jmt_ev2_top6 = eigen_vals_norm_top6[2] ;

        //---- variables related to leptons.

        TLorentzVector lep1_tlv ;
        TLorentzVector lep2_tlv ;
        bool lep1_is_mu(false) ;
        bool lep1_is_e(false) ;
        bool lep2_is_mu(false) ;
        bool lep2_is_e(false) ;
        int lep1_q(0) ;
        int lep2_q(0) ;
        nleptons = 0 ;

        for ( unsigned int mi=0; mi<Muons.size(); mi++ ) 
        {
            if ( Muons.at(mi).Pt() > 0 ) nleptons ++ ;
            if ( Muons.at(mi).Pt() > lep1_tlv.Pt() ) 
            {
                if ( lep1_tlv.Pt() > 0 ) 
                {
                    lep2_tlv = lep1_tlv ;
                    lep2_q = lep1_q ;
                    lep2_is_mu = lep1_is_mu ;
                    lep2_is_e  = lep1_is_e  ;
                }
                lep1_tlv =  Muons.at(mi) ;
                lep1_q = Muons_charge.at(mi) ;
                lep1_is_mu = true ;
            } 
            else if ( Muons.at(mi).Pt() > lep2_tlv.Pt() ) 
            {
                lep2_tlv =  Muons.at(mi) ;
                lep2_q = Muons_charge.at(mi) ;
                lep2_is_mu = true ;
            }
        }
        for ( unsigned int ei=0; ei<Electrons.size(); ei++ ) 
        {
            if ( Electrons.at(ei).Pt() > 0 ) nleptons ++ ;
            if ( Electrons.at(ei).Pt() > lep1_tlv.Pt() ) 
            {
                if ( lep1_tlv.Pt() > 0 ) 
                {
                    lep2_tlv = lep1_tlv ;
                    lep2_q = lep1_q ;
                    lep2_is_mu = lep1_is_mu ;
                    lep2_is_e  = lep1_is_e ;
                }
                lep1_tlv =  Electrons.at(ei) ;
                lep1_q = Electrons_charge.at(ei) ;
                lep1_is_e = true ;
            } 
            else if ( Electrons.at(ei).Pt() > lep2_tlv.Pt() ) 
            {
                lep2_tlv =  Electrons.at(ei) ;
                lep2_q = Electrons_charge.at(ei) ;
                lep2_is_e = true ;
            }
        }

        leppt1 = lep1_tlv.Pt() ;
        leppt2 = lep2_tlv.Pt() ;

        int lowest_m_lep1_b_bjet_rji(-1) ;
        int lowest_m_lep2_b_bjet_rji(-1) ;

        TLorentzVector lowest_m_lep1_b_tlv ;
        TLorentzVector lowest_m_lep2_b_tlv ;

        if ( nleptons > 0 ) 
        {
            double lowest_m_lep1_b(999.) ;
            double lowest_m_lep2_b(999.) ;

            for ( unsigned int rji=0; rji<Jets.size(); rji++ ) 
            {
                TLorentzVector j_tlv = Jets.at(rji) ;

                if ( lep1_tlv.Pt() > 5 && lep1_tlv.DeltaR( j_tlv ) < 0.05 ) 
                {
                    if ( verb_  ) 
                    {
                        printf( "\n  jet %2u is a DR match for lep1 %2s : %7.3f  p ratio = %7.3f  mu Efr = %5.3f  EMfr = %5.3f  \n", rji,
                                ( lep1_is_mu ? "mu" : "e " ),
                                lep1_tlv.DeltaR( j_tlv ),
                                j_tlv.P() / lep1_tlv.P(),
                                Jets_muonEnergyFraction.at(rji),
                                Jets_neutralEmEnergyFraction.at(rji) ) ;
                        printf("   lepton  Pt = %7.1f, Eta = %7.3f, Phi = %7.3f\n", lep1_tlv.Pt(), lep1_tlv.Eta(), lep1_tlv.Phi() ) ;
                        printf("   jet     Pt = %7.1f, Eta = %7.3f, Phi = %7.3f\n", j_tlv.Pt(), j_tlv.Eta(), j_tlv.Phi() ) ;
                    }
                    continue ;
                }
                if ( lep2_tlv.Pt() > 5 && lep2_tlv.DeltaR( j_tlv ) < 0.05 ) 
                {
                    if ( verb_  ) 
                    {
                        printf( "\n  jet %2u is a DR match for lep2 %2s : %7.3f  p ratio = %7.3f  mu Efr = %5.3f  EMfr = %5.3f  \n", rji,
                                ( lep2_is_mu ? "mu" : "e " ),
                                lep2_tlv.DeltaR( j_tlv ),
                                j_tlv.P() / lep2_tlv.P(),
                                Jets_muonEnergyFraction.at(rji),
                                Jets_neutralEmEnergyFraction.at(rji) ) ;
                        printf("   lepton  Pt = %7.1f, Eta = %7.3f, Phi = %7.3f\n", lep2_tlv.Pt(), lep2_tlv.Eta(), lep2_tlv.Phi() ) ;
                        printf("   jet     Pt = %7.1f, Eta = %7.3f, Phi = %7.3f\n", j_tlv.Pt(), j_tlv.Eta(), j_tlv.Phi() ) ;
                    }
                    continue ;
                }

                if ( Jets_bDiscriminatorCSV.at(rji) > 0.85 ) 
                {
                    TLorentzVector lep1_b_tlv = lep1_tlv + j_tlv ;
                    TLorentzVector lep2_b_tlv = lep2_tlv + j_tlv ;

                    if ( lep1_b_tlv.M() < lowest_m_lep1_b ) 
                    {
                        lowest_m_lep1_b = lep1_b_tlv.M() ;
                        m_lep1_b = lep1_b_tlv.M() ;
                        lowest_m_lep1_b_bjet_rji = rji ;
                        lowest_m_lep1_b_tlv = lep1_b_tlv ;
                    }
                    if ( lep2_tlv.Pt() > 0 && lep2_b_tlv.M() < lowest_m_lep2_b ) 
                    {
                        lowest_m_lep2_b = lep2_b_tlv.M() ;
                        m_lep2_b = lep2_b_tlv.M() ;
                        lowest_m_lep2_b_bjet_rji = rji ;
                        lowest_m_lep2_b_tlv = lep2_b_tlv ;
                    }

                }

            } // rji

            //--- If have two lepton-b pairs, make sure they don't both use the same b jet.

            if ( nleptons >= 2 && m_lep1_b > 0 && m_lep2_b > 0 ) 
            {
                if ( lowest_m_lep1_b_bjet_rji == lowest_m_lep2_b_bjet_rji ) 
                {
                    if ( verb_ ) 
                    {
                        printf("  both leptons pair with the same b jet.  Choosing one with mass in (30,180).\n" ) ;
                        printf("     rji for lep1 = %d, mass = %7.1f, rji for b with lep2 = %d, mass = %7.1f\n",
                               lowest_m_lep1_b_bjet_rji, m_lep1_b, lowest_m_lep2_b_bjet_rji, m_lep2_b ) ;
                    }

                    bool m_lep1_b_in_window(false) ;
                    bool m_lep2_b_in_window(false) ;
                    if ( m_lep1_b > 30 && m_lep1_b < 180 ) m_lep1_b_in_window = true ;
                    if ( m_lep2_b > 30 && m_lep2_b < 180 ) m_lep2_b_in_window = true ;

                    if ( m_lep1_b_in_window && !m_lep2_b_in_window ) m_lep2_b = 0 ;
                    if ( m_lep2_b_in_window && !m_lep1_b_in_window ) m_lep1_b = 0 ;
                    if ( m_lep1_b_in_window && m_lep2_b_in_window ) 
                    {
                        if ( m_lep1_b < m_lep2_b ) 
                        {
                            m_lep2_b = 0. ;
                        } 
                        else 
                        {                            
                            m_lep1_b = 0. ;
                        }
                    }
                    if ( !m_lep1_b_in_window && !m_lep2_b_in_window ) 
                    {
                        if ( m_lep1_b < m_lep2_b ) 
                        {
                            m_lep2_b = 0. ;
                        } 
                        else 
                        {
                            m_lep1_b = 0. ;
                        }
                    }

                } // same b jet?

            } // two lepton-b pairs?

        } // any leptons?

        if ( verb_ ) 
        {
            printf("\n  Lepton info:   nleptons = %d\n", nleptons) ;
            if ( leppt1 > 0 ) 
            {
                printf("  Lepton 1 :  Pt = %7.1f,  %2s,  m(b,lep) = %7.1f\n", leppt1, (lep1_is_mu?"mu":" e"), m_lep1_b ) ;
            }
            if ( leppt2 > 0 ) 
            {
                printf("  Lepton 2 :  Pt = %7.1f,  %2s,  m(b,lep) = %7.1f\n", leppt2, (lep2_is_mu?"mu":" e"), m_lep2_b ) ;
            }
            printf("\n") ;
            printf("  Njets(pt>20,|eta|<2.4) : %2d\n", njets_pt20_eta24 ) ;
            printf("  Njets(pt>30,|eta|<2.4) : %2d\n", njets_pt30_eta24 ) ;
            printf("  Njets(pt>45,|eta|<2.4) : %2d\n", njets_pt45_eta24 ) ;
            printf("  Njets(pt>20,|eta|<5.0) : %2d\n", njets_pt20_eta50 ) ;
            printf("  Njets(pt>30,|eta|<5.0) : %2d\n", njets_pt30_eta50 ) ;
            printf("  Njets(pt>45,|eta|<5.0) : %2d\n", njets_pt45_eta50 ) ;
            printf("  Number of CVS>0.85 b tags : %2d\n", nbtag_csv85_pt30_eta24 ) ;
            printf("\n") ;
            printf("  PFHT = %7.1f  \n", pfht_pt40_eta24 ) ;
            printf("\n") ;
        }

        tr.registerDerivedVec("cm_jets", cm_jets);
        tr.registerDerivedVar("fwm2_top6", fwm2_top6);
        tr.registerDerivedVar("fwm3_top6", fwm3_top6);
        tr.registerDerivedVar("fwm4_top6", fwm4_top6);
        tr.registerDerivedVar("fwm5_top6", fwm5_top6);
        tr.registerDerivedVar("fwm6_top6", fwm6_top6);
        tr.registerDerivedVar("jmt_ev0_top6", jmt_ev0_top6);
        tr.registerDerivedVar("jmt_ev1_top6", jmt_ev1_top6);
        tr.registerDerivedVar("jmt_ev2_top6", jmt_ev2_top6);
        tr.registerDerivedVar("fwm2_top6_tr_v3pt30", fwm2_top6_tr_v3pt30);
        tr.registerDerivedVar("fwm3_top6_tr_v3pt30", fwm3_top6_tr_v3pt30);
        tr.registerDerivedVar("fwm4_top6_tr_v3pt30", fwm4_top6_tr_v3pt30);
        tr.registerDerivedVar("fwm5_top6_tr_v3pt30", fwm5_top6_tr_v3pt30);
        tr.registerDerivedVar("fwm6_top6_tr_v3pt30", fwm6_top6_tr_v3pt30);
        tr.registerDerivedVar("jmt_ev0_top6_tr_v3pt30", jmt_ev0_top6_tr_v3pt30);
        tr.registerDerivedVar("jmt_ev1_top6_tr_v3pt30", jmt_ev1_top6_tr_v3pt30);
        tr.registerDerivedVar("jmt_ev2_top6_tr_v3pt30", jmt_ev2_top6_tr_v3pt30);
        tr.registerDerivedVar("event_beta_z", event_beta_z);
        tr.registerDerivedVar("njets_pt45_eta24", njets_pt45_eta24);
        tr.registerDerivedVar("njets_pt30_eta24", njets_pt30_eta24);
        tr.registerDerivedVar("njets_pt20_eta24", njets_pt20_eta24);
        tr.registerDerivedVar("njets_pt45_eta50", njets_pt45_eta50);
        tr.registerDerivedVar("njets_pt30_eta50", njets_pt30_eta50);
        tr.registerDerivedVar("njets_pt20_eta50", njets_pt20_eta50);
        tr.registerDerivedVar("pfht_pt40_eta24", pfht_pt40_eta24);
        tr.registerDerivedVar("pfht_pt45_eta24", pfht_pt45_eta24);
        tr.registerDerivedVar("nleptons", nleptons);
        tr.registerDerivedVar("nbtag_csv85_pt30_eta24", nbtag_csv85_pt30_eta24); 
        tr.registerDerivedVar("leppt1", leppt1);
        tr.registerDerivedVar("m_lep1_b", m_lep1_b);
        tr.registerDerivedVar("leppt2", leppt2);
        tr.registerDerivedVar("m_lep2_b", m_lep2_b);
        tr.registerDerivedVar("evt_count", evt_count);
        tr.registerDerivedVar("run", run);
        tr.registerDerivedVar("lumi", lumi);
        tr.registerDerivedVar("event", event);
    }

public:
    MakeMVAVariables(const bool verb = false) : verb_(verb), ei_(0)
    {
    }

    void operator()(NTupleReader& tr)
    {
        makeMVAVariables(tr);
    }
};

#endif
