#ifndef MAKEMVAVARIABLES_H
#define MAKEMVAVARIABLES_H 

#include "Framework/Framework/src/EventShapeVariables.cc"
#include "Framework/Framework/src/get_cmframe_jets.c"

class MakeMVAVariables
{
private:
    bool verb_;

    void makeMVAVariables(NTupleReader& tr)
    {
        const auto& Jets = tr.getVec<TLorentzVector>("Jets");

        //--- Initialize all derived ntuple variables.
        double fwm2_top6 = 0.0;
        double fwm3_top6 = 0.0;
        double fwm4_top6 = 0.0;
        double fwm5_top6 = 0.0;
        double fwm6_top6 = 0.0;
        double jmt_ev0_top6 = 0.0;
        double jmt_ev1_top6 = 0.0;
        double jmt_ev2_top6 = 0.0;
        double event_beta_z = -9.0;

        TLorentzVector rlv_all_jets ;

        for ( unsigned int rji=0; rji < Jets.size() ; rji++ ) 
        {
            TLorentzVector jlv( Jets.at(rji) ) ;
            rlv_all_jets += jlv ;

        } // rji

        double reco_jets_beta = rlv_all_jets.Pz() / rlv_all_jets.E() ;
        event_beta_z = reco_jets_beta ;

        TVector3 rec_boost_beta_vec( 0.0, 0.0, -1.0*reco_jets_beta ) ;

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

        // register Variables
        tr.registerDerivedVec("cm_jets", cm_jets);
        tr.registerDerivedVar("fwm2_top6", fwm2_top6);
        tr.registerDerivedVar("fwm3_top6", fwm3_top6);
        tr.registerDerivedVar("fwm4_top6", fwm4_top6);
        tr.registerDerivedVar("fwm5_top6", fwm5_top6);
        tr.registerDerivedVar("fwm6_top6", fwm6_top6);
        tr.registerDerivedVar("jmt_ev0_top6", jmt_ev0_top6);
        tr.registerDerivedVar("jmt_ev1_top6", jmt_ev1_top6);
        tr.registerDerivedVar("jmt_ev2_top6", jmt_ev2_top6);
        tr.registerDerivedVar("event_beta_z", event_beta_z);
    }

public:
    MakeMVAVariables(const bool verb = false) : verb_(verb)
    {
    }

    void operator()(NTupleReader& tr)
    {
        makeMVAVariables(tr);
    }
};

#endif
