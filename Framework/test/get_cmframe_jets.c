#ifndef get_cmframe_jets_c
#define get_cmframe_jets_c

#include "TMath.h"
#include "Vector3D.h"
#include "TLorentzVector.h"
#include "TSystem.h"


#include <vector>
using std::vector ;

static bool compare_p( math::RThetaPhiVector v1, math::RThetaPhiVector v2 ) { return (v1.R() > v2.R() ) ; }

static void get_cmframe_jets( std::vector<TLorentzVector>* lab_frame_jets,
                              std::vector<math::RThetaPhiVector>& cm_frame_jets,
                              int max_number_of_jets = -1 ) {

      if ( lab_frame_jets == 0x0 ) {
         printf("\n\n *** get_cmframe_jets : null pointer for input vector of lab frame jets (1st argument).\n\n" ) ;
         gSystem -> Exit(-1) ;
      }

       //--- Find boost to CM frame.  Use all jets to get best estimate.

      TLorentzVector rlv_all_jets ;

      for ( unsigned int rji=0; rji < lab_frame_jets->size() ; rji++ ) {

         rlv_all_jets += lab_frame_jets->at(rji) ;

      } // rji



      double reco_jets_beta = rlv_all_jets.Pz() / rlv_all_jets.E() ;

      TVector3 rec_boost_beta_vec( 0., 0., -1.*reco_jets_beta ) ;



     //--- Fill vector of jet momenta in CM frame.

      std::vector<math::RThetaPhiVector> cm_jets ;

      for ( unsigned int rji=0; rji < lab_frame_jets->size() ; rji++ ) {

         TLorentzVector jlvcm( lab_frame_jets->at(rji) ) ;
         jlvcm.Boost( rec_boost_beta_vec ) ;

         math::RThetaPhiVector cmvec( jlvcm.P(), jlvcm.Theta(), jlvcm.Phi() ) ;
         cm_jets.push_back( cmvec ) ;

      } // rji

      if ( max_number_of_jets < 0 ) {
         cm_frame_jets = cm_jets ;
         return ;
      }


     //--- If requesting only the top N jets, figure out which those are.
     //--- This sorting is by the CM frame momentum (vector magnitude), not Pt.

      std::vector<math::RThetaPhiVector> cm_jets_psort = cm_jets ;
      std::sort( cm_jets_psort.begin(), cm_jets_psort.end(), compare_p ) ;
      std::vector<math::RThetaPhiVector> cm_jets_topn ;
      for ( unsigned int ji=0; ji<cm_jets.size(); ji++ ) {
         int ji_int = ji ;  // get rid of compiler warning in if statement below
         if ( ji_int < max_number_of_jets ) {
            cm_jets_topn.push_back( cm_jets_psort.at(ji) ) ;
         }
      } // ji

      cm_frame_jets = cm_jets_topn ;



   } // get_cmframe_jets



#endif


