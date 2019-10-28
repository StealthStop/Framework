#include "Framework/Framework/include/Utility.h"
#include <math.h>
#include "TSystem.h"

namespace utility
{
    double calcDPhi(const double phi1, const double phi2)
    {
        double dphi = phi1 - phi2 ;
        if ( dphi >  M_PI ) dphi -= 2*M_PI ;
        if ( dphi < -M_PI ) dphi += 2*M_PI ;
        return dphi;
    }
    
    double calcDR(const double eta1, const double eta2, const double phi1, const double phi2)
    {
        const double deta = fabs( eta1 - eta2 ) ;
        double dphi = phi1 - phi2 ;
        if ( dphi > M_PI ) dphi -= 2*M_PI ;
        if ( dphi <-M_PI ) dphi += 2*M_PI ;        
        return sqrt( dphi*dphi + deta*deta ) ;
    }

    double calcMT(const TLorentzVector& lepton, const TLorentzVector& met)
    {
        // Assuming that both lepton and met are massless
        const double mt_sq = 2 * lepton.Pt() * met.Pt() * ( 1-cos(met.Phi()-lepton.Phi()) );
        return sqrt(mt_sq);
    }    

    const std::string color(const std::string& text, const std::string& color)
    {
        std::string c;
        if(color=="red") c = "31";
        else if(color=="green") c = "32";
        else if(color=="yellow") c = "33";
        else if(color=="blue") c = "34";
        else if(color=="white") c = "37";       
        return "\033[1;"+c+"m"+ text +"\033[0m";
    }    

    bool compare_p( math::RThetaPhiVector v1, math::RThetaPhiVector v2 ) 
    { 
        return ( v1.R() > v2.R() ); 
    }

    void get_cmframe_jets(const std::vector<TLorentzVector>* lab_frame_jets, std::vector<math::RThetaPhiVector>& cm_frame_jets, int max_number_of_jets ) 
    {
      if( lab_frame_jets == nullptr ) 
      {
         printf("\n\n *** get_cmframe_jets : null pointer for input vector of lab frame jets (1st argument).\n\n" );
         gSystem -> Exit(-1);
      }

       //--- Find boost to CM frame.  Use all jets to get best estimate.
      TLorentzVector rlv_all_jets;
      for( unsigned int rji=0; rji < lab_frame_jets->size() ; rji++ ) 
      {
         rlv_all_jets += lab_frame_jets->at(rji);
      }

      double reco_jets_beta = rlv_all_jets.Pz() / rlv_all_jets.E();
      TVector3 rec_boost_beta_vec( 0., 0., -1.*reco_jets_beta );

      //--- Fill vector of jet momenta in CM frame.
      std::vector<math::RThetaPhiVector> cm_jets;
      for( unsigned int rji=0; rji < lab_frame_jets->size() ; rji++ ) 
      {
         TLorentzVector jlvcm( lab_frame_jets->at(rji) );
         jlvcm.Boost( rec_boost_beta_vec );

         math::RThetaPhiVector cmvec( jlvcm.P(), jlvcm.Theta(), jlvcm.Phi() );
         cm_jets.push_back( cmvec );
      }

      if( max_number_of_jets < 0 ) 
      {
         cm_frame_jets = cm_jets;
         return;
      }

      //--- If requesting only the top N jets, figure out which those are.
      //--- This sorting is by the CM frame momentum (vector magnitude), not Pt.
      std::vector<math::RThetaPhiVector> cm_jets_psort = cm_jets;
      std::sort( cm_jets_psort.begin(), cm_jets_psort.end(), compare_p );
      std::vector<math::RThetaPhiVector> cm_jets_topn;
      for( unsigned int ji=0; ji<cm_jets.size(); ji++ ) 
      {
          if( int(ji) < max_number_of_jets ) 
         {
            cm_jets_topn.push_back( cm_jets_psort.at(ji) );
         }
      }
      cm_frame_jets = cm_jets_topn;
   } // get_cmframe_jets
}
