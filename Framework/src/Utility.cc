#include "Framework/Framework/include/Utility.h"
#include <math.h>
#include "TSystem.h"

namespace utility
{

    double DeltaR(const LorentzVector& v1, const LorentzVector& v2)
    {
        return ROOT::Math::VectorUtil::DeltaR(v1, v2);
    }

    double DeltaPhi(const LorentzVector& v1, const LorentzVector& v2)
    {
        return ROOT::Math::VectorUtil::DeltaPhi(v1, v2);
    }

    LorentzVector Boost(const LorentzVector& v, const BoostVector& b)
    {
        return ROOT::Math::VectorUtil::boost(v, b.BetaVector());
    }

    LorentzVector RotateZ(const LorentzVector& v, const double r)
    {

        auto v3D = v.Vect();
        double t = v.T();

        auto v3Drot = ROOT::Math::VectorUtil::RotateZ(v3D, r);

        LorentzVector lv;
        lv.SetXYZT(v3Drot.X(), v3Drot.Y(), v3Drot.Z(), t);

        return lv;
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

    std::string split(const std::string& half, const std::string& s, const std::string& h)
    {
        std::string token;
        if      ("first"==half) token = s.substr(0, s.find(h));
        else if ("last" ==half) token = s.substr(s.find(h) + h.length(), std::string::npos);
        return token;
    } 

    bool compare_p(const math::RThetaPhiVector& v1, const math::RThetaPhiVector& v2 ) 
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

    const std::vector<std::vector<LorentzVector>> nestVecOfVec(const std::vector<LorentzVector>& vec, const std::vector<int>& counts)
    {
    
        std::vector<std::vector<LorentzVector>> nestedVecs(counts.size());

        unsigned int iSubjet = 0;
        unsigned int iSubjetsVec = 0;
        for (unsigned int j = 0; j < counts.size(); j++)
        {

            unsigned nSubjets = counts.at(j);

            for (unsigned int k = 0; k < nSubjets; k++)
            {
                nestedVecs.at(iSubjetsVec).push_back(vec.at(iSubjet));
                iSubjet++;
            }

            iSubjetsVec++;
        }

        return nestedVecs;
    }

    std::vector<TLorentzVector> convertVectorOfLV(const std::vector<utility::LorentzVector>& vec)
    {
        std::vector<TLorentzVector> newvec(vec.size());
        for(unsigned int i = 0; i < vec.size(); i++)
        {
            auto& lv = vec[i];
            newvec.at(i).SetPtEtaPhiM(lv.Pt(), lv.Eta(), lv.Phi(), lv.M());
        }
        return newvec;
    }

    std::vector<LorentzVector> convertVectorOfTLV(const std::vector<TLorentzVector>& vec)
    {
        std::vector<LorentzVector> newvec(vec.size());
        for(unsigned int i = 0; i < vec.size(); i++)
        {
            auto& tlv = vec[i];
            newvec.at(i).SetPt(tlv.Pt()); newvec.at(i).SetEta(tlv.Eta()); newvec.at(i).SetPhi(tlv.Phi()); newvec.at(i).SetE(tlv.E());
        }
        return newvec;
    }

    std::vector<std::vector<TLorentzVector> > convertVecVecOfLV(const std::vector<std::vector<utility::LorentzVector> >& vecvec)
    {
        std::vector<std::vector<TLorentzVector> > newvecvec(vecvec.size());
        for(unsigned int i = 0; i < vecvec.size(); i++)
        {
            newvecvec.at(i) = convertVectorOfLV(vecvec.at(i));
        }
        return newvecvec;
    }

    TLorentzVector convertLV(const utility::LorentzVector& lv)
    {
        TLorentzVector newTLV;
        newTLV.SetPtEtaPhiM(lv.Pt(), lv.Eta(), lv.Phi(), lv.M());

        return newTLV;
    }

    LorentzVector convertTLV(const TLorentzVector& tlv)
    {
        LorentzVector newLV;
        newLV.SetPt(tlv.Pt()); newLV.SetEta(tlv.Eta()); newLV.SetPhi(tlv.Phi()); newLV.SetE(tlv.E());

        return newLV;
    }

    LorentzVector convertTLV(const TLorentzVector* tlv)
    {
        LorentzVector newLV;
        newLV.SetPt(tlv->Pt()); newLV.SetEta(tlv->Eta()); newLV.SetPhi(tlv->Phi()); newLV.SetE(tlv->E());

        return newLV;
    }
}
