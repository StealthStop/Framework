#ifndef Utility_h
#define Utility_h

#include "Framework/Framework/include/Vector3D.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "Math/Boost.h"
#include <cmath>

namespace utility
{
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>> LorentzVector;
    typedef ROOT::Math::Boost BoostVector;

    double DeltaR(const LorentzVector& v1, const LorentzVector& v2);
    double DeltaPhi(const LorentzVector& v1, const LorentzVector& v2);
    LorentzVector RotateZ(const LorentzVector& v, const double r);
    LorentzVector Boost(const LorentzVector& v, const BoostVector& b);
    const std::string color(const std::string& text, const std::string& color);
    std::string split(const std::string& half, const std::string& s, const std::string& h);
    bool compare_p(const math::RThetaPhiVector& v1, const math::RThetaPhiVector& v2);
    void get_cmframe_jets(const std::vector<TLorentzVector>* lab_frame_jets, std::vector<math::RThetaPhiVector>& cm_frame_jets, int max_number_of_jets = -1 );
    std::vector<std::vector<TLorentzVector> > convertVecVecOfLV(const std::vector<std::vector<LorentzVector> >& vecvec);

    std::vector<std::string> splitString(const std::string& string, const char delim = ',');

    template<typename T> T sum2(T v) { return v*v; }
    template<typename T, typename... Args> T sum2(T v, Args... args) { return v*v + sum2(args...); }
    template<typename T, typename... Args> T addInQuad(T v, Args... args) { return sqrt(sum2(v, args...)); }

    template<typename LV1, typename LV2>
    double calcMT(const LV1& lepton, const LV2& met)
    {
        // Assuming that both lepton and met are massless
        const double mt_sq = 2 * lepton.Pt() * met.Pt() * ( 1-cos(met.Phi()-lepton.Phi()) );
        return sqrt(mt_sq);
    }    

    template<typename LV1, typename LV2>
    bool compare_pt_TLV(const LV1& v1, const LV1& v2 )
    {
        return ( v1.Pt() > v2.Pt() );
    }

    template<typename LV1, typename LV2>
    LV1 convertLV(const LV2& lv)
    {
        LV1 newLV;
        newLV.SetPxPyPzE(lv.Px(), lv.Py(), lv.Pz(), lv.E());

        return newLV;
    }

    template<typename LV1, typename LV2>
    LV1 convertLV(const LV2* lv)
    {
        LV1 newLV;
        newLV.SetPxPyPzE(lv->Px(), lv->Py(), lv->Pz(), lv->E());

        return newLV;
    }

    template<typename LV1, typename LV2>
    std::vector<LV1> convertVectorOfLV(const std::vector<LV2>& vec)
    {
        std::vector<LV1> newvec(vec.size());
        for(unsigned int i = 0; i < vec.size(); i++)
        {
            auto& tlv = vec[i];
            newvec.at(i).SetPxPyPzE(tlv.Px(), tlv.Py(), tlv.Pz(), tlv.E());
        }
        return newvec;
    }

    template<typename LV1, typename LV2>
    std::vector<std::vector<LV1> > convertVecVecOfLV(const std::vector<std::vector<LV2> >& vecvec)
    {
        std::vector<std::vector<LV1> > newvecvec(vecvec.size());
        for(unsigned int i = 0; i < vecvec.size(); i++)
        {
            newvecvec.at(i) = convertVectorOfLV<LV1, LV2>(vecvec.at(i));
        }
        return newvecvec;
    }

    template<typename LV1, typename LV2>
    std::vector<std::vector<LV1>> nestVecOfVec(const std::vector<LV2>& vec, const std::vector<int>& counts)
    {

        std::vector<std::vector<LV1>> nestedVecs(counts.size());

        unsigned int iSubjet = 0;
        unsigned int iSubjetsVec = 0;
        for (unsigned int j = 0; j < counts.size(); j++)
        {

            unsigned nSubjets = counts.at(j);

            for (unsigned int k = 0; k < nSubjets; k++)
            {
                nestedVecs.at(iSubjetsVec).push_back(convertLV<LV1, LV2>(vec.at(iSubjet)));
                iSubjet++;
            }

            iSubjetsVec++;
        }

        return nestedVecs;
    }
}

#endif
