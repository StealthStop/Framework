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
    std::vector<TLorentzVector> convertVectorOfLV(const std::vector<utility::LorentzVector>& vec);
    std::vector<LorentzVector> convertVectorOfTLV(const std::vector<TLorentzVector>& vec);
    std::vector<std::vector<TLorentzVector> > convertVecVecOfLV(const std::vector<std::vector<LorentzVector> >& vecvec);

    TLorentzVector convertLV(const LorentzVector& lv);
    LorentzVector  convertTLV(const TLorentzVector& tlv);
    LorentzVector  convertTLV(const TLorentzVector* tlv);

    const std::vector<std::vector<LorentzVector>> nestVecOfVec(const std::vector<LorentzVector>& vec, const std::vector<int>& counts);

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
}

#endif
