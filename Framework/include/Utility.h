#ifndef Utility_h
#define Utility_h

#include "Framework/Framework/include/Vector3D.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include <cmath>

namespace utility
{
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>> LorentzVector;

    double calcDPhi(const double phi1, const double phi2);
    double calcDR(const double eta1, const double eta2, const double phi1, const double phi2);
    double DeltaR(const LorentzVector& v1, const LorentzVector& v2);
    double DeltaPhi(const LorentzVector& v1, const LorentzVector& v2);
    const std::string color(const std::string& text, const std::string& color);
    std::string split(const std::string& half, const std::string& s, const std::string& h);
    bool compare_p(const math::RThetaPhiVector& v1, const math::RThetaPhiVector& v2);
    bool compare_pt_TLV(const TLorentzVector& v1, const TLorentzVector& v2);
    void get_cmframe_jets(const std::vector<TLorentzVector>* lab_frame_jets, std::vector<math::RThetaPhiVector>& cm_frame_jets, int max_number_of_jets = -1 );
    std::vector<TLorentzVector> convertVectorOfLV(const std::vector<utility::LorentzVector>& vec);

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
}

#endif
