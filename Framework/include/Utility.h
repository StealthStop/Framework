#ifndef Utility_h
#define Utility_h

#include "Framework/Framework/include/Vector3D.h"
#include "TLorentzVector.h"
#include <cmath>

namespace utility
{
    double calcDPhi(const double phi1, const double phi2);
    double calcDR(const double eta1, const double eta2, const double phi1, const double phi2);
    double calcMT(const TLorentzVector& lepton, const TLorentzVector& met);
    const std::string color(const std::string& text, const std::string& color);
    bool compare_p( math::RThetaPhiVector v1, math::RThetaPhiVector v2 );
    void get_cmframe_jets(const std::vector<TLorentzVector>* lab_frame_jets, std::vector<math::RThetaPhiVector>& cm_frame_jets, int max_number_of_jets = -1 );

    template<typename T> T sum2(T v) { return v*v; }
    template<typename T, typename... Args> T sum2(T v, Args... args) { return v*v + sum2(args...); }
    template<typename T, typename... Args> T addInQuad(T v, Args... args) { return sqrt(sum2(v, args...)); }
}

#endif
