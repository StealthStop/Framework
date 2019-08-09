#ifndef Utility_h
#define Utility_h

#include "TLorentzVector.h"
#include <cmath>

namespace utility
{
    double calcDPhi(double phi1, double phi2);
    double calcDR(double eta1, double eta2, double phi1, double phi2);
    double calcMT(const TLorentzVector& lepton, const TLorentzVector& met);

    template<typename T> T sum2(T v) { return v*v; }
    template<typename T, typename... Args> T sum2(T v, Args... args) { return v*v + sum2(args...); }
    template<typename T, typename... Args> T addInQuad(T v, Args... args) { return sqrt(sum2(v, args...)); }
}

#endif
