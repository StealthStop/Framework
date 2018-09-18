#ifndef Utility_h
#define Utility_h

#include "TLorentzVector.h"

namespace utility
{
    double calcDPhi(double phi1, double phi2);
    double calcDR(double eta1, double eta2, double phi1, double phi2);
    double calcMT(const TLorentzVector& lepton, const TLorentzVector& met);
}

#endif
