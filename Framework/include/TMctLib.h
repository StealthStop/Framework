#ifndef TMctLib_H
#define TMctLib_H

#include <TLorentzVector.h>
#include <TVector2.h>
#include <Framework/Framework/include/mctlib.h>

class TMctLib : public mctlib {

 public:

  TMctLib();
  ~TMctLib();

  double mctcorr(const TLorentzVector v1, const TLorentzVector v2, const TLorentzVector vds, const TVector2 ptm, const double ecm=14000.0, const double mxlo=0.0);
  double mct(const TLorentzVector v1, const TLorentzVector v2);
  double mt2(const TLorentzVector v1, const TLorentzVector v2, const TLorentzVector vds, const TVector2 ptm, const double ecm=14000.0, const double mxlo=0.0);
  double mt2neg(const TLorentzVector v1, const TLorentzVector v2, const TVector2 ptm, const double mxlo=0.0);
  double mcy(const TLorentzVector v1, const TLorentzVector v2, const TLorentzVector vds, const TVector2 ptm);
  double mcx(const TLorentzVector v1, const TLorentzVector v2, const TLorentzVector vds, const TVector2 ptm);

  using mctlib::mctcorr;
  using mctlib::mct;
  using mctlib::mt2;
  using mctlib::mt2neg;
  using mctlib::mcy;
  using mctlib::mcx;

};

#endif
