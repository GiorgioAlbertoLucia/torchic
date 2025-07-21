#pragma once

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooSillPdf : public RooAbsPdf {
public:
  RooSillPdf() {} // For serialization only
  RooSillPdf(const char *name, const char *title,
             RooAbsReal& _x,
             RooAbsReal& _mass,
             RooAbsReal& _gamma,
             RooAbsReal& _eth);
  RooSillPdf(const RooSillPdf& other, const char* name = nullptr);
  virtual TObject* clone(const char* newname) const override { return new RooSillPdf(*this, newname); }
  inline virtual ~RooSillPdf() {}

protected:
  RooRealProxy x;       // Observable (E)
  RooRealProxy mass;    // Mass (M)
  RooRealProxy gamma;   // Width (Gamma)
  RooRealProxy eth;     // Threshold (Eth)

  Double_t evaluate() const override;

private:
  ClassDefOverride(RooSillPdf, 1)
};

///////////////////////////// Class Implementation /////////////////////////////

#include <RooRealVar.h>
#include <RooMath.h>
#include <cmath>

ClassImp(RooSillPdf)

RooSillPdf::RooSillPdf(const char *name, const char *title,
                       RooAbsReal& _x,
                       RooAbsReal& _mass,
                       RooAbsReal& _gamma,
                       RooAbsReal& _eth)
  : RooAbsPdf(name, title),
    x("x", "Observable", this, _x),
    mass("mass", "Mass", this, _mass),
    gamma("gamma", "Width", this, _gamma),
    eth("eth", "Threshold", this, _eth)
{}

RooSillPdf::RooSillPdf(const RooSillPdf& other, const char* name)
  : RooAbsPdf(other, name),
    x("x", this, other.x),
    mass("mass", this, other.mass),
    gamma("gamma", this, other.gamma),
    eth("eth", this, other.eth)
{}

double RooSillPdf::evaluate() const {
  double E = x;
  double M = mass;
  double G = gamma;
  double Eth = eth;

  if (E <= Eth) return 0.0;

  double E2 = E * E;
  double M2 = M * M;
  double Eth2 = Eth * Eth;

  double denom_sqrt = std::sqrt(M2 - Eth2);
  if (denom_sqrt == 0.0) return 0.0;

  double gamma_tilde = G * M / denom_sqrt;
  double root_term = std::sqrt(E2 - Eth2);

  double numerator = root_term * gamma_tilde;
  double denominator = (E2 - M2)*(E2 - M2) + (root_term * gamma_tilde)*(root_term * gamma_tilde);

  return (2.0 * E / M_PI) * (numerator / denominator);
}
