#pragma once

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooSillGeneralizedPdf : public RooAbsPdf {
public:
  RooSillGeneralizedPdf() {} // For serialization only
  RooSillGeneralizedPdf(const char *name, const char *title,
             RooAbsReal& _x,
             RooAbsReal& _mass,
             RooAbsReal& _gamma,
             RooAbsReal& _eth,
             RooAbsReal& _l);
  RooSillGeneralizedPdf(const RooSillGeneralizedPdf& other, const char* name = nullptr);
  virtual TObject* clone(const char* newname) const override { return new RooSillGeneralizedPdf(*this, newname); }
  inline virtual ~RooSillGeneralizedPdf() {}

protected:
  RooRealProxy x;       // Observable (E)
  RooRealProxy mass;    // Mass (M)
  RooRealProxy gamma;   // Width (Gamma)
  RooRealProxy eth; // Threshold energy (Eth)
  RooRealProxy l;  // Tota angular momentum (l)

  Double_t evaluate() const override;

private:
  ClassDefOverride(RooSillGeneralizedPdf, 1)
};

///////////////////////////// Class Implementation /////////////////////////////

#include <RooRealVar.h>
#include <RooMath.h>
#include <cmath>

ClassImp(RooSillGeneralizedPdf)

RooSillGeneralizedPdf::RooSillGeneralizedPdf(const char *name, const char *title,
                       RooAbsReal& _x,
                       RooAbsReal& _mass,
                       RooAbsReal& _gamma,
                       RooAbsReal& _eth,
                       RooAbsReal& _l)
  : RooAbsPdf(name, title),
    x("x", "E", this, _x),
    mass("mass", "Mass", this, _mass),
    gamma("gamma", "Width", this, _gamma),
    eth("eth", "Threshold energy", this, _eth),
    l("l", "Total angular momentum", this, _l)
{
}

RooSillGeneralizedPdf::RooSillGeneralizedPdf(const RooSillGeneralizedPdf& other, const char* name)
  : RooAbsPdf(other, name),
    x("x", this, other.x),
    mass("mass", this, other.mass),
    gamma("gamma", this, other.gamma),
    eth("eth", this, other.eth),
    l("l", this, other.l)
{}

double RooSillGeneralizedPdf::evaluate() const {
  
  double E = x;
  double M = mass;
  double G = gamma;
  double Eth = eth;

  if (E <= Eth) return 0.0;

  double E2 = E * E;
  double M2 = M * M;
  double Eth2 = Eth * Eth;
  
  double gamma_tilde_denom = (std::pow(M2 - Eth2, l + 0.5));
  if (gamma_tilde_denom == 0.0) return 0.0;
  
  
  double gamma_tilde = G * (std::pow(M, 2. * l + 1)) / gamma_tilde_denom;

  double numerator_fraction = std::pow(E2 - Eth2, l + 0.5) / std::pow(E, 2. * l);
  double numerator = gamma_tilde * numerator_fraction;
  double denominator = (E2 - M2)*(E2 - M2) + (gamma_tilde * numerator_fraction)*(gamma_tilde * numerator_fraction);

  return (2.0 * E / M_PI) * (numerator / denominator);
}
