#pragma once

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooSillKstarPdf : public RooAbsPdf {
public:
  RooSillKstarPdf() {} // For serialization only
  RooSillKstarPdf(const char *name, const char *title,
             RooAbsReal& _x,
             RooAbsReal& _mass,
             RooAbsReal& _gamma,
             RooAbsReal& _mass_daughter1,
             RooAbsReal& _mass_daughter2,
             RooAbsReal& _l);
  RooSillKstarPdf(const RooSillKstarPdf& other, const char* name = nullptr);
  virtual TObject* clone(const char* newname) const override { return new RooSillKstarPdf(*this, newname); }
  inline virtual ~RooSillKstarPdf() {}

protected:
  RooRealProxy x;       // Observable (E)
  RooRealProxy mass;    // Mass (M)
  RooRealProxy gamma;   // Width (Gamma)
  RooRealProxy mass_daughter1; // Mass of daughter particle 1
  RooRealProxy mass_daughter2; // Mass of daughter particle 2
  RooRealProxy l;  // Tota angular momentum (l)

  Double_t evaluate() const override;

private:
  ClassDefOverride(RooSillKstarPdf, 1)
};

///////////////////////////// Class Implementation /////////////////////////////

#include <RooRealVar.h>
#include <RooMath.h>
#include <cmath>

ClassImp(RooSillKstarPdf)

RooSillKstarPdf::RooSillKstarPdf(const char *name, const char *title,
                       RooAbsReal& _x,
                       RooAbsReal& _mass,
                       RooAbsReal& _gamma,
                       RooAbsReal& _mass_daughter1,
                       RooAbsReal& _mass_daughter2,
                       RooAbsReal& _l)
  : RooAbsPdf(name, title),
    x("x", "kstar", this, _x),
    mass("mass", "Mass", this, _mass),
    gamma("gamma", "Width", this, _gamma),
    mass_daughter1("mass_daughter1", "Mass of daughter particle 1", this, _mass_daughter1),
    mass_daughter2("mass_daughter2", "Mass of daughter particle 2", this, _mass_daughter2),
    l("l", "Total angular momentum", this, _l)
{
}

RooSillKstarPdf::RooSillKstarPdf(const RooSillKstarPdf& other, const char* name)
  : RooAbsPdf(other, name),
    x("x", this, other.x),
    mass("mass", this, other.mass),
    gamma("gamma", this, other.gamma),
    mass_daughter1("mass_daughter1", this, other.mass_daughter1),
    mass_daughter2("mass_daughter2", this, other.mass_daughter2),
    l("l", this, other.l)
{}

double RooSillKstarPdf::evaluate() const {
  
  double kstar = x;
  double E = std::sqrt(kstar * kstar + mass_daughter1 * mass_daughter1) + std::sqrt(kstar * kstar + mass_daughter2 * mass_daughter2);
  double M = mass;
  double G = gamma;
  double Eth = mass_daughter1 + mass_daughter2;

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
