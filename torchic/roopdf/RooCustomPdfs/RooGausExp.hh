#pragma once

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

//#define _AN_INT_

class RooGausExp : public RooAbsPdf {
public:
  RooGausExp() {} ;
  RooGausExp(const char *name, const char *title,
          RooAbsReal& _x,
          RooAbsReal& _mu,
          RooAbsReal& _sig,
          RooAbsReal& _tau);
  RooGausExp(const RooGausExp& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooGausExp(*this,newname); }
  inline virtual ~RooGausExp() { }
  
#ifdef _AN_INT_
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* r=0) const;
  virtual Double_t analyticalIntegral(Int_t code,const char* rangeName=0) const;
#endif
protected:
  double IntExp(double x,double tau) const;
  double IntGaus(double x) const;
  
  RooRealProxy x ;
  RooRealProxy mu ;
  RooRealProxy sig ;
  RooRealProxy tau ;
  
  Double_t evaluate() const ;
  
private:
  
  ClassDef(RooGausExp,1) // Your description goes here...
};

///////////////////////////// Class Implementation /////////////////////////////

#include <Riostream.h>
#include <TMath.h>

ClassImp(RooGausExp)

RooGausExp::RooGausExp(const char *name, const char *title,
                 RooAbsReal& _x,
                 RooAbsReal& _mu,
                 RooAbsReal& _sig,
                 RooAbsReal& _tau)
: RooAbsPdf(name,title)
, x("x","x",this,_x)
, mu("mu","mu",this,_mu)
, sig("sig","sig",this,_sig)
, tau("tau","tau",this,_tau) {
}


RooGausExp::RooGausExp(const RooGausExp& other, const char* name)
: RooAbsPdf(other,name)
, x("x",this,other.x)
, mu("mu",this,other.mu)
, sig("sig",this,other.sig)
, tau("tau",this,other.tau) {
}

Double_t RooGausExp::evaluate() const {
  double u = (x - mu) / sig;
  double abstau = tau;
  if (tau < 0) {
    u = -u;
    abstau = -tau;
  }
  if (u <= abstau)
    return TMath::Exp(-u * u * 0.5);
  else
    return TMath::Exp(-abstau * (u - 0.5 * abstau));
}

#ifdef _AN_INT_
Int_t RooGausExp::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char*) const {
  if (matchArgs(allVars,analVars,x)) return 1 ;
  return 0 ;
}

Double_t RooGausExp::analyticalIntegral(Int_t code, const char* r) const
{
  double umin = (x.min(r) - mu) / sig;
  double umax = (x.max(r) - mu) / sig;
  R__ASSERT(code==1);
  double integral = 0.;
  
  double abstau = tau;
  if (tau < 0) {
    double a = umin;
    umin = -umax;
    umax = -a;
    abstau = -tau;
  }
  
  if (umin < abstau) {
    integral -= IntGaus(umin);
    if (umax <= abstau)
      integral += IntGaus(umax);
    else
      integral += IntGaus(abstau) - IntExp(abstau,abstau) + IntExp(umax, abstau);
  } else
    integral = - IntExp(umin,abstau) + IntExp(umax, abstau);
  
  return sig * integral;
}
#endif

double RooGausExp::IntExp(double x,double tau) const {
  return -1. * TMath::Exp(tau * (0.5 * tau - x)) / tau;
}

double RooGausExp::IntGaus(double x) const {
  static const double rootPiBy2 = TMath::Sqrt(TMath::PiOver2());
  return rootPiBy2 * (TMath::Erf(x / TMath::Sqrt2()));
}
