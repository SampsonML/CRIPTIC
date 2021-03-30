// This class defines cosmic ray sources. It is built on top of the
// AMReX Particle class, and provides methods for injecting CR
// packets. CR sources emit CRs are a rate per unit time per unit
// momentum given by
// d^2n / dt dp = A p^q,
// where p = momentum in units of mp * c, mp = proton mass, and A and
// q are free parameters; A has units of 1/time. The spectrum applies
// over a momentum range from p0 to p1 (again in units of mp * c).

#ifndef _CRSOURCE_H_
#define _CRSOURCE_H_

#include <AMReX_Particles.H>
#include <cmath>
#include "Constants.H"
#include "CRPacket.H"
#include "SR.H"
#include <gsl/gsl_sf.h>

// Holder for indices of data
namespace CRSourceIdx {
  const int m = 0;      // mass / m_p of CRs produced by this source
  const int A = 1;      // coefficient of CR production rate
  const int p0 = 2;     // minimum CR momentum in units of m_p * c
  const int p1 = 3;     // maximum CR momentum in units of m_p * c
  const int q = 4;      // index of momentum distribution, dn/dp ~ p^q
  const int Z = 0;      // charge of CRs produced by this source
  const int nReal = 5;  // number of real attributes
  const int nInt = 1;   // number of integer attributes
}

class CRSource :
  public amrex::Particle<CRSourceIdx::nReal, CRSourceIdx::nInt> {

public:

  // Constructor; here x = position, lum = luminosity (defined as
  // energy emitted per unit time in erg/s or J/s), T0_GeV and T1_GeV
  // = minimum and maximum CR energy in units of GeV, q = index of
  // momentum distribution, m = mass of CR particles in proton masses,
  // Z = charge of CR particles in units of elementary charge
  CRSource(amrex::Vector<amrex::Real> x,
	   amrex::Real lum_metric,
	   amrex::Real T0_GeV,
	   amrex::Real T1_GeV,
	   amrex::Real q,
	   amrex::Real m = 1.0,
	   int Z = 1);

  // Constructor from particle; provided so that we can static_cast
  // particles to CR sources
  CRSource(amrex::Particle<CRSourceIdx::nReal, CRSourceIdx::nInt> p) :
    amrex::Particle<CRSourceIdx::nReal, CRSourceIdx::nInt>(p) { }

  ~CRSource() {}

  // nFac returns 1 \int_{p0}^{p1} p^q dp. Thus the
  // number of CRs emitter per unit time is A * nFac(). This function
  // comes in three versions -- one that takes p0, p1, and q as
  // arguments, one that takes only q and uses the p0 and p1
  // associated with this object, and one that takes no arguments and
  // uses the p0, p1, and q for this object
  amrex::Real nFac(amrex::Real p0, amrex::Real p1, amrex::Real q) const {
    if (q != -1.0) {
      return (std::pow(p1, q+1) - std::pow(p0, q+1)) / (q + 1);
    } else {
      return std::log(p1 / p0);
    }
  }
  amrex::Real nFac(amrex::Real q) const {
    return nFac(rdata(CRSourceIdx::p0), rdata(CRSourceIdx::p1), q);
  }
  amrex::Real nFac() const {
    return nFac(rdata(CRSourceIdx::p0), rdata(CRSourceIdx::p1),
		rdata(CRSourceIdx::q));
  }

  // pBar returns \int_{p0}^{p1} p (p^q dp, which is the
  // mean CR momentum. Thus the total scalar momentum injected into
  // CRs per unit time is A * pBar(). As with nFac, this comes in
  // three flavors.
  amrex::Real pBar(amrex::Real p0, amrex::Real p1, amrex::Real q) const {
    if (q != -2.0) {
      return (std::pow(p1, q+2) - std::pow(p0, q+2)) / (q + 2);
    } else {
      return std::log(p1 / p0);
    }
  }
  amrex::Real pBar(amrex::Real q) const {
    return pBar(rdata(CRSourceIdx::p0), rdata(CRSourceIdx::p1), q);
  }
  amrex::Real pBar() const {
    return pBar(rdata(CRSourceIdx::p0), rdata(CRSourceIdx::p1),
		rdata(CRSourceIdx::q));
  }

  // TBar returns \int_{p0}^{p1} T(p) p^q dp, which is
  // the mean CR kinetic energy. Thus the CR energy injected per unit
  // time is A * TBar(). This has the same three versions as
  // nFac and pBar, but note that even the version that takes p0, p1,
  // and q as arguments still makes use of the rest mass m.
  amrex::Real TBar(amrex::Real p0, amrex::Real p1, amrex::Real q) const {
    amrex::Real m = rdata(CRSourceIdx::m);
    amrex::Real m2 = m*m, p02 = p0*p0, p12 = p1*p1;
    if (q != -2.0) {
      amrex::Real a = -1-q/2;
      amrex::Real b = (1+q)/2;
      amrex::Real term1 = m/(1+q) * (std::pow(p0, 1+q) - std::pow(p1, 1+q));
      amrex::Real term2 = 0.5 * std::pow(m, 2+q) * gsl_sf_beta(a, b) *
	( gsl_sf_beta_inc( a, b, 1/(1+p02/m2) ) -
	  gsl_sf_beta_inc( a, b, 1/(1+p12/m2) ) );
      return term1 + term2;
    } else {
      return m * (1/p1 - 1/p0) +
	std::sqrt(1 + m2/p02) - std::sqrt(1 + m2/p12) +
	std::asinh(p1/m) - std::asinh(p0/m);
    } 
  }
  amrex::Real TBar(amrex::Real q) const {
    return TBar(rdata(CRSourceIdx::p0), rdata(CRSourceIdx::p1), q);
  }
  amrex::Real TBar() const {
    return TBar(rdata(CRSourceIdx::p0), rdata(CRSourceIdx::p1),
		rdata(CRSourceIdx::q));
  }
  
  // Convenience methods to get the rate per unit time at which CRs
  // are produced, and the linear momentum and energy per unit time
  // they carry
  amrex::Real nLum() const { return rdata(CRSourceIdx::A) * nFac(); }
  amrex::Real pLum() const { return rdata(CRSourceIdx::A) * pBar(); }
  amrex::Real lum() const { return rdata(CRSourceIdx::A) * TBar(); }

  // Method to draw packets; here n = number of packets to draw, t =
  // start time of injection, dt = time interval of injection, qSamp =
  // index of the sampling distribution (packets distributed as dn/dp
  // ~ p^qSamp
  amrex::Vector<CRPacket> drawPackets(const int n,
				      const amrex::Real t,
				      const amrex::Real dt,
				      const amrex::Real qSamp) const;

};
#endif
// _CRSOURCE_H_
