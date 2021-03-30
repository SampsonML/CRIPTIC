// This header file contains convenience functions to do special
// relativistic transformations between mass, momentum, energy,
// etc. In the following functions, variable names correspond to
// physical quantities as follows:
//    m = rest mass in units of m_p
//    p = momentum in units of m_p c
//    T = kinetic energy in units of m_p c^2
//    v = velocity in units of c
//    gamma = Lorentz factor
// All these functions use appropriate series expansions to ensure
// reasonable behavior in both the non-relativistic and
// ultra-relativistic limits. Errors should be <~ 10^-6 regardless of
// the value of gamma.

#ifndef _SR_H_
#define _SR_H_

#include <cmath>
#include <AMReX_Vector.H>

namespace sr {

  // Scalar versions
  constexpr amrex::Real gamma_from_p(amrex::Real m, amrex::Real p) {
    amrex::Real pm = p / m;
    if (pm < 1.0e-6)
      return 1.0 + pm*pm / 2.0;
    else
      return std::sqrt(1.0 + pm*pm);
  }
  
  constexpr amrex::Real gamma_from_T(amrex::Real m, amrex::Real T) {
    return 1.0 + T/m;
  }
  
  constexpr amrex::Real p_from_gamma(amrex::Real m, amrex::Real gamma) {
    amrex::Real g = gamma-1.0;
    if (g < 1.0e-6)
      return m * (std::sqrt(2.0)*g + (g/2.0)*std::sqrt(g/2.0));
    else
      return m * std::sqrt(gamma*gamma - 1.0);
  }
  
  constexpr amrex::Real p_from_T(amrex::Real m, amrex::Real T) {
    amrex::Real Tm = T / m;
    if (Tm < 1.0e-6)
      return std::sqrt(2.0*Tm) + m * 
	(Tm/2.0) * std::sqrt(Tm/2.0);
    else
      return m * std::sqrt(Tm * (2.0+Tm));
  }
  
  constexpr amrex::Real T_from_gamma(amrex::Real m, amrex::Real gamma) {
    return m * (gamma-1.0);
  }
  
  constexpr amrex::Real T_from_p(amrex::Real m, amrex::Real p) {
    amrex::Real pm = p / m;
    if (pm < 1.0e-6)
      return p*p / (2.0*m) * (1.0 - pm*pm/4.0);
    else
      return m * (std::sqrt(1.0 + pm*pm) - 1.0);
  }
  
  constexpr amrex::Real dT_dp(amrex::Real m, amrex::Real p) {
    amrex::Real pm = p / m;
    if (pm < 1.0e-6)
      return pm * (1.0 - pm*pm/2.0);
    else
      return pm / std::sqrt(1.0 + pm*pm);
  }
  
  constexpr amrex::Real v_from_p(amrex::Real m, amrex::Real p) {
    amrex::Real pm = p / m;
    if (pm < 1.0e-6)
      return pm - 0.5 / (pm*pm*pm);
    else
      return 1.0 / std::sqrt(1.0 + 1.0/(pm*pm));
  }
}

#endif
// _SR_H_
