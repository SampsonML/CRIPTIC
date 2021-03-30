// This header file contains physical constants, numerical constants,
// and unit conversions

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <AMReX_REAL.H>
#include "Definitions.H"

// Unit-dependent stuff
#if (UNITS == MKS)
// MKS units
#include <gsl/gsl_const_mksa.h>
namespace constants {
  static constexpr amrex::Real c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
  static constexpr amrex::Real mp = GSL_CONST_MKSA_MASS_PROTON;
  static constexpr amrex::Real me = GSL_CONST_MKSA_MASS_ELECTRON;
}
namespace units {
  static constexpr amrex::Real amu = GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS;
  static constexpr amrex::Real eV = GSL_CONST_MKSA_ELECTRON_VOLT;
  static constexpr amrex::Real pc = GSL_CONST_MKSA_PARSEC;
  static constexpr amrex::Real lyr = GSL_CONST_MKSA_LIGHT_YEAR;
  static constexpr amrex::Real mb = 1e-31;   // 1 millibarn in MKS units
}
#elif (UNITS == CGS)
// CGS units
#include <gsl/gsl_const_cgsm.h>
namespace constants {
  static constexpr amrex::Real c = GSL_CONST_CGSM_SPEED_OF_LIGHT;
  static constexpr amrex::Real mp = GSL_CONST_CGSM_MASS_PROTON;
  static constexpr amrex::Real me = GSL_CONST_CGSM_MASS_ELECTRON;
}
namespace units {
  static constexpr amrex::Real amu = GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS;
  static constexpr amrex::Real eV = GSL_CONST_CGSM_ELECTRON_VOLT;
  static constexpr amrex::Real pc = GSL_CONST_CGSM_PARSEC;
  static constexpr amrex::Real lyr = GSL_CONST_CGSM_LIGHT_YEAR;
  static constexpr amrex::Real mb = 1e-27;   // 1 millibarn in CGS units
}
#endif

// Stuff that is the same in CGS or MKS
namespace units {
  static constexpr amrex::Real GeV = 1e9*eV;
  static constexpr amrex::Real TeV = 1e12*eV;
  static constexpr amrex::Real PeV = 1e15*eV;
  static constexpr amrex::Real yr = lyr / constants::c;
  static constexpr amrex::Real kyr = 1e3*yr;
  static constexpr amrex::Real Myr = 1e6*yr;
}
namespace constants {
  static constexpr amrex::Real c2 = c*c;
  static constexpr amrex::Real mH = mp + me;
  static constexpr amrex::Real mp_c = mp * c;
  static constexpr amrex::Real mp_c2 = mp * c2;
  static constexpr amrex::Real MW_nuc_per_p =
    1.17;    // Nucleons/proton for MW ISM composition
  static constexpr amrex::Real MW_mass_per_H =
    1.4;     // MW ISM mass per H, expressed in units of m_H
  static constexpr amrex::Real MW_nuc_per_mass =
    MW_nuc_per_p * MW_mass_per_H / mH; // Nucleons / unit mass for MW composition
  static constexpr amrex::Real Tpth =
    0.2792 * units::GeV / mp_c2;  // Pion production threshold in units of m_p c^2
}

// Runge-Kutta-Fehlberg coefficients, copied from GSL
namespace rkf {
  static constexpr amrex::Real ah[] =
    { 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0 };
  static constexpr amrex::Real b3[] = { 3.0 / 32.0, 9.0 / 32.0 };
  static constexpr amrex::Real b4[] =
    { 1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0 };
  static constexpr amrex::Real b5[] =
    { 8341.0 / 4104.0, -32832.0 / 4104.0, 29440.0 / 4104.0, -845.0 / 4104.0 };
  static constexpr amrex::Real b6[] =
    { -6080.0 / 20520.0, 41040.0 / 20520.0, -28352.0 / 20520.0,
      9295.0 / 20520.0, -5643.0 / 20520.0
    };
  static constexpr amrex::Real c1 = 902880.0 / 7618050.0;
  static constexpr amrex::Real c3 = 3953664.0 / 7618050.0;
  static constexpr amrex::Real c4 = 3855735.0 / 7618050.0;
  static constexpr amrex::Real c5 = -1371249.0 / 7618050.0;
  static constexpr amrex::Real c6 = 277020.0 / 7618050.0;
  static constexpr amrex::Real ec[] = { 0.0,
    1.0 / 360.0,
    0.0,
    -128.0 / 4275.0,
    -2197.0 / 75240.0,
    1.0 / 50.0,
    2.0 / 55.0
  };
}

#endif
// _CONSTANTS_H_
