// This class defines the data structure that holds CR packets. This
// is just a light wrapper over the AMReX Particle class that adds
// convenience functions for manipulating packets. It adds no extra
// data, so typecasting between Particles and CRPackets is safe.

#ifndef _CRPACKET_H_
#define _CRPACKET_H_

#include <AMReX_Particles.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParallelDescriptor.H>
#include "Constants.H"
#include "SR.H"

// Holder for indices of data
namespace CRPacketIdx {

  // Real attributes
  const int m = 0;      // mass in units of m_p
  const int p = 1;      // momentum in units of m_p c
  const int w = 2;      // statistical weight
  const int t = 3;      // current time
  const int tInj = 4;   // time of injection
  const int wInj = 5;   // weight at time of injection
  const int gr = 6;     // grammage traversed
  const int vT = 7;     // effective velocity in the T direction
  const int vN = 8;     // effective velocity in the N direction
  const int vB = 9;     // effective velocity in the B direction
  const int tStep = 10; // time at which current step should stop

  // Integer attributes
  const int src = 0;    // index of source
  const int Z = 1;      // charge (in units of e)

  // Number of integer and real attributes
  const int nReal = 11; // number of Real attributes
  const int nInt = 2;   // number of integer attributes
}

class CRPacket :
  public amrex::Particle<CRPacketIdx::nReal, CRPacketIdx::nInt> {

public:

  CRPacket() {
    id()  = NextID();
    cpu() = amrex::ParallelDescriptor::MyProc();
  }

  // Constructor from particle; provided so that we can static_cast
  // particles to CRPackets
  CRPacket(amrex::Particle<CRPacketIdx::nReal, CRPacketIdx::nInt> p) :
    amrex::Particle<CRPacketIdx::nReal, CRPacketIdx::nInt>(p) { }

  ~CRPacket() {}

  // Methods to derive quantities from particle data; here gamma =
  // Lorentz factor, T = kinetic energy, v = velocity
  amrex::Real gamma() const {
    amrex::Real m = rdata(CRPacketIdx::m);
    amrex::Real p = rdata(CRPacketIdx::p);
    return sr::gamma_from_p(m, p);
  }
  amrex::Real T() const {
    amrex::Real m = rdata(CRPacketIdx::m);
    amrex::Real p = rdata(CRPacketIdx::p);
    return sr::T_from_p(m, p);
  }
  amrex::Real v() const {
    amrex::Real m = rdata(CRPacketIdx::m);
    amrex::Real p = rdata(CRPacketIdx::p);
    return sr::v_from_p(m, p);
  }
};

#endif
// _CRPACKET_H_
