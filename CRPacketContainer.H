// This class is derived from AMReX ParticleContainer; it contains CR packets

#ifndef _CRPACKETCONTAINER_H_
#define _CRPACKETCONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include "Accum.H"
#include "Gas.H"
#include "CRPacket.H"
#include "PropModel.H"
#include <cstring>
#if _cplusplus == 201703L
#   include <filesystem>
#endif
#ifdef _OPENMP
#   include <omp.h>
#endif

class CRPacketContainer :
  public amrex::ParticleContainer<CRPacketIdx::nReal, CRPacketIdx::nInt> {

public:

  // Constructor; this just records some trivial information
  CRPacketContainer(amrex::ParmParse &pp,
		    const amrex::Geometry &geom,
		    const amrex::DistributionMapping &dm,
		    const amrex::BoxArray &ba,
		    const PropModel &prop_) :
    ParticleContainer<CRPacketIdx::nReal, CRPacketIdx::nInt>(geom, dm, ba),
    dxInv(m_gdb->Geom(0).InvCellSize()),
    pLo(m_gdb->Geom(0).ProbLo()),
    pDomain(m_gdb->Geom(0).ProbDomain()),
    volInv(dxInv.product()),
    dxInvMax(dxInv[dxInv.maxDir(false)]),
    domVol(m_gdb->Geom(0).ProbSize()),
    prop(prop_)
  {
    pp.query("cr_step_size", cStep);
    pp.query("err_tol", errTol);
    pp.query("c_glob", cGlob);
    pp.query("c_frac_min", cFracMin);
    pp.query("min_cr_wgt", minWgtFrac);
    pp.query("max_dt_incr", maxDtIncr);
#ifdef _OPENMP
    localAcc.resize(omp_get_max_threads());
    localNAvg.resize(omp_get_max_threads());
#endif
  };

  // Method to advance packets
  amrex::Real advancePackets(const amrex::Real &tStart,
			     const amrex::Real &tStop,
			     Gas &gas,
			     Accum &acc);

  // Methods to write out CR packets; note that internally we work in
  // units where CR masses are in proton masses, CR momenta are in
  // m_p * c, and CR energies are in m_p * c**2. However, for IO it is
  // more convenient to use more conventional units. We therefore
  // convert all CR masses to amu, and all momenta to amu * cm/s (or m/s
  // if we are using MKS units) before writing; we then convert back
  // to our internal units after writing
  void writePackets(const std::string &dir) {
    int lev = 0;
    // Convert mass to amu, momentum to amu * cm/s or amu * m/s
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::ParIter<CRPacketIdx::nReal, CRPacketIdx::nInt>
	   pti(*this, lev); pti.isValid(); ++pti) {
      auto& packets = pti.GetArrayOfStructs();
      const size_t np = packets.numParticles();
      ParticleType* pstruct = packets().dataPtr();
      amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (size_t n)
      {
	ParticleType &packet = pstruct[n];
	packet.rdata(CRPacketIdx::m) *= constants::mp / units::amu;
	packet.rdata(CRPacketIdx::p) *= constants::mp * constants::c
	  / units::amu;
      });
    }
    // Write data
#if 1
    WritePlotFile(dir, "CRPackets", writeRealComp, writeIntComp,
		  realVarNames, intVarNames);
#else
    WriteAsciiFile(dir + "/CRPackets");
#endif
    // Convert back
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::ParIter<CRPacketIdx::nReal, CRPacketIdx::nInt>
	   pti(*this, lev); pti.isValid(); ++pti) {
      auto& packets = pti.GetArrayOfStructs();
      const size_t np = packets.numParticles();
      ParticleType* pstruct = packets().dataPtr();
      amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (size_t n)
      {
	ParticleType &packet = pstruct[n];
	packet.rdata(CRPacketIdx::m) /= constants::mp / units::amu;
	packet.rdata(CRPacketIdx::p) /= constants::mp * constants::c
	  / units::amu;
      });
    } 
  }

private:

  // Internal data
  const amrex::RealVect dxInv, pLo;
  const amrex::RealBox pDomain;
  const amrex::Real volInv, dxInvMax, domVol;

  // Control parameters, used to set time steps
  amrex::Real cStep = 0.5;
  amrex::Real errTol = 1.0e-3;
  amrex::Real minWgtFrac = 1.0e-2;
  amrex::Real cGlob = 0.1;
  amrex::Real cFracMin = 1.0e-3;
  amrex::Real maxDtIncr = 1.1;

  // Pointer to the propagation model
  const PropModel &prop;

  // Names of CR packet entries
  const amrex::Vector<std::string> realVarNames = {
    "mass", "momentum", "weight", "time", "injection_time",
    "injection_weight", "grammage", "velocity_T",
    "velocity_N", "velocity_B", "tStep" };
  const amrex::Vector<std::string> intVarNames = { "source", "Z" };

  // Flags for which entries to write
  const amrex::Vector<int> writeRealComp = {
    1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0 };
  const amrex::Vector<int> writeIntComp = { 1, 1 };

  // Holders for thread-local quantites, required if running with openMP
#ifdef _OPENMP
  amrex::Vector<amrex::FArrayBox> localAcc;
  amrex::Vector<amrex::Vector<amrex::Real> > localNAvg;
#endif
};

#endif
// _CRPACKETCONTAINER_H_
