// This class is derived from AMReX ParticleContainer; its purpose is
// to hold CR sources.

#ifndef _CRSOURCECONTAINER_H_
#define _CRSOURCECONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include "CRSource.H"
#include "CRPacketContainer.H"
#include <cstring>
#if _cplusplus == 201703L
#   include <filesystem>
#endif

class CRSourceContainer :
  public amrex::ParticleContainer<CRSourceIdx::nReal, CRSourceIdx::nInt> {

public:

  // Constructor
  CRSourceContainer(amrex::ParmParse &pp,
		    const amrex::Geometry &geom,
		    const amrex::DistributionMapping &dm,
		    const amrex::BoxArray &ba) :
    ParticleContainer<5, 1>(geom, dm, ba) {
    pp.query("q_samp", q);
    pp.get("packet_rate", packetRate);
  };

  // Method to get souces ready to inject CRs; must be called whenever
  // a new source is added, or when the sampling index is changed
  void init();

  // Method to inject CR packets; here n = number of packets to
  // inject, t = start time for injection, dt = time interval over
  // which packets are injected, packets = the CR packet container to
  // which the Packets will be added
  void injectPackets(const amrex::Real t,
		     const amrex::Real dt,
		     CRPacketContainer &packets);

  // Methods to write out data on CR sources
  void writeSources(const std::string &dir) {
    WritePlotFile(dir, "CRSources", realVarNames, intVarNames);
  }
  void writeSourcesASCII(const std::string &dir) {
#if _cplusplus == 201703L
    std::path d(dir);
    WriteAsciiFile(d / "CRSources");
#else
    std::string d = dir + "/CRSources";
    WriteAsciiFile(d);
#endif
  }

private:

  // Names of CR source entries
  const amrex::Vector<std::string> realVarNames = {
    "nLum", "p0", "p1", "q", "m" };
  const amrex::Vector<std::string> intVarNames = { "Z" };

  // Private data
  // q = index used to sample momentum distributions
  // packetRate = # packets injected / time across all sources
  // wTot = total weight of all sources -- used to divide particles
  // between sources
  amrex::Real q = -1.0;
  amrex::Real packetRate;
  amrex::Real wTot;
};

#endif
// _CRSOURCECONTAINER_H_
