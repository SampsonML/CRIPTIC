// This class defines the background thermal gas

#ifndef _GAS_H_
#define _GAS_H_

#include <AMReX_Vector.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <cstring>
#include "BCFill.H"

// Indexing for gas quantities
namespace GasIdx {

  const int TotalDensity = 0;         // Total mass density
  const int IonDensity = 1;           // Ionized mass density
  const int Vx = 2;                   // X-velocity
#if AMREX_SPACEDIM > 1
  const int Vy = 3;                   // Y-velocity
#endif
#if AMREX_SPACEDIM > 2
  const int Vz = 4;                   // Z-velocity
#endif
  const int Bx = 2+AMREX_SPACEDIM;    // X-magnetic field
#if AMREX_SPACEDIM > 1
  const int By = 3+AMREX_SPACEDIM;    // Y-magnetic field
#endif
#if AMREX_SPACEDIM > 2
  const int Bz = 4+AMREX_SPACEDIM;    // Z-magnetic field
#endif
  
  // Total number of gas components
  const int nComp = 2 + 2*AMREX_SPACEDIM;

  // Variable names
  const amrex::Vector<std::string> varNames = { "total_density",
						"ion_density"
						, "x-vel"
#if AMREX_SPACEDIM > 1
						, "y-vel"
#endif
#if AMREX_SPACEDIM > 2
						, "z-vel"
#endif
						, "x-mag-field"
#if AMREX_SPACEDIM > 1
						, "y-mag-field"
#endif
#if AMREX_SPACEDIM > 2
						, "z-mag-field"
#endif
  };
}


// The gas class
class Gas : public amrex::MultiFab {

public:

  // Constructor
  Gas(const amrex::Geometry& geom_,
      const amrex::BoxArray& ba,
      const amrex::DistributionMapping &dm) :
    MultiFab(ba, dm, GasIdx::nComp, 2),
    geom(geom_)
  {
    // Set BCs
    bcs.resize(GasIdx::nComp);
    for (int n=0; n<GasIdx::nComp; n++) {
      for (int i=0; i<AMREX_SPACEDIM; i++) {
	if (geom.isPeriodic(i)) {
	  bcs[n].setLo(i, amrex::BCType::int_dir);
	  bcs[n].setHi(i, amrex::BCType::int_dir);
	} else {
	  bcs[n].setLo(i, amrex::BCType::foextrap);
	  bcs[n].setHi(i, amrex::BCType::foextrap);
	}
      }
    }
  }

  // Boundary condition fill routine
  void BCFill() {
    FillBoundary(geom.periodicity());  // Periodic BCs
    if (!geom.isAllPeriodic()) {
      amrex::CpuBndryFuncFab bf(nullptr);
      amrex::PhysBCFunct<amrex::CpuBndryFuncFab> bcFiller(geom, bcs, bf);
      bcFiller(*this, 0, this->nComp(), this->nGrowVect(), 0.0, 0);
    }
  }

private:
  
  const amrex::Geometry& geom;            // Geometry
  amrex::Vector<amrex::BCRec> bcs;        // Boundary condition filler
  const int nGhost = 2;                   // Number of ghost cells

};

#endif
// _GAS_H_
