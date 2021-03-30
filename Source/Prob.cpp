// This problem setup file initializes a uniform domain

#include <AMReX_Vector.H>
#include <AMReX_Particles.H>
#include "Prob.H"
#include "Constants.H"
#include "Gas.H"
#include "SR.H"
#include <math.h> 



using namespace amrex;
namespace GIdx = GasIdx;

//////////////////////////////////////////////////////////////////////
// Routine to initialize the gas state
//////////////////////////////////////////////////////////////////////
void init_gas(Gas& gas, 
	      Geometry const& geom,
	      ParmParse& pp) {

  // No-op to prevent warning
  (void) geom;

  // Extract the density, velocity, and magnetic field from input deck
  Real rho, chi, B_Phi;
  Vector<Real>  v(AMREX_SPACEDIM);
  pp.get("rho", rho);
  pp.get("chi", chi);
  pp.queryarr("v", v);
  pp.get("B_Phi", B_Phi);
  

   // Extract physical coordinates from AMREX Geometry
  const Real*dx = geom.CellSize();
  const Real*plo = geom.ProbLo();
  
  // Loop over boxes
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(gas,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    // Grab the box on which to operate
    const Box& bx = mfi.tilebox();

    // Get pointer to data
    Array4<Real> const& data = gas[mfi].array();

    // Initialize data
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

    // Getting coords to convert to polars
    #define PI 3.1415926 
    Real x = plo[0] + (i + 0.5)*dx[0];
    Real y = plo[1] + (j + 0.5)*dx[1];
    Real z = plo[2] + (k + 0.5)*dx[2];
   
    data(i,j,k,GIdx::TotalDensity) = rho;
    data(i,j,k,GIdx::IonDensity) = rho * chi;
    data(i,j,k,GIdx::Vx) = v[0];
    data(i,j,k,GIdx::Bx) = -B_Phi*y/sqrt(x*x + y*y);

 


#if AMREX_SPACEDIM > 1
	data(i,j,k,GIdx::Vy) = v[1];
	data(i,j,k,GIdx::By) = B_Phi*x/sqrt(x*x + y*y);
#endif
#if AMREX_SPACEDIM > 2
	data(i,j,k,GIdx::Vz) = v[2];
	data(i,j,k,GIdx::Bz) = 0;
#endif
    });

  }

}


//////////////////////////////////////////////////////////////////////
// Routine to initialize CR sources
//////////////////////////////////////////////////////////////////////
void init_sources(CRSourceContainer &sources,
		  ParmParse &pp,
		  ParmParse &pp_cr) {

  // Initialize CR sources; note that we initialize only one one
  // processor, then call redistribute to scatter the particles to the
  // right processor
  if (ParallelDescriptor::MyProc() == 
      ParallelDescriptor::IOProcessorNumber()) {

    // Read list of sources from parameter file
    Vector<Real> x, y, z, L, q, eMin, eMax;
    pp.queryarr("source_x", x);
    pp.queryarr("source_y", y);
    pp.queryarr("source_z", z);
    pp.queryarr("source_L", L);
    pp.queryarr("source_q", q);
    pp.queryarr("source_e_min", eMin);
    pp.queryarr("source_e_max", eMax);

    // Consistency check
    if (x.size() != y.size() ||
	x.size() != z.size() ||
	x.size() != L.size() ||
	x.size() != q.size() ||
	x.size() != eMin.size() ||
	x.size() != eMax.size()) {
      Abort("Found inconsistent number of sources!");
    }

    // Loop over sources
    for (Vector<Real>::size_type i=0; i<L.size(); i++) {

      Vector<Real> pos;
      pos.push_back(x[i]);
#if AMREX_SPACEDIM > 1
      pos.push_back(y[i]);
#endif
#if AMREX_SPACEDIM > 2
      pos.push_back(z[i]);
#endif

      // Create a source
      CRSource s(pos, L[i], eMin[i], eMax[i], q[i]);

      // Add to container at level 0, grid 0, tile 0
      std::pair<int,int> key {0,0};
      auto& particle_tile = sources.GetParticles(0)[key];
      particle_tile.push_back(s);
    }

  }

  // Distribute particles to processors
  sources.Redistribute();

}
