// Main function

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>
#include <AMReX_RealBox.H>
#include <AMReX_CoordSys.H>
#include <AMReX_PlotFileUtil.H>
#include <time.h>
#include "Accum.H"
#include "CRSourceContainer.H"
#include "CRPacketContainer.H"
#include "Gas.H"
#include "Prob.H"

using namespace amrex;

int main(int argc, char *argv[]) {

  // Start up amrex
  Initialize(argc, argv);
  amrex::InitRandom(clock());

  // Parse the input parameters file
  ParmParse pp;

  // Set up top-level grid data
  BoxArray ba;
  Geometry geom;
  {

    // Define the top-level logical grid
    Vector<int> ncells;
    pp.getarr("ncells", ncells);
    IntVect dom_lo(AMREX_D_DECL(0, 0, 0));
    IntVect dom_hi(AMREX_D_DECL(ncells[0]-1, ncells[1]-1, ncells[2]-1));
    Box domain(dom_lo, dom_hi);

    // Create the top-level box array from the domain
    ba.define(domain);
    int max_grid_size;
    if (!pp.query("max_grid_size", max_grid_size)) {
      max_grid_size = ncells[0];
      for (int i=1; i<AMREX_SPACEDIM; i++) 
	max_grid_size = max_grid_size < ncells[i] ? ncells[i] : max_grid_size;
    }

    // Define the physical domain
    Vector<Real> xlo, xhi;
    pp.getarr("xlo", xlo);
    pp.getarr("xhi", xhi);
    RealBox real_box(AMREX_D_DECL(xlo[0], xlo[1], xlo[2]),
		     AMREX_D_DECL(xhi[0], xhi[1], xhi[2]));

    // Set periodicity of domain; defaults to non-periodic
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    Vector<int> periodic;
    pp.queryarr("periodic", periodic);
    for (int i=0; i<AMREX_SPACEDIM; i++) {
      if (i >= periodic.size()) break;
      if (periodic[i] != 0) is_periodic[i] = 1;
    }

    // Set up the geometry object
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

  }

  // Distribute boxes across MPI ranks
  DistributionMapping dm(ba);

  // Build a particle container to hold CR sources
  ParmParse pp_cr("cr");
  CRSourceContainer sources(pp_cr, geom, dm, ba);

  // Build the accumulation object, which will store CR quantities
  Accum acc(pp_cr, geom, ba, dm);

  // Initialize the background gas state
  ParmParse pp_prob("prob");
  Gas gas(geom, ba, dm);
  init_gas(gas, geom, pp_prob);
  gas.BCFill();

  // Initialize CR sources
  init_sources(sources, pp_prob, pp_cr);
  sources.init();

  // Write initial gas and source state
  WriteSingleLevelPlotfile("init",
			   gas,
			   GasIdx::varNames,
			   geom,
			   0.0,
			   0);
  sources.writeSourcesASCII("init");

  // Initialize the object to hold CR propagation parameters
  ParmParse pp_crprop("crprop");
  PropModel prop(pp_crprop);
  
  // Build a container to hold CR packets
  CRPacketContainer packets(pp, geom, dm, ba, prop);

  // Write an initial plot file containing the CR energy density and
  // initial packet list
  {
    const std::string& pfname = amrex::Concatenate("plt", 0, 5);
    WriteSingleLevelPlotfile(pfname, acc.newData(), acc.varNames(),
			     geom, 0.0, 0);
    packets.WritePlotFile(pfname, "CRPackets");
  }
  
  // Get program control parameters
  Real dt;
  pp.get("dt_init", dt);
  Real maxDtIncr = 1.1;
  pp.query("max_dt_incr", maxDtIncr);
  int maxStep = 0;
  pp.query("max_step", maxStep);
  if (maxStep <= 0) maxStep = std::numeric_limits<int>::max();
  int plotInt = 0;
  pp.query("plot_int", plotInt);
  Real plotTime = 0.0;
  pp.query("plot_time", plotTime);
  Real maxTime = 0.0;
  pp.query("max_time", maxTime);

  // Main loop
  Real t = 0.0, setNextTime = 0.0, saveDt = 0.0;
  bool writePlot = false;
  Real nextPlotTime;
  if (plotTime != 0) nextPlotTime = plotTime;
  else nextPlotTime = 0.0;
  int plotCtr = 1;
  for (int step = 0; step < maxStep; ++step) {
    
    // Shift new data to old
    acc.shiftData();
    
    // Inject new packets
    sources.injectPackets(t, dt, packets);

    // Advance packets, getting back estimate of next time step
    Real dtNew = packets.advancePackets(t, dt, gas, acc);

    // Update times and time step, taking into account max simulation
    // time and plot file writing times
    if (setNextTime == 0.0) {
      t += dt;
    } else {
      // Ensure that we hit target time with no roundoff error
      t = setNextTime;
      setNextTime = 0.0;
    }
    if (dtNew > saveDt * maxDtIncr) dt = dt * maxDtIncr;
    else dt = dtNew;
    saveDt = dt;
    if (t + dt > maxTime && maxTime > 0.0) {
      dt = maxTime - t;
      setNextTime = maxTime;
    }
    if (t != nextPlotTime && t + dt > nextPlotTime) {
      dt = nextPlotTime - t;
      setNextTime = nextPlotTime;
    }
    if (t == maxTime) break;

    // Print status
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Step " << step+1 << ", t = " << t
		<< ", dt = " << dt
		<< ", total packets = "
		<< packets.NumberOfParticlesAtLevel(0)
		<< ", mean CR energy in volume during last step = "
		<< acc.totEn()
#if defined(MKS_UNITS)
		<< " J"
#else
		<< " erg"
#endif
		<< std::endl;
    }
    
    // Write an output if we have reached a step number or time when
    // we are supposed to do so
    writePlot = false;
    if (plotInt > 0) {
      if (step+1 % plotInt == 0) writePlot = true;
    }
    if (t == nextPlotTime) {
      writePlot = true;
      nextPlotTime += plotTime;
    }
    if (writePlot) {
      const std::string& pfname = amrex::Concatenate("plt", plotCtr, 5);
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Writing plot file " << pfname << std::endl;
      }
      WriteSingleLevelPlotfile(pfname, acc.newData(), acc.varNames(),
			       geom, t, step+1);
      packets.writePackets(pfname);
      plotCtr++;
    }
  }

  // Print final status
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Calculation done at t = " << t
	      << ", total packets = "
	      << packets.NumberOfParticlesAtLevel(0)
	      << ", mean CR energy in volume during last step = "
	      << acc.totEn()
#if defined(MKS_UNITS)
	      << " J"
#else
	      << " erg"
#endif
	      << std::endl;
  }
    
  // Write final plot file if we didn't just write one
  if (!writePlot) {
    const std::string& pfname = amrex::Concatenate("plt", plotCtr, 5);
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Writing plot file " << pfname << std::endl;
    }
    WriteSingleLevelPlotfile(pfname, acc.newData(), acc.varNames(),
			     geom, t, maxStep);
    packets.writePackets(pfname);
  }

  // Close down amrex
  Finalize();

  // End
  return 0;
}
