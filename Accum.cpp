#include "Accum.H"
#include "BCFill.H"
#include "Constants.H"
#include <cmath>
#include <sstream>
#include <cstring>

using namespace amrex;
using namespace std;

// Constructor
Accum::Accum(ParmParse &pp,
	     const amrex::Geometry &geom_,
	     const BoxArray &ba,
	     const DistributionMapping &dm) : geom(geom_)
{

  // Set up CR energy bins and associated names; note that internally
  // we use m_p c^2 as our energy unit, but we report take input and
  // output in units of GeV
  Real eMin, eMax;
  int nEnergy;
  pp.get("n_energy", nEnergy);
  pp.get("e_min", eMin);
  pp.get("e_max", eMax);
  TBins.resize(nEnergy+1);
  TBins[0] = eMin * units::GeV / constants::mp_c2;
  TBins[nEnergy] = eMax * units::GeV / constants::mp_c2;
  dlogT = log(eMax / eMin) / nEnergy;
  for (Vector<Real>::size_type i=1; i<TBins.size()-1; i++)
    TBins[i] = TBins[0] * exp(i * dlogT);

  // Allocate multifabs to hold accumulated quantities
  accOld = new MultiFab(ba, dm, 2*nEnergy, nGhost);
  accNew = new MultiFab(ba, dm, 2*nEnergy, nGhost);
  (*accOld) = 0.0;
  (*accNew) = 0.0;

  // Set up vector of strings to hold field names
  accNames.resize(2*nEnergy);
  for (Vector<Real>::size_type i=0; i<nEnergy; i++) {
    std::stringstream ss;
    ss << "nCR_" 
       << eMin * exp(i * dlogT)
       << "_to_";
    if (i != nEnergy - 1) ss << eMin * exp((i+1) * dlogT);
    else ss << eMax;
    accNames[i] = ss.str();
  }
  for (Vector<Real>::size_type i=0; i<nEnergy; i++) {
    std::stringstream ss;
    ss << "ECR_" 
       << eMin * exp(i * dlogT)
       << "_to_";
    if (i != nEnergy - 1) ss << eMin * exp((i+1) * dlogT);
    else ss << eMax;
    accNames[i+nEnergy] = ss.str();
  }

  // Set up our boundary conditions
  bcs.resize(nComp());
  for (int n=0; n<nComp(); n++) {
    for (int i=0; i<AMREX_SPACEDIM; i++) {
      if (geom.isPeriodic(i)) {
	bcs[n].setLo(i, BCType::int_dir);
	bcs[n].setHi(i, BCType::int_dir);
      } else {
	bcs[n].setLo(i, BCType::foextrap);
	bcs[n].setHi(i, BCType::foextrap);
      }
    }
  }
}


// BC fill routine
void Accum::BCFill() {

  // Fill interior and periodic ghost cells
  (*accOld).FillBoundary(geom.periodicity());

  // Fill physical boundary conditions if any exist
  if (!geom.isAllPeriodic()) {
    CpuBndryFuncFab bf(nullptr);
    PhysBCFunct<CpuBndryFuncFab> bcFiller(geom, bcs, bf);
    bcFiller(*accOld, 0, accOld->nComp(), accOld->nGrowVect(), 0.0, 0);
  }
}
