// This defines the Accum class. The purpose of this class is to
// accumulate CR quantities, for example CR energy densities, CR
// ionization rates, gamma ray luminosities, etc. It contains a
// multifab to store these quantities, along with auxiliary data
// describing how data are binned in CR energy, photon energy, etc.

#ifndef _ACCUM_H_
#define _ACCUM_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>
#include <AMReX_RealBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <cstring>

// Indexing for Accum quantities
namespace AccumIdx {
  const int num = 0;         // number density
  const int en = 1;          // energy density
}

class Accum {

public:

  // The constructor
  Accum(amrex::ParmParse &pp,
	const amrex::Geometry &geom_,
	const amrex::BoxArray &ba,
	const amrex::DistributionMapping &dm);

  // The destructor
  ~Accum() {
    delete accNew;
    delete accOld;
  }

  // Access to accumulated quantities
  amrex::MultiFab &newData() { return *accNew; }
  amrex::MultiFab &oldData() { return *accOld; }
  const amrex::MultiFab &newData() const { return *accNew; }
  const amrex::MultiFab &oldData() const { return *accOld; }

  // Boundary condition fill routine; this fills the boundary values
  // for the old data only
  void BCFill();

  // Shift new data into old data, and clear new data
  void shiftData() {
    amrex::MultiFab *tmp = accOld;
    accOld = accNew;
    accNew = tmp;
    (*accNew) = 0.0;
  }

  // Return volume-averaged and summed quantities
  amrex::Real nAvg(int c) const {
    return accNew->norm1(c) / accNew->boxArray().numPts();
  }
  amrex::Vector<amrex::Real> nAvg() const {
    amrex::Vector<int> comps(nT());
    for (int i=0; i<nT(); i++) comps[i] = i;
    amrex::Vector<amrex::Real> n = accNew->norm1(comps);
    for (int i=0; i<nT(); i++) n[i] /= accNew->boxArray().numPts();
    return n;
  }
  amrex::Real enAvg(int c) const {
    return accNew->norm1(c+nT()) / accNew->boxArray().numPts();
  }
  amrex::Vector<amrex::Real> enAvg() const {
    amrex::Vector<int> comps(nT());
    for (int i=0; i<nT(); i++) comps[i] = i + nT();
    amrex::Vector<amrex::Real> en = accNew->norm1(comps);
    for (int i=0; i<nT(); i++) en[i] /= accNew->boxArray().numPts();
    return en;
  }
  // Total energy in the computational volume, either in a single
  // energy bin (for c != -1) or summed over all energy bins (for c =
  // -1)
  amrex::Real totEn(int c = -1) const {
    if (c != -1) return enAvg(c) * geom.ProbSize();
    else {
      amrex::Real en = 0.0;
      for (int comp=0; comp<nT(); comp++) en += enAvg(comp);
      return en * geom.ProbSize();
    }
  }

  // Access to metadata
  int ngrow() const { return nGhost; }
  amrex::Real get_dlogT() const { return dlogT; }
  amrex::Real get_Tmin() const { return TBins[0]; }
  amrex::Real get_Tmax() const { return TBins[TBins.size()-1]; }
  const amrex::Vector<amrex::Real> &get_TBins() const { return TBins; }
  amrex::Long nT() const { return TBins.size()-1; }
  amrex::Long nComp() const { return accOld->nComp(); }
  const amrex::Vector<std::string> &varNames() const { return accNames; }

private:

  // Private data
  const amrex::Geometry& geom;          // Geometry
  amrex::Vector<amrex::BCRec> bcs;      // Boundary condition filler
  amrex::Real dlogT;                    // Energy bin size
  amrex::Vector<amrex::Real> TBins;     // CR energy bin edges
  amrex::Vector<std::string> accNames;  // Names of accumulated quantities
  amrex::MultiFab *accNew;              // Holder for quantities at new time
  amrex::MultiFab *accOld;              // Holder for quantities at old time
  const int nGhost = 2;                 // Number of ghost cells
};

#endif
// _ACCUM_H_
