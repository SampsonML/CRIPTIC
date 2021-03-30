// This class provides a mechanism to get runtime parameters for the
// microphysical propagation model into the update kernel. Its purpose
// in life is to read parameters from the criptic.in file, then
// package them in a GpuArray for export to the kernel. Which
// parameters are read depends on the microphysical propagration
// model, which is set at compile time, so there are lots of ifdef's
// here.

#ifndef _PROPMODEL_H_
#define _PROPMODEL_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Array.H>
#include "Definitions.H"

#if (PROPMODEL == POWERLAW)
#   define NPROP_PARAM   8
namespace PropModelIdx {
  const int kPar0 = 0;    // kPar = kPar0 * (p / m_p c)^kParIdx
  const int kParIdx = 1;  // kPar = kPar0 * (p / m_p c)^kParIdx
  const int kPerp0 = 2;   // kPerp = kPerp0 * (p / m_p c)^kPerpIdx
  const int kPerpIdx = 3; // kPerp = kPerp0 * (p / m_p c)^kPerpIdx
  const int kPP0 = 4;     // kPP = kPP0 * (p / m_p c)^kPPIdx
  const int kPPIdx = 5;   // kPP = kPP0 * (p / m_p c)^kPPIdx
  const int vStr0 = 6;    // vStr = vStr0 * (p / m_p c)^vStrIdx
  const int vStrIdx = 7;  // vStr = vStr0 * (p / m_p c)^vStrIdx
}
#endif

class PropModel {

public:

  // The constructor and destructor
  PropModel(amrex::ParmParse &pp);
  ~PropModel() {}

  // The export function
  const amrex::GpuArray<amrex::Real,NPROP_PARAM>&
  getParams() const { return params; }
	     
private:

#if (PROPMODEL == POWERLAW)
  amrex::GpuArray<amrex::Real,NPROP_PARAM> params;
#endif

};

#endif
// _PROPMODEL_H_
