// Implementation of the PropModel class

#include "PropModel.H"
#include "Constants.H"
#include "SR.H"
#include <cmath>

using namespace std;
using namespace amrex;
namespace PMIdx = PropModelIdx;

// Constructor
PropModel::PropModel(ParmParse &pp) {

#if (PROPMODEL == POWERLAW)
  // This simple model (inteded for code testing) just sets all the
  // coefficients to simple powerlaw functions of the CR momentum, of
  // the form
  //    kPar = kPar0 * (p/m_p c)^kParIdx
  // and similar for all other quantities
  pp.get("kPar0", params[PMIdx::kPar0]);
  params[PMIdx::kParIdx] = 0.0;
  pp.query("kParIdx", params[PMIdx::kParIdx]);
  pp.get("kPerp0", params[PMIdx::kPerp0]);
  params[PMIdx::kPerpIdx] = 0.0;
  pp.query("kPerpIdx", params[PMIdx::kPerpIdx]);
  pp.get("kPP0", params[PMIdx::kPP0]);
  params[PMIdx::kPPIdx] = 0.0;
  pp.query("kPPIdx", params[PMIdx::kPPIdx]);
  pp.get("vStr0", params[PMIdx::vStr0]);
  params[PMIdx::vStrIdx] = 0.0;
  pp.query("vStrIdx", params[PMIdx::vStrIdx]);
#endif

}
