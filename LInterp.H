// This header module defines some convenience functions for doing
// linear interpolation on data arrays. In the functions defined here,
// we use the following naming conventions for arguments:
//
// a = Array4, array to be interpolated
// i = IntVect, index of lower left corner of cell containing
//     interpolation position
// f = pointer to Real[3], giving position within cell
// c = component number, for functions that expect a single component
// scomp = integer, starting component, for functions that expect
//     multiple compoments
// ncomp = integer, number of components, for functions that expect
//     multiple compoments
// dir = direction, for functions that evaluate the gradient along a
//     particular direction
//
// Function names follow the convention TYPE or TYPE_vec, where TYPE
// indicates the type of quantity returned
//
// lint = function returns a linear interpolation of a at the input
//        point
// partial = function returns the partial derivative of a at the input
//        point
// grad = function returns the gradient of a at the input point
//
// Functions that end in _vec operate on multiple components, while
// those that do not end in _vec operate on a single component.

#ifndef _LINTERP_H_
#define _LINTERP_H_

#include <AMReX_Array4.H>
#include <AMReX_IntVect.H>

// A function that turns a position expressed as 3 real numbers into a
// position expressed as an index within a cell defined by the centers
// of grid cells and a position relative to that cell
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void convertPos(amrex::RealVect const& p,
		amrex::RealVect const& pLo,
		amrex::RealVect const& dxInv,
		amrex::IntVect &i,
		amrex::RealVect &f) {
  for (int n=0; n<AMREX_SPACEDIM; n++) {
    amrex::Real tmp = (p[n] - pLo[n]) * dxInv[n] - 0.5;
    i[n] = (int) floor(tmp);
    f[n] = tmp - i[n];
  }
}


// Linear interpolation, one component version
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
lint(amrex::Array4<amrex::Real const> const& a,
     amrex::IntVect const& i,
     amrex::RealVect const& f,
     int const c = 0) {
  amrex::Real fx, fy, fz;
  amrex::Real q = 0.0;
  for (int k_=0; k_<=1; k_++) {
    if (k_ == 0) fz = 1.0 - f[2]; else fz = f[2];
    for (int j_=0; j_<=1; j_++) {
      if (j_ == 0) fy = 1.0 - f[1]; else fy = f[1];
      for (int i_=0; i_<=1; i_++) {
	if (i_ == 0) fx = 1.0 - f[0]; else fx = f[0];
	q += fx * fy * fz * a(i[0] + i_,
			      i[1] + j_,
			      i[2] + k_,
			      c);
      }
    }
  }
  return q;
}

// Linear interoplation on a vector quantity
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::RealVect
lint_vec(amrex::Array4<amrex::Real const> const& a,
	 amrex::IntVect const& i,
	 amrex::RealVect const& f,
	 int const scomp = 0) {
  amrex::Real fx, fy, fz;
  amrex::RealVect q;
  for (int k_=0; k_<=1; k_++) {
    if (k_ == 0) fz = 1.0 - f[2]; else fz = f[2];
    for (int j_=0; j_<=1; j_++) {
      if (j_ == 0) fy = 1.0 - f[1]; else fy = f[1];
      for (int i_=0; i_<=1; i_++) {
	if (i_ == 0) fx = 1.0 - f[0]; else fx = f[0];
	for (int c=0; c<3; c++) {
	  q[c] += fx * fy * fz * a(i[0] + i_,
				   i[1] + j_,
				   i[2] + k_,
				   c + scomp);
	}
      }
    }
  }
  return q;
}

// Partial derivative of the linearly-interpolated quantity --
// single component along a single direction
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
partial(amrex::Array4<amrex::Real const> const& a,
	amrex::IntVect const& i,
	amrex::RealVect const& f,
	int const dir,
	int const c = 0) {
  amrex::Real fx, fy, fz;
  amrex::Real q = 0.0;
  for (int k_=0; k_<=1; k_++) {
    if (dir != 2) {
      if (k_ == 0) fz = 1.0 - f[2]; else fz = f[2];
    } else {
      if (k_ == 0) fz = -1.0; else fz = 1.0;
    }
    for (int j_=0; j_<=1; j_++) {
      if (dir != 1) {
	if (j_ == 0) fy = 1.0 - f[1]; else fy = f[1];
      } else {
	if (j_ == 0) fy = -1.0; else fy = 1.0;
      }
      for (int i_=0; i_<=1; i_++) {
	if (dir != 0) {
	  if (i_ == 0) fx = 1.0 - f[0]; else fx = f[0];
	} else {
	  if (i_ == 0) fx = -1.0; else fx = 1.0;
	}
	q += fx * fy * fz * a(i[0] + i_,
			      i[1] + j_,
			      i[2] + k_,
			      c);
      }
    }
  }
  return q;
}

// Partial derivative of the linearly-interpolated quantity --
// multiple components along a single direction
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::RealVect
partial_vec(amrex::Array4<amrex::Real const> const& a,
	    amrex::IntVect const& i,
	    amrex::RealVect const& f,
	    int const dir,
	    int const scomp = 0) {
  amrex::Real fx, fy, fz;
  amrex::RealVect q;
  for (int k_=0; k_<=1; k_++) {
    if (dir != 2) {
      if (k_ == 0) fz = 1.0 - f[2]; else fz = f[2];
    } else {
      if (k_ == 0) fz = -1.0; else fz = 1.0;
    }
    for (int j_=0; j_<=1; j_++) {
      if (dir != 1) {
	if (j_ == 0) fy = 1.0 - f[1]; else fy = f[1];
      } else {
	if (j_ == 0) fy = -1.0; else fy = 1.0;
      }
      for (int i_=0; i_<=1; i_++) {
	if (dir != 0) {
	  if (i_ == 0) fx = 1.0 - f[0]; else fx = f[0];
	} else {
	  if (i_ == 0) fx = -1.0; else fx = 1.0;
	}
	for (int c=0; c<3; c++) {
	  q[c] += fx * fy * fz * a(i[0] + i_,
				   i[1] + j_,
				   i[2] + k_,
				   c + scomp);
	}
      }
    }
  }
  return q;
}

// Gradient -- partial derivative along each direction, for one
// component
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::RealVect
grad(amrex::Array4<amrex::Real const> const& a,
     amrex::IntVect const& i,
     amrex::RealVect const& f,
     int const c = 0) {
  amrex::Real fx, fy, fz;
  amrex::RealVect q;
  for (int d=0; d<3; d++) {
    for (int k_=0; k_<=1; k_++) {
      if (d != 2) {
	if (k_ == 0) fz = 1.0 - f[2]; else fz = f[2];
      } else {
	if (k_ == 0) fz = -1.0; else fz = 1.0;
      }
      for (int j_=0; j_<=1; j_++) {
	if (d != 1) {
	  if (j_ == 0) fy = 1.0 - f[1]; else fy = f[1];
	} else {
	  if (j_ == 0) fy = -1.0; else fy = 1.0;
	}
	for (int i_=0; i_<=1; i_++) {
	  if (d != 0) {
	    if (i_ == 0) fx = 1.0 - f[0]; else fx = f[0];
	  } else {
	    if (i_ == 0) fx = -1.0; else fx = 1.0;
	  }
	  q[d] += fx * fy * fz * a(i[0] + i_,
				   i[1] + j_,
				   i[2] + k_,
				   c);
	}
      }
    }
  }
  return q;
}

// Gradient over multiple components. The results are returned in a
// vector indexed so that d/dx_j q_i is stored in element j + 3*i;
// here i is the index of the component in the return vector, not in the
// underlying Array4, so, for example, if one sets scomp = 3 and
// ncomp = 5, then i runs from 0 to 2.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,9>
grad_vec(amrex::Array4<amrex::Real const> const& a,
	 amrex::IntVect const& i,
	 amrex::RealVect const& f,
	 int const scomp = 0) {
  amrex::Real fx, fy, fz;
  amrex::GpuArray<amrex::Real,9> q;
  for (int d=0; d<3; d++) {
    for (int k_=0; k_<=1; k_++) {
      if (d != 2) {
	if (k_ == 0) fz = 1.0 - f[2]; else fz = f[2];
      } else {
	if (k_ == 0) fz = -1.0; else fz = 1.0;
      }
      for (int j_=0; j_<=1; j_++) {
	if (d != 1) {
	  if (j_ == 0) fy = 1.0 - f[1]; else fy = f[1];
	} else {
	  if (j_ == 0) fy = -1.0; else fy = 1.0;
	}
	for (int i_=0; i_<=1; i_++) {
	  if (d != 0) {
	    if (i_ == 0) fx = 1.0 - f[0]; else fx = f[0];
	  } else {
	    if (i_ == 0) fx = -1.0; else fx = 1.0;
	  }
	  for (int c=0; c<3; c++) {
	    q[3*c+d] += fx * fy * fz * a(i[0] + i_,
					 i[1] + j_,
					 i[2] + k_,
					 c + scomp);
	  }
	}
      }
    }
  }
  return q;
}

#endif
// _LINTERP_H_
