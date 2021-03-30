// This defines the interface to the problem-specific initialization
// routine. Which implementation of this initialization is used is
// determined by the prob=... flag set at compile time.

#ifndef PROB_H_
#define PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include "CRSourceContainer.H"
#include "Gas.H"

void init_gas(Gas& gas, 
	      amrex::Geometry const& geom,
	      amrex::ParmParse& pp);

void init_sources(CRSourceContainer &sources,
		  amrex::ParmParse& pp,
		  amrex::ParmParse& pp_cr);

#endif
// _PROB_H_
