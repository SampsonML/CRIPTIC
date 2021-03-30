// This module defines a no-op boundary condition filling function;
// this is required by amrex, but never gets called because we do not
// allow Dirichlet BC's

#ifndef BCFILL_H
#define BCFILL_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

struct ExtBCFill
{
    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& /*iv*/, amrex::Array4<amrex::Real> const& /*data*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     amrex::GeometryData const& /*geom*/, const amrex::Real /*time*/,
                     const amrex::BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
        {
            // no op
        }
};

#endif
