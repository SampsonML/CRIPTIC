// This module defines standalone functions that act as a kernel for
// propagating CRs. They are defined here as standalone functions so
// that they can put on GPU if desired.

#ifndef _PROPKERNEL_H_
#define _PROPKERNEL_H_

#include <AMReX_Particles.H>
#include <AMReX_IntVect.H>
#include <AMReX_RealVect.H>
#include <AMReX_FArrayBox.H>
#include "Constants.H"
#include "CRPacket.H"
#include "Gas.H"
#include "PropModel.H"
#include "SR.H"
#include "LInterp.H"
#include "Definitions.H"

// A function that takes as input a position described by an i, f
// pair, together with a magnetic field array, and returns the
// tangent, normal, and binormal vectors at the specified point
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void computeTNB(amrex::RealVect const& p,
		amrex::RealVect const& pLo,
		amrex::RealVect const& dxInv,
		amrex::Array4<amrex::Real const> const& gas,
		amrex::RealVect& eT,
		amrex::RealVect& eN,
		amrex::RealVect& eB) {

  using amrex::Real;
  using amrex::IntVect;
  using amrex::RealVect;
  namespace GIdx = GasIdx;

  // Step 1: interpolate B field target point, and get tangent vector
  // there
  IntVect ix;
  RealVect f;
  convertPos(p, pLo, dxInv, ix, f);
  RealVect B = lint_vec(gas, ix, f, GIdx::Bx);
  Real BMag = B.vectorLength();
  eT = B / BMag;

  // Step 2: compute gradient of B; note that, for the object we
  // construct, element 3*i + j contains the quantity dB_i / dx_j
  amrex::GpuArray<Real,9> gradB = grad_vec(gas, ix, f, GIdx::Bx);
  Real gradBMax = 0.0;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      gradB[3*i+j] *= dxInv[j];
      if (std::fabs(gradB[3*i+j]) > gradBMax)
	gradBMax = std::fabs(gradB[3*i+j]);
    }
  }

  // Step 3: construct e_N, using the relation
  //    e_N ~ e_T,j d(e_T,j) / dx_i
  //        ~ e_T,j ( |B| dB_i/dx_j - B_i B_k dB_k/dx_j )
  if (gradBMax > 1.0e3 * std::numeric_limits<double>::epsilon()) {

    // Normal case where gradB is non-zero

    // Compute B_k dB_k / dx_j = 2 grad |B|
    RealVect twoGradB;
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	twoGradB[j] += B[k] * gradB[j + 3*k];
      }
    }
    
    // Compute the direction of e_N
    for (int i=0; i<3; i++) {
      eN[i] = 0.0;
      for (int j=0; j<3; j++) {
	eN[i] += eT[j] * (BMag * gradB[3*i+j] - B[i]*twoGradB[j]);
      }
    }

    // Normalize
    eN /= eN.vectorLength();

  } else {

    // Special gas where grad B = 0, indicating that the field is
    // straight. In this case the normal and binormal vectors are
    // degenerate, and can be rotated arbitrarily around the tangent
    // vector. We break the degeneracy by making the normal vector lie
    // in the xy plane (or in the yz plane in case e_T = \hat{z}
    if (eT[2] < 1.0) {
      eN[0] = eT[1];
      eN[1] = -eT[0];
      eN[2] = 0.0;
    } else {
      eN[0] = 1.0;
      eN[1] = 0.0;
      eN[2] = 0.0;
    }
  }

  // Step 5: construct the binormal vector as eB = eT x eN
  eB[0] = eT[1] * eN[2] - eT[2] * eN[1];
  eB[1] = eT[2] * eN[0] - eT[0] * eN[2];
  eB[2] = eT[0] * eN[1] - eT[1] * eN[0];
}


// A function that advances a CR packet through a time step of fixed
// size, using error control on the deviation from field line
// following, and domain checking to prevent going outside the region
// covered by the gas data. On return, the particle position and time
// (but NOT other properties) will be updated. We also return the time
// step completed.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
movePacket(amrex::Real const dt,
	   amrex::Real const errTol,
	   amrex::RealVect const& pLo,
	   amrex::RealVect const& dxInv,
	   amrex::Real const dxInvMax,
	   amrex::RealBox const& rbx,
	   amrex::Array4<amrex::Real const> const& gas,
	   amrex::RealVect const &va,
	   amrex::Particle<CRPacketIdx::nReal, CRPacketIdx::nInt> &packet) {
  
  using amrex::Real;
  using amrex::RealVect;
  using amrex::IntVect;
  namespace PIdx = CRPacketIdx;

  // Aliases
  Real& t = packet.rdata(PIdx::t);
  const Real &vT = packet.rdata(PIdx::vT);
  const Real &vN = packet.rdata(PIdx::vN);
  const Real &vB = packet.rdata(PIdx::vB);

  // Temporary storage
  RealVect x, eT, eN, eB;
  
  // Save starting position
  const RealVect xInit = packet.pos();

  // Set number of sub-steps, enforcing the condition that each
  // substep can move the particle no more than one cell size -- this
  // ensures that, if we start inside the valid box, and we have at
  // least one ghost cell of gas data, we won't leave the region
  // covered by the gas data
  int nStep = 1;
  Real vaMag = va.vectorLength();
  Real vTNBMag = std::sqrt(vT*vT + vN*vN + vB*vB);
  while (dt * dxInvMax * (vaMag + vTNBMag) / nStep > 1.0) nStep *= 2;
#if 0
  std::cout << "movePacket called with dt = " << dt
	    << ", nStep = " << nStep
	    << std::endl;
#endif
  
  // Try to update position
  bool success = false;
  while (!success) {

    // Set substep size for this attempt
    Real dtStep = dt / nStep;
    success = true;

    // Advance through substeps; at each substep, we recompute the TNB
    // basis at that position, then move the particle
    Real err2 = 0.0;
    for (int s=0; s<nStep; s++) {

      // Check if we have moved outside the valid box, and, if so,
      // update the time to as far as we advanced, then return
      if (s > 0) {
	if (!rbx.contains(packet.pos())) {
	  t += s*dtStep;
	  return s*dtStep;
	}
      }
      
      // Stage 1
      computeTNB(xInit, pLo, dxInv, gas, eT, eN, eB);
      RealVect k1 = va + eT*vT + eB*vB + eN*vN;
      x = xInit + dtStep * rkf::ah[0] * k1;

      // Stage 2
      computeTNB(x, pLo, dxInv, gas, eT, eN, eB);
      RealVect k2 = va + eT*vT + eB*vB + eN*vN;
      x = xInit + dtStep * (rkf::b3[0] * k1 + rkf::b3[1] * k2);

      // Stage 3
      computeTNB(x, pLo, dxInv, gas, eT, eN, eB);
      RealVect k3 = va + eT*vT + eB*vB + eN*vN;
      x = xInit + dtStep *
	(rkf::b4[0] * k1 + rkf::b4[1] * k2 + rkf::b4[2] * k3);

      // Stage 4
      computeTNB(x, pLo, dxInv, gas, eT, eN, eB);
      RealVect k4 = va + eT*vT + eB*vB + eN*vN;
      x = xInit + dtStep *
	(rkf::b5[0] * k1 + rkf::b5[1] * k2 + rkf::b5[2] * k3 +
	 rkf::b5[3] * k4);

      // Stage 5
      computeTNB(x, pLo, dxInv, gas, eT, eN, eB);
      RealVect k5 = va + eT*vT + eB*vB + eN*vN;
      x = xInit + dtStep *
	(rkf::b6[0] * k1 + rkf::b6[1] * k2 + rkf::b6[2] * k3 +
	 rkf::b6[3] * k4 + rkf::b6[4] * k5);

      // Stage 6
      computeTNB(x, pLo, dxInv, gas, eT, eN, eB);
      RealVect k6 = va + eT*vT + eB*vB + eN*vN;
      for (int n=0; n<3; n++)
	packet.pos(n) += dtStep *
	  (rkf::c1 * k1[n] + rkf::c3 * k3[n] + rkf::c4 * k4[n] +
	   rkf::c5 * k5[n] + rkf::c6 * k6[n]);

      // Compute error estimate on this step, and add to total
      RealVect err = dtStep * dxInvMax *
	(rkf::ec[1] * k1 + rkf::ec[3] * k3 +
	 rkf::ec[4] * k4 + rkf::ec[5] * k5 +
	 rkf::ec[6] * k6);
      err2 += err.radSquared();

      // Check if we have exceeded the error budget; if so, halve
      // the step, restore the original position, and try again
      if (err2 > errTol*errTol) {
	for (int n=0; n<3; n++) packet.pos(n) = xInit[n];
	nStep *= 2;
	success = false;
	break;
      }
    }
  }

  // If we are here, we completed the full time step, so update the
  // time and return the full time step
  t += dt;
  return dt;
}

// A function that advances a single packet to a specified time,
// using time step control, and accumulating data on CR number
// density, energy density, etc. On return, the packet is updated to
// the target time, or to a time less than that if it exits the valid
// region. We indicate a return due to leaving the valid region by a
// return value of 1, and a return for any other reason by a return
// value of 0
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int
advancePacket(amrex::Real const& tStop,
	      amrex::Real const& dtGlob,
	      amrex::Real &dtPartMin,
	      amrex::Real const& cStep,
	      amrex::Real const& errTol,
	      amrex::RealVect const& pLo,
	      amrex::RealVect const& dxInv,
	      amrex::Real const& dxInvMax,
	      amrex::Real const& volInv,
	      amrex::GpuArray<amrex::Real,NPROP_PARAM> const& propParam,
	      amrex::Array4<amrex::Real const> const& gas,
	      int const& nT,
	      amrex::Real const& Tmin,
	      amrex::Real const& dlogT,
	      amrex::Real const& minWgtFrac,
	      amrex::RealBox const& rbx,
	      amrex::Array4<amrex::Real const> const& accOld,
	      amrex::Array4<amrex::Real> const& accNew,
	      amrex::Particle<CRPacketIdx::nReal, CRPacketIdx::nInt>& packet) {
  
  using amrex::Real;
  using amrex::RealVect;
  using amrex::IntVect;
  namespace GIdx = GasIdx;
  namespace PIdx = CRPacketIdx;
  namespace PMIdx = PropModelIdx;
  namespace AIdx = AccumIdx;

#if 0
  std::cout
#ifdef _OPENMP
    << "Thread " << omp_get_thread_num() << ": "
#endif
    << "PK start for packet " << packet.id()
    << " from time " << packet.rdata(PIdx::t)
    << " to time " << tStop
    << ", dtGlob = " << dtGlob
    << ", starting position = "
    << packet.pos()
    << ", tile boundaries = " << rbx
    << std::endl;
#endif
  
  // Convenient aliases
  Real& t = packet.rdata(PIdx::t);
  Real& p = packet.rdata(PIdx::p);
  Real& w = packet.rdata(PIdx::w);
  Real& gr = packet.rdata(PIdx::gr);
  Real& vT = packet.rdata(PIdx::vT);
  Real& vN = packet.rdata(PIdx::vN);
  Real& vB = packet.rdata(PIdx::vB);
  Real& tStep = packet.rdata(PIdx::tStep);
  const Real& m = packet.rdata(PIdx::m);
  const Real& wInj = packet.rdata(PIdx::wInj);
  const Real eps = std::numeric_limits<Real>::epsilon();
  const Real Tmax = Tmin * std::exp(nT * dlogT);

  // Loop until advance is complete
  while (true) {

#if 0
    if (!rbx.contains(packet.pos(), 1.0e-3)) {
      std::cout << "containment error for packet " << packet.id()
		<< " at time " << packet.rdata(PIdx::t)
		<< ", position = "
		<< packet.pos()
		<< ", tile boundaries = " << rbx
		<< std::endl;
    }
#endif

    // Step 1: get packet index and location in position and energy
    // space
    IntVect iAcc;
    {
      RealVect iTmp = (packet.pos() - pLo) * dxInv;
      for (int n=0; n<3; n++) iAcc[n] = (int) iTmp[n];
    }
    Real T = sr::T_from_p(m, p);
    int accIdx;
    Real accWgt;
    {
      Real tmp = std::log(T / Tmin) / dlogT;
      accIdx = (int) floor(tmp);
      accWgt = tmp - accIdx;
    }

    // Step 2: use the propagation model to compute the diffusion
    // coefficients and streaming speed for the current packet position
    // and properties; note: vStr > 0 corresponds to streaming in the
    // direciton aligned with the magnetic field vStr < corresponds to
    // streaming against the magnetic field direction
    Real kPar, kPerp, kPP, vStr;
#if (PROPMODEL == POWERLAW)
    kPar = propParam[PMIdx::kPar0] * std::pow(p, propParam[PMIdx::kParIdx]);
    kPerp = propParam[PMIdx::kPerp0] * std::pow(p, propParam[PMIdx::kPerpIdx]);
    kPP = propParam[PMIdx::kPP0] * std::pow(p, propParam[PMIdx::kPPIdx]);
    vStr = propParam[PMIdx::vStr0] * std::pow(p, propParam[PMIdx::vStrIdx]);
#endif
    IntVect i;
    RealVect f;
    convertPos(packet.pos(), pLo, dxInv, i, f);
    RealVect gradNumCR = (1.0 - accWgt) *
      grad(accOld, i, f, AIdx::num + accIdx);
    for (int n = accIdx+1; n < nT; n++)
      gradNumCR += grad(accOld, i, f, AIdx::num + n);
    RealVect B = lint_vec(gas, i, f, GIdx::Bx);
    if (gradNumCR.dotProduct(B) > 0) vStr *= -1.0;

    // Step 3: compute the continuous and catastrophic loss rates per
    // distance travelled
    Real rho = lint(gas, i, f, GIdx::TotalDensity);
    Real opacity, dp_dl;
#if (LOSSMODEL == NONE)
    opacity = dp_dl = 0.0;
#else
    if (T <= constants::Tpth) {
      opacity = 0.0;
    } else {
      {
	// PP cross section from Kafexhiu et al. (2014, PRD, 90 123014)
	Real fac1 = std::log(T/constants::Tpth);
	Real fac2 = 1.0 - std::pow(constants::Tpth/T, 1.9);
	Real sigmaPP = (30.7 - 0.96 * fac1 + 0.18 * fac1*fac1) *
	  fac2*fac2*fac2 * units::mb;
	Real etaPP = 0.5;
	Real nucDen = rho / constants::MW_nuc_per_mass;
	opacity = nucDen * sigmaPP * etaPP;
      }
      // TODO: put ionisation loss terms here
      dp_dl = 0.0;
    }
#endif

    // Step 4: compute advection velocity at current position
    RealVect va = lint_vec(gas, i, f, GIdx::Vx);
    Real vaMag = va.vectorLength();

    // Step 5: compute time step
    Real kMag = std::sqrt(kPar*kPar + 2*kPerp*kPerp);
    Real dtAdv = cStep / (dxInvMax * (vaMag + std::fabs(vStr) + eps));
    Real dtDiff = cStep / (dxInvMax * dxInvMax * (kMag + eps));
    Real dtPP = cStep * p * p / (kPP + eps);
    Real dt = 1.0 / (1.0/dtAdv + 1.0/dtDiff + 1.0/dtPP);
#if AMREX_USE_GPU
    dtPartMin = amrex::Gpu::Atomic::Min(dtPartMin, dt);
#else
    if (dtPartMin > dt) dtPartMin = dt;
#endif
    bool lastStep;
    if (t + dt >= tStop) {
      dt = tStop - t;
      lastStep = true;
    } else {
      lastStep = false;
    }
#if 0
    if (!(dt > 0.0)) {
      std::cout << "dtAdv = " << dtAdv
		<< ", dtDiff = " << dtDiff
		<< ", dtPP = " << dtPP
		<< ", dt = " << dt
		<< ", t = " << t
		<< ", tStop = " << tStop
		<< std::endl;
    }
    if (!(dtAdv > 0.0)) {
      std::cout << "vaMag = " << vaMag
		<< ", vStrMag = " << std::fabs(vStr)
		<< ", eps = " << eps
		<< std::endl;
    }
#endif

    // Step 6: add to accumulator
#ifdef AMREX_USE_GPU
    amrex::Gpu::Atomic::
      AddNoRet( &accNew(iAcc[0], iAcc[1], iAcc[2], accIdx),
		w * dt * volInv / dtGlob );
    amrex::Gpu::Atomic::
      AddNoRet( &accNew(iAcc[0], iAcc[1], iAcc[2], accIdx+nT),
		w * T * constants::mp_c2 * dt * volInv / dtGlob );
#else
    // For openMP case, the accNew object we have been passed is a
    // copy local to this thread, so we can safely add to it;
    // reduction will occur at the end
    accNew(iAcc[0], iAcc[1], iAcc[2], accIdx)
      += w * dt * volInv / dtGlob;
    accNew(iAcc[0], iAcc[1], iAcc[2], accIdx+nT)
      += w * T * constants::mp_c2 * dt * volInv / dtGlob;
#endif
    
    // Step 7: update grammage, momentum, and statistical weight
    Real etaP = amrex::RandomNormal(0, 1);
    Real pathlength = dt * constants::c * sr::v_from_p(m, T);
    gr += rho * pathlength;
    //gr += dt;
#if 0
    std::cout << "increasing dt by " << dt
	      << ", new time will be " << t + dt << std::endl;
#endif
    p += dt * (kPP / p - dp_dl * pathlength) +
      std::sqrt(2 * kPP * dt) * etaP;
    w *= std::exp(-opacity * pathlength);
    if (w / wInj < minWgtFrac) {
#if 0
      std::cout << "Propagation kernel exit due to packet weight "
		<< "falling below minimum: w = " << w
		<< ", wInj = " << wInj << ", w/wInj = " << w/wInj
		<< ", min = " << minWgtFrac
		<< std::endl;
#endif
      packet.id() = -1;  // Flag for deletion
      return 0;
    }
    T = sr::T_from_p(m, p);
    if (T < Tmin || T >= Tmax) {
#if 0
      std::cout
#ifdef _OPENMP
	<< "Thread " << omp_get_thread_num() << ": "
#endif
	<< "Propagation kernel exit due to energy out of "
	<< "range: T = " << T
	<< ", Tmin = " << Tmin
	<< ", Tmax = " << Tmax
	<< std::endl;
#endif
      packet.id() = -1;  // Flag for deletion
      return 0;
    }

    // Step 8: compute effective velocities for this step
    Real etaT = amrex::RandomNormal(0, 1);
    Real etaN = amrex::RandomNormal(0, 1);
    Real etaB = amrex::RandomNormal(0, 1);
    vT = vStr + sqrt(2*kPar/dt) * etaT;
    vN = sqrt(2*kPerp/dt) * etaN;
    vB = sqrt(2*kPerp/dt) * etaB;
    
    // Step 9: update position
#if 0
    std::cout << "calling movePacket with dt = " << dt
	      << ", packet position = " << packet.pos()
	      << std::endl;
#endif
    Real dtDone = movePacket(dt, errTol, pLo, dxInv, dxInvMax,
			     rbx, gas, va, packet);
#if 0
    std::cout << "movePacket returned dtDone = " << dtDone
	      << ", packet position = " << packet.pos()
	      << std::endl;
#endif

    // Step 10: check for termination of update, either because we
    // left the box, or because the step is done
    if (dtDone != dt) {

      // Case 1: if we couldn't complete the step because we left the
      // valid box during it, we need to flag that the step is incomplete,
      // so we can finish it in the next cycle
      if (!lastStep) tStep = t + dt - dtDone;
      else tStep = tStop;
#if 0
      std::cout
#ifdef _OPENMP
	<< "Thread " << omp_get_thread_num() << ": "
#endif
	<< "PK done for packet " << packet.id()
	<< "; packet left tile -- "
	<< "final position = "
	<< packet.pos()
	<< ", time = "
	<< t
	<< ", target time for this step = "
	<< tStep
	<< std::endl << std::endl;
#endif
      return 1;

    } else {

      // Case 2: if this was the last step, set all times to tStop to
      // indicate that this packet is done; if it was not the last
      // step, check that the packet is still inside the valid box,
      // since it may have left the valid box during the last part of
      // its update
      if (lastStep) {
	t = tStep = tStop;
#if 0
	std::cout
#ifdef _OPENMP
	  << "Thread " << omp_get_thread_num() << ": "
#endif
	  << "PK done for packet " << packet.id()
	  << "; complete step -- "
	  << "final position = "
	  << packet.pos()
	  << ", time = "
	  << t
	  << std::endl << std::endl;
#endif
	return 0;
      } else {
	if (!rbx.contains(packet.pos())) return 1;
      }
    }
  }
}

#endif
// _PROPKERNEL_H_
