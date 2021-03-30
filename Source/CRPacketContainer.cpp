// Implementation of the CRPacketContainer class

#include "CRPacketContainer.H"
#include "Constants.H"
#include "Gas.H"
#include "PropKernel.H"
#include <cmath>
#include <limits>

using namespace amrex;
using namespace std;
namespace PIdx = CRPacketIdx;
namespace GIdx = GasIdx;
namespace AIdx = AccumIdx;
using PIter = ParIter<PIdx::nReal, PIdx::nInt>;


// The advance function
Real CRPacketContainer::
advancePackets(const Real &tStart, const Real &dtGlob,
	       Gas &gas, Accum &acc) {

  // Single level for now
  const int lev = 0;
  
  // Fill ghost cells in accumulator
  acc.BCFill();

  // Grab information on energy gridding
  int nT = acc.nT();
  Real Tmin = acc.get_Tmin();
  Real dlogT = acc.get_dlogT();

  // Grab parameters of the propagation model
  GpuArray<Real,NPROP_PARAM> propParams = prop.getParams();

  // Loop until all packets are done
  Real tStop = tStart + dtGlob;
  Real dtPartMin = numeric_limits<Real>::max();
  while (true) {

    // Pointers to things we need to accumulate across threads on CPU
    // or GPU
    int npRemain = 0;
#if AMREX_USE_GPU
    AsyncArray<int> npRemain_aa(&npRemain, 1);
    int *npRemain_ptr = npRemain_aa.data();
    AsyncArray<Real> dtPartMin_aa(&dtPartMin, 1);
    Real *dtPartMin_ptr = dtPartMin_aa.data();
#endif

    // Loop over packet tiles; this is threaded
#ifdef _OPENMP
#pragma omp parallel reduction(+:npRemain) reduction(min:dtPartMin) if (Gpu::notInLaunchRegion())
#endif
    for (PIter pti(*this, lev); pti.isValid(); ++pti) {

      // Grab the packet data on this tile
      auto& packets = pti.GetArrayOfStructs();
      const size_t np = packets.numParticles();
      ParticleType* pstruct = packets().dataPtr();

      // Grab the gas and accumulator data on this tile
      const FArrayBox &gasFab = gas[pti];
      const FArrayBox &accOldFab = acc.oldData()[pti];
      const Array4<Real const> &gasData = gasFab.const_array();
      const Array4<Real const> &accOldData = accOldFab.const_array();
      FArrayBox &accNewFab = acc.newData()[pti];

      // Grab the limits of the tile box
      const Box &bx = pti.tilebox();
      RealVect xLo = bx.smallEnd() / dxInv + pLo;
      RealVect xHi = (bx.bigEnd() + IntVect::TheUnitVector()) / dxInv + pLo;
      const RealBox &rbx = RealBox(xLo.dataPtr(), xHi.dataPtr());

      // If we're using openMP, we need to make dummy copies of the
      // accumulator fabs and all counters for each thread, which
      // we can then pass to the propagation kernel; we will then
      // sum these up after propagation is complete.
#ifdef _OPENMP
      FArrayBox &accNew = localAcc[omp_get_thread_num()];
      Box growbox(bx);
      growbox.grow(acc.ngrow());
      accNew.resize(growbox, accNewFab.nComp());
      int npLeft = 0;
      int *npLeft_ptr = &npLeft;
      Real dtPartMinLoc = dtPartMin;
      Real *dtPartMinLoc_ptr = &dtPartMinLoc;
#else
      FArrayBox &accNew = accNewFab;
      int *npLeft_ptr = npRemain_ptr;
      Real *dtPartMinLoc_ptr = dtPartMin_ptr;
#endif
      accNew.setVal(0.0);
      Array4<Real> const& accNewData = accNew.array();

#if 0
      std::cout
#ifdef _OPENMP
	<< "Thread " << omp_get_thread_num() << ": "
#endif
	<< "Advancing " << np << " packets on tile " << bx
	<< std::endl;
#endif

	
      // Loop over packets
      amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (size_t n)
      {
	// Grab the packet we're working on
	ParticleType &packet = pstruct[n];

	// Flag for if the packet has left the currently active box
	int leftBox = 0;

	// If this packet failed to complete its previous step
	// because it left its previous active box, try to complete
	// the step now
	if (packet.rdata(PIdx::t) < packet.rdata(PIdx::tStep)) {
	  Real dt = packet.rdata(PIdx::tStep) - packet.rdata(PIdx::t);
#if 0
	  std::cout << "Packet " << packet.id()
		    << " did not complete last propagation step;"
		    << " trying to finish remaining time = "
		    << dt
		    << std::endl;
#endif
	  IntVect i;
	  RealVect f;
	  convertPos(packet.pos(), pLo, dxInv, i, f);
	  RealVect va = lint_vec(gasData, i, f, GIdx::Vx);
	  Real dtDone = movePacket(dt, errTol, pLo, dxInv, dxInvMax,
				   rbx, gasData, va, packet);
	  if (dtDone != dt) leftBox = 1;
	}

	// If the packet has not reached the target time, and has
	// not left our box, call the propagation kernel
	if (packet.rdata(PIdx::t) < tStop && !leftBox) {
	  leftBox =
	    advancePacket(tStop, dtGlob, *dtPartMinLoc_ptr,
			  cStep, errTol,
			  pLo, dxInv, dxInvMax, volInv,
			  propParams, gasData,
			  nT, Tmin, dlogT,
			  minWgtFrac, rbx, accOldData, accNewData,
			  packet);
	}

	// If we stopped because the particle left the active box,
	// check if it left the computational domain entirely; if
	// so, flag it for deletion, and, if not, increment the
	// number of incomplete particles
	if (leftBox) {
	  if (pDomain.contains(packet.pos())) {
#if AMREX_USE_GPU
	    amrex::Gpu::Atomic::AddNoRet(npLeft_ptr, 1);
#else
	    (*npLeft_ptr)++;
#endif
	  } else {
	    packet.id() = -1;
	  }
	}
      });

      // If running in openMP mode, reduce quantities held by the
      // individual threads
#ifdef _OPENMP
      accNewFab.atomicAdd(accNew, growbox, growbox, 0, 0, accNew.nComp());
      npRemain += npLeft;
      if (dtPartMin > dtPartMinLoc) dtPartMin = dtPartMinLoc;
#endif
    }

    // If running in GPU mode, copy data back to the CPU from the GPU
#if AMREX_USE_GPU
    npRemain_aa.copyToHost(&npRemain, 1);
    dtPartMin_aa.copyToHost(dtPartMin_ptr, 1);
#endif

    // Redistribute particles
    Redistribute();

    // Sum number of remaining particles across all MPI ranks
    ParallelDescriptor::ReduceIntSum(npRemain);

#if 0
    std::cout << "****************************************"
	      << std::endl
	      << "Done with packet propagation all boxes"
	      << "; number of incomplete packets = "
	      << npRemain
	      << std::endl;
#endif
    
    // If there are no packets remaining that have not reached their
    // final time, sum guard cells across processors and then exit;
    // otherwise start another pass
    if (npRemain == 0) {
      acc.newData().SumBoundary(m_gdb->Geom(0).periodicity());
      break;
    }
  }

  // Final step: now that we have filled in the new accumulator values
  // everywhere, use them to estimate the next time step; we set this
  // to the global minimum of:
  //
  // dt = cGlob * dtPrev * (n_i,new + cFracMin * nAvg_i) / abs(n_i,new-n_i,old)
  //
  // where dtPrev is the previous time step, n_i is the number density
  // in energy bin i, and nAvg_i is the mean number density in that
  // energy bin over the whole domain, computed at the new time.
  Real dtNext = std::numeric_limits<Real>::max();
  const MultiFab& newAcc = acc.newData();
  const MultiFab& oldAcc = acc.oldData();
  Vector<Real> nAvg = acc.nAvg();
#if AMREX_USE_GPU
  {
    // GPU version
    AsyncArray<Real> dtNext_aa(&dtNext, 1);
    AsyncArray<Real> nAvg_aa(nAvg.data(), nAvg.size());
    Real *dtNext_data = dtNext_aa.data();
    Real *nAvg_data = nAvg_aa.data();
    for (MFIter mfi(newAcc,false); mfi.isValid(); ++mfi) {

      // Grab old and new data on this tile
      Array4<Real const> const& n = newAcc[mfi].const_array();
      Array4<Real const> const& o = oldAcc[mfi].const_array();
      const Box& bx = mfi.tilebox();

      // Compute dt for every component in every cell
      ParallelFor(bx, acc.nT(),
		  [=] AMREX_GPU_DEVICE (int i, int j, int k, int c)
      {
	Real dtCell =
	  cGlob * dtGlob * (n(i,j,k,c) + cFracMin * nAvg_data[c]) /
	  amrex::Math::abs(n(i,j,k,c)-o(i,j,k,c));
	Gpu::Atomic::Min(dtNext_data, dtCell);
      });
    }

    // Copy time step data back to CPU
    dtNext_aa.copyToHost(&dtNext, 1);

  }
#else
  {
    // CPU / openMP version
#ifdef _OPENMP
#pragma omp parallel reduction(min:dtNext)
#endif
    for (MFIter mfi(newAcc,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      // Grab old and new data on this tile
      Array4<Real const> const& n = newAcc[mfi].const_array();
      Array4<Real const> const& o = oldAcc[mfi].const_array();
      const Box& bx = mfi.tilebox();

      // Compute dt for every cell
      Real *dtNextDat = &dtNext;
      ParallelFor(bx, acc.nT(),
		  [=] AMREX_GPU_DEVICE (int i, int j, int k, int c)
      {
	Real dtCell =
	  cGlob * dtGlob * (n(i,j,k,c) + cFracMin * nAvg[c]) /
	  amrex::Math::abs(n(i,j,k,c)-o(i,j,k,c));
#if 0
#ifdef VERBOSEDEBUG
	if (*dtNextDat > dtCell)
	  cout << "old dtNext = " << *dtNextDat
	       << ", new dtNext = " << dtCell
	       << ", n = " << n(i,j,k,c)
	       << ", o = " << o(i,j,k,c)
	       << ", ijkc = " << i << " " << j << " " << k << " " << c
	       << endl;
#endif
#endif
	(*dtNextDat) = std::min(*dtNextDat, dtCell);
      });
    }
  }
#endif

#if 0
#ifdef VERBOSEDEBUG
  cout << "dtNext = " << dtNext
       << ", dtPartMin = " << dtPartMin
       << endl;
#endif
#endif
  
  // Set minimum of next time step to be equal to smallest individual
  // particle time step
  if (dtNext < 100*dtPartMin) dtNext = 100*dtPartMin;

  // Reduce dtNext across MPI ranks, then return
  ParallelDescriptor::ReduceRealMin(dtNext);
  return dtNext;
}
		 

		 
