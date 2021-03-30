// Implementation of the CRSourceContainer class

#include "CRSourceContainer.H"
#include <AMReX_ParIter.H>

using namespace amrex;
using namespace std;
namespace SIdx = CRSourceIdx;
namespace PIdx = CRPacketIdx;
using PIter = ParIter<SIdx::nReal, SIdx::nInt>;

// A trivial little functor for use with ReduceSum. Its purpose is
// to store the sampling index, and use it to compute the relative
// weight of each source -- which in turn is used to compute the
// number of packets that should be drawn for that source.
class srcWgt {
public:
  srcWgt(Real q_) : q(q_) { }
  Real operator()(const CRSourceContainer::SuperParticleType& p) const {
    Real b0 = p.rdata(SIdx::p0) / constants::mp_c;
    Real b1 = p.rdata(SIdx::p1) / constants::mp_c;
    if (q != -1.0) {
      return (std::pow(b1, q+1) - std::pow(b0, q+1)) / (q + 1);
    } else {
      return std::log(b1 / b0);
    }
  }
private:
  Real q;
};

// Method to get the CR sources ready; this method copies the sampling
// index to every CR source, uses that to compute the weight, and then
// sums the weight over all sources
void CRSourceContainer::init() {

  // Get the local sum on this processor
  srcWgt wFunc(q);
  wTot = amrex::ReduceSum(*this, wFunc);
  
  // Add up weights over all processors
  ParallelDescriptor::ReduceRealSum(wTot);
}

// Method to inject packets
void CRSourceContainer::
injectPackets(const Real t,
	      const Real dt,
	      CRPacketContainer &packets) {

#if 0
  std::cout << "starting injectPackets" << std::endl;
#endif
  
  // Iterate over sources
  const int lev = 0;
  {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (PIter src_iter(*this, lev); src_iter.isValid(); ++src_iter) {

#if 0
      std::cout << "inject packets called on tile "
		<< src_iter.index() << ", "
		<< src_iter.LocalTileIndex() << std::endl;
#endif
      
      // Grab source list
      const auto& sources = src_iter.GetArrayOfStructs();

      // Grab the corresponding packet list
      pair<int,int> idx(src_iter.index(), src_iter.LocalTileIndex());
      auto &packetTile = packets.GetParticles(lev)[idx];

      // Set up functor to calculate weights
      srcWgt wFunc(q);

      // Iterate through sources
      for (const auto& s : sources) {

	// Set number of packets to draw for this source
	int n = std::lround(dt * packetRate * wFunc(s) / wTot);

	// Get new packets for this source
	Vector<CRPacket> new_packets =
	  static_cast<const CRSource>(s).drawPackets(n, t, dt, q);

	// Add new packets to packet list
	for (const auto& np : new_packets)
	  packetTile.push_back(np);

      }
    }
  }

#if 0
  std::cout << "injectPackets done" << std::endl;
#endif
}
