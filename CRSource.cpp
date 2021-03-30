// Implementation of CRSource

#include "CRSource.H"

using namespace std;
using namespace amrex;
namespace SIdx = CRSourceIdx;
namespace PIdx = CRPacketIdx;

// Constructor
CRSource::CRSource(Vector<Real> x,
		   Real lum_metric,
		   Real T0_GeV,
		   Real T1_GeV,
		   Real q,
		   Real m,
		   int Z) {

  // Set ID number and CPU number
  id()  = NextID();
  cpu() = ParallelDescriptor::MyProc();

  // Store input data in the particle structure
  pos(0) = x[0];
#if AMREX_SPACEDIM > 1
  pos(1) = x[1];
#endif
#if AMREX_SPACEDIM > 2
  pos(2) = x[2];
#endif

  // Store mass, charge and spectral index
  rdata(SIdx::q) = q;
  rdata(SIdx::m) = m;
  idata(SIdx::Z) = Z;

  // Derive momentum limits and normalization factor
  rdata(SIdx::p0) = sr::p_from_T(m, T0_GeV * units::GeV / constants::mp_c2);
  rdata(SIdx::p1) = sr::p_from_T(m, T1_GeV * units::GeV / constants::mp_c2);
  rdata(SIdx::A) = lum_metric / (TBar() * constants::mp_c2);
}

// Packet drawing method
Vector<CRPacket> CRSource::drawPackets(const int n,
				       const Real t,
				       const Real dt,
				       const Real qSamp) const {

#if 0
  std::cout << "drawPackets called for source " << id()
	    << " with n = " << n
	    << ", t = " << t
	    << ", dt = " << dt
	    << ", qSamp = " << qSamp
	    << std::endl;
#endif
  
  // Allocate holder for output
  Vector<CRPacket> packets(n);

  // Loop over packets to be drawn
  const Real q = rdata(SIdx::q);
  Real w = 0.0;
  for (int i=0; i<n; i++) {

    // Set packet properties that don't depend on this source
    packets[i].rdata(PIdx::gr) = 0.0;
    packets[i].rdata(PIdx::vT) = 0.0;
    packets[i].rdata(PIdx::vN) = 0.0;
    packets[i].rdata(PIdx::vB) = 0.0;

    // Set CR packet properties that do depend on this source
    packets[i].pos(0) = pos(0);
#if AMREX_SPACEDIM > 1
    packets[i].pos(1) = pos(1);
#endif
#if AMREX_SPACEDIM > 2
    packets[i].pos(2) = pos(2);
#endif
    packets[i].rdata(PIdx::m) = rdata(SIdx::m);
    packets[i].idata(PIdx::Z) = idata(SIdx::Z);
    packets[i].idata(PIdx::src) = id();

    // Set time counters
    packets[i].rdata(PIdx::t) = t + Random()*dt;
    packets[i].rdata(PIdx::tInj) = packets[i].rdata(PIdx::t);
    packets[i].rdata(PIdx::tStep) = packets[i].rdata(PIdx::t);

    // Set momentum
    Real y = Random();
    if (qSamp == -1.0) {
      packets[i].rdata(PIdx::p) = rdata(SIdx::p0) *
	exp(y * log(rdata(SIdx::p1) / rdata(SIdx::p0)));
    } else {
      Real p0q = pow(rdata(SIdx::p0), 1+qSamp);
      Real p1q = pow(rdata(SIdx::p1), 1+qSamp);
      packets[i].rdata(PIdx::p) =
	pow(y*p0q + (1.0-y)*p1q, 1.0/(1.0+qSamp));
    }

    // Set weight scaling
    if (qSamp != q)
      packets[i].rdata(PIdx::w) = pow(packets[i].rdata(PIdx::p), q-qSamp);
    w += packets[i].rdata(PIdx::w);
  }

  // Normalize weights to add up to correct sum
  Real scal = nLum() * dt / w;
  for (int i=0; i<n; i++) {
    packets[i].rdata(PIdx::w) *= scal;
    packets[i].rdata(PIdx::wInj) = packets[i].rdata(PIdx::w);
  }

#if 0
  Real nTot = 0.0, Ttot = 0.0;
  for (int i=0; i<n; i++) {
    nTot += packets[i].rdata(PIdx::w);
    Ttot += packets[i].rdata(PIdx::w) *
      sr::T_from_p(packets[i].rdata(PIdx::m),
		   packets[i].rdata(PIdx::p)) * constants::mp_c2;
  }
  std::cout << "Source " << id() << " injected " << n
	    << " packets, total w = "
	    << nTot << ", T = " << Ttot
	    << std::endl;
#endif

  // Return the packets we drew
  return packets;
}
