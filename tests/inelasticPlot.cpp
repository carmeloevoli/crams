#include "cgs.h"
#include "inelastic.h"
#include "pid.h"

int main() {
  const auto xC12 = CRAMS::InXsecTripathi99(CRAMS::C12);
  const auto xO16 = CRAMS::InXsecTripathi99(CRAMS::O16);
  const auto xFe56 = CRAMS::InXsecTripathi99(CRAMS::Fe56);

  for (double E = 0.1 * CRAMS::CGS::GeV; E < 3e3 * CRAMS::CGS::GeV; E *= 1.1) {
    std::cout << std::scientific;
    std::cout << E / CRAMS::CGS::GeV << " ";
    std::cout << xC12.getXsecOnHtarget(E) / CRAMS::CGS::mbarn << " ";
    std::cout << xO16.getXsecOnHtarget(E) / CRAMS::CGS::mbarn << " ";
    std::cout << xFe56.getXsecOnHtarget(E) / CRAMS::CGS::mbarn << " ";
    std::cout << std::endl;
  }
}
