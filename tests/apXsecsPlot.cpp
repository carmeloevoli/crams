#include <iostream>

#include "cgs.h"
#include "xsecs/Korsmeier2018.h"

using namespace CRAMS;

int main() {
  const auto xs = Korsmeier2018SecAp();
  const auto T_p = 100. * CGS::GeV;
  const auto units = CGS::mbarn / CGS::GeV;
  for (double T_ap = 0.1 * CGS::GeV; T_ap < T_p; T_ap *= 1.1) {
    std::cout << std::scientific;
    std::cout << T_ap / CGS::GeV << " ";
    std::cout << xs.get(PbarChannel::pp, T_p, T_ap) / units << " ";
    std::cout << xs.get(PbarChannel::pHe, T_p, T_ap) / units << " ";
    std::cout << xs.get(PbarChannel::Hep, T_p, T_ap) / units << " ";
    std::cout << xs.get(PbarChannel::HeHe, T_p, T_ap) / units << " ";
    std::cout << std::endl;
  }
}