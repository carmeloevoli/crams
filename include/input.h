#ifndef INCLUDE_INPUT_H_
#define INCLUDE_INPUT_H_

#include <map>

#include "cgs.h"
#include "particlelist.h"

namespace CRAMS {

class Input {
 private:
  double m_TSimMin = 0.1 * CGS::GeV;
  double m_TSimMax = 10. * CGS::TeV;
  size_t m_TSimSize = 100;

  double m_ROutputMin = 2 * CGS::GeV;
  double m_ROutputMax = 10. * CGS::TeV;
  size_t m_ROutputSize = 100;

  bool m_doSecondary = true;

  double m_H = 7. * CGS::kpc;
  double m_mu = 2.3 * CGS::mgram / CGS::cm2;
  double m_v_A = 4.40940 * CGS::km / CGS::sec;
  double m_R_b = 312. * CGS::GeV;  // Evoli, PRD, 2019
  double m_delta = 5.65132e-01;
  double m_ddelta = 0.2;      // Evoli, PRD, 2019
  double m_smoothness = 0.1;  // Evoli, PRD, 2019
  double m_D_0 = 2.48255e28 * CGS::cm2 / CGS::sec;
  double m_X_s = -1;  // CGS::gram / CGS::cm2;
  double m_modulationPotential = 4.87754e-01 * CGS::GeV;
  double m_xsecsFudge = 1;
  size_t m_id = 0;

 public:
  Input() {}

  virtual ~Input();
  void print() const;
  void readParamsFromFile(const std::string& filename);
  void setParam(const std::string& key, const double& value);

  const double& TSimMin = m_TSimMin;
  const double& TSimMax = m_TSimMax;
  const size_t& TSimSize = m_TSimSize;

  const double& ROutputMin = m_ROutputMin;
  const double& ROutputMax = m_ROutputMax;
  const size_t& ROutputSize = m_ROutputSize;

  const bool& doSecondary = m_doSecondary;

  const double& H = m_H;
  const double& mu = m_mu;
  const double& v_A = m_v_A;
  const double& R_b = m_R_b;
  const double& delta = m_delta;
  const double& ddelta = m_ddelta;
  const double& smoothness = m_smoothness;
  const double& D_0 = m_D_0;
  const double& X_s = m_X_s;
  const double& modulationPotential = m_modulationPotential;
  const double& xsecsFudge = m_xsecsFudge;
  const size_t& id = m_id;
};

}  // namespace CRAMS

#endif  // INCLUDE_INPUT_H_