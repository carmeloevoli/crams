#ifndef INCLUDE_PARTICLE_H_
#define INCLUDE_PARTICLE_H_

#include <algorithm>
#include <memory>
#include <vector>

#include "grammage.h"
#include "inelastic.h"
#include "losses.h"
#include "pid.h"
#include "primary.h"
#include "secondary.h"
#include "xsecs/Evoli2019.h"

namespace CRAMS {
class Particle {
 public:
  Particle(const PID& pid);
  Particle(const PID& pid, const NucleusParameters& nucleusParameters);
  virtual ~Particle();

  bool operator==(const Particle& other) const { return m_pid == other.m_pid; }
  const PID getPid() const { return m_pid; }
  const bool isDone() const { return m_isDone; }
  const bool isChargeZ(const int& Z) const { return (m_pid.getZ() == Z) ? true : false; }
  const bool isStable() const { return m_decayTime < 0.; }
  const double getDecayTime() const { return m_decayTime; }
  const double getAbundance() const { return m_abundance; }
  const double getSlope() const { return m_slope; }
  const std::vector<double>& getEnergyVector() const { return m_T; }
  const std::vector<double>& getIntensityVector() const { return m_I_T; }
  void setDone() { m_isDone = true; }
  void unsetDone() { m_isDone = false; }

  void buildVectors(const Input& input);
  void buildGrammage(const Input& input);
  void buildPrimarySource(const Input& input);
  void buildLosses(const Input& input);
  void buildInelasticXsecs(const Input& input);
  void buildSecondarySource(const Input& input, const std::vector<Particle>& particles);
  void buildGrammageAtSource(const Input& input, const std::vector<Particle>& particles);
  void buildTertiarySource(const std::vector<Particle>& particles);
  void buildAntiprotonSource(const std::vector<Particle>& particles);
  void reset();
  void computeIntensity();
  void dump() const;
  double I_T_interpol(const double& T) const;
  double I_R_TOA(const double& R, const double& modulationPotential) const;

 public:
  // publicsolver.cpp
  double Q_total(const double& T) const;
  double Lambda_1(const double& T) const;
  double Lambda_2(const double& T) const;
  double externalIntegrand(const double& T_prime, const double& T);
  double internalIntegrand(const double& T_second);
  double ExpIntegral(const double& T, const double& T_prime);
  double computeFluxAtEnergy(const double& T);

 protected:
  double productionProfileFromUnstable(const Input& input, const double& T, const double& decayTimeAtRest);

 protected:
  PID m_pid;
  bool m_isDone = false;
  bool m_doGrammageAtSource = false;
  bool m_doSecondary = false;
  double m_abundance = 0;
  double m_slope = 0;
  double m_decayTime = -1;
  std::vector<double> m_T;
  std::vector<double> m_I_T;
  std::shared_ptr<Grammage> m_X;
  std::shared_ptr<PrimarySource> m_Q_p;
  std::shared_ptr<SecondarySource> m_Q_sec;
  std::shared_ptr<SecondarySource> m_Q_ter;
  std::shared_ptr<SecondarySource> m_Q_ap;
  std::shared_ptr<SecondarySource> m_Q_Xs;
  std::shared_ptr<InelasticXsec> m_sigmaIn;
  std::shared_ptr<Losses> m_dEdX;
};

typedef std::vector<Particle> Particles;
typedef std::pair<bool, std::vector<Particle>::iterator> itParticle;

}  // namespace CRAMS

#endif  // INCLUDE_PARTICLE_H_
