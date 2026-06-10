#ifndef CRAMS_PARTICLE_H_
#define CRAMS_PARTICLE_H_

#include <memory>
#include <utility>
#include <vector>

#include "crams/core/pid.h"

namespace CRAMS {

class Grammage;
class InelasticXsec;
class Input;
class Losses;
class PrimarySource;
class SecondarySource;
struct NucleusParameters;

class Particle {
 public:
  explicit Particle(const PID& pid);
  Particle(const PID& pid, const NucleusParameters& nucleusParameters);
  Particle(const Particle& other) = delete;
  Particle& operator=(const Particle& other) = delete;
  Particle(Particle&& other) noexcept;
  Particle& operator=(Particle&& other) noexcept;
  ~Particle();

  bool operator==(const Particle& other) const { return m_pid == other.m_pid; }
  const PID& getPid() const { return m_pid; }
  bool isDone() const { return m_isDone; }
  bool isChargeZ(int Z) const { return m_pid.getZ() == Z; }
  bool isStable() const { return m_decayTime < 0.; }
  double getDecayTime() const { return m_decayTime; }
  double getAbundance() const { return m_abundance; }
  double getSlope() const { return m_slope; }
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
  void reset();
  void computeIntensity(const Input& input);
  void dump() const;
  void computeFluxAtEnergy_num();
  double I_T_interpol(double T) const;
  double I_T_TOA(double T, double modulationPotential) const;
  double I_R_TOA(double R, double modulationPotential) const;

 public:
  double Q_total(double T) const;
  double Lambda_1(double T) const;
  double Lambda_2(double T) const;
  double externalIntegrand(double T_prime, double T) const;
  double internalIntegrand(double T_second) const;
  double ExpIntegral(double T, double T_prime) const;
  double computeFluxAtEnergy(double T) const;

 protected:
  double productionProfileFromUnstable(const Input& input, double T, double decayTimeAtRest) const;

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
  std::unique_ptr<Grammage> m_X;
  std::unique_ptr<PrimarySource> m_Q_p;
  std::unique_ptr<SecondarySource> m_Q_sec;
  std::unique_ptr<SecondarySource> m_Q_ter;
  std::unique_ptr<SecondarySource> m_Q_Xs;
  std::unique_ptr<InelasticXsec> m_sigmaIn;
  std::unique_ptr<Losses> m_dEdX;
};

using Particles = std::vector<Particle>;
using itParticle = std::pair<bool, Particles::iterator>;

}  // namespace CRAMS

#endif  // CRAMS_PARTICLE_H_
