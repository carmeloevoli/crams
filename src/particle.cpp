#include "crams/particle.h"

#include <plog/Log.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

#include "crams/core/cgs.h"
#include "crams/core/input.h"
#include "crams/fragmentation.h"
#include "crams/inelastic.h"
#include "crams/particlelist.h"
#include "crams/physics/grammage.h"
#include "crams/physics/losses.h"
#include "crams/physics/primary.h"
#include "crams/secondary.h"
#include "crams/utils/numeric.h"
#include "crams/utils/utilities.h"

namespace CRAMS {
namespace {

using Utilities::pow2;

constexpr size_t kSourceGridSize = 400;
constexpr double kUnavailable = std::numeric_limits<double>::quiet_NaN();

struct DecayContribution {
  PID child;
  PID parent;
};

const DecayContribution kDecayContributions[] = {
    {B10, Be10}, {N14, C14}, {Mg26, Al26}, {Ar36, Cl36}, {Fe54, Mn54},
    {Ni60, Fe60},
};

double coth(double x) { return 1. / std::tanh(x); }

std::vector<double> makeSourceEnergyGrid() {
  return Utilities::LogAxis(0.1 * CGS::GeV, 10. * CGS::TeV, kSourceGridSize);
}

double safeReciprocal(double value) { return (value != 0.) ? 1. / value : kUnavailable; }

void writeParticleDumpColumns(std::ostream& out) {
  out << "# Columns\n";
  out << "# T [GeV] -> 1\n";
  out << "# R [GV] -> 2\n";
  out << "# Q_pri -> 3\n";
  out << "# Q_sec -> 4\n";
  out << "# X [g/cm2] -> 5\n";
  out << "# tau_esc [Myr] -> 6\n";
  out << "# tau_adv [Myr] -> 7\n";
  out << "# X_cr [g/cm2] -> 8\n";
  out << "# dEdX -> 9\n";
  out << "# tau_ion [Myr] -> 10\n";
  out << "# tau_inelastic [Myr] -> 11\n";
}

const Particle* findParticle(const Particles& particles, const PID& pid) {
  for (const auto& particle : particles) {
    if (particle.getPid() == pid) return &particle;
  }
  return nullptr;
}

const Particle& findParticleOrThrow(const Particles& particles, const PID& pid) {
  const auto particle = findParticle(particles, pid);
  if (particle == nullptr) throw std::runtime_error("Particle: required particle " + pid.toString() + " not found");
  return *particle;
}

}  // namespace

Particle::Particle(const PID& pid, const NucleusParameters& nucleusParameters)
    : m_pid(pid),
      m_abundance(nucleusParameters.abundance),
      m_slope(nucleusParameters.slope),
      m_decayTime(nucleusParameters.decayTime) {}

Particle::Particle(const PID& pid) : m_pid(pid) {}

Particle::Particle(Particle&& other) noexcept = default;

Particle& Particle::operator=(Particle&& other) noexcept = default;

Particle::~Particle() { LOGD << "released memory of particle " << m_pid; }

void Particle::reset() {
  m_X.reset();
  m_Q_p.reset();
  m_Q_sec.reset();
  m_Q_ter.reset();
  m_Q_Xs.reset();
  m_sigmaIn = nullptr;
  m_dEdX.reset();
}

void Particle::buildVectors(const Input& input) {
  m_T = Utilities::LogAxis(input.TSimMin(), input.TSimMax(), input.TSimSize());
  m_I_T.assign(input.TSimSize(), 0.);
}

void Particle::buildGrammage(const Input& input) {
  if (isStable())
    m_X = std::make_unique<Grammage>(m_pid, input);
  else
    m_X = std::make_unique<Grammage>(m_pid, input, m_decayTime);
}

void Particle::buildPrimarySource(const Input& input) {
  m_Q_p = std::make_unique<PrimarySource>(m_pid, m_abundance, m_slope, input.mu());
}

void Particle::buildLosses(const Input& input) { m_dEdX = std::make_unique<Losses>(m_pid, input); }

void Particle::buildInelasticXsecs(const InelasticXsec& sigmaIn) { m_sigmaIn = &sigmaIn; }

double Particle::productionProfileFromUnstable(const Input& input, double T, double decayTimeAtRest) const {
  const double v = Utilities::T2beta(T) * CGS::cLight;
  const double u = input.v_A();
  const double H = input.H();
  const double D = m_X->D(T);
  const double value = u / input.mu() / v;
  const double decayTimeAtT = Utilities::T2gamma(T) * decayTimeAtRest;
  const double Delta = std::sqrt(1. + 4. * D / (pow2(u) * decayTimeAtT));
  const double profile = Delta * coth(u * H * Delta / 2. / D) - coth(u * H / 2. / D);
  return value * profile;
}

void Particle::buildSecondarySource(const Input& input, const std::vector<Particle>& particles,
                                    const NucFragXsec& nucfrag) {
  m_doSecondary = input.doSecondary();
  const auto T_s = makeSourceEnergyGrid();
  std::vector<double> Q_s;
  Q_s.reserve(T_s.size());

  for (const auto T : T_s) {
    double value = 0.;
    for (const auto& particle : particles) {
      const auto& parentPid = particle.getPid();
      if (parentPid.getA() <= m_pid.getA() || !particle.isDone()) continue;
      value += nucfrag.getXsecOnISM(parentPid, m_pid, T) * particle.I_T_interpol(T);
    }
    Q_s.push_back(value / CGS::meanISMmass);
  }

  for (const auto& contribution : kDecayContributions) {
    if (m_pid != contribution.child) continue;

    const auto& parent = findParticleOrThrow(particles, contribution.parent);
    const double parentDecayTime = parent.getDecayTime();
    for (size_t i = 0; i < T_s.size(); ++i) {
      Q_s[i] += productionProfileFromUnstable(input, T_s[i], parentDecayTime) * parent.I_T_interpol(T_s[i]);
    }
    break;
  }

  m_Q_sec = std::make_unique<SecondarySource>(m_pid, T_s, Q_s);
}

void Particle::buildTertiarySource(const std::vector<Particle>& particles) {
  const auto T_t = makeSourceEnergyGrid();
  const double mp = CGS::protonMassC2;
  const auto proton = findParticle(particles, H1);
  const bool useProtonFlux = proton != nullptr && proton->isDone();
  std::vector<double> Q_t;
  Q_t.reserve(T_t.size());

  for (const auto T : T_t) {
    const double T_prime = T / CGS::inelasticity;
    double sigma_ISM = sigma_pp(T_prime);
    sigma_ISM *= (1. + CGS::K_He * CGS::f_He) / (1. + CGS::f_He);
    double value = sigma_ISM / CGS::inelasticity;
    value *= (T_prime + mp) / (T + mp);
    value *= std::pow(T * (T + 2. * mp), 1.5) / std::pow(T_prime * (T_prime + 2. * mp), 1.5);
    if (useProtonFlux) value *= proton->I_T_interpol(T_prime);
    value /= CGS::meanISMmass;
    Q_t.push_back(value);
  }
  m_Q_ter = std::make_unique<SecondarySource>(m_pid, T_t, Q_t);
}

void Particle::buildGrammageAtSource(const Input& input, const std::vector<Particle>& particles,
                                     const NucFragXsec& nucfrag) {
  m_doGrammageAtSource = true;
  const auto T_X = makeSourceEnergyGrid();
  std::vector<double> Q_X(T_X.size(), 0.);

  for (const auto& particle : particles) {
    const auto& parentPid = particle.getPid();
    if (parentPid.getA() <= m_pid.getA() || !particle.isDone()) continue;

    const auto Q_p = PrimarySource(parentPid, particle.getAbundance(), particle.getSlope(), input.mu());
    for (size_t i = 0; i < T_X.size(); ++i) {
      const double T = T_X[i];
      const double r = input.X_s() / CGS::meanISMmass * nucfrag.getXsecOnISM(parentPid, m_pid, T);
      Q_X[i] += r * Q_p.get(T);
    }
  }
  m_Q_Xs = std::make_unique<SecondarySource>(m_pid, T_X, Q_X);
}

double Particle::I_T_interpol(double T) const {
  if (m_T.empty() || T <= m_T.front() || T >= m_T.back()) return 0.;
  return Numeric::LinearInterpolatorLog<double>(m_T, m_I_T, T);
}

double Particle::I_T_TOA(double T, double modulationPotential) const {
  // see arXiv:1511.08790
  const double Phi = m_pid.getZoverA() * modulationPotential;
  const double T_ISM = T + Phi;
  double factor = T * (T + 2. * CGS::protonMassC2);
  factor /= (T + Phi) * (T + Phi + 2. * CGS::protonMassC2);
  return factor * I_T_interpol(T_ISM);
}

double Particle::I_R_TOA(double R, double modulationPotential) const {
  // see arXiv:1511.08790
  constexpr double mpSquared = pow2(CGS::protonMassC2);
  const double ZOverASquared = pow2(m_pid.getZoverA());
  const double ESquared = pow2(R) * ZOverASquared + mpSquared;
  const double T = std::sqrt(ESquared) - CGS::protonMassC2;
  const double Phi = m_pid.getZoverA() * modulationPotential;
  const double T_ISM = T + Phi;
  double dTdR = R * ZOverASquared;
  dTdR /= std::sqrt(ZOverASquared * pow2(R) + mpSquared);
  double factor = T * (T + 2. * CGS::protonMassC2);
  factor /= (T + Phi) * (T + Phi + 2. * CGS::protonMassC2);
  return factor * I_T_interpol(T_ISM) * dTdR;
}

double Particle::Q_total(double T) const {
  const double Q_ter = (m_pid == H1_ter && m_Q_ter) ? m_Q_ter->get(T) : 0.;
  const double Q_sec = (m_doSecondary && m_Q_sec) ? m_Q_sec->get(T) : 0.;
  const double Q_sec_source = (m_doGrammageAtSource && m_Q_Xs) ? m_Q_Xs->get(T) : 0.;
  const double Q_p = (m_abundance > 0. && m_Q_p) ? m_Q_p->get(T) : 0.;
  return Q_p + Q_sec + Q_sec_source + Q_ter;
}

void Particle::computeIntensity(const Input& input) {
  switch (input.fluxSolver()) {
    case FluxSolver::Analytical:
      for (size_t i = 0; i < m_T.size(); ++i) {
        m_I_T[i] = computeFluxAtEnergy(m_T[i]);
      }
      break;
    case FluxSolver::CrankNicolson:
      computeFluxAtEnergyCrankNicolson();
      break;
    case FluxSolver::Exponential:
      computeFluxAtEnergyExponential();
      break;
  }

  if (Utilities::isGoodAndPositive(m_I_T))
    setDone();
  else
    throw std::runtime_error("Houston, we've had a problem here.");
}

std::string makeParticleFilename(const PID& pid) {
  std::string filename = "output/crams_particle_dump";
  filename += "_" + std::to_string(pid.getZ());
  filename += "_" + std::to_string(pid.getA());
  if (pid.isTertiary()) filename += "_tertiary";
  filename += ".txt";
  return filename;
}

void Particle::dump() const {
  std::ofstream outfile(makeParticleFilename(m_pid));
  if (!outfile.is_open()) throw std::runtime_error("cannot open for writing: " + makeParticleFilename(m_pid));

  writeParticleDumpColumns(outfile);
  outfile << std::scientific;
  for (const auto T : m_T) {
    const double R = Utilities::T2pc(T, m_pid) / std::fabs((double)m_pid.getZ());
    const double Q_p = (m_Q_p) ? m_Q_p->get(T) : 0.;
    const double Q_sec = (m_Q_sec) ? m_Q_sec->get(T) : 0.;
    const double X = (m_X) ? m_X->get(T) / (CGS::gram / CGS::cm2) : kUnavailable;
    const double tauDiff = (m_X) ? m_X->diffusionTimescale(T) / CGS::Myr : kUnavailable;
    const double tauAdv = (m_X) ? m_X->advectionTimescale() / CGS::Myr : kUnavailable;
    const double sigmaISM = (m_sigmaIn) ? m_sigmaIn->getXsecOnISM(m_pid, T) : kUnavailable;
    const double Xcr = (std::isfinite(sigmaISM) && sigmaISM > 0.) ? CGS::meanISMmass / sigmaISM / (CGS::gram / CGS::cm2)
                                                                  : kUnavailable;
    const double dEdX = (m_dEdX) ? m_dEdX->get(T) : kUnavailable;
    const double tauIon = (m_dEdX) ? T / m_dEdX->dTdt_ionization(T, 1. / CGS::cm3) / CGS::Myr : kUnavailable;
    const double tauInelastic =
        (std::isfinite(sigmaISM) && sigmaISM > 0.)
            ? safeReciprocal(Utilities::T2beta(T) * sigmaISM * CGS::cLight / CGS::cm3) / CGS::Myr
            : kUnavailable;

    outfile << T / CGS::GeV << "\t";
    outfile << R / CGS::GeV << "\t";
    outfile << Q_p << "\t";
    outfile << Q_sec << "\t";
    outfile << X << "\t";
    outfile << tauDiff << "\t";
    outfile << tauAdv << "\t";
    outfile << Xcr << "\t";
    outfile << dEdX << "\t";
    outfile << tauIon << "\t";
    outfile << tauInelastic << "\t";
    outfile << "\n";
  }
  outfile.close();
  LOGD << "dumped " << m_pid << " to file " << makeParticleFilename(m_pid);
}

}  // namespace CRAMS
