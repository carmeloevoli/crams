#include "particle.h"

#include <plog/Log.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "cgs.h"
#include "utilities.h"
#include "xsecs/Evoli2019.h"

namespace CRAMS {

#define ARRAYSIZE 400

Particle::Particle(const PID& pid, const NucleusParameters& nucleusParameters)
    : m_pid(pid),
      m_abundance(nucleusParameters.abundance),
      m_slope(nucleusParameters.slope),
      m_decayTime(nucleusParameters.decayTime) {}

Particle::Particle(const PID& pid) : m_pid(pid) {}

Particle::~Particle() { LOGD << "released memory of particle " << m_pid; }

void Particle::reset() {
  m_X.reset();
  m_Q_p.reset();
  m_Q_sec.reset();
  m_Q_Xs.reset();
  m_sigmaIn.reset();
  m_dEdX.reset();
}

void Particle::buildVectors(const Input& input) {
  m_T = Utilities::LogAxis(input.TSimMin, input.TSimMax, input.TSimSize);
  m_I_T.resize(input.TSimSize);
}

void Particle::buildGrammage(const Input& input) {
  if (isStable())
    m_X = std::make_shared<Grammage>(m_pid, input);
  else
    m_X = std::make_shared<Grammage>(m_pid, input, m_decayTime);
}

void Particle::buildPrimarySource(const Input& input) {
  m_Q_p = std::make_shared<PrimarySource>(m_pid, m_abundance, m_slope, input.mu);
}

void Particle::buildLosses(const Input& input) { m_dEdX = std::make_shared<Losses>(m_pid, input); }

void Particle::buildInelasticXsecs(const Input& input) {
  m_sigmaIn = (input.id == 0) ? std::make_shared<InXsecTripathi99>(m_pid, false)
                              : std::make_shared<InXsecTripathi99>(m_pid, true);
}

#define COTH(A) (1. / std::tanh(A))

double Particle::productionProfileFromUnstable(const Input& input, const double& T, const double& decayTimeAtRest) {
  const double v = Utilities::T2beta(T) * CGS::cLight;
  const double u = input.v_A;
  const double H = input.H;
  const double D = m_X->D(T);
  const double value = 2. * u / input.mu / v;
  const double decayTimeAtT = Utilities::T2gamma(T) * decayTimeAtRest;
  const double Delta = std::sqrt(1. + 4. * D / (pow2(u) * decayTimeAtT));
  const double profile = Delta * COTH(u * H * Delta / 2. / D) - COTH(u * H / 2. / D);
  return value * profile;
}

void Particle::buildSecondarySource(const Input& input, const std::vector<Particle>& particles) {
  m_doSecondary = input.doSecondary;
  const auto id = input.id;
  const auto fudge = input.xsecsFudge;
  const auto xsecs = (id == 0) ? SpallationXsecs(m_pid, fudge) : SpallationXsecs(m_pid, fudge, true);
  const auto T_s = Utilities::LogAxis(0.1 * CGS::GeV, 10. * CGS::TeV, ARRAYSIZE);
  std::vector<double> Q_s;
  for (auto& T : T_s) {
    double value = 0;
    for (auto& particle : particles) {
      if (particle.getPid().getA() > m_pid.getA() && particle.isDone()) {
        value += xsecs.getXsecOnISM(particle.getPid(), T) * particle.I_T_interpol(T);
      }
    }
    value /= CGS::meanISMmass;
    Q_s.push_back(value);
  }

  if (m_pid == B10) {
    const auto ptr_Be10 = std::find(particles.begin(), particles.end(), Particle(Be10));
    const auto tau_Be10 = ptr_Be10->getDecayTime();
    size_t counter = 0;
    for (auto& T : T_s) {
      const double Q_Be10 = productionProfileFromUnstable(input, T, tau_Be10) * ptr_Be10->I_T_interpol(T);
      Q_s.at(counter) += Q_Be10;
      counter++;
    }
  }

  if (m_pid == N14) {
    const auto ptr_C14 = std::find(particles.begin(), particles.end(), Particle(C14));
    const auto tau_C14 = ptr_C14->getDecayTime();
    size_t counter = 0;
    for (auto& T : T_s) {
      const double Q_C14 = productionProfileFromUnstable(input, T, tau_C14) * ptr_C14->I_T_interpol(T);
      Q_s.at(counter) += Q_C14;
      counter++;
    }
  }

  if (m_pid == Mg26) {
    const auto ptr_Al26 = std::find(particles.begin(), particles.end(), Particle(Al26));
    const auto tau_Al26 = ptr_Al26->getDecayTime();
    size_t counter = 0;
    for (auto& T : T_s) {
      const double Q_Al26 = productionProfileFromUnstable(input, T, tau_Al26) * ptr_Al26->I_T_interpol(T);
      Q_s.at(counter) += Q_Al26;
      counter++;
    }
  }

  if (m_pid == Ar36) {
    const auto ptr_Cl36 = std::find(particles.begin(), particles.end(), Particle(Cl36));
    const auto tau_Cl36 = ptr_Cl36->getDecayTime();
    size_t counter = 0;
    for (auto& T : T_s) {
      const double Q_Cl36 = productionProfileFromUnstable(input, T, tau_Cl36) * ptr_Cl36->I_T_interpol(T);
      Q_s.at(counter) += Q_Cl36;
      counter++;
    }
  }

  if (m_pid == Fe54) {
    const auto ptr_Mn54 = std::find(particles.begin(), particles.end(), Particle(Mn54));
    const auto tau_Mn54 = ptr_Mn54->getDecayTime();
    size_t counter = 0;
    for (auto& T : T_s) {
      const double Q_Mn54 = productionProfileFromUnstable(input, T, tau_Mn54) * ptr_Mn54->I_T_interpol(T);
      Q_s.at(counter) += Q_Mn54;
      counter++;
    }
  }

  m_Q_sec = std::make_shared<SecondarySource>(m_pid, T_s, Q_s);
}

// void Particle::build_tertiary_source(const std::vector<Particle>& particles) {
//   auto T_t = LogAxis(0.1 * cgs::GeV, 10. * cgs::TeV, ARRAYSIZE);
//   std::vector<double> Q_t;
//   for (auto& T : T_t) {
//     const double T_prime = T / cgs::inelasticity;
//     double sigma_ISM = sigma_pp(T_prime);
//     sigma_ISM *= (1. + cgs::K_He * cgs::f_He) / (1. + cgs::f_He);
//     double value = sigma_ISM / cgs::inelasticity;
//     value *= (T_prime + cgs::proton_mass_c2) / (T + cgs::proton_mass_c2);
//     value *= std::pow(T * (T + 2. * cgs::proton_mass_c2), 1.5) /
//              std::pow(T_prime * (T_prime + 2. * cgs::proton_mass_c2), 1.5);
//     for (auto& particle : particles) {
//       if (particle.get_pid() == H1 && particle.isDone()) {
//         value *= particle.I_T_interpol(T_prime);
//       }
//     }
//     value /= cgs::mean_ism_mass;
//     Q_t.push_back(value);
//   }
//   _Q_ter = new SourceTerm(T_t, Q_t);
// }

void Particle::buildGrammageAtSource(const Input& input, const std::vector<Particle>& particles) {
  m_doGrammageAtSource = true;
  const auto id = input.id;
  const auto fudge = input.xsecsFudge;
  const auto xsecs = (id == 0) ? SpallationXsecs(m_pid, fudge) : SpallationXsecs(m_pid, fudge, true);
  const auto T_X = Utilities::LogAxis(0.1 * CGS::GeV, 10. * CGS::TeV, ARRAYSIZE);
  std::vector<double> Q_X;
  for (const auto& T : T_X) {
    double value = 0;
    for (const auto& particle : particles) {
      if (particle.getPid().getA() > m_pid.getA() && particle.isDone()) {
        const double r = input.X_s / CGS::meanISMmass * xsecs.getXsecOnISM(particle.getPid(), T);
        const auto Q_p = PrimarySource(particle.getPid(), particle.getAbundance(), particle.getSlope(), input.mu);
        value += r * Q_p.get(T);
      }
    }
    Q_X.push_back(value);
  }
  m_Q_Xs = std::make_shared<SecondarySource>(m_pid, T_X, Q_X);
}

double Particle::I_T_interpol(const double& T) const {
  double value = 0;
  if (T > m_T.front() && T < m_T.back()) {
    value = Utilities::LinearInterpolatorLog(m_T, m_I_T, T);
  }
  return value;
}

double Particle::Q_total(const double& T) const {
  //  if (_pid == H1_ter) value = _Q_ter->get(T); TODO do this
  const double Q_sec_total = m_Q_sec->get(T) + ((m_doGrammageAtSource) ? m_Q_Xs->get(T) : 0.);
  return m_Q_p->get(T) + ((m_doSecondary) ? Q_sec_total : 0.);
}

void Particle::computeIntensity() {
#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
  for (size_t i = 0; i < m_T.size(); ++i) {
    m_I_T[i] = computeFluxAtEnergy(m_T[i]);
  }

  if (Utilities::isGoodAndPositive(m_I_T))
    setDone();
  else
    throw std::runtime_error("Houston, we've had a problem here.");
}

std::string makeParticleFilename(const PID& pid) {
  std::string filename = "output/crams_particle_dump";
  filename += "_" + std::to_string(pid.getZ());
  filename += "_" + std::to_string(pid.getA()) + ".txt";
  return filename;
}

void Particle::dump() const {
  std::ofstream outfile(makeParticleFilename(m_pid));
  outfile << "# T [GeV] - R [GV] - Q_pri - Q_sec - X [gr/cm2] - tau_esc [yr] - tau_adv [yr] - X_cr [gr/cm2] - dEdX\n";
  outfile << std::scientific;
  for (double T = CGS::GeV; T < 1.1 * CGS::TeV; T *= 1.1) {
    outfile << T / CGS::GeV << "\t";
    double R = Utilities::T2pc(T, m_pid) / (double)m_pid.getZ();
    outfile << R / CGS::GeV << "\t";
    outfile << m_Q_p->get(T) << "\t";
    outfile << m_Q_sec->get(T) << "\t";
    outfile << m_X->get(T) / (CGS::gram / CGS::cm2) << "\t";
    outfile << m_X->diffusionTimescale(T) / CGS::year << "\t";
    outfile << m_X->advectionTimescale() / CGS::year << "\t";
    outfile << CGS::meanISMmass / m_sigmaIn->getXsecOnISM(T) / (CGS::gram / CGS::cm2) << "\t";
    outfile << m_dEdX->get(T) << "\t";  // TODO how to transform in time?
    outfile << "\n";
  }
  outfile.close();
}

}  // namespace CRAMS