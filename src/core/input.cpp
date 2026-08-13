#include "crams/core/input.h"

#include <plog/Log.h>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "crams/utils/utilities.h"

namespace {

void eraseExtension(std::string& s, const std::string& ext) {
  const auto pos = s.rfind(ext);
  if (pos != std::string::npos && pos == s.length() - ext.length())
    s.erase(pos);
  else
    throw std::runtime_error("Input filename must end with '" + ext + "'");
}

std::string fluxSolverName(CRAMS::FluxSolver solver) {
  switch (solver) {
    case CRAMS::FluxSolver::Analytical:
      return "analytical";
    case CRAMS::FluxSolver::CrankNicolson:
      return "crank_nicolson";
    case CRAMS::FluxSolver::Exponential:
      return "exponential";
  }

  return "unknown";
}

std::string inelasticModelName(CRAMS::InelasticModel model) {
  switch (model) {
    case CRAMS::InelasticModel::Tripathi99:
      return "tripathi99";
    case CRAMS::InelasticModel::Glauber:
      return "glauber";
    case CRAMS::InelasticModel::Crosec:
      return "crosec";
  }

  return "unknown";
}

std::string fragmentationModelName(CRAMS::FragmentationModel model) {
  switch (model) {
    case CRAMS::FragmentationModel::Fluka4Dragon:
      return "fluka4dragon";
    case CRAMS::FragmentationModel::UsineGalprop17Opt12:
      return "usine_galprop17_opt12";
    case CRAMS::FragmentationModel::UsineGalprop17Opt22:
      return "usine_galprop17_opt22";
    case CRAMS::FragmentationModel::UsineWebber03Coste12:
      return "usine_webber03_coste12";
    case CRAMS::FragmentationModel::Evoli2019:
      return "evoli2019";
    case CRAMS::FragmentationModel::Evoli2026W93:
      return "evoli2026_w93";
    case CRAMS::FragmentationModel::Evoli2026St99:
      return "evoli2026_st99";
  }

  return "unknown";
}

}  // namespace

namespace CRAMS {

FluxSolver parseFluxSolver(const std::string& value) {
  const auto solver = CRAMS::Utilities::simplifyKey(value);
  if (solver == "analytical") return CRAMS::FluxSolver::Analytical;
  if (solver == "cranknicolson") return CRAMS::FluxSolver::CrankNicolson;
  if (solver == "exponential") return CRAMS::FluxSolver::Exponential;

  throw std::runtime_error("Input: unknown flux solver '" + value + "'");
}

InelasticModel parseInelasticModel(const std::string& value) {
  const auto model = CRAMS::Utilities::simplifyKey(value);
  if (model == "tripathi99" || model == "tripathi1999") return CRAMS::InelasticModel::Tripathi99;
  if (model == "glauber") return CRAMS::InelasticModel::Glauber;
  if (model == "crosec") return CRAMS::InelasticModel::Crosec;

  throw std::runtime_error("Input: unknown inelastic model '" + value + "'");
}

FragmentationModel parseFragmentationModel(const std::string& value) {
  const auto model = CRAMS::Utilities::simplifyKey(value);
  if (model == "fluka4dragon") return CRAMS::FragmentationModel::Fluka4Dragon;
  if (model == "usinegalprop17opt12") return CRAMS::FragmentationModel::UsineGalprop17Opt12;
  if (model == "usinegalprop17opt22") return CRAMS::FragmentationModel::UsineGalprop17Opt22;
  if (model == "usinewebber03coste12") return CRAMS::FragmentationModel::UsineWebber03Coste12;
  if (model == "evoli2019") return CRAMS::FragmentationModel::Evoli2019;
  if (model == "evoli2026w93") return CRAMS::FragmentationModel::Evoli2026W93;
  if (model == "evoli2026st99") return CRAMS::FragmentationModel::Evoli2026St99;

  throw std::runtime_error("Input: unknown fragmentation model '" + value + "'");
}

std::string Input::fluxSolverName() const { return ::fluxSolverName(m_fluxSolver); }

std::string Input::inelasticModelName() const { return ::inelasticModelName(m_inelasticModel); }

std::string Input::fragmentationModelName() const { return ::fragmentationModelName(m_fragmentationModel); }

void Input::setParam(const std::string& KEY, double value) {
  const auto key = Utilities::simplifyKey(KEY);
  if (key == "d0") {
    m_D_0 = value * 1e28 * CGS::cm2 / CGS::sec;
    LOGD << "changed D_0 to " << m_D_0 / (CGS::cm2 / CGS::sec) << " cm2/s";
  } else if (key == "xs") {
    m_X_s = value * CGS::gram / CGS::cm2;
    LOGD << "changed X_s to " << m_X_s / (CGS::gram / CGS::cm2) << " g/cm2";
  } else if (key == "h") {
    m_H = value * CGS::kpc;
    LOGD << "changed H to " << m_H / CGS::kpc << " kpc";
  } else if (key == "delta") {
    m_delta = value;
    LOGD << "changed delta to " << m_delta;
  } else if (key == "ddelta") {
    m_ddelta = value;
    LOGD << "changed ddelta to " << m_ddelta;
  } else if (key == "rb") {
    m_R_b = value * CGS::GeV;
    LOGD << "changed R_b to " << m_R_b / CGS::GeV << " GV";
  } else if (key == "va") {
    m_v_A = value * CGS::km / CGS::sec;
    LOGD << "changed v_A to " << m_v_A / (CGS::km / CGS::sec) << " km/s";
  } else if (key == "phi") {
    m_modulationPotential = value * CGS::GeV;
    LOGD << "changed phi to " << m_modulationPotential / CGS::GeV << " GV";
  } else if (key == "fudgebe7") {  // .ini key "fudge_be7" (simplifyKey strips '_')
    m_fudgeBe7 = value;
    LOGD << "changed fudge_Be7 to " << m_fudgeBe7;
  } else if (key == "fudgebe9") {
    m_fudgeBe9 = value;
    LOGD << "changed fudge_Be9 to " << m_fudgeBe9;
  } else if (key == "fudgebe10") {
    m_fudgeBe10 = value;
    LOGD << "changed fudge_Be10 to " << m_fudgeBe10;
  } else if (key == "sourcefeaturer") {
    m_sourceSpectrumFeatureR = value * CGS::GeV;
    LOGD << "changed source feature R to " << value << " GV";
  } else if (key == "sourcebreakdslope") {
    m_sourceSpectrumBreakDeltaSlope = value;
    LOGD << "changed source break delta slope to " << value;
  } else if (key == "sourcebreakomega") {
    m_sourceSpectrumBreakOmega = value;
    LOGD << "changed source break omega to " << value;
  } else if (key == "id") {
    m_id = static_cast<size_t>(value);
  }
}

void Input::readParamsFromFile(const std::string& filename) {
  std::ifstream infile(filename);
  if (!infile.is_open()) throw std::runtime_error("Input: cannot open file '" + filename + "'");
  std::string line;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::string key;
    std::string valueToken;
    if (!(iss >> key >> valueToken)) continue;

    if (Utilities::simplifyKey(key) == "solver") {
      m_fluxSolver = parseFluxSolver(valueToken);
      LOGD << "changed flux solver to " << fluxSolverName();
      continue;
    }

    if (Utilities::simplifyKey(key) == "inelasticmodel") {
      m_inelasticModel = parseInelasticModel(valueToken);
      LOGD << "changed inelastic model to " << inelasticModelName();
      continue;
    }

    if (Utilities::simplifyKey(key) == "fragmentationmodel") {
      m_fragmentationModel = parseFragmentationModel(valueToken);
      LOGD << "changed fragmentation model to " << fragmentationModelName();
      continue;
    }

    double value = 0.;
    try {
      value = std::stod(valueToken);
    } catch (const std::exception&) {
      continue;
    }
    setParam(key, value);
  }
}

void Input::setSimname(const std::string& inifilename) {
  m_simname = inifilename;
  // Drop any leading directory components so the output path
  // (output/<simname>_...) stays valid no matter where the .ini lives.
  const auto slash = m_simname.find_last_of("/\\");
  if (slash != std::string::npos) m_simname.erase(0, slash + 1);
  eraseExtension(m_simname, ".ini");
}

std::string Input::describe() const {
  std::stringstream out;
  out << "H      [kpc]        : " << std::setprecision(4) << m_H / CGS::kpc << std::endl;
  out << "mu     [mg/cm2]     : " << std::setprecision(4) << m_mu / (CGS::mgram / CGS::cm2) << std::endl;
  out << "v_A    [km/s]       : " << std::setprecision(4) << m_v_A / (CGS::km / CGS::sec) << std::endl;
  out << "D_0    [1e28 cm2/s] : " << std::setprecision(4) << m_D_0 / (1e28 * CGS::cm2 / CGS::sec) << std::endl;
  out << "delta  []           : " << std::setprecision(4) << m_delta << std::endl;
  out << "ddelta []           : " << std::setprecision(4) << m_ddelta << std::endl;
  if (m_X_s > 0.)
    out << "X_s    [g/cm2]      : " << m_X_s / (CGS::gram / CGS::cm2) << std::endl;
  else
    out << "X_s    [g/cm2]      : none" << std::endl;
  out << "R_b    [GV]         : " << m_R_b / CGS::GeV << std::endl;
  out << "s      []           : " << m_smoothness << std::endl;
  out << "phi    [GeV]        : " << m_modulationPotential / CGS::GeV << std::endl;
  if (m_fudgeBe7 != 1. || m_fudgeBe9 != 1. || m_fudgeBe10 != 1.)
    out << "fudge Be7/9/10      : " << m_fudgeBe7 << " / " << m_fudgeBe9 << " / " << m_fudgeBe10 << std::endl;
  out << "E_min  [GeV]        : " << m_TSimMin / CGS::GeV << std::endl;
  out << "E_max  [GeV]        : " << m_TSimMax / CGS::GeV << std::endl;
  out << "E_size []           : " << m_TSimSize << std::endl;
  out << "doSecondary         : " << std::boolalpha << m_doSecondary << std::endl;
  out << "flux solver         : " << fluxSolverName() << std::endl;
  out << "inelastic model     : " << inelasticModelName() << std::endl;
  out << "fragmentation model : " << fragmentationModelName() << std::endl;

  return out.str();
}

void Input::print() const { LOGD << describe(); }

}  // namespace CRAMS
