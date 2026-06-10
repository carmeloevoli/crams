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

}  // namespace

namespace CRAMS {

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
    LOGD << "changed R_b to " << m_R_b / CGS::GeV << " GeV";
  } else if (key == "va") {
    m_v_A = value * CGS::km / CGS::sec;
    LOGD << "changed v_A to " << m_v_A / (CGS::km / CGS::sec) << " km/s";
  } else if (key == "phi") {
    m_modulationPotential = value * CGS::GeV;
    LOGD << "changed phi to " << m_modulationPotential / CGS::GeV << " GV";
  } else if (key == "xsecsfudge") {
    m_xsecsFudge = value;
    LOGD << "changed xsecsFudge to " << m_xsecsFudge;
  } else if (key == "num") {
    m_num = static_cast<bool>(value);
    LOGD << "changed num to " << m_num;
  } else if (key == "id") {
    m_id = static_cast<size_t>(value);
  }
}

void Input::readParamsFromFile(const std::string& filename) {
  std::ifstream infile(filename);
  if (!infile.is_open())
    throw std::runtime_error("Input: cannot open file '" + filename + "'");
  std::string line;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::string key;
    double value;
    if (!(iss >> key >> value)) continue;
    setParam(key, value);
  }
}

void Input::setSimname(const std::string& inifilename) {
  m_simname = inifilename;
  eraseExtension(m_simname, ".ini");
}

void Input::print() const {
  LOGD << "H      [kpc]        : " << std::setprecision(4) << m_H / CGS::kpc;
  LOGD << "mu     [mg/cm2]     : " << std::setprecision(4) << m_mu / (CGS::mgram / CGS::cm2);
  LOGD << "v_A    [km/s]       : " << std::setprecision(4) << m_v_A / (CGS::km / CGS::sec);
  LOGD << "D_0    [1e28 cm2/s] : " << std::setprecision(4) << m_D_0 / (1e28 * CGS::cm2 / CGS::sec);
  LOGD << "delta  []           : " << std::setprecision(4) << m_delta;
  LOGD << "ddelta []           : " << std::setprecision(4) << m_ddelta;
  if (m_X_s > 0.)
    LOGD << "X_s    [g/cm2]      : " << m_X_s / (CGS::gram / CGS::cm2);
  else
    LOGD << "X_s    [g/cm2]      : none";
  LOGD << "R_b    [GV]         : " << m_R_b / CGS::GeV;
  LOGD << "s      []           : " << m_smoothness;
  LOGD << "phi    [GeV]        : " << m_modulationPotential / CGS::GeV;
  LOGD << "xsecs_f[]           : " << m_xsecsFudge;
  LOGD << "E_min  [GeV]        : " << m_TSimMin / CGS::GeV;
  LOGD << "E_max  [GeV]        : " << m_TSimMax / CGS::GeV;
  LOGD << "E_size []           : " << m_TSimSize;
  LOGD << "doSecondary         : " << std::boolalpha << m_doSecondary;
  LOGD << "numerical method    : " << std::boolalpha << m_num;
}

}  // namespace CRAMS
