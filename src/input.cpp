#include "input.h"

#include <plog/Log.h>

#include <fstream>
#include <iostream>
#include <sstream>

#include "utilities.h"

namespace CRAMS {

Input::~Input() { LOGD << "released memory from Input"; }

void Input::setParam(const std::string& key, const double& value) {
  const auto simpleKey = Utilities::simplifyKey(key);
  if (key == "d0") {
    m_D_0 = value * 1e28 * CGS::cm2 / CGS::sec;
    LOGD << "changed D_0 value to " << m_D_0 / CGS::cm2 * CGS::sec << " cm2/s";
  } else if (key == "xs") {
    m_X_s = value * CGS::gram / CGS::cm2;
    LOGD << "changed X_s value to " << m_X_s / CGS::gram * CGS::cm2 << " gr/cm2";
  } else if (key == "h") {
    m_H = value * CGS::kpc;
    LOGD << "changed H value to " << m_H / CGS::kpc << " kpc";
  } else if (key == "delta") {
    m_delta = value;
    LOGD << "chanded delta value to " << m_delta;
  } else if (key == "ddelta") {
    m_ddelta = value;
    LOGD << "chanded ddelta value to " << m_ddelta;
  } else if (key == "rb") {
    m_R_b = value * CGS::GeV;
    LOGD << "chanded R_b value to " << m_R_b / CGS::GeV << " GeV";
  } else if (key == "va") {
    m_v_A = value * CGS::km / CGS::sec;
    LOGD << "changed v_A value to " << m_v_A / CGS::km * CGS::sec << " km/s";
  } else if (key == "phi") {
    m_modulationPotential = value * CGS::GeV;
    LOGD << "changed modulationPotential value to " << m_modulationPotential / CGS::GeV << " GV";
  } else if (key == "xsecsfudge") {
    m_xsecsFudge = value;
    LOGD << "changed xsecs fudge value to " << m_xsecsFudge;
  } else if (key == "id")
    m_id = (int)value;
}

void Input::readParamsFromFile(const std::string& filename) {
  std::ifstream infile(filename.c_str());
  std::string line;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::string key;
    double value;
    if (!(iss >> key >> value)) {
      break;
    }  // error
    setParam(key, value);
  }
}

void eraseSubStr(std::string* main_string, const std::string& sub_string) {
  // Search for the substring in string
  size_t pos = main_string->find(sub_string);
  if (pos != std::string::npos) {
    main_string->erase(pos, sub_string.length());
  }
}

void Input::setSimname(const std::string& params_filename) {
  m_simname = params_filename;
  const auto size = m_simname.length();
  eraseSubStr(&m_simname, ".ini");
  if (m_simname.length() != size - 4) throw std::runtime_error("Problem with the input filename!");
}

void Input::print() const {
  LOGD << "H      [kpc]        : " << std::setprecision(4) << m_H / CGS::kpc;
  LOGD << "mu     [mg/cm2]     : " << std::setprecision(4) << m_mu / (CGS::mgram / CGS::cm2);
  LOGD << "v_A    [km/s]       : " << std::setprecision(4) << m_v_A / (CGS::km / CGS::sec);
  LOGD << "D_0    [1e28 cm2/s] : " << std::setprecision(4) << m_D_0 / (1e28 * CGS::cm2 / CGS::sec);
  LOGD << "delta  []           : " << std::setprecision(4) << delta;
  LOGD << "ddelta []           : " << std::setprecision(4) << m_ddelta;
  if (m_X_s > 0.)
    LOGD << "X_s    [gr/cm2]" << m_X_s / (CGS::gram / CGS::cm2);
  else
    LOGD << "X_s    [gr/cm2]     : none";
  LOGD << "R_b    [GV]         : " << m_R_b / CGS::GeV;
  LOGD << "s      []           : " << m_smoothness;
  LOGD << "phi    [GeV]        : " << m_modulationPotential / CGS::GeV;
  LOGD << "xsecs_f[]           : " << m_xsecsFudge;
  LOGD << "E_min  [GeV]        : " << m_TSimMin / CGS::GeV;
  LOGD << "E_max  [GeV]        : " << m_TSimMax / CGS::GeV;
  LOGD << "E_size []           : " << m_TSimSize;
  LOGD << "do Secondary        : " << std::boolalpha << m_doSecondary;
}

}  // namespace CRAMS