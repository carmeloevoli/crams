#include "xsecs/Korsmeier2018.h"

#include "cgs.h"
#include "gsl.h"
#include "utilities.h"

#define NHEADERLINES 44

namespace CRAMS {

Korsmeier2018SecAp::Korsmeier2018SecAp() {
  if (!Utilities::fileExists(datafile))
    throw std::runtime_error("error in reading an input file for Korsmeier2018 model.");
  readDataFiles();
}

void Korsmeier2018SecAp::readDataFiles() {
  double TprojMin = 1.0 * CGS::GeV;
  double TprojMax = 1e7 * CGS::GeV;
  size_t TprojSize = 30 * 7 + 1;
  m_lgTprojAxis = Utilities::LinAxis(std::log(TprojMin), std::log(TprojMax), TprojSize);

  double TsecMin = 0.1 * CGS::GeV;
  double TsecMax = 1e4 * CGS::GeV;
  size_t TsecSize = 30 * 5 + 1;
  m_lgTapAxis = Utilities::LinAxis(std::log(TsecMin), std::log(TsecMax), TsecSize);

  m_sigma_pp = Grid<double>(TprojSize, TsecSize);
  m_sigma_pHe = Grid<double>(TprojSize, TsecSize);
  m_sigma_dp = Grid<double>(TprojSize, TsecSize);
  m_sigma_dHe = Grid<double>(TprojSize, TsecSize);
  m_sigma_He3p = Grid<double>(TprojSize, TsecSize);
  m_sigma_He3He = Grid<double>(TprojSize, TsecSize);
  m_sigma_He4p = Grid<double>(TprojSize, TsecSize);
  m_sigma_He4He = Grid<double>(TprojSize, TsecSize);
  m_sigma_C12p = Grid<double>(TprojSize, TsecSize);
  m_sigma_C12He = Grid<double>(TprojSize, TsecSize);
  m_sigma_O16p = Grid<double>(TprojSize, TsecSize);
  m_sigma_O16He = Grid<double>(TprojSize, TsecSize);

  double units = CGS::m2 / CGS::GeV;
  {
    const auto sigma = Utilities::loadColumn(datafile, 2, NHEADERLINES);
    m_sigma_pp.copy(sigma);
    m_sigma_pp.for_each([units](double& s) { s *= units; });
  }
  {
    const auto sigma = Utilities::loadColumn(datafile, 3, NHEADERLINES);
    m_sigma_pHe.copy(sigma);
    m_sigma_pHe.for_each([units](double& s) { s *= units; });
  }
  {
    const auto sigma = Utilities::loadColumn(datafile, 4, NHEADERLINES);
    m_sigma_dp.copy(sigma);
    m_sigma_dp.for_each([units](double& s) { s *= units; });
  }
  {
    const auto sigma = Utilities::loadColumn(datafile, 5, NHEADERLINES);
    m_sigma_dHe.copy(sigma);
    m_sigma_dHe.for_each([units](double& s) { s *= units; });
  }
  {
    const auto sigma = Utilities::loadColumn(datafile, 6, NHEADERLINES);
    m_sigma_He3p.copy(sigma);
    m_sigma_He3p.for_each([units](double& s) { s *= units; });
  }
  {
    const auto sigma = Utilities::loadColumn(datafile, 7, NHEADERLINES);
    m_sigma_He3He.copy(sigma);
    m_sigma_He3He.for_each([units](double& s) { s *= units; });
  }
  {
    const auto sigma = Utilities::loadColumn(datafile, 8, NHEADERLINES);
    m_sigma_He4p.copy(sigma);
    m_sigma_He4p.for_each([units](double& s) { s *= units; });
  }
  {
    const auto sigma = Utilities::loadColumn(datafile, 9, NHEADERLINES);
    m_sigma_He4He.copy(sigma);
    m_sigma_He4He.for_each([units](double& s) { s *= units; });
  }
  {
    const auto sigma = Utilities::loadColumn(datafile, 10, NHEADERLINES);
    m_sigma_C12p.copy(sigma);
    m_sigma_C12p.for_each([units](double& s) { s *= units; });
  }
  {
    const auto sigma = Utilities::loadColumn(datafile, 11, NHEADERLINES);
    m_sigma_C12He.copy(sigma);
    m_sigma_C12He.for_each([units](double& s) { s *= units; });
  }
  {
    const auto sigma = Utilities::loadColumn(datafile, 18, NHEADERLINES);
    m_sigma_O16p.copy(sigma);
    m_sigma_O16p.for_each([units](double& s) { s *= units; });
  }
  {
    const auto sigma = Utilities::loadColumn(datafile, 19, NHEADERLINES);
    m_sigma_O16He.copy(sigma);
    m_sigma_O16He.for_each([units](double& s) { s *= units; });
  }
}

double Korsmeier2018SecAp::get(PbarChannel ch, const double& T_proj, const double& T_ap) const {
  using std::log;
  double value = 0;

  const auto lgTproj = log(T_proj);
  const auto lgTprojRange = std::make_pair(m_lgTprojAxis.front(), m_lgTprojAxis.back());

  const auto lgTap = log(T_ap);
  const auto lgTapRange = std::make_pair(m_lgTapAxis.front(), m_lgTapAxis.back());

  if (Utilities::inRange(lgTproj, lgTprojRange) && Utilities::inRange(lgTap, lgTapRange)) {
    if (ch == PbarChannel::pp) {
      auto z = m_sigma_pp.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else if (ch == PbarChannel::pHe) {
      auto z = m_sigma_pHe.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else if (ch == PbarChannel::dp) {
      auto z = m_sigma_dp.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else if (ch == PbarChannel::dHe) {
      auto z = m_sigma_dHe.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else if (ch == PbarChannel::He3p) {
      auto z = m_sigma_He3p.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else if (ch == PbarChannel::He3He) {
      auto z = m_sigma_He3He.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else if (ch == PbarChannel::He4p) {
      auto z = m_sigma_He4p.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else if (ch == PbarChannel::He4He) {
      auto z = m_sigma_He4He.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else if (ch == PbarChannel::C12p) {
      auto z = m_sigma_C12p.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else if (ch == PbarChannel::C12He) {
      auto z = m_sigma_C12He.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else if (ch == PbarChannel::O16p) {
      auto z = m_sigma_O16p.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else if (ch == PbarChannel::O16He) {
      auto z = m_sigma_O16He.get();
      value = GSL::interpolate2d<double>(m_lgTprojAxis, m_lgTapAxis, z, lgTproj, lgTap);
    } else {
      throw std::runtime_error("channel not implemented in Korsmeier2018 model");
    }
  }
  return std::max(value, 0.);
}

}  // namespace CRAMS
