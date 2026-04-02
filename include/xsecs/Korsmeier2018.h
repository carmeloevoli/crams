#ifndef INCLUDE_SECAP_H_
#define INCLUDE_SECAP_H_

#include <string>

#include "include/grid.h"

namespace CRAMS {

enum class PbarChannel { pp, pHe, dp, dHe, He3p, He3He, He4p, He4He, C12p, C12He, O16p, O16He };

class Korsmeier2018SecAp {
 public:
  Korsmeier2018SecAp();

  double get(PbarChannel ch, const double &T_n, const double &T_ap) const;

 private:
  void init();
  void checkDatafilesExist();
  void readDataFiles();

 private:
  const std::string datafile = "data/Korsmeier2018/supplementary__XS_table_Param_II_B.txt";

  std::vector<double> m_lgTprojAxis;
  std::vector<double> m_lgTapAxis;

  Grid<double> m_sigma_pp;
  Grid<double> m_sigma_pHe;
  Grid<double> m_sigma_dp;
  Grid<double> m_sigma_dHe;
  Grid<double> m_sigma_He3p;
  Grid<double> m_sigma_He3He;
  Grid<double> m_sigma_He4p;
  Grid<double> m_sigma_He4He;
  Grid<double> m_sigma_C12p;
  Grid<double> m_sigma_C12He;
  Grid<double> m_sigma_O16p;
  Grid<double> m_sigma_O16He;
};

}  // namespace CRAMS

#endif