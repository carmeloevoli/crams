#ifndef INCLUDE_SECAP_H_
#define INCLUDE_SECAP_H_

#include <string>

#include "include/grid.h"

namespace CRAMS {

enum class PbarChannel { pp, pHe, Hep, HeHe };

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
  Grid<double> m_sigma_Hep;
  Grid<double> m_sigma_pHe;
  Grid<double> m_sigma_HeHe;
};

}  // namespace CRAMS

#endif