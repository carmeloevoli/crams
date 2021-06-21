#include "chi2.h"

#include <plog/Log.h>

#define max_num_of_char_in_a_line 512
#define num_of_header_lines 7

namespace CRAMS {

Chi2::Chi2() {}

Chi2::Chi2(const std::string& name, const Particles& particles, const double& phi)
    : m_name(name), m_particles(particles), m_phi(phi) {}

void Chi2::readKISSfile(const std::string& filename, const double& xunits, const double& yunits) {
  LOGI << "reading data from " << filename << "... ";
  std::ifstream file_to_read(filename.c_str());
  if (file_to_read.is_open()) {
    for (int i = 0; i < num_of_header_lines; ++i) file_to_read.ignore(max_num_of_char_in_a_line, '\n');
    file_to_read.ignore(max_num_of_char_in_a_line, '\n');
    double values[6];
    while (!file_to_read.eof()) {
      file_to_read >> values[0] >> values[1] >> values[2] >> values[3] >> values[4] >> values[5];
      dataPoint point;
      point.R = values[0] * xunits;
      point.I = values[1] * yunits;
      point.I_err = std::make_pair<double, double>(values[2] * yunits, values[3] * yunits);
      if (file_to_read.good()) m_data.push_back(point);
    }
  } else {
    throw std::runtime_error("data file cannot be open!");
  }
  LOGI << " with size : " << m_data.size();
  file_to_read.close();
}

double Chi2::computeChi2(const double& R_min, const double& R_max) const {
  double chi2 = 0.0;
  size_t ndata = 0;
  for (auto idata = m_data.begin(); idata != m_data.end(); idata++) {
    if (idata->R > R_min && idata->R < R_max) {
      const double I_R_TOA = getModel(idata->R, m_phi);
      double delta_chi2 = pow2(I_R_TOA - idata->I);
      delta_chi2 /= (I_R_TOA < idata->I) ? pow2(idata->I_err.first) : pow2(idata->I_err.second);
      chi2 += delta_chi2;
      ndata++;
    }
  }
  return chi2 / (double)ndata;
}

double Chi2IH::getModel(const double& R, const double& phi) const {
  double value = (itH1.first) ? itH1.second->I_R_TOA(R, phi) : 0.;
  value += (itH2.first) ? itH2.second->I_R_TOA(R, phi) : 0.;
  value += (itH1_ter.first) ? itH1_ter.second->I_R_TOA(R, phi) : 0.;
  return value;
}

double Chi2IHe::getModel(const double& R, const double& phi) const {
  double value = (itHe3.first) ? itHe3.second->I_R_TOA(R, phi) : 0.;
  value += (itHe4.first) ? itHe4.second->I_R_TOA(R, phi) : 0.;
  return value;
}

double Chi2BC::getModel(const double& R, const double& phi) const {
  double B = (itB10.first) ? itB10.second->I_R_TOA(R, phi) : 0.;
  B += (itB11.first) ? itB11.second->I_R_TOA(R, phi) : 0.;
  double C = (itC12.first) ? itC12.second->I_R_TOA(R, phi) : 0.;
  C += (itC13.first) ? itC13.second->I_R_TOA(R, phi) : 0.;
  C += (itC14.first) ? itC14.second->I_R_TOA(R, phi) : 0.;
  return B / C;
}

}  // namespace CRAMS

// double Chi2_B::get_model(const double& R, const double& phi) const {
// 	double value = (ptr_B10.isPresent) ? ptr_B10.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_B11.isPresent) ? ptr_B11.it->I_R_TOA(R, phi) : 0.;
// 	return value;
// }

// double Chi2_C::get_model(const double& R, const double& phi) const {
// 	double value = (ptr_C12.isPresent) ? ptr_C12.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_C13.isPresent) ? ptr_C13.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_C14.isPresent) ? ptr_C14.it->I_R_TOA(R, phi) : 0.;
// 	return value;
// }

// double Chi2_N::get_model(const double& R, const double& phi) const {
// 	double value = (ptr_N14.isPresent) ? ptr_N14.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_N15.isPresent) ? ptr_N15.it->I_R_TOA(R, phi) : 0.;
// 	return value;
// }

// double Chi2_O::get_model(const double& R, const double& phi) const {
// 	double value = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, phi) : 0.;
// 	return value;
// }

// double Chi2_HeO::get_model(const double& R, const double& phi) const {
// 	double O = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, phi) : 0.;
// 	double He = (ptr_He3.isPresent) ? ptr_He3.it->I_R_TOA(R, phi) : 0.;
// 	He += (ptr_He4.isPresent) ? ptr_He4.it->I_R_TOA(R, phi) : 0.;
// 	return He / O;
// }

// double Chi2_BeB::get_model(const double& R, const double& phi) const {
// 	double B = (ptr_B10.isPresent) ? ptr_B10.it->I_R_TOA(R, phi) : 0.;
// 	B += (ptr_B11.isPresent) ? ptr_B11.it->I_R_TOA(R, phi) : 0.;
// 	double Be = (ptr_Be7.isPresent) ? ptr_Be7.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be9.isPresent) ? ptr_Be9.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be10.isPresent) ? ptr_Be10.it->I_R_TOA(R, phi) : 0.;
// 	return Be / B;
// }

// double Chi2_BeB_statsonly::get_model(const double& R, const double& phi) const {
// 	double B = (ptr_B10.isPresent) ? ptr_B10.it->I_R_TOA(R, phi) : 0.;
// 	B += (ptr_B11.isPresent) ? ptr_B11.it->I_R_TOA(R, phi) : 0.;
// 	double Be = (ptr_Be7.isPresent) ? ptr_Be7.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be9.isPresent) ? ptr_Be9.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be10.isPresent) ? ptr_Be10.it->I_R_TOA(R, phi) : 0.;
// 	return Be / B;
// }

// double Chi2_BeC::get_model(const double& R, const double& phi) const {
// 	double C = (ptr_C12.isPresent) ? ptr_C12.it->I_R_TOA(R, phi) : 0.;
// 	C += (ptr_C13.isPresent) ? ptr_C13.it->I_R_TOA(R, phi) : 0.;
// 	C += (ptr_C14.isPresent) ? ptr_C14.it->I_R_TOA(R, phi) : 0.;
// 	double Be = (ptr_Be7.isPresent) ? ptr_Be7.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be9.isPresent) ? ptr_Be9.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be10.isPresent) ? ptr_Be10.it->I_R_TOA(R, phi) : 0.;
// 	return Be / C;
// }

// double Chi2_BeO::get_model(const double& R, const double& phi) const {
// 	double O = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, phi) : 0.;
// 	double Be = (ptr_Be7.isPresent) ? ptr_Be7.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be9.isPresent) ? ptr_Be9.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be10.isPresent) ? ptr_Be10.it->I_R_TOA(R, phi) : 0.;
// 	return Be / O;
// }

// double Chi2_BO::get_model(const double& R, const double& phi) const {
// 	double B = (ptr_B10.isPresent) ? ptr_B10.it->I_R_TOA(R, phi) : 0.;
// 	B += (ptr_B11.isPresent) ? ptr_B11.it->I_R_TOA(R, phi) : 0.;
// 	double O = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, phi) : 0.;
// 	return B / O;
// }

// double Chi2_CO::get_model(const double& R, const double& phi) const {
// 	double O = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, phi) : 0.;
// 	double C = (ptr_C12.isPresent) ? ptr_C12.it->I_R_TOA(R, phi) : 0.;
// 	C += (ptr_C13.isPresent) ? ptr_C13.it->I_R_TOA(R, phi) : 0.;
// 	C += (ptr_C14.isPresent) ? ptr_C14.it->I_R_TOA(R, phi) : 0.;
// 	return C / O;
// }
