#include "particle.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

//double Particle::f(double T_, double Y) {
//	return Q->get(T_) - Y / -dEdx->get(T_) * (1. / X->get(T_) + sigma->get(T_));
//}
//
//int Particle::run(const LogAxis& T_) {
//	_T = T_.get_axis();
//	_I_T.resize(_T.size());
//	double Y_n = 0;
//	for (size_t i = _T.size() - 1; i > 1; --i) {
//		for (double T_n = _T.at(i); T_n > _T.at(i - 1); T_n /= 1.01) {
//			double h = T_n * (1. - 1. / 1.01);
//			double k_1 = h * f(T_n, Y_n);
//			double k_2 = h * f(T_n - .5 * h, Y_n + .5 * k_1);
//			double k_3 = h * f(T_n - .5 * h, Y_n + .5 * k_2);
//			double k_4 = h * f(T_n - h, Y_n + k_3);
//			Y_n = Y_n + 1. / 6. * (k_1 + 2. * k_2 + 2. * k_3 + k_4);
//			std::cout << std::scientific << T_n / cgs::GeV << " " << h / cgs::GeV << " "
//					<< Y_n / -dEdx->get(T_n) << "\n";
//		}
//	}
//	return 0;
//}
//
//struct gsl_f_params {
//	Particle * pt_Particle;
//};
//
//int func(double T, const double Y[], double f[], void *params) {
//	auto p = *(gsl_f_params *) params;
//	Particle * particle = p.pt_Particle;
//	T = std::fabs(T);
//	double Q = particle->get_Q(T);
//	double dEdx = std::fabs(particle->get_dEdx(T));
//	double X = particle->get_X(T);
//	double sigma_m = particle->get_sigma_m(T);
//	f[0] = Q - Y[0] / dEdx * (1. / X + sigma_m);
//	return GSL_SUCCESS;
//}
//
//int jac(double T, const double Y[], double *dfdY, double dfdT[], void *params) {
//	auto p = *(gsl_f_params *) params;
//	Particle * particle = p.pt_Particle;
//	T = std::fabs(T);
//	double dEdx = std::fabs(particle->get_dEdx(T));
//	double X = particle->get_X(T);
//	double sigma_m = particle->get_sigma_m(T);
//	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdY, 1, 1);
//	gsl_matrix * m = &dfdy_mat.matrix;
//	gsl_matrix_set(m, 0, 0, -1. / dEdx * (1. / X + sigma_m));
//	dfdT[0] = 0.0;
//	return GSL_SUCCESS;
//}
//
//int Particle::run_gsl(const LogAxis& T_) {
//	_T = T_.get_axis();
//	gsl_f_params params = { this };
//	gsl_odeiv2_system sys = { func, jac, 1, &params };
//	double Y[1] = { 0.0 };
//	_I_T.reserve(_T.size());
//	int status = 0;
//	_I_T.emplace_back(0);
//	for (size_t i = _T.size() - 1; i > 0; --i) {
//		double T_0 = -_T.at(i);
//		double T_1 = -_T.at(i - 1);
//		double h_start = std::fabs(T_0 - T_1);
//		//auto d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, h_start, 1e-7, 0.);
//		auto d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4imp, h_start, 1e-7, 0.);
//		status = gsl_odeiv2_driver_apply(d, &T_0, T_1, Y);
//		_I_T.emplace_back(Y[0] / -dEdx->get(-T_0));
//		gsl_odeiv2_driver_free(d);
//	}
//	std::reverse(std::begin(_I_T), std::end(_I_T));
//	assert(_T.size() == _I_T.size());
//	return status;
//}

//int Particle::run(const LogAxis& T) {
//	I_T.resize(T.get_size());
//	double Y_n = 0;
//	double r = T.get_ratio();
//	for (double T_n = cgs::TeV; T_n > cgs::GeV; T_n /= 1.01) {
//		//for (auto it = T.get_axis().rbegin(); it != T.get_axis().rend(); ++it) {
//		//	double h_big = *it * (1. - 1. / r);
//		//	double T_n = *it;
//		double h = T_n * (1. - 1. / 1.01);
//		double k_1 = h * f(T_n, Y_n);
//		double k_2 = h * f(T_n - .5 * h, Y_n + .5 * k_1);
//		double k_3 = h * f(T_n - .5 * h, Y_n + .5 * k_2);
//		double k_4 = h * f(T_n - h, Y_n + k_3);
//		Y_n = Y_n + 1. / 6. * (k_1 + 2. * k_2 + 2. * k_3 + k_4);
//		std::cout << T_n / cgs::GeV << " " << Y_n / -dEdx->get(T_n) << "\n";
//
//		//std::cout << std::scientific << *it / cgs::GeV << " " << (*it - h_big) / cgs::GeV << " "
//		//		<< Y_n / -dEdx->get(*it) << "\n";
//	}
//
//	exit(1);
//	return 0;
//}
