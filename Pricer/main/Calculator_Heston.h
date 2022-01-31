#ifndef CALCULATOR_HESTON_h
#define CALCULATOR_HESTON_h

#include <cmath>
#include <complex>

class CalculatorHeston{
    public:
		double r, v, T, s0, alpha, k, theta, sigma, p, v0, eta;
		int N;
		CalculatorHeston(double r, double v, double T, double s0, double alpha, double k, double theta, double sigma, double p, double v0, double eta, int N);
		double pricer(std::complex<double> *x_out, std::complex<double> *k_u) const;
};

#endif