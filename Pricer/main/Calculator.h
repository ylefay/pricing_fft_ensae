#ifndef CALCULATOR_h
#define CALCULATOR_h

#include <cmath>
#include <complex>

class Calculator{
    public:
		double r, v, T, s0, alpha, eta;
		int N;
		Calculator(double r, double v, double T, double s0, double alpha, double eta, int N);
		double pricer(std::complex<double> *x_out, std::complex<double> *k_u) const;
};

#endif