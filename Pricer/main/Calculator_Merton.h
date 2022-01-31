#ifndef CALCULATOR_MERTON_h
#define CALCULATOR_MERTON_h

#include <cmath>
#include <complex>

class CalculatorMerton{
    public:
		double r, v, T, s0, alpha, lambda, mu, delta, eta;
		int N;
		CalculatorMerton(double r, double v, double T, double s0, double alpha, double lambda, double mu, double delta, double eta, int N);
		double pricer(std::complex<double> *x_out, std::complex<double> *k_u) const;
};

#endif