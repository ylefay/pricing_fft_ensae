#define _USE_MATH_DEFINES

#ifndef CF_h
#define CF_h

#include <cmath>
#include <complex>

class CF{
    public:
        CF(double r, double T, double alpha, int N, double eta, double s0);
        complex<double> psi(complex<double> u) const;
        virtual complex<double> phi(complex<double> u) const = 0;
	    double pricer(std::complex<double> *x_out, std::complex<double> *k_u) const;

    protected:
        double r, T, alpha, eta, s0;
        int N;
};

#endif