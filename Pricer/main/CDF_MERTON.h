#define _USE_MATH_DEFINES

#ifndef CDF_MERTON_h
#define CDF_MERTON_h

#include <cmath>
#include <complex>


class CdfMerton {
    public:
        double r, v, T, s0, alpha, lambda, mu, delta;
    
        CdfMerton(double r, double v, double T, double s0, double alpha, double lambda, double mu, double delta);
        std::complex<double> cdfBsMerton(const std::complex<double> z) const;
        std::complex<double> psi(const std::complex<double> t) const;
};

#endif

