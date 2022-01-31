#define _USE_MATH_DEFINES

#ifndef CDF_h
#define CDF_h

#include <cmath>
#include <complex>


class Cdf {
    public:
        double r, v, T, s0, alpha;    
    
        Cdf(double r, double v, double T, double s0, double alpha);
        std::complex<double> cdfBs(const std::complex<double> t) const;
        std::complex<double> psi(const std::complex<double> t) const;
};

#endif