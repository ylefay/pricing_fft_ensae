#define _USE_MATH_DEFINES

#ifndef CDF_HESTON_h
#define CDF_HESTON_h

#include <cmath>
#include <complex>


class CdfHeston {
    public:
        double r, v, T, s0, alpha;
        double k; // Speed of mean reversion
        double theta; // avg level of volatility
        double sigma; // Volatility of volatility
        double p; // Correlation between the 2 Brownian components
        double v0; // Initial volatility
    
        CdfHeston(double r, double v, double T, double s0, double alpha, double k, double theta, double sigma, double p, double v0);
        std::complex<double> cdfHeston(const std::complex<double> z) const;
        std::complex<double> psi(const std::complex<double> t) const;
};

#endif