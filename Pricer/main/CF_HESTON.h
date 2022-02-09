#define _USE_MATH_DEFINES

#ifndef CF_HESTON_h
#define CF_HESTON_h

#include <cmath>
#include <complex>
#include "CF.h"


class CfHeston:public CF{
    protected:
        double v; //volatility
        double k; // Speed of mean reversion
        double theta; // avg level of volatility
        double sigma; // Volatility of volatility
        double p; // Correlation between the 2 Brownian components
        double v0; // Initial volatility
    
    public:
        virtual complex<double> phi(complex<double> u) const;
        CfHeston(double r, double T, double alpha, int N, double eta, double s0, double v, double k, double theta, double sigma, double p, double v0);
};

#endif