#define _USE_MATH_DEFINES

#ifndef CF_MERTON_h
#define CF_MERTON_h

#include <cmath>
#include <complex>
#include "CF.h"

class CfMerton:public CF{
    protected:
        double v, lambda, mu, delta;
    
    public:
        virtual complex<double> phi(complex<double> u) const;
        CfMerton(double r, double T, double alpha, int N, double eta, double s0, double v, double lambda, double mu, double delta);
};

#endif

