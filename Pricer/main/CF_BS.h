#define _USE_MATH_DEFINES

#ifndef CF_BS_h
#define CF_BS_h

#include <cmath>
#include <complex>
#include "CF.h"

class CF_BS:public CF{
    protected:
        double v; //volatility    
    
    public:
        virtual complex<double> phi(complex<double> u) const;
        CF_BS(double r, double T, double alpha, int N, double eta, double s0, double v);
       
};

#endif