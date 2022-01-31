#define _USE_MATH_DEFINES

#ifndef ClosedForm_h
#define ClosedForm_h

#include <cmath>
#include <complex>


class ClosedForm {
    public:
        double r, v, T, S;
    
        ClosedForm(double r, double v, double T, double S);
        double cdf_norm(double x);
        double d_j(int j, double K);
        double call_price(double K);
};

#endif