#include "CDF.h"

using namespace std;

Cdf::Cdf(double r, double v, double T, double s0, double alpha) : r(r), v(v), T(T), s0(s0), alpha(alpha){
}

// Cumulative Distribution Function for BS model
complex<double> Cdf::cdfBs(const std::complex<double> t) const {
    complex<double> i(0.0, 1.0);
    complex<double> phi = exp(i * t * (s0 + (r - 0.5 * v * v) * T) - 0.5 * v * v * t * t * T);
    return phi;
}


// Fourier tranform of c_T
complex<double> Cdf::psi(const std::complex<double> t) const {
    complex<double> i(0.0, 1.0);
    return exp(-r * T) * cdfBs(t - i * (alpha+1)) / (alpha * alpha + alpha - t * t + i * t * (2 * alpha + 1));
}