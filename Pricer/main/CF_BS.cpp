#include "CF_BS.h"

using namespace std;

CF_BS::CF_BS(double r, double T, double alpha, int N, double eta, double s0, double v) : CF(r, T, alpha, N, eta, s0), v(v) {
}

// Characteristic Function for BS model
complex<double> CF_BS::phi(const std::complex<double> t) const {
    complex<double> i(0.0, 1.0);
    complex<double> phi = exp(i * t * (s0 + (r - 0.5 * v * v) * T) - 0.5 * v * v * t * t * T);
    return phi;
}