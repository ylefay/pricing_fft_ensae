#include "CDF_MERTON.h"

using namespace std;

CdfMerton::CdfMerton(double r, double v, double T, double s0, double alpha, double lambda, double mu, double delta) : r(r), v(v), T(T), s0(s0), alpha(alpha), lambda(lambda), mu(mu), delta(delta) {
}

// Cumulative Distribution Function for BS Merton model
complex<double> CdfMerton::cdfBsMerton(const complex<double> z) const {
    complex<double> i(0.0, 1.0);
    complex<double> mu_M = r - v * v - lambda * (exp(mu + 0.5 * delta * delta) - 1);
    complex<double> phi = exp(T * (-0.5 * pow(v * z, 2) + i * mu_M * z + lambda * exp(-pow(delta * z, 2) / (complex<double>(2, 0)) + i * mu * z - complex<double> (1.0, 0.0))));
    return phi * exp(i * z * s0);
}

// Fourier tranform of c_T
complex<double> CdfMerton::psi(const complex<double> t) const {
    complex<double> i(0.0, 1.0);
    return (exp(-r * T) * cdfBsMerton(t - i * (alpha + 1))) / (alpha * alpha + alpha - t * t + i * t * (2 * alpha + 1));
}