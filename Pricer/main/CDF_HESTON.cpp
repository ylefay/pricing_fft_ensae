#include "CDF_HESTON.h"

using namespace std;

CdfHeston::CdfHeston(double r, double v, double T, double s0, double alpha, double k, double theta, double sigma, double p, double v0) : r(r), v(v), T(T), s0(s0), alpha(alpha), k(k), theta(theta), sigma(sigma), p(p), v0(v0) {
}

// Cumulative Distribution Function for Heston model
complex<double> CdfHeston::cdfHeston(const complex<double> z) const {
    complex<double> i(0.0, 1.0);
    complex<double> gamma = sqrt(sigma * sigma * (z * z + i * z) + pow(k - i * p * sigma * z, 2));
    complex<double> up = exp(k * theta * T * (k - i * p * sigma * z) / (sigma * sigma) + i * z * T * r);
    complex<double> down = pow(cosh(0.5 * gamma * T) + sinh(0.5 * gamma * T) * (k - i * p * sigma * z) / gamma, 2 * k * theta / (sigma * sigma));
    complex<double> sec_term = exp(-(z * z + i * z) * v0 / (gamma * cosh(0.5 * gamma * T) / sinh(0.5 * gamma * T) + k - i * p * sigma * z));
    return up * sec_term * exp(i * z * s0) / down;
}

// Fourier tranform of c_T
complex<double> CdfHeston::psi(const complex<double> t) const {
    complex<double> i(0.0, 1.0);
    return exp(-r * T) * cdfHeston(t - i * (alpha + 1)) / (alpha * alpha + alpha - t * t + i * t * (2 * alpha + 1));
}