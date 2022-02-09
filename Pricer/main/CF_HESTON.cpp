#include "CF_HESTON.h"

using namespace std;

CfHeston::CfHeston(double r, double T, double alpha, int N, double eta, double s0, double v, double k, double theta, double sigma, double p, double v0) : CF(r, T, alpha, N, eta, s0), v(v), k(k), theta(theta), sigma(sigma), p(p), v0(v0) {
}

// Characteristic Function for Heston model
complex<double> CfHeston::phi(const complex<double> z) const {
    complex<double> i(0.0, 1.0);
    complex<double> gamma = sqrt(sigma * sigma * (z * z + i * z) + pow(k - i * p * sigma * z, 2));
    complex<double> up = exp(k * theta * T * (k - i * p * sigma * z) / (sigma * sigma) + i * z * T * r);
    complex<double> down = pow(cosh(0.5 * gamma * T) + sinh(0.5 * gamma * T) * (k - i * p * sigma * z) / gamma, 2 * k * theta / (sigma * sigma));
    complex<double> sec_term = exp(-(z * z + i * z) * v0 / (gamma * cosh(0.5 * gamma * T) / sinh(0.5 * gamma * T) + k - i * p * sigma * z));
    return up * sec_term * exp(i * z * s0) / down;
}
