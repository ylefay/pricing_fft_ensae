#include "CF.h"

using namespace std;

CF::CF(double r, double T, double alpha, int N, double eta, double s0):r(r), T(T), alpha(alpha), N(N), eta(eta), s0(s0){}

// Fourier tranform of c_T
complex<double> CF::psi(const complex<double> t) const {
    complex<double> i(0.0, 1.0);
    return exp(-r * T) * phi(t - i * (alpha + 1)) / (alpha * alpha + alpha - t * t + i * t * (2 * alpha + 1));
}

double CF::pricer(complex<double> *x_out, complex<double> *k_u) const{
    double zeta = 2 * M_PI / (N * eta);
    complex<double> i(0.0, 1.0);

    // Calculate the array of v_j values
    complex<double> v[N];
    for (int j = 0; j < N; j++){
        v[j] = eta * j;
    }

    // Calculate the array of x_j values
    complex<double> x_in[N];
    for (int j = 0; j < N; j++){
    x_in[j] = exp(i * (0.5 * N * zeta - s0) * v[j]) * psi(v[j]) ;
    }

    // Calculate the array of k_u values
    for (int j = 0; j < N; j++){
        k_u[j] = (-0.5 * N * zeta + zeta * j) + s0;
    }

    //Apply the FFT
    fft(x_in, x_out, N);

    for (int j = 0; j < N; j++){
        x_out[j] *= eta * exp(-alpha * k_u[j]) / M_PI;
    }

    return 0;
}