#include "CDF_HESTON.h"
#include "FFT.h"
#include "Calculator_Heston.h"

using namespace std;


CalculatorHeston::CalculatorHeston(double r, double v, double T, double s0, double alpha, double k, double theta, double sigma, double p, double v0, double eta, int N) : r(r), v(v), T(T), s0(s0), alpha(alpha), k(k), theta(theta), sigma(sigma), p(p), v0(v0), eta(eta), N(N) {
}


double CalculatorHeston::pricer(complex<double> *x_out, complex<double> *k_u) const {
    double zeta = 2 * M_PI / (N * eta);
    CdfHeston cdf(r, v, T, s0, alpha, k ,theta, sigma, p, v0);
    complex<double> i(0.0, 1.0);


    // Create the array of v_j values
    complex<double> v[N];
    for (int j = 0; j < N; j++){
        v[j] = eta * j;
    }

    // Calculate the array of x_j values
    complex<double> x_in[N];
    for (int j = 0; j < N; j++){
		x_in[j] = exp(i * (0.5 * N * zeta - s0) * v[j]) * cdf.psi(v[j]);
    }

    // Calculate the array of k_u values
    for (int j = 0; j < N; j++){
        k_u[j] = -0.5 * N * zeta + zeta * j + s0;
    }
    
    //Apply the FFT
    fft(x_in, x_out, N);

    for (int j = 0; j < N; j++){
        x_out[j] *= eta * exp(-alpha * k_u[j]) / M_PI;
    }

    return 0;
}