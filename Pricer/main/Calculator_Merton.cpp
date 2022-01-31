#include "CDF_MERTON.h"
#include "FFT.h"
#include "Calculator_Merton.h"

using namespace std;

CalculatorMerton::CalculatorMerton(double r, double v, double T, double s0, double alpha, double lambda, double mu, double delta, double eta, int N) : r(r), v(v), T(T), s0(s0), alpha(alpha), lambda(lambda), mu(mu), delta(delta), eta(eta), N(N){
}


double CalculatorMerton::pricer(complex<double> *x_out, complex<double> *k_u) const {
    double zeta = 2 * M_PI / (N * eta);
    complex<double> i(0.0, 1.0);
    CdfMerton cdf(r, v, T, s0, alpha, lambda, mu, delta);
    

    // Create the array of v_j values
    complex<double> v[N];
    for (int j = 0; j < N; j++){
        v[j] = eta * j;
    }

    // Calculate the array of x_j values
    complex<double> x_in[N];
    for (int j = 0; j < N; j++){
    x_in[j] = exp(i * (0.5 * N * zeta - s0) * v[j]) * cdf.psi(v[j]) ;
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