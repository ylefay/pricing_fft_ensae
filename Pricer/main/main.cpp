#include <iostream>
#include<cstdio>
#include<cmath>
#include "FFT.cpp"
#include "CF.cpp"
// BS Model
#include "closed_form.cpp"
#include "CF_BS.cpp"
// BS Merton Model
#include "CF_MERTON.cpp"
// Heston Model
#include "CF_HESTON.cpp"


using namespace std;

void extract_results_fft(const char *file_name, complex<double> *x_out, complex<double> *k_u, int N){
	FILE * f;
	f = fopen (file_name, "wt");

	for (int i = 0; i < N; i++){
		fprintf (f, "{ %f , %f } , ", exp(real(k_u[i])), real(x_out[i]));
	}
	
	fclose(f);
}

void extract_results_closed_form(const char *file_name, double *call, double *K_u, int N){
	FILE * f;
	f = fopen (file_name, "wt");

	for (int i = 0; i < N; i++){
		fprintf (f, "{ %f , %f } , ", K_u[i], call[i]);
	}
	
	fclose(f);
}

int main(){
	// Parameters used in BS, Heston and BSM
	int N = 8192;
	double r = 0.05;
	double v = 0.2;
	double T = 1;
	double S = 100;
	double s0 = log(S);
	double alpha = 0.75;
	double eta = 0.0125;
	// Parameters for Heston
	double k = 10;
    double theta = 0.2;
    double sigma = 0.7;
    double p = -0.5;
    double v0 = 0.2;
	// Parameters for BSMerton
	double lambda = 0.13;
	double mu = 0;
	double delta = 0.0004;

	ClosedForm pricer_cf(r, v, T, S);
	CF_BS cf_bs(r, T, alpha, N, eta, s0, v);
	//CF_MERTON cf_merton(r, T, alpha, N, eta, s0, v, lambda, mu, delta);
	//CF_HESTON cf_heston(r, T, alpha, N, eta, s0, v, k, theta, sigma, p, v0);

	complex<double> x_out[N];
	complex<double> k_u[N]; // Strikes for the FFT

	double call[N];
	double K_u[N]; //Strikes for the closed Form


	double zeta = 2 * M_PI / (N * eta);
    for (int j = 0; j < N; j++){
        K_u[j] = exp(-0.5 * N * zeta + zeta * j + s0);
    }
	for (int j=0; j<N; j++){
		call[j] = pricer_cf.call_price(K_u[j]);
	}
	cout << "Closed Form : " << call[N/2] << endl;
	extract_results_closed_form("values_closed_form.txt", call, K_u, N);

	cf_bs.pricer(x_out, k_u);
	cout << "BS : " << real(x_out[N / 2]) << endl;
	extract_results_fft("values_BS.txt", x_out, k_u, N);

	//cf_merton.pricer(x_out, k_u);
	//cout << "BS MERTON : " << real(x_out[N / 2]) << endl;
	//extract_results_fft("values_BSM.txt", x_out, k_u, N);

	//cf_heston.pricer(x_out, k_u);
	//cout << "HESTON : " << real(x_out[N / 2]) << endl;
	//extract_results_fft("values_H.txt", x_out, k_u, N);
}