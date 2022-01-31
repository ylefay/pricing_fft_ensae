#include "FFT.h"

using namespace std;

void fft_recurrence(complex<double> *x, int N) {
	if (N <= 1) {
		return;
	}

	complex<double> odd[N / 2];
	complex<double> even[N / 2];
	
	for (int i = 0; i < N / 2; i++) {
		even[i] = x[i * 2];
		odd[i] = x[i * 2 + 1];
	}

	fft_recurrence(even, N / 2);
	fft_recurrence(odd, N / 2);


	// Apply the DFT to even and odd
	for (int k = 0; k < N / 2; k++) {
		complex<double> t = exp(complex<double>(0, -2 * M_PI * k / N)) * odd[k];
		x[k] = even[k] + t;
		x[N / 2 + k] = even[k] - t;
	}
}

void fft(complex<double> *x_in, complex<double> *x_out, int N) {

	for (int i = 0; i < N; i++) {
		x_out[i] = x_in[i];
	}

	fft_recurrence(x_out, N);
}