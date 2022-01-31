#ifndef FFT_h
#define FFT_h

#include <cmath>
#include <complex>

void fft_recurrence(std::complex<double> *x, int N);
void fft(std::complex<double> *x_in, std::complex<double> *x_out, int N);

#endif