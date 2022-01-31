#include <iostream>
#include <cmath>
#include "closed_form.h"

using namespace std;

ClosedForm::ClosedForm(double r, double v, double T, double S) : r(r), v(v), T(T), S(S) {}


double ClosedForm::cdf_norm(double x) {
    double k = 1.0 / (0.2316419 * x + 1.0);
    double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

    if (x >= 0.0) {
        return 1.0 - 1.0 / pow(2 * M_PI, 0.5) * exp(-0.5 * x * x) * k_sum;
    } else {
        return 1.0 - cdf_norm(-x);
    }
}

double ClosedForm::d_j(int j, double K) {
    return (log(S / K) + (r + pow(-1, j - 1) * 0.5 * v * v) * T) / (v * sqrt(T));
}

double ClosedForm::call_price(double K) {
    return S * cdf_norm(d_j(1, K)) - K * exp(-r * T) * cdf_norm(d_j(2, K));
}