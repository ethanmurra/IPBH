#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <boost/math/special_functions/lambert_w.hpp>  // Include Boost for Lambert W
#include<Rcpp.h>

// [[Rcpp::export]]
std::vector<double> qIPBH(const std::vector<double>& p, const std::vector<double>& alpha, const std::vector<double>& eta, bool lower_tail = true, bool log_p = false) {
  // Validate inputs: Check if all vectors are of the same length
  if (p.size() != alpha.size() || alpha.size() != eta.size()) {
    throw std::invalid_argument("All input vectors must have the same length.");
  }

  std::vector<double> results;
  results.reserve(p.size());  // Reserve space for results to optimize memory allocation

  for (size_t i = 0; i < p.size(); ++i) {
    // Validate each set of parameters
    if (alpha[i] <= 0 || eta[i] <= 0) {
      throw std::invalid_argument("alpha and eta must be greater than 0 for all elements.");
    }
    if (p[i] <= 0 || p[i] >= 1) {
      throw std::invalid_argument("p must be strictly between 0 and 1 for all elements.");
    }

    // Adjust p if necessary
    double adjusted_p = p[i];
    if (log_p) {
      adjusted_p = exp(p[i]);
    }
    if (!lower_tail) {
      adjusted_p = 1 - adjusted_p;
    }

    // Calculation
    double z = alpha[i] * exp(alpha[i]) / adjusted_p;
    double Wz = boost::math::lambert_w0(z);  // Using Boost for Lambert W function
    double Q = pow((1 / alpha[i] * Wz - 1), -1 / eta[i]);

    results.push_back(Q);
  }

  return results;
}
