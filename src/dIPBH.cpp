#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>
#include<Rcpp.h>

// [[Rcpp::export]]
std::vector<double> dIPBH(const std::vector<double>& x, const std::vector<double>& alpha, const std::vector<double>& eta, bool log_transform = false) {
  // Validate inputs: Check if all vectors are of the same length
  if (x.size() != alpha.size() || alpha.size() != eta.size()) {
    throw std::invalid_argument("All input vectors must have the same length.");
  }

  std::vector<double> results;
  results.reserve(x.size()); // Reserve space for results to optimize memory allocation

  for (size_t i = 0; i < x.size(); ++i) {
    // Validate each set of parameters
    if (x[i] <= 0 || alpha[i] <= 0 || eta[i] <= 0) {
      throw std::invalid_argument("x, alpha, and eta must be greater than 0 for all elements.");
    }

    // Calculate density for each set
    double log_f = log(eta[i]) - alpha[i] * pow(x[i], -eta[i]) + log(alpha[i] + (alpha[i] + 1) * pow(x[i], eta[i])) - log(x[i]) - 2 * log(pow(x[i], eta[i]) + 1);

    if (log_transform) {
      results.push_back(log_f);
    } else {
      results.push_back(exp(log_f));
    }
  }

  return results;
}
