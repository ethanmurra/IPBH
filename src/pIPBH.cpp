#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>
#include<Rcpp.h>

// [[Rcpp::export]]
std::vector<double> pIPBH(const std::vector<double>& x, const std::vector<double>& alpha, const std::vector<double>& eta, bool lower_tail = true, bool log_p = false) {
  size_t n = x.size(); // Use size once and store it for repeated use

  // Validate inputs: Check if all vectors are of the same length
  if (n != alpha.size() || n != eta.size()) {
    throw std::invalid_argument("All input vectors must have the same length.");
  }

  std::vector<double> results;
  results.reserve(n); // Pre-allocate memory for results to optimize memory allocation

  for (size_t i = 0; i < n; ++i) {
    // Validate each set of parameters
    if (x[i] <= 0 || alpha[i] <= 0 || eta[i] <= 0) {
      throw std::invalid_argument("x, alpha, and eta must be greater than 0 for all elements.");
    }

    // Calculate the CDF more efficiently
    double x_neg_eta = pow(x[i], -eta[i]);
    double log_F = -alpha[i] * x_neg_eta - log(x_neg_eta + 1); // Avoid recalculating pow(x[i], -eta[i])

    // Adjust for tail
    if (!lower_tail) {
      log_F = log1p(-exp(log_F)); // Use log1p for better numerical stability when computing log(1 + value)
    }

    if (log_p) {
      results.push_back(log_F);
    } else {
      results.push_back(exp(log_F)); // Use direct exp of log_F for non-log probabilities
    }
  }

  return results;
}
