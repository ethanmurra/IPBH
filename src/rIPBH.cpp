#include <vector>
#include <stdexcept>
#include <random>
#include <cmath>
#include <boost/math/special_functions/lambert_w.hpp>
#include<Rcpp.h>

// [[Rcpp::export]]
std::vector<std::vector<double>> rIPBH(int n, const std::vector<double>& alpha, const std::vector<double>& eta) {
  if (n <= 0) {
    throw std::invalid_argument("N must be a positive integer");
  }
  if (alpha.size() != eta.size()) {
    throw std::invalid_argument("Alpha and eta vectors must have the same length.");
  }

  std::vector<std::vector<double>> results(alpha.size(), std::vector<double>(n));
  static std::mt19937 gen(std::random_device{}()); // Static to initialize once
  std::uniform_real_distribution<double> dis(0.0, 1.0);

  for (size_t i = 0; i < alpha.size(); ++i) {
    if (alpha[i] <= 0 || eta[i] <= 0) {
      throw std::invalid_argument("Alpha and eta must be greater than 0.");
    }
    long double alphaExpAlpha = alpha[i] * exp(alpha[i]); // Precompute outside the loop, with higher precision
    for (int j = 0; j < n; ++j) {
      double p = dis(gen); // Use double for the random number generation
      long double z = alphaExpAlpha / p;
      long double Wz = boost::math::lambert_w0(z);
      long double Q = powl((1 / alpha[i] * Wz - 1), -1 / eta[i]);
      results[i][j] = static_cast<double>(Q); // Cast back to double to store in results
    }
  }

  return results;
}
