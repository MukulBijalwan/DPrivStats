// Fast DP primitives for DPrivStats.
//
// - Inverse-CDF Laplace sampler: u ~ U(-0.5, 0.5),
//   x = -scale * sign(u) * log(1 - 2|u|).
// - Row-wise gradient clipping: g_i * min(1, C / ||g_i||_2) for DP-SGD.

#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rlaplace_cpp(int n, double scale) {
  NumericVector u = runif(n, -0.5, 0.5);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double ui = u[i];
    out[i] = -scale * ((ui > 0) - (ui < 0)) *
             std::log(1.0 - 2.0 * std::fabs(ui));
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix clip_gradients_cpp(NumericMatrix grads, double max_grad_norm) {
  int nr = grads.nrow(), nc = grads.ncol();
  NumericMatrix out(nr, nc);
  for (int i = 0; i < nr; ++i) {
    double norm2 = 0.0;
    for (int j = 0; j < nc; ++j) norm2 += grads(i, j) * grads(i, j);
    double norm = std::sqrt(norm2);
    double factor = (norm > max_grad_norm && norm > 0.0)
                      ? max_grad_norm / norm : 1.0;
    for (int j = 0; j < nc; ++j) out(i, j) = grads(i, j) * factor;
  }
  return out;
}
