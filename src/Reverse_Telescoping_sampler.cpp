// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <random>
#include <limits>
using namespace arma;
using namespace Rcpp;

// ===== Gamma sampler using std::gamma_distribution with local RNG ===== //
inline double rgamma_std(double shape, double scale, std::mt19937& local_gen) {
  std::gamma_distribution<double> gamma_dist(shape, scale);
  return gamma_dist(local_gen);
}

// Inline helper: recover Omega_tilde from Omega
inline arma::mat recover_Omegatilde_from_Omega(const arma::mat& Omega) {
  int p = Omega.n_rows;
  arma::mat Omega_tilde = Omega;
  for (int j = p - 1; j >= 1; --j) {
    arma::vec wtilde = Omega_tilde.submat(0, j, j - 1, j);
    double wjj = Omega_tilde(j, j);
    double denom = std::max(wjj, 1e-8);
    arma::mat Gammajj = wtilde * wtilde.t() / denom;
    Omega_tilde.submat(0, 0, j - 1, j - 1) -= Gammajj;
  }
  return Omega_tilde;
}


// ===== Exact corrected RT Gibbs helpers ===== //
// The corrected sampler retains the full RT state.  When updating column j,
// theta_minus stores the contribution of every RT block except j to the
// current upper-left j x j block of the original precision matrix.
inline double runif01(std::mt19937& local_gen) {
  static thread_local std::uniform_real_distribution<double> unif(0.0, 1.0);
  double u = unif(local_gen);
  // Guard only against an exact endpoint from an implementation-specific RNG.
  if (u <= 0.0) u = std::numeric_limits<double>::min();
  if (u >= 1.0) u = std::nextafter(1.0, 0.0);
  return u;
}

inline double log_pivot_density_z(double z,
                                  int n,
                                  double sjj,
                                  double chi,
                                  const arma::vec& omega_col,
                                  const arma::mat& theta_minus,
                                  const arma::mat& Lambda2_mat,
                                  double tau2) {
  if (!std::isfinite(z) || z > 700.0 || z < -700.0) {
    return -std::numeric_limits<double>::infinity();
  }
  const double d = std::exp(z);
  if (!(d > 0.0) || !std::isfinite(d)) {
    return -std::numeric_limits<double>::infinity();
  }

  // Jacobian for d = exp(z) changes d^{n/2} to exp{(n/2+1) z}.
  double out = (n / 2.0 + 1.0) * z - 0.5 * sjj * d - 0.5 * chi / d;

  const int j = static_cast<int>(omega_col.n_elem);
  for (int r = 0; r < j; ++r) {
    for (int k = r + 1; k < j; ++k) {
      const double lam2 = Lambda2_mat(r, k);
      if (!(lam2 > 0.0) || !(tau2 > 0.0)) {
        Rcpp::stop("Non-positive local/global scale encountered in RT pivot update.");
      }
      const double theta_rk = theta_minus(r, k) + omega_col(r) * omega_col(k) / d;
      out -= 0.5 * theta_rk * theta_rk / (tau2 * lam2);
    }
  }
  return out;
}

// Univariate slice update on z = log(tilde_theta_jj).  The density evaluation
// uses the unexpanded sum-of-squares correction for numerical stability.
inline double sample_log_pivot_slice(double current_d,
                                     int n,
                                     double sjj,
                                     double chi,
                                     const arma::vec& omega_col,
                                     const arma::mat& theta_minus,
                                     const arma::mat& Lambda2_mat,
                                     double tau2,
                                     std::mt19937& local_gen,
                                     double width = 1.0) {
  if (!(current_d > 0.0) || !std::isfinite(current_d)) {
    Rcpp::stop("Current RT pivot must be finite and positive.");
  }

  const double z0 = std::log(current_d);
  const double f0 = log_pivot_density_z(z0, n, sjj, chi, omega_col,
                                        theta_minus, Lambda2_mat, tau2);
  if (!std::isfinite(f0)) {
    Rcpp::stop("Non-finite log density at current RT pivot.");
  }
  const double logy = f0 + std::log(runif01(local_gen));

  double left = z0 - width * runif01(local_gen);
  double right = left + width;

  int guard = 0;
  while (log_pivot_density_z(left, n, sjj, chi, omega_col,
                             theta_minus, Lambda2_mat, tau2) > logy) {
    left -= width;
    if (++guard > 10000) Rcpp::stop("Slice stepping-out failed on left side.");
  }
  guard = 0;
  while (log_pivot_density_z(right, n, sjj, chi, omega_col,
                             theta_minus, Lambda2_mat, tau2) > logy) {
    right += width;
    if (++guard > 10000) Rcpp::stop("Slice stepping-out failed on right side.");
  }

  guard = 0;
  while (true) {
    const double z = left + (right - left) * runif01(local_gen);
    const double fz = log_pivot_density_z(z, n, sjj, chi, omega_col,
                                          theta_minus, Lambda2_mat, tau2);
    if (fz >= logy) return std::exp(z);
    if (z < z0) left = z; else right = z;
    if (++guard > 100000) Rcpp::stop("Slice shrinkage failed for RT pivot.");
  }
}

// Inline helper: safe matrix inversion
// this is used for generate Omega_init = safe_inverse(S /n + delta I_n) //S=YY^T
inline arma::mat safe_inverse(const arma::mat& A) {
  arma::mat A_inv;
  bool success = false;

  try {
    A_inv = arma::inv_sympd(A);
    success = true;
  } catch (...) {
    // inv_sympd failed
  }

  if (!success) {
    try {
      A_inv = arma::inv(A);
      success = true;
    } catch (...) {
      // inv failed
    }
  }

  if (!success) {
    Rcpp::Rcerr << "Matrix inversion failed: possibly singular or ill-conditioned." << std::endl;
    return arma::mat(A.n_rows, A.n_cols, arma::fill::none);
  }

  return A_inv;
}

/**
 * @title Graphical Horseshoe MCMC Sampler (Tilde fast Sampler)
 *
 * Key Features:
 *   - Exact scalar Gaussian Gibbs updates for RT off-diagonal coordinates
 *   - Exact log-scale slice update for each positive RT pivot
 *   - Local (lambda) and global (tau) shrinkage updated using standard and
 *     inverse gamma distributions
 *   - Parallel-safe random number generation with std::mt19937 and arma::arma_rng
 *
 * @param p         Number of variables (dimensionality of precision matrix)
 * @param Y         n x p data matrix
 * @param M         Total number of MCMC iterations
 * @param burnin    Number of burn-in iterations to discard
 * @param seed      RNG seed for reproducibility (sets both C++ and Armadillo RNGs)
 *
 * @return A list containing:
 *   - Omega_save      : p x p x (M - burnin) posterior samples of precision matrix
 *   - OmegaTilde_save : intermediate conditional precision components
 *   - Lambda2_save    : local shrinkage lambda^2 values
 *   - tau2_save       : global shrinkage parameter tau^2
 *
 * This implementation is designed for integration with R via Rcpp and
 * supports safe parallel execution using per-call RNG seeding.
 */



////////////////main function: GHS posterior tilde sampler ////////////////////////////////////////////////

// [[Rcpp::export]]
List RT_sampler_cpp(const arma::mat& Y, int M, int burnin, int seed, const std::string& prior) {
  // validate prior
  if (prior != "GHS" && prior != "GHSL") {
    Rcpp::stop("`prior` must be either \"GHS\" or \"GHSL\".");
  }

  // dim of Y
  int n = Y.n_rows;
  int p = Y.n_cols;

  // Scale Y and set up RNG
  double k = p;
  arma::mat Y_scaled = Y / std::sqrt(k);

  std::mt19937 gen(seed);
  arma::arma_rng::set_seed(seed);
  int save_dim = M - burnin;

  // ==== initialize Omega using ridge-like inverse ====
  arma::mat S = Y_scaled.t() * Y_scaled;
  double delta = 0.01;
  arma::mat Sigma_hat = S / n + delta * arma::eye<arma::mat>(p, p);
  arma::mat Omega = safe_inverse(Sigma_hat);

  // ==== initialize Omega_tilde via helper ====
  arma::mat Omega_tilde = recover_Omegatilde_from_Omega(Omega);

  // Initialize other parameters
  arma::mat Lambda2_mat(p, p, fill::ones);
  arma::mat Nu_mat(p, p, fill::ones);
  double tau2 = 1.0;
  double xi = 1.0;

  // Storage
  arma::cube Omega_save(p, p, save_dim, fill::zeros);
  arma::cube Lambda2_save(p, p, save_dim, fill::zeros);
  arma::vec tau2_save(save_dim, fill::zeros);

  // ==== MCMC loop ====
  for (int m = 0; m < M; ++m) {
    // if (m % 100 == 0) {
    //   Rcpp::Rcout << "m = " << m << std::endl;
    // }
    Rcpp::Rcout << "\rIteration m = " << m << std::flush;

    arma::mat C(p, p, fill::zeros);

    // Deterministic-scan Gibbs over RT columns, from p down to 1.
    // C contains the rank-one contributions from already-updated future columns.
    for (int j = p-1; j >= 0; --j) {
      if (j == 0) {
        // No earlier off-diagonal prior factors exist for the first pivot.
        double shape = n/2.0 + 1.0;
        double rate = S(0,0) / 2.0;
        Omega_tilde(0,0) = rgamma_std(shape, 1.0/rate, gen);
        Omega(0,0) = Omega_tilde(0,0) + C(0,0);
      } else {
        // Current (old) RT column and its fixed future shift.
        arma::vec omega_sub = Omega_tilde.submat(0, j, j-1, j);
        const double old_d = Omega_tilde(j,j);
        arma::vec gamma_sub = C.submat(0, j, j-1, j);
        const double gamma_jj = C(j,j);

        if (!(old_d > 0.0)) {
          Rcpp::stop("Encountered a non-positive RT pivot before column update.");
        }

        // Remove the old rank-one contribution of RT block j.  Conditional on
        // all other RT blocks this matrix is fixed throughout the column update.
        arma::mat old_outer = omega_sub * omega_sub.t();
        arma::mat theta_minus = Omega.submat(0, 0, j-1, j-1) - old_outer / old_d;
        theta_minus = 0.5 * (theta_minus + theta_minus.t());

        // ---- exact diagonal full conditional: slice on log(tilde_theta_jj) ----
        const arma::mat Ssub = S.submat(0, 0, j-1, j-1);
        const double chi = as_scalar(omega_sub.t() * Ssub * omega_sub);
        const double sjj = S(j,j);
        const double new_d = sample_log_pivot_slice(old_d, n, sjj, chi,
                                                     omega_sub, theta_minus,
                                                     Lambda2_mat, tau2, gen);
        Omega_tilde(j,j) = new_d;

        // ---- exact scalar Gaussian full conditionals for tilde_theta_{rj} ----
        // Use a Gauss-Seidel sweep: omega_sub always contains the most recent
        // value of every coordinate already visited in this column.
        for (int r = 0; r < j; ++r) {
          const double w_rj = 1.0 / (tau2 * Lambda2_mat(r,j));
          double precision = S(r,r) / new_d + w_rj;
          double natural = -S(r,j) - w_rj * gamma_sub(r);

          for (int k2 = 0; k2 < j; ++k2) {
            if (k2 == r) continue;
            const double w_rk = 1.0 / (tau2 * Lambda2_mat(r,k2));
            precision += w_rk * omega_sub(k2) * omega_sub(k2) / (new_d * new_d);
            natural -= S(r,k2) * omega_sub(k2) / new_d;
            natural -= w_rk * theta_minus(r,k2) * omega_sub(k2) / new_d;
          }

          if (!(precision > 0.0) || !std::isfinite(precision) || !std::isfinite(natural)) {
            Rcpp::stop("Invalid scalar Gaussian full conditional in corrected RT sampler.");
          }
          const double mean = natural / precision;
          const double sd = 1.0 / std::sqrt(precision);
          std::normal_distribution<double> normal_dist(mean, sd);
          omega_sub(r) = normal_dist(gen);
        }

        Omega_tilde.submat(0, j, j-1, j) = omega_sub;
        Omega_tilde.submat(j, 0, j, j-1) = omega_sub.t();

        // Replace the old contribution by the newly sampled rank-one term,
        // keeping the original precision matrix synchronized with the RT state.
        arma::mat new_outer = omega_sub * omega_sub.t();
        Omega.submat(0, 0, j-1, j-1) = theta_minus + new_outer / new_d;
        Omega.submat(0, j, j-1, j) = omega_sub + gamma_sub;
        Omega.submat(j, 0, j, j-1) = Omega.submat(0, j, j-1, j).t();
        Omega(j,j) = new_d + gamma_jj;

        // The updated block now contributes to all earlier future shifts.
        C.submat(0, 0, j-1, j-1) += new_outer / new_d;

        // Prior-specific local-scale updates are unchanged and use the
        // reconstructed original precision entries theta_{ij} = Omega(i,j).
        if (prior == "GHS") {
          for (int i = 0; i < j; ++i) {
            double rate_lambda = 1.0/Nu_mat(i,j) + std::pow(Omega(i,j),2)/(2.0*tau2);
            Lambda2_mat(i,j) = 1.0/rgamma_std(1.0, 1.0/rate_lambda, gen);
            Lambda2_mat(j,i) = Lambda2_mat(i,j);

            double rate_nu = 1.0 + 1.0/Lambda2_mat(i,j);
            Nu_mat(i,j) = 1.0/rgamma_std(1.0, 1.0/rate_nu, gen);
            Nu_mat(j,i) = Nu_mat(i,j);
          }
        } else { /* GHSL */
          for (int i = 0; i < j; ++i) {
            double rate_lambda = Nu_mat(i,j)/2.0 + std::pow(Omega(i,j),2)/(2.0*tau2);
            double Lambda2_sample = 1.0/rgamma_std(1.0, 1.0/rate_lambda, gen);
            Lambda2_mat(i,j) = Lambda2_sample;
            Lambda2_mat(j,i) = Lambda2_mat(i,j);

            double rand_val = arma::randu<double>();
            double Nu_sample = -2.0 * Lambda2_sample * log(1.0 - rand_val * (1.0 - exp(-0.5 / Lambda2_sample)));
            Nu_mat(i,j) = Nu_sample;
            Nu_mat(j,i) = Nu_mat(i,j);
          }
        }
      }
    }

    // update tau & xi
    arma::uvec lower_indices = find(trimatl(ones<arma::mat>(p,p), -1));
    arma::vec omega_vec = Omega.elem(lower_indices);
    arma::vec lambda_sq_vec = Lambda2_mat.elem(lower_indices);

    double shape_tau = (p*(p-1)/2.0 + 1)/2.0;
    double rate_tau = 1.0/xi + accu(pow(omega_vec,2)/(2.0*lambda_sq_vec));
    tau2 = 1.0 / rgamma_std(shape_tau, 1.0 / rate_tau, gen);
    xi = 1.0 / rgamma_std(1.0, 1.0 / (1.0 + 1.0 / tau2), gen);

    if (m >= burnin) {
      int idx = m - burnin;
      Omega_save.slice(idx) = Omega / k;
      Lambda2_save.slice(idx) = Lambda2_mat;
      tau2_save(idx) = tau2;
    }
  }
  Rcpp::Rcout << std::endl;

  return List::create(
    Named("Omega_save") = Omega_save,
    Named("Lambda2_save") = Lambda2_save,
    Named("tau2_save") = tau2_save
  );
}
