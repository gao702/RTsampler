// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <random>
using namespace arma;
using namespace Rcpp;

// ===== Gamma sampler using std::gamma_distribution with local RNG ===== //
inline double rgamma_std(double shape, double scale, std::mt19937& local_gen) {
  std::gamma_distribution<double> gamma_dist(shape, scale);
  return gamma_dist(local_gen);
}

// ===== Bhattacharya Sampler ===== //
arma::vec Bhattacharya_cpp(int n, int p,
                           const arma::mat& X,
                           const arma::vec& B,
                           double sig,
                           const arma::vec& y) {
  // make sure the input dim is correct
  if (X.n_rows != n || X.n_cols != p) {
    stop("X dimensions must match n and p");
  }

  // --- Step 1: build Phi_mat ---
  arma::mat Phi_mat = X / sig;

  // --- Step 2: build D mat ---
  arma::mat D;
  if (p == 1) {
    D = B(0) * pow(sig, 2) * eye<mat>(1, 1);  // vector to a 1 by 1 matrix
  } else {
    D = diagmat(B) * pow(sig, 2);  // diagnal matrix
  }

  // --- Step 3: generate random number ---
  arma::vec u = randn<vec>(p);
  u %= sqrt(D.diag());  // %= means time element by element

  arma::vec delta = randn<vec>(n);
  arma::vec v = Phi_mat * u + delta;

  // --- Step 4: core matrix computation ---
  arma::vec alp = y / sig;
  arma::mat mat_temp = Phi_mat * D * Phi_mat.t() + eye<mat>(n, n);
  arma::vec vec_temp = alp - v;

  // --- Step 5: linear matrix system Ax=b get x=A^{-1}b ---
  arma::vec w_temp = solve(mat_temp, vec_temp, solve_opts::fast);  // refine/fast solution

  // --- Step 6: compute beta ---
  arma::vec beta = u + D * Phi_mat.t() * w_temp;

  // --- return the result---
  return beta;
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

// Inline helper: sample GMRF N(Q^{-1}b,Q^{-1}) with nearPD
/*** -------- nearPD helpers (minimal) -------- **/
// Access Matrix::nearPD
inline Function get_nearPD() {
  return Environment::namespace_env("Matrix")["nearPD"];
}

// Call nearPD safely: convert to base matrix in/out to avoid S4 class issues
inline arma::mat call_nearPD(const arma::mat& A) {
  Function nearPD   = get_nearPD();
  Function asMatrix = Environment::base_env()["as.matrix"];

  // Ensure input is a base R matrix
  SEXP Ain = asMatrix(wrap(A));

  // Run nearPD (let it symmetrize internally)
  SEXP resSEXP = nearPD(_["x"] = Ain, _["doSym"] = true);
  List res(resSEXP);

  // Extract repaired matrix (usually in element "x")
  SEXP X = res.containsElementNamed("x") ? (SEXP)res["x"] : (SEXP)res[0];

  // Convert back to base matrix, then Armadillo
  SEXP Xmat = asMatrix(X);
  arma::mat out = as<arma::mat>(Xmat);

  // Enforce numeric symmetry
  return 0.5 * (out + out.t());
}
/*** -------- end helpers -------- ***/


// Inline helper: sample GMRF x ~ N(Q^{-1} b, Q^{-1})
// Strategy: try Cholesky on symmetrized Q; if it fails, repair with nearPD; if still fails, stop.
inline arma::vec sample_GMRF(const arma::mat& Q,
                             const arma::vec& b,
                             bool verbose = false)
{
  const arma::uword n = Q.n_rows;
  if (Q.n_rows != Q.n_cols) Rcpp::stop("Q must be square.");
  if (b.n_elem != n)        Rcpp::stop("Dimension mismatch: length(b) != nrow(Q).");

  // 1) Symmetrize Q to reduce numerical asymmetry
  arma::mat Q_sym = 0.5 * (Q + Q.t());
  if (!Q_sym.is_finite() || !b.is_finite()) Rcpp::stop("Q or b has NaN/Inf.");

  // 2) Try Cholesky: Q = L L^T (lower-triangular L)
  arma::mat L;
  bool ok = arma::chol(L, Q_sym, "lower");


  // 3) If Cholesky fails, project to nearest PD via Matrix::nearPD
  if (!ok) {
    if (verbose) Rcpp::Rcout << "[sample_GMRF] Q not SPD; trying nearPD..." << std::endl;
    arma::mat Q_spd = call_nearPD(Q_sym);
    ok = arma::chol(L, Q_spd, "lower");
    if (!ok) Rcpp::stop("sample_GMRF failed: chol() still fails after nearPD.");
    if (verbose) Rcpp::Rcout << "[sample_GMRF] nearPD succeeded." << std::endl;
  }

  // 4) Solve for the mean: L w = b  -> w;  then  L^T mu = w  -> mu
  arma::vec w;
  if (!arma::solve(w, arma::trimatl(L), b)) {
    if (verbose) Rcpp::Rcout << "min diag(L) = " << L.diag().min() << "\n";
    Rcpp::stop("sample_GMRF: triangular solve (L w = b) failed.");
  }

  arma::vec mu;
  if (!arma::solve(mu, arma::trimatu(L.t()), w)) {
    if (verbose) Rcpp::Rcout << "min diag(L) = " << L.diag().min() << "\n";
    Rcpp::stop("sample_GMRF: triangular solve (L^T mu = w) failed.");
  }

  // 5) Draw standard normal z and solve L^T v = z  -> v
  arma::vec z = arma::randn<arma::vec>(n);

  arma::vec v;
  if (!arma::solve(v, arma::trimatu(L.t()), z)) {
    Rcpp::stop("sample_GMRF: triangular solve (L^T v = z) failed.");
  }

  // 6) Return x = mu + v
  return mu + v;
}


/**
 * @title Graphical Horseshoe MCMC Sampler (Tilde fast Sampler)
 *
 * Key Features:
 *   - Bhattacharya et al. fast sampling algorithm for high-dimensional regression
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

  // call rgig
  Function rgig("rgig", Environment::namespace_env("GIGrvg"));

  // ==== MCMC loop ====
  for (int m = 0; m < M; ++m) {
    // if (m % 100 == 0) {
    //   std::cout << "m = " << m << std::endl;
    // }
    std::cout << "\rIteration m = " << m << std::flush;

    arma::mat C(p, p, fill::zeros);

    for (int j = p-1; j >= 0; --j) {
      if (j == 0) {
        double shape = n/2.0 + 1;
        double rate = accu(Y_scaled.col(0) % Y_scaled.col(0)) / 2.0;
        Omega_tilde(0,0) = rgamma_std(shape, 1.0/rate, gen);
        Omega(0,0) = Omega_tilde(0,0) + C(0,0);
      } else {
        arma::vec Yj = Y_scaled.col(j);
        arma::mat Ysub = Y_scaled.cols(0, j-1);
        arma::vec omega_sub = Omega_tilde.submat(0, j, j-1, j);

        double a = dot(Yj, Yj);
        double b = dot(Ysub * omega_sub, Ysub * omega_sub);

        // using rgig from GIGrvg package
        NumericVector val = rgig(1, n/2.0 + 1, b, a);
        double sampled = val[0];
        Omega_tilde(j,j) = sampled;

        double sig = 1.0 / std::sqrt(Omega_tilde(j,j));
        arma::mat X = -Ysub * sig;
        arma::vec y = Yj - Ysub * C.submat(0, j, j-1, j) * (sig*sig);
        arma::vec B = tau2 * Lambda2_mat.submat(0, j, j-1, j);

        arma::vec beta;
        if (j > n) {
          beta = Bhattacharya_cpp(n, j, X, B, sig, y);
        } else {
          arma::mat A = X.t() * X + arma::diagmat(1.0 / B);
          beta = sample_GMRF(A / (sig * sig), (X.t() * y) / (sig * sig),true);
        }

        Omega_tilde.submat(0, j, j-1, j) = beta / sig - C.submat(0, j, j-1, j);
        Omega_tilde.submat(j, 0, j, j-1) = Omega_tilde.submat(0, j, j-1, j).t();

        arma::mat outer_prod = Omega_tilde.submat(0, j, j-1, j) * Omega_tilde.submat(0, j, j-1, j).t();
        C.submat(0, 0, j-1, j-1) += outer_prod / Omega_tilde(j,j);

        Omega.submat(0, j, j, j) = Omega_tilde.submat(0, j, j, j) + C.submat(0, j, j, j);
        Omega.submat(j, 0, j, j-1) = Omega.submat(0, j, j-1, j).t();

        // need to switch between GHS and GHSL for local parameter
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

            double rand_val = arma::randu<double>();  // Uniform random value between 0 and 1
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
  std::cout << std::endl;

  return List::create(
    Named("Omega_save") = Omega_save,
    Named("Lambda2_save") = Lambda2_save,
    Named("tau2_save") = tau2_save
  );
}
