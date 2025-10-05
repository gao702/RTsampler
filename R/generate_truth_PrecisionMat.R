#' Generate a ground-truth precision matrix under common sparse structures
#'
#' @description
#' Construct a \eqn{p \times p} positive-definite precision matrix \eqn{\Theta_0}
#' (the inverse covariance) with one of four toy graph structures:
#' \itemize{
#'   \item \code{"tridiagonal"} - unit diagonal with band-1 off-diagonals (0.25);
#'   \item \code{"hubs"} - nodes partitioned into equal-size groups with a hub node per group connected to all others (weight 0.25);
#'   \item \code{"cliques neg"} - random 3-node cliques with **negative** partial correlations (off-diagonals -0.45);
#'   \item \code{"cliques pos"} - random 3-node cliques with **positive** partial correlations (off-diagonals 0.75).
#' }
#'
#'
#' @param p Integer. Dimension of the precision matrix.
#' @param structure Character. One of
#'   \code{c("tridiagonal","hubs","cliques neg","cliques pos")}.
#'
#' @return A numeric \eqn{p \times p} symmetric positive-definite matrix
#'   representing the precision matrix.
#'
#'
#' @examples
#' p <- 10
#' Om1 <- generate_truth_PrecisionMat(p, "tridiagonal")
#' Om2 <- generate_truth_PrecisionMat(p, "hubs")
#' Om3 <- generate_truth_PrecisionMat(p, "cliques neg")
#' Om4 <- generate_truth_PrecisionMat(p, "cliques pos")
#'
#' @export

generate_truth_PrecisionMat <- function(p, structure = c("tridiagonal", "hubs", "cliques neg", "cliques pos")) {
  structure <- match.arg(structure)

  if (structure == "tridiagonal") {
    sigma_inv <- diag(1, p)
    offdiag_val = 0.25
    for (i in 1:(p-1)) {
      sigma_inv[i, i+1] <- offdiag_val
      sigma_inv[i+1, i] <- offdiag_val
    }
  }

  if (structure == "hubs") {
    sigma_inv <- diag(1, p, p)
    alpha <- 0.25
    num_hub <- p / 10
    for (i in 1:num_hub) {
      start <- 1 + p / num_hub * (i - 1)
      end <- p / num_hub + p / num_hub * (i - 1)
      sigma_inv[start, (start + 1):end] <- alpha
      sigma_inv[(start + 1):end, start] <- alpha
    }
  }

  if (structure == "cliques neg") {
    k=3    #number of features in a group; ten features make a group
    # alpha=-0.45 for positive correlation; alpha=0.75 for negative correlation
    alpha = rep(0.75,k*(k-1)/2)    #for magnitude of nonzero partial covariance
    # alpha = rep(-0.45,k*(k-1)/2)
    set.seed(2016)
    sigma_inv <- matrix(0, p, p)
    for (i in 1:(p / 10)) {
      rind <- sample((1 + 10 * (i - 1)):(10 + 10 * (i - 1)), size = k)
      submat <- matrix(0, k, k)
      submat[upper.tri(submat)] <- alpha
      submat <- submat + t(submat)
      sigma_inv[rind, rind] <- sigma_inv[rind, rind] + submat
    }
    diag(sigma_inv) <- 1
  }

  if (structure == "cliques pos") {
    k=3    #number of features in a group; ten features make a group
    # alpha=-0.45 for positive correlation; alpha=0.75 for negative correlation
    # alpha = rep(0.75,k*(k-1)/2)    #for magnitude of nonzero partial covariance
    alpha = rep(-0.45,k*(k-1)/2)
    set.seed(2016)
    sigma_inv <- matrix(0, p, p)
    for (i in 1:(p / 10)) {
      rind <- sample((1 + 10 * (i - 1)):(10 + 10 * (i - 1)), size = k)
      submat <- matrix(0, k, k)
      submat[upper.tri(submat)] <- alpha
      submat <- submat + t(submat)
      sigma_inv[rind, rind] <- sigma_inv[rind, rind] + submat
    }
    diag(sigma_inv) <- 1
  }

  # Check if positive definite
  eigvals <- eigen(sigma_inv, only.values = TRUE)$values
  if (any(eigvals <= 0)) stop("Generated sigma_inv is not positive definite.")

  return(sigma_inv)
}
