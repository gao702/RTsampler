#' Cyclical MCMC sampler (column-wise, data-augmented Gibbs)
#'
#' \code{Cyclical_sampler} returns a list with:
#' \itemize{
#'   \item \code{omega}: array \eqn{p \times p \times M-burnin} of saved precision matrices;
#'   \item \code{lambda_sq}: array \eqn{p \times p \times M-burnin} of local parameters；
#'   \item \code{tau_sq}: numeric vector (\eqn{M-burnin}) of the global parameter.
#' }
#'
#' @param S sample covariance matrix
#' @param n sample size
#' @param M Total number of MCMC iterations
#' @param burnin number of MCMC burnins
#' @param prior character, one of \code{"GHS"} or \code{"GHSL"}.
#' @return A list \code{list(omega, lambda_sq, tau_sq)} as described above.
#' @export
#' @importFrom stats rgamma runif

Cyclical_sampler <- function(S,n,M, burnin,prior)
{
  if (!(prior %in% c("GHS", "GHSL"))) {
    stop("`prior` must be either \"GHS\" or \"GHSL\".")
  }
  p <- nrow(S)
  k = p
  S = S / k
  nmc = M - burnin
  omega_save <- array(0, dim=c(p,p,nmc))
  lambda_sq_save <- array(0,dim=c(p,p,nmc)) # save matrix instead of vector
  tau_sq_save <- array(0,dim=c(1,nmc))

  ind_all = array(0,dim=c(p-1,p))
  for (i in seq(1,p,1))
  {
    if (i==1)
    {
      ind <- t(array(seq(2,p,1),dim=c(1,p-1)))
    }else if (i==p)
    {
      ind <- t(array(seq(1,p-1,1),dim=c(1,p-1)))
    }else
    {
      ind <- t(array(c(seq(1,i-1,1),seq(i+1,p,1)),dim=c(1,p-1)))
    }
    ind_all[,i] <- ind
  }

  # set initial values
  Omega <- diag(p)
  Sigma <- diag(p)
  tau_sq <- 1
  xi <- 1
  Lambda_sq <- array(1,dim=c(p,p))
  Nu <- array(1,dim=c(p,p))

  for (iter in seq(1,M,1))
  {
    if (iter %% 1000 == 0) {
      cat(sprintf("iter = %d\n", iter))
    }
    for (i in seq(1,p,1))
    {
      ind <- ind_all[,i]
      Sigma_11 <- Sigma[ind,ind]
      sigma_12 <- Sigma[ind,i]
      sigma_22 <- Sigma[i,i]
      s_21 <- S[ind,i]
      s_22 <- S[i,i];
      lambda_sq_12 <- Lambda_sq[ind,i]
      nu_12 <- Nu[ind,i];

      # sample gamma and beta, gamma in the code is gamma_sample here and beta in code is beta_def here
      gamma_sample <- rgamma(1,(n/2)+1,s_22/2) # random gamma with shape=n/2+1, rate=2/s_22
      inv_Omega_11 <- Sigma_11 - (sigma_12%*%t(sigma_12))/sigma_22
      inv_C <- s_22*inv_Omega_11 + diag(1/(as.vector(lambda_sq_12)*tau_sq))
      inv_C_chol <- chol(inv_C)
      mu_i <- matrix(solve(-inv_C,s_21,tol = 1e-25))
      beta_def <- mu_i+ matrix(solve(inv_C_chol,rnorm(p-1,0,1)))
      omega_12 <- beta_def
      omega_22 <- gamma_sample + (t(beta_def)%*%inv_Omega_11%*%beta_def)
      # sample local parameter
      if(prior == "GHS" ){
        rate <- omega_12^2/(2*tau_sq) + 1/nu_12
        lambda_sq_12 <- matrix(1/rgamma(rep(1,dim(rate)[1]),rep(1,dim(rate)[1]),rate=rate),dim(rate)[1],1)  #' random inv gamma with shape=1, rate=rate
        nu_12 <- matrix(1/rgamma(rep(1,dim(lambda_sq_12)[1]),rep(1,dim(lambda_sq_12)[1]),rate=1+1/lambda_sq_12),dim(lambda_sq_12)[1],1)#' random inv gamma with shape=1, rate=1+1/lambda_sq_12
      }else{ # GHSL
        rate <- omega_12^2/(2*tau_sq) + nu_12/2
        lambda_sq_12 <- matrix(1/rgamma(rep(1,dim(rate)[1]),rep(1,dim(rate)[1]),rate=rate),dim(rate)[1],1)  #' random inv gamma with shape=1, rate=rate
        nu_12 <- matrix(-2 * lambda_sq_12 * log(1 - runif(p-1) * (1 - exp(-0.5 / lambda_sq_12))),dim(rate)[1],1) # m in [0,1]
      }


      # update Omega, Sigma, Lambda_sq, Nu
      Omega[i,ind] <- omega_12
      Omega[ind,i] <- omega_12
      Omega[i,i] <- omega_22
      temp <- inv_Omega_11%*%beta_def
      Sigma_11 <- inv_Omega_11 + temp%*%t(temp)/gamma_sample
      sigma_12 <- -temp/gamma_sample
      sigma_22 <- 1/gamma_sample
      Sigma[ind,ind] <- Sigma_11
      Sigma[i,i] <- sigma_22
      Sigma[i,ind] <- sigma_12
      Sigma[ind,i] <- sigma_12
      Lambda_sq[i,ind] <- lambda_sq_12
      Lambda_sq[ind,i] <- lambda_sq_12
      Nu[i,ind] <- nu_12
      Nu[ind,i] <- nu_12
    }
    # Sample tau_sq and xi
    omega_vector <- matrix(Omega[lower.tri(Omega)])
    lambda_sq_vector <- matrix(Lambda_sq[lower.tri(Omega)])
    rate <- (1/xi) + sum((omega_vector^2)/(2*lambda_sq_vector));
    tau_sq = 1/rgamma(1,((p*(p-1)/2)+1)/2,rate)
    xi = 1/rgamma(1,1,1+1/tau_sq);    #' inv gamma w/ shape=1, rate=1+1/tau_sq

    # save Omega, lambda_sq, tau_sq
    if (iter >burnin)
    {
      omega_save[,,iter-burnin] <- Omega / k
      lambda_sq_save[,,iter-burnin] <- Lambda_sq
      tau_sq_save[iter-burnin] <- tau_sq
    }
  }
  output <- list(omega=omega_save,lambda_sq=lambda_sq_save,tau_sq=tau_sq_save)
  return(output)
}
