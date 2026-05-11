#' Numerical gradient of the joint log-likelihood
#'
#' Computes the numerical gradient of the joint log-likelihood with respect
#' to the divergence time parameters \code{tau}. This function serves as a
#' wrapper around the gradient calculation for
#' \code{jointloglikelihood_tau_4param()}.
#'
#' @param tau Numeric vector of length 4 containing the divergence time
#' parameters.
#' @param theta Numeric mutation parameter used in the likelihood model.
#' @param r_matrix Matrix of allele counts with three columns corresponding
#' to populations A, B, and C.
#' @param n Integer sample size per population.
#'
#' @return Numeric vector containing the gradient of the joint log-likelihood
#' with respect to the four \code{tau} parameters.
#'
#' @import rootSolve
#' @noRd

gradient_formula <- function(tau, theta, r_matrix, n) {
  initial_guess <- c(tau[1], tau[2], tau[3], tau[4])
  gradientvalue <- rootSolve::gradient(
    f = jointloglikelihood_tau_4param,
    x = initial_guess,
    theta = theta,
    r_matrix = r_matrix,
    n = n
  )
  return(gradientvalue)
}

#' Estimate internal branch lengths between population pairs
#'
#' Estimates the internal branch length parameter (\code{tau_x}) for all
#' population pairs using allele count data across loci. For each pair of
#' populations \code{(f, g)}, the function constructs a three-population
#' configuration consisting of populations \code{f}, \code{g}, and the
#' outgroup (assumed to be the \code{(P + 1)}-th row of
#' \code{mat_allele_count}). The internal branch length is then estimated
#' by maximizing the joint log-likelihood under the four-parameter model
#' using numerical optimization.
#'
#' The optimization is performed using the \code{"L-BFGS-B"} method with
#' numerical gradients, and the computation is parallelized across
#' population pairs.
#'
#' @param mat_allele_count Matrix of allele counts where rows correspond to
#' populations and columns correspond to loci.
#' @param P Integer. Number of populations (excluding the outgroup).
#' @param theta Numeric mutation parameter used in the likelihood model.
#' @param fixed.n.at.tips Integer sample size per population.
#' @param n.cores Optional integer specifying the number of CPU cores used
#' for parallel computation. Defaults to \code{detectCores() - 1}.
#'
#' @return A symmetric \code{P x P} matrix where element \code{(i, j)}
#' represents the estimated internal branch length \code{tau_x} between
#' populations \code{i} and \code{j}.
#'
#' @details
#' For each pair of populations, the function constructs a submatrix of
#' allele counts for populations \code{(g, f)} and the outgroup, and
#' estimates the parameters of the four-parameter coalescent model using
#' \code{optim()}. Only the estimated internal branch length
#' (\code{tau_x}) is stored.
#'
#' Parallel computation is implemented using \code{pbapply::pblapply()},
#' with explicit cluster support on Windows systems.
#'
#' @import pbapply
#' @import parallel
#' @import rootSolve
#' @import stats
#'
#' @noRd

estimate_tau_x_paralell <- function(mat_allele_count,
                           P,
                           theta,
                           fixed.n.at.tips,
                           n.cores = NULL) {

  # Generate all pairs where f > g
  pairs <- which(outer(1:P, 1:P, ">"), arr.ind = TRUE)

  # Set up cores
  if (is.null(n.cores)) {
    n.cores <- max(1, parallel::detectCores() - 1)
  }

  # Grab references to internal functions once, in the parent environment.
  # This ensures worker processes (both fork and PSOCK) can always find them,
  # regardless of whether the code is run as a script or as an installed package.
  .jll  <- jointloglikelihood_tau_4param
  .grad <- gradient_formula

  # Worker function for a single pair
  compute_pair <- function(idx) {
    f <- pairs[idx, 1]
    g <- pairs[idx, 2]
    tip.allele.counts.submatrix <- cbind(
      mat_allele_count[g, ],
      mat_allele_count[f, ],
      mat_allele_count[P + 1, ]
    )
    eps <- 1e-10
    result <- optim(
      par = c(0, 0, 0, 0),
      fn = .jll,
      gr = .grad,
      method = "L-BFGS-B",
      lower = rep(eps, 4),
      upper = rep(Inf, 4),
      control = list(fnscale = -1),
      theta = theta,
      r_matrix = tip.allele.counts.submatrix,
      n = fixed.n.at.tips
    )
    return(list(f = f, g = g, x_est = result$par[4]))
  }

  message(sprintf("Running %d pairs on %d cores...", nrow(pairs), n.cores))

  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl))

    # Export all variables needed by compute_pair, including the resolved
    # function references. envir = environment() captures .jll and .grad
    # along with the data variables, so workers get everything they need
    # without relying on package namespace lookup.
    parallel::clusterExport(
      cl,
      varlist = c(
        "pairs", "mat_allele_count", "P", "theta", "fixed.n.at.tips",
        ".jll", ".grad"
      ),
      envir = environment()
    )

    # rootSolve must be loaded on each worker because gradient_formula
    # calls rootSolve::gradient() internally.
    parallel::clusterEvalQ(cl, library(rootSolve))

  } else {
    # Mac / Linux: mclapply-style fork — child processes inherit the full
    # parent environment, so no export is needed.
    cl <- n.cores
  }

  # pblapply works with both a cluster object (Windows) and an integer (Mac/Linux)
  results <- pbapply::pblapply(seq_len(nrow(pairs)), compute_pair, cl = cl)

  # Fill symmetric matrix
  mat_est_tau_x <- matrix(nrow = P, ncol = P)
  for (res in results) {
    mat_est_tau_x[res$f, res$g] <- res$x_est
    mat_est_tau_x[res$g, res$f] <- res$x_est
  }

  return(mat_est_tau_x)
}
