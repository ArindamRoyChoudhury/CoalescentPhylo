#' Probability of lineage reduction under the coalescent
#'
#' Computes the probability that the number of ancestral lineages decreases
#' from \code{i.current.n} to \code{i.dash.previous.n} after time \code{tau}
#' under the standard neutral coalescent process. The calculation uses the
#' eigenvalue formulation of the coalescent transition probabilities, where
#' the coalescence rates are \eqn{\lambda_k = k(k-1)/2}.
#'
#' This function is used internally when computing likelihoods over
#' population trees by integrating over possible ancestral lineage counts.
#'
#' @param i.dash.previous.n Integer. Number of ancestral lineages after time \code{tau}.
#' @param i.current.n Integer. Number of lineages at the start of the time interval.
#' @param tau Numeric. Time interval (coalescent time units).
#'
#' @return Numeric value giving the probability of transitioning from
#' \code{i.current.n} lineages to \code{i.dash.previous.n} lineages after
#' time \code{tau}.
#'
#' @details
#' The probability is computed using the closed-form expression for the
#' coalescent transition probabilities based on exponential waiting times
#' between coalescent events.
#'
#' @noRd

prob.n.coal <- function(i.dash.previous.n, i.current.n, tau )
{

  one.to.current.n <- (1:i.current.n)

  lambda <- (one.to.current.n)*(one.to.current.n - 1)/2


  prod.lambda <- 1

  if ((i.dash.previous.n+1) <= (i.current.n))
    prod.lambda <- prod( lambda[ (i.dash.previous.n+1):(i.current.n) ] )


  sum.of.terms <- 0

  for (j in i.dash.previous.n:i.current.n)
  {
    numerator <- exp(-lambda[j]*tau)

    index_to_remove <- (j - i.dash.previous.n + 1)
    indices_to_multiple <- ((i.dash.previous.n):(i.current.n))[-index_to_remove]

    denominator <- prod( lambda[indices_to_multiple] - lambda[j] )

    sum.of.terms <- sum.of.terms + (numerator/denominator)
  }

  final.prob <- prod.lambda*sum.of.terms

  return(final.prob)
}

#' Joint log-likelihood under a four-parameter coalescent model
#'
#' Computes the joint log-likelihood of allele count data across loci under
#' a four-parameter coalescent model with divergence times \code{tau_a},
#' \code{tau_b}, \code{tau_c}, and \code{tau_x}. The likelihood is calculated
#' by integrating over possible ancestral lineage counts using dynamic
#' programming based on the formulation in RoyChoudhury et al. (2008).
#'
#' For each locus, the function evaluates the probability of observing the
#' allele counts in three populations (A, B, C) given the coalescent process
#' along the species tree. Intermediate likelihood values are cached to avoid
#' redundant computation when identical allele configurations appear across
#' loci.
#'
#' @param tau Numeric vector of length 4 representing divergence time parameters:
#' \describe{
#'   \item{`tau[1]`}{Coalescent time for population A}
#'   \item{`tau[2]`}{Coalescent time for population B}
#'   \item{`tau[3]`}{Coalescent time for population C}
#'   \item{`tau[4]`}{Coalescent time for the internal ancestral population}
#' }
#' @param theta Mutation rate parameter of the Beta-Binomial model.
#' @param r_matrix Matrix of allele counts with three columns corresponding
#' to populations A, B, and C. Each row represents a locus.
#' @param n Integer sample size per population.
#'
#' @return A numeric value representing the joint log-likelihood across loci.
#' @noRd
#'
#' @details
#' The likelihood calculation follows the coalescent framework described in
#' RoyChoudhury et al. (2008), where probabilities of ancestral lineage counts
#' are computed using recursive dynamic programming. The method integrates
#' over all possible ancestral configurations and uses hypergeometric and
#' Beta-Binomial distributions to model allele sampling and mutation processes.
#'
#' Intermediate likelihoods for identical allele count configurations
#' \code{(r_A, r_B, r_C)} are stored to reduce repeated computation across loci.
#'
#' @references
#' RoyChoudhury, A., Felsenstein, J., & Thompson, E. A. (2008).
#' A two-stage pruning algorithm for likelihood computation for a population
#' tree. \emph{Genetics}.
#' @import stats

jointloglikelihood_tau_4param <- function(tau, theta, r_matrix, n)
{
  eps <- 1e-8

  tau_a <- tau[1]
  tau_b <- tau[2]
  tau_c <- tau[3]
  tau_x <- tau[4]


  joint_loglikelihood <- 0

  num_locus <- dim(r_matrix)[1]

  store_loglikelihood <- array(dim=c(n+1,n+1,n+1))

  for (i_locus in 1:num_locus)
  {

    r_A <- r_matrix[i_locus,1]
    r_B <- r_matrix[i_locus,2]
    r_C <- r_matrix[i_locus,3]

    if( (!is.na(r_A)) && (!is.na(r_B)) && (!is.na(r_C)) )
    {

      if (is.na(store_loglikelihood[r_A+1,r_B+1,r_C+1]))
      {

        ArrayA_A <- integer(n)
        ArrayA_A[n] <- 1

        ArrayB_A <- array(dim=c(n,n + 1))
        for (i in 1:n)
          ArrayB_A[i,(1:(i+1))] <- 0
        ArrayB_A[n,(r_A + 1)] <- 1



        ArrayA_B <- integer(n)
        ArrayA_B[n] <- 1

        ArrayB_B <- array(dim=c(n,n + 1))
        for (i in 1:n)
          ArrayB_B[i,(1:(i+1))] <- 0
        ArrayB_B[n,(r_B + 1)] <- 1



        ArrayA_C <- integer(n)
        ArrayA_C[n] <- 1


        ArrayB_C <- array(dim=c(n,n + 1))
        for (i in 1:n)
          ArrayB_C[i,(1:(i+1))] <- 0
        ArrayB_C[n,(r_C + 1)] <- 1



        if(tau_a == 0){

          ArrayA_XL <- ArrayA_A
          ArrayB_XL <- ArrayB_A

        }else if (tau_a != 0){

          ArrayA_XL <- array(dim=n)
          ArrayB_XL <- array(0,dim=c(n, n + 1))

          for (n_X_left in 1:n)
          {ArrayA_XL[n_X_left] <- prob.n.coal(n_X_left, n, tau_a)}

          if(r_A == 0){
            ArrayB_XL[,1] <- 1
          } else if (r_A == n){
            for (i in 1:n){
              ArrayB_XL[i, i + 1] <- 1
            }
          } else {


            for (n_X_left in 1:n)
            {
              max_r_X_left <- min(r_A,n_X_left)
              min_r_X_left <- max(0,(n_X_left - (n - r_A)))

              for (i in min_r_X_left:max_r_X_left){
                {ArrayB_XL[n_X_left, (i + 1)] <- (beta( r_A , n - r_A )/ beta( i , n_X_left - i )) * choose(n - n_X_left , r_A - i)}

              }
            }
          }
        }



        if (tau_b == 0){
          ArrayA_XR <- ArrayA_B
          ArrayB_XR <- ArrayB_B

        } else if (tau_b != 0){


          ArrayA_XR <- array(dim=n)
          ArrayB_XR <- array(0,dim=c(n, n + 1))

          for (n_X_right in 1:n){
            ArrayA_XR[n_X_right] <- prob.n.coal(n_X_right, n, tau_b)
          }


          if (r_B == 0){
            ArrayB_XR[,1] <- 1
          } else if (r_B == n){
            for (i in 1:n){
              ArrayB_XR[i, i + 1] <- 1
            }
          } else {

            for (n_X_right in 1:n)
            {
              max_r_X_right <- min(r_B,n_X_right)
              min_r_X_right <- max(0,(n_X_right - (n - r_B)))

              for (i in min_r_X_right:max_r_X_right){
                {ArrayB_XR[n_X_right, (i + 1)] <- (beta( r_B , n - r_B )/ beta( i , n_X_right - i )) * choose(n - n_X_right , r_B - i)}

              }
            }
          }
        }
        #r_x_right is i now

        ArrayA_X <- integer(2*n)

        for (n_X_left in 1:n){
          for (n_X_right in 1:n)
          {
            n_X <- n_X_left + n_X_right
            ArrayA_X[n_X] <- ArrayA_X[n_X] + ArrayA_XL[n_X_left]*ArrayA_XR[n_X_right]
          }
        }


        ArrayB_X <- array(0,dim=c(2*n,(2*n + 1)))
        ## Eq (7) at RoyChoudhury et al. 2008


        for (i in 2:(2*n))
          if(ArrayA_X[i] > 0)
          {
            for (j in 0:i){
              for (i_dash in 1:(i - 1)){
                for (j_dash in 0:i_dash){
                  for (j_dash_dash in 0:(i - i_dash)){
                    if ((i_dash <= n) && ( (i-i_dash) <= n ))
                      if ((j_dash + j_dash_dash) == j)
                        ArrayB_X[i,(j + 1)] <- ArrayB_X[i, (j + 1)] + ArrayB_XL[i_dash,(j_dash + 1)] * ArrayB_XR[(i - i_dash),(j_dash_dash + 1)] * dhyper(x = j_dash, m = j, n = i - j, k = i_dash) * ArrayA_XL[i_dash] * ArrayA_XR[(i - i_dash)] / ArrayA_X[i]
                  }
                }
              }
            }
          }

        #Here n represents n and i. r_x_left's definition seems like the same as j_dash_dash? nope, j_dash_dash is r_x_right j_dash is r_x_left.n_x_left is i_dash



        if (tau_x == 0){
          ArrayA_rootL <- ArrayA_X
          ArrayB_rootL <- ArrayB_X

        } else if (tau_x != 0){

          ArrayA_rootL <- array(0,dim=c(2*n))
          ## Eq (4) at RoyChoudhury et al. 2008

          for (i_dash in 1:(2*n))
          {
            for (i in i_dash:(2*n))
            {
              ArrayA_rootL[i_dash] <- ArrayA_rootL[i_dash] + ArrayA_X[i] * prob.n.coal(i_dash, i, tau_x)
            }
          }

          ArrayB_rootL <- array(0,dim=c(2*n,(2*n + 1)))
          ## Eq (5) at RoyChoudhury et al. 2008
          for(i in 1:(2*n))
          {
            if ((r_A + r_B) == 0)
            {ArrayB_rootL[i,1] <- 1}
            else if ((r_A + r_B) == (2*n))
            {ArrayB_rootL[i,(i+1)] <- 1}
            else if (i >= 2)
            {
              for(j in 1:(i-1))
              {
                for(i_dash in i:(2*n))
                {
                  for(j_dash in 1:(i_dash-1))
                    if ((j_dash >= j) && ((i_dash-j_dash) >= (i-j)))
                      ArrayB_rootL[i,j + 1] <- ArrayB_rootL[i,j + 1] + ArrayB_X[i_dash, j_dash + 1] * (prob.n.coal(i, i_dash, tau_x) * ArrayA_X[i_dash]/ArrayA_rootL[i]) * (beta( j_dash , i_dash - j_dash )/beta( j ,  i - j)) * choose(i_dash - i, j_dash - j )
                  #reverse the prob of n here.
                }
              }
            }
          }
        }



        if (tau_c == 0){
          ArrayA_rootR <- ArrayA_C
          ArrayB_rootR <- ArrayB_C
        } else if ( tau_c != 0){

          ArrayA_rootR <- array(dim=n)

          for (n_root_right in 1:n){
            ArrayA_rootR[n_root_right] <- prob.n.coal(n_root_right, n, tau_c)
          }

          ArrayB_rootR <- array(0,dim=c(n,n + 1))

          if (r_C == 0){
            ArrayB_rootR[,1] <- 1
          } else if (r_C == n){
            for (i in 1:n){
              ArrayB_rootR[i, i + 1] <- 1
            }
          } else {

            for (n_root_right in 1:n)
            {
              max_r_root_right <- min(r_C,n_root_right)
              min_r_root_right <- max(0,(n_root_right - (n - r_C)))

              for (r_root_right in min_r_root_right:max_r_root_right){
                ArrayB_rootR[n_root_right, (r_root_right + 1)] <- (beta( r_C , n - r_C )/ beta( r_root_right , n_root_right - r_root_right )) * choose(n - n_root_right , r_C - r_root_right)
              }
            }
          }
        }

        ArrayA_root <- integer(3*n)

        for (n_root_left in 1:(2*n)){
          for (n_root_right in 1:n)
          {
            n_root <- n_root_left + n_root_right
            ArrayA_root[n_root] <- ArrayA_root[n_root] + ArrayA_rootL[n_root_left]*ArrayA_rootR[n_root_right]
          }
        }



        ArrayB_root <- array(0, dim = c(3*n, 3*n + 1))
        ## Eq (7) at RoyChoudhury et al. 2008

        for (i in 2:(3*n))
          if(ArrayA_root[i] > 0)
          {
            for (j in 0:i)
            {
              for (i_dash in 1:(i - 1))
              {
                for (j_dash in 0:i_dash)
                {
                  for (j_dash_dash in 0:(i - i_dash))
                  {

                    if ((i_dash <= (2*n)) && ( (i-i_dash) <= n ))
                      if ((j_dash + j_dash_dash) == j)
                      {
                        add_term <- ArrayB_rootL[i_dash, j_dash + 1] * ArrayB_rootR[i-i_dash, j_dash_dash + 1] * dhyper(x = j_dash, m = j, n = i - j, k = i_dash) * ArrayA_rootL[i_dash] * ArrayA_rootR[(i - i_dash)] / ArrayA_root[i]
                        ArrayB_root[i, j + 1] <- ArrayB_root[i, j + 1] + add_term
                      }
                  }
                }
              }
            }
          }


        likelihood_m_zero <- 0
        for(i in 1:(3*n)){
          for(j in 0:i){
            likelihood_m_zero <- likelihood_m_zero + ArrayB_root[i,j+1] * ArrayA_root[i] * choose(i, j) * (beta(j + theta, i - j + theta)/beta(theta, theta))
          }
        }

        loglikelihood <- log(likelihood_m_zero)

        store_loglikelihood[r_A+1,r_B+1,r_C+1] <- loglikelihood
      }
      else
        loglikelihood <- store_loglikelihood[r_A+1,r_B+1,r_C+1]

      joint_loglikelihood <- joint_loglikelihood + loglikelihood
    }
  }

  return(joint_loglikelihood)
}
