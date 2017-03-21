# Loris Michel, EPFL
# locus horseshoe based on an exponential family reduction
locus_core_horseshoeExp <- function(Y, X, d, n, p, list_hyper, b_vb, c_vb, sigma2_bv, mu_beta_vb,
                        sig2_beta_vb, tau_vb, tol, maxit, batch, verbose, scheme = "noPrec",
                        full_output = FALSE) {

  # Y must have been centered, and X, standardized.

  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects

    # second moment of the \beta's
    m2_beta <- (mu_beta_vb ^ 2)  +  sig2_beta_vb
    mat_x_m1 <-  X %*% mu_beta_vb

    a_inv_vb <- A^{2} / 2
    converged <- FALSE
    lb_old <- -Inf
    it <- 1

    while ((!converged) & (it <= maxit)) {

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % # update of sigma^{-2}

      eta_vb <- (p*d+1)/2
      kappa_vb <- update_kappa_vb_horseshoeExp(a_inv_vb, m2_beta, b_vb, tau_vb, scheme)
      sig2_inv_vb <- eta_vb / kappa_vb

      # % # update of a^{-1}
      a_inv_vb <- 1 / (sig2_inv_vb + A^{-2})

      # % # # update of tau_{t}
      if(scheme == "noPrec") {
        lambda_vb <- lambda + (n/2)
      } else {
        lambda_vb <- lambda + ((n+p)/2)
      }

      nu_vb <- update_nu_vb_horseshoeExp(Y, X,  mat_x_m1, b_vb, d, n, p, mu_beta_vb, m2_beta, nu, sig2_inv_vb, scheme)
      tau_vb <- lambda_vb / nu_vb

      # % # update of the variance of the \beta's (inefficient replication of the tau value)
      if(scheme == "noPrec") {
        sig2_beta_vb <- 1 / ((n-1)*matrix(tau_vb,ncol=d, nrow=p, byrow = T) + (sig2_inv_vb * b_vb))
      } else {
        sig2_beta_vb <- 1 / (matrix(tau_vb,ncol=d, nrow=p, byrow = T)*((n-1) + (sig2_inv_vb * b_vb)))
      }




      if (batch) { # some updates are made batch-wise

        for (j in 1:p) {

          mat_x_m1 <- mat_x_m1 - tcrossprod(X[, j], mu_beta_vb[j, ])

          # % # update of the \mu_beta
          mu_beta_vb[j, ] <- sig2_beta_vb[j,] * (tau_vb *
                                               crossprod(Y - mat_x_m1, X[, j]))

          mat_x_m1 <- mat_x_m1 + tcrossprod(X[, j], mu_beta_vb[j, ])

          m2_beta <- (mu_beta_vb ^ 2)  +  sig2_beta_vb

        # % # update of the G values
          if(scheme == "noPrec") {
            G_vb[j,] <- (1/2)*sig2_inv_vb * m2_beta[j,]
          } else {
            G_vb[j,] <- (1/2)* sig2_inv_vb * tau_vb * m2_beta[j,]
          }

          # % # update of the b values



          #b_vb[j,] <- 1/(G_vb[j,] * Q_approx(G_vb[j,])) - 1
           b_vb[j,] <- 1/(G_vb[j,] + c_vb[j,])
           c_vb[j,] <- 1/(b_vb[j,]+1)



        }


      } else {

        for (k in 1:d) {

          vec_x_j_k <-  X %*% mu_beta_vb[, k]

          for (j in 1:p) {

            vec_x_j_k <- vec_x_j_k - X[, j] * mu_beta_vb[j, k]

            # % # update of the \mu values
            mu_beta_vb[j, k] <- sig2_beta_vb[k] * tau_vb[k] *
              crossprod(X[, j], Y[,k] - vec_x_j_k)

            vec_x_j_k <- vec_x_j_k + X[, j] * mu_beta_vb[j, k]

            m2_beta[j,k] <- (mu_beta_vb[j,k])^2 + sig2_beta_vb[j,k]

            # % # update of the G values
            G_vb[j,k] <- (1/1)*sig2_inv_vb * m2_beta[j,k]

            # % # update of the b values
            b_vb[j,k] <- (1 / (G[j,k] * exp(G[j,k]) * expint_E1(G[j,k]))) - 1
          }

        }

      }



    # % # computation of the lower bound
    lb_new <- lower_bound_horseshoeExp(Y, X, d, n, p, sig2_beta_vb, sig2_inv_vb, tau_vb,
                           eta, kappa, lambda, nu, c_vb,b_vb, mat_x_m1, mu_beta_vb,
                           m2_beta, G_vb,  a_inv_vb, scheme)


      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Lower bound = ", format(lb_new), "\n\n", sep = ""))

        converged <- (abs(lb_new - lb_old) < tol)


        lb_old <- lb_new
        it <- it + 1
    }



    if (verbose) {
      if (converged) {
        cat(paste("Convergence obtained after ", format(it),
                  " iterations with variational lower bound = ",
                  format(lb_new), ". \n\n",
                  sep = ""))

      } else {
        cat("Maximal number of iterations reached before convergence. Exit.")
      }
    }

    lb_opt <- lb_new





    if (full_output) { # for internal use only
      create_named_list_(mu_beta_vb, sig2_beta_vb, sig2_inv_vb, tau_vb,c_vb, b_vb,
                         eta_vb, kappa_vb, lambda_vb, nu_vb,
                         lambda, nu, a_inv_vb, A, m2_beta)
    } else {
      names_x <- colnames(X)
      names_y <- colnames(Y)

      rownames(mu_beta_vb) <- names_x
      colnames(mu_beta_vb) <- names_y
      rownames(b_vb) <- names_x
      colnames(b_vb) <- name_y
      #names(x_prpnst) <- names_x
      #names(y_prpnst) <- names_y

      create_named_list_(lb_opt, mu_beta_vb, b_vb, sig2_beta_vb, tau_vb)
    }
  })

}


update_kappa_vb_horseshoeExp <- function(a_inv, m2_beta, b_vb, tau_vb, scheme) {

  # % # updating the \kappa
  if(scheme == "noPrec") {
    as.numeric(a_inv + (1/2)*sum(m2_beta * b_vb))
  } else {
    as.numeric(a_inv + (1/2)*sum(tau_vb * colSums(m2_beta * b_vb)))
  }

}

update_nu_vb_horseshoeExp <- function(Y_mat, X_mat, mat_x_m1, b_vb, d, n, p, m1_beta,
                             m2_beta, nu, sig2_inv_vb, scheme) {
  # put X_mat and Y_mat instead of X and Y to avoid conflicts with the function sapply,
  # which has also an "X" argument with different meaning...

  # X must be standardized as we use (X \hadamard X)^T \one_n = (n-1)\one_p
  if(scheme == "noPrec") {
    nu_vb <- nu + (colSums(Y_mat ^ 2) - 2*colSums(Y_mat * (mat_x_m1)) +
                +  colSums((n-1) * m2_beta) +
                +  colSums(mat_x_m1^2) - (n - 1) * colSums(m1_beta^2))/ 2
  } else {

    nu_vb <- nu + (colSums(Y_mat ^ 2) - 2*colSums(Y_mat * (mat_x_m1)) +
                +  colSums(((n-1) + sig2_inv_vb * b_vb) * m2_beta) +
                +  colSums(mat_x_m1^2) - (n - 1) * colSums(m1_beta^2))/ 2

  }

}


# this function should be changed adequately
lower_bound_horseshoeExp <- function(Y, X, d, n, p, sig2_beta_vb, sig2_inv_vb, tau_vb,
                         eta, kappa, lambda, nu,c_vb,  b_vb, mat_x_m1, m1_beta,
                         m2_beta, G_vb, a_inv_vb, scheme) {

  # update for \tau_{t}
  if(scheme == "noPrec") {
    lambda_vb <- lambda + (n/2)
  } else {
    lambda_vb <- lambda + ((n+p)/2)
  }
  nu_vb <- update_nu_vb_horseshoeExp(Y, X, mat_x_m1, b_vb, d, n, p, m1_beta,
                         m2_beta, nu, sig2_inv_vb, scheme)

  # updates for \sigma^{-2}
  eta_vb <- (p*d+1)/2
  kappa_vb <- update_kappa_vb_horseshoeExp(a_inv_vb, m2_beta, b_vb, tau_vb, scheme)

  # mean of a^{-1}
  a_inv_vb <- 1 / (sig2_inv_vb + A^{-2})

  # log values for |tau_{t} and |sig^{-2}

  log_tau_vb <- digamma(lambda_vb) - log(nu_vb)
  log_sig2_inv_vb <- digamma(eta_vb) - log(kappa_vb)
  log_b_vb <- digamma(1) - log(G_vb + c_vb)
  log_c_vb <- digamma(1) - log(b_vb + 1)


 # do the computation for the beta's 
  L_A <- -(1/2)*(n*d)*log(2*pi) + (n/2)*sum(log_tau_vb) -  crossprod(tau_vb/2,(nu_vb - nu))
  
  if(scheme == "noPrec") {
    L_B <- -(1/2)*(p*d)*log(2*pi) + (1/2)*(p*d)*log_sig2_inv_vb - (1/2) * sum(b_vb * sig2_inv_vb * m2_beta)
  } else {
     L_B <- -(1/2)*(p*d)*log(2*pi) + (1/2)*(p*d)*log_sig2_inv_vb +(1/2)*(p*d)*log_tau_vb - (1/2) *  sig2_inv_vb * sum(tau_vb * colSums(b_vb * m2_beta))
  }
  H_beta <-  (1/2)*(sum(-log(sig2_beta_vb))- p*d*(log(2*pi)+1))
 
# do the computation for the tau's
 
   L_H_tau <- sum((lambda-lambda_vb) * log_tau_vb + lambda_vb * (1 - (nu/nu_vb)) - lgamma(lambda) + lgamma(lambda_vb) +
       lambda*log(nu) - lambda_vb*log(nu_vb))
# do computation for sigma

 L_sigma <- -(1/2)*log_sig2_inv_vb - a_inv_vb*sig2_inv_vb -lgamma(1/2) # + (1/2)*log_a_inv compensate
 H_sigma <- (eta_vb -1)*log_sig2_inv_vb + eta_vb -lgamma(eta_vb) + eta_vb*log(kappa_vb)

# do the computation for a^{-1}

 L_a <- -A^{-2}*a_inv_vb -lgamma(1/2) -log(A)
 H_a <- -1 + log(A^{-2}+sig2_inv_vb)-lgamma(1)

# do the computation for b 

H_b <- sum(-1 + log(G_vb + c_vb) - lgamma(1))
L_b <- sum(-(1/2)*log_b_vb - c_vb*b_vb - lgamma(1/2) + (1/2)*log_c_vb) 

# do the computation for c
H_c <- sum(-1 + log(b_vb+1) - lgamma(1))
L_c <- sum(-(1/2)*log_c_vb - c_vb - lgamma(1/2) + (1/2)*log(1))


l <- L_A + L_B + L_sigma + L_H_tau + L_a + L_b + L_c - H_beta - H_sigma - H_a - H_b - H_c
 # lower bound as in the paper
#  if(scheme == "noPrec") {
#    l <- (-1)*(n*d/2)* log(2*pi) + p*d/2 + (n/2)*sum(log_tau_vb) - (1/2)*sum(log(sig2_beta_vb)) -
#       crossprod(tau_vb/2,(nu_vb - nu)) +
#       +(p*d/2)*log_sig2_inv_vb - (1/2) * sum(b_vb * sig2_inv_vb * m2_beta) +
#       sum((lambda-lambda_vb) * log_tau_vb + lambda_vb * (1 - (nu/nu_vb)) - lgamma(lambda) + lgamma(lambda_vb) +
#       lambda*log(nu) - lambda_vb*log(nu_vb)) +
#       ((1/2)-eta_vb)*log_sig2_inv_vb - a_inv_vb*(sig2_inv_vb + A^{-2}) -2*lgamma(1/2) + (1/2)*log(A^{-2}) +
#       1 - log(sig2_inv_vb + A^{-2}) + eta_vb + lgamma(eta_vb) - eta_vb*log(kappa_vb) - (p*d*log(pi)) +
#       sum(log(Q_approx(G_vb))) + p*d - sum(G_vb*Q_approx(G_vb))
#  } else {
#    l <- (-1)*(n*d/2)* log(2*pi) + p*d/2 + ((n+p)/2)*sum(log_tau_vb) - (1/2)*sum(log(sig2_beta_vb)) -
#        crossprod(tau_vb/2,(nu_vb - nu)) +
#       + (p*d/2)*log_sig2_inv_vb - (1/2) *  sig2_inv_vb * sum(tau_vb * colSums(b_vb * m2_beta)) +
#       sum((lambda-lambda_vb) * log_tau_vb + lambda_vb * (1 - (nu/nu_vb)) - lgamma(lambda) + lgamma(lambda_vb) +
#       lambda*log(nu) - lambda_vb*log(nu_vb)) +
#       ((1/2)-eta_vb)*log_sig2_inv_vb - a_inv_vb*(sig2_inv_vb + A^{-2}) -2*lgamma(1/2) + (1/2)*log(A^{-2}) +
#       1 - log(sig2_inv_vb + A^{-2}) + eta_vb + lgamma(eta_vb) - eta_vb*log(kappa_vb) - (p*d*log(pi)) +
#      sum(log(Q_approx(G_vb))) + p*d - sum(G_vb*Q_approx(G_vb))
#  }

  return(l)

}

