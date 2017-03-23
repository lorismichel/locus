#' @title HS locus Exp
#' @description HS locus Exp
#' @export
#'
locus_core_horseshoeExp <- function(Y, X, d, n, p, list_hyper, b_vb, c_vb, sigma2_vb, mu_beta_vb,
                        sig2_beta_vb, tau_vb, tol, maxit, batch, verbose, scheme = "noPrec",loop="c++",
                        full_output = FALSE) {

  # Y must have been centered, and X, standardized.

  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects

    # second moment of the \beta's
    m2_beta <- (mu_beta_vb ^ 2)  +  sig2_beta_vb
    mat_x_m1 <-  X %*% mu_beta_vb
    a_inv_vb <- 1 / 2*(A^{-2})

    converged <- FALSE
    lb_old <- -Inf
    it <- 1
    ELBO <- c()

    while ((!converged) & (it <= maxit)) {

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % # update of sigma^{-2}
      eta_vb <- (p*d+1)/2
      kappa_vb <- update_kappa_vb_horseshoeExp(a_inv_vb, m2_beta,
                                               b_vb, tau_vb, scheme)
      sig2_inv_vb <- eta_vb / kappa_vb

      # % # update of a^{-1}
      a_inv_vb <- 1 / (sig2_inv_vb + A^{-2})

      # % # # update of tau_{t}
      if(scheme == "noPrec") {
        lambda_vb <- lambda + (n/2)
      } else {
        lambda_vb <- lambda + ((n+p)/2)
      }

      nu_vb <- update_nu_vb_horseshoeExp(Y, X,  mat_x_m1, b_vb, d, n, p,
                                         mu_beta_vb, m2_beta, nu,
                                         sig2_inv_vb, scheme)
      tau_vb <- lambda_vb / nu_vb

      # % # update of the variance of the \beta's (inefficient replication of the tau value)
      if(scheme == "noPrec") {
        sig2_beta_vb <- 1 / sweep(sig2_inv_vb * b_vb, 2, (n-1)*tau_vb,`+`)
      } else {
        sig2_beta_vb <- 1 / sweep((n-1) + (sig2_inv_vb * b_vb), 2, tau_vb,`*`)
      }




      if (batch) { # some updates are made batch-wise

        if(loop != "c++") {
          for (j in 1:p) {

            mat_x_m1 <- mat_x_m1 - tcrossprod(X[, j], mu_beta_vb[j, ])

            # % # update of the \mu_beta
            mu_beta_vb[j, ] <- sig2_beta_vb[j,] * (tau_vb *
                                                     crossprod(Y - mat_x_m1, X[, j]))

            mat_x_m1 <- mat_x_m1 + tcrossprod(X[, j], mu_beta_vb[j, ])
          }
        } else {
          coreHorseShoeLoop(X, Y, mat_x_m1, mu_beta_vb, mu_beta_vb, sig2_beta_vb, tau_vb)
        }

        # % # update of the G values
        if(scheme == "noPrec") {
          G_vb <- (1/2)* sig2_inv_vb * m2_beta
        } else {
          G_vb <- (1/2)* sig2_inv_vb * tau_vb * m2_beta
        }

        # % # update of the b values and c values
        b_vb <- 1/(G_vb + c_vb)

        c_vb <- 1 / (b_vb + 1)
        m2_beta <- (mu_beta_vb ^ 2)  +  sig2_beta_vb


      }



    # % # computation of the lower bound
    lb_new <- lower_bound_horseshoeExp(Y, X, d, n, p, sig2_beta_vb, sig2_inv_vb, tau_vb,
                           eta, kappa, lambda, nu, c_vb, b_vb, mat_x_m1, mu_beta_vb,
                           m2_beta, G_vb,  a_inv_vb,  nu_vb, kappa_vb, scheme)
    ELBO <- c(ELBO, lb_new)

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
      create_named_list_(mu_beta_vb, sig2_beta_vb, sig2_inv_vb, tau_vb,c_vb, b_vb, c_vb,
                         eta_vb, kappa_vb, lambda_vb, nu_vb,
                         lambda, nu, a_inv_vb, A, m2_beta, ELBO)
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
                         eta, kappa, lambda, nu, c_vb,  b_vb, mat_x_m1, m1_beta,
                         m2_beta, G_vb, a_inv_vb, nu_vb, kappa_vb, scheme) {

  # update for \tau_{t}
  if(scheme == "noPrec") {
    lambda_vb <- lambda + (n/2)
  } else {
    lambda_vb <- lambda + ((n+p)/2)
  }
 # nu_vb <- update_nu_vb_horseshoeExp(Y, X, mat_x_m1, b_vb, d, n, p, m1_beta,
  #                       m2_beta, nu, sig2_inv_vb, scheme)

  # ESS <-  update_nu_vb_horseshoeExp(Y, X, mat_x_m1, b_vb, d, n, p, m1_beta,
  #                                   m2_beta, nu, sig2_inv_vb, scheme = "noPrec")

  # updates for \sigma^{-2}
  eta_vb <- (p*d+1)/2
  # kappa_vb <- update_kappa_vb_horseshoeExp(a_inv_vb, m2_beta, b_vb, tau_vb, scheme)

  # mean of a^{-1}
  a_inv_vb <- 1 / (sig2_inv_vb + A^{-2})

  # log values for |tau_{t} and |sig^{-2}

  log_tau_vb <- digamma(lambda_vb) - log(nu_vb)
  log_sig2_inv_vb <- digamma(eta_vb) - log(kappa_vb)
  log_b_vb <- digamma(1) - log(G_vb + c_vb)
  log_c_vb <- digamma(1) - log(b_vb + 1)
  log_a_inv <- digamma(1) - log(A^{-2} + sig2_inv_vb)

 # do the computation for the beta's
  L_1 <- sum(-(n/2) * log(2*pi) +
         (n/2)*log_tau_vb - tau_vb * (nu_vb - colSums(m2_beta * b_vb) * sig2_inv_vb / 2 - nu))

  if(scheme == "noPrec") {
    L_2 <- -(1/2)*(p*d)*log(2*pi) + (1/2)*(p*d)*log_sig2_inv_vb +
      (1/2) * sum(log_b_vb)  -
      (1/2) * sig2_inv_vb * sum(b_vb * m2_beta)
  } else {
    L_2 <- -(1/2)*(p*d)*log(2*pi) + (1/2)*(p*d)*log_sig2_inv_vb +
      (1/2)* sum(log_b_vb) + (1/2)*p*sum(log_tau_vb) -
      (1/2) * sig2_inv_vb * sum(tau_vb * colSums(b_vb * m2_beta))
  }
  H_2 <-  (1/2)*(-sum(log(sig2_beta_vb)) - p*d*(log(2*pi)+1))

# do the computation for the tau's
 L_H_3 <- sum((lambda-lambda_vb) * log_tau_vb + lambda_vb * (1 - (nu/nu_vb)) - lgamma(lambda) + lgamma(lambda_vb) +
           lambda*log(nu) - lambda_vb*log(nu_vb))

# do computation for sigma
 L_4 <- -(1/2)*log_sig2_inv_vb - a_inv_vb*sig2_inv_vb - lgamma(1/2) + (1/2)*log_a_inv
 H_4 <- (eta_vb-1)*log_sig2_inv_vb - eta_vb - lgamma(eta_vb) + eta_vb*log(kappa_vb)

# do the computation for a^{-1}
L_5 <- -(1/2)*log_a_inv - (A^{-2}*a_inv_vb) - lgamma(1/2) + (1/2)*log(A^{-2})
H_5 <- -1 + log((A^{-2}) + sig2_inv_vb)

# do the computation for b
L_6 <- -(1/2)*sum(log_b_vb) - sum(c_vb*b_vb) - (p*d)*lgamma(1/2) + (1/2)*sum(log_c_vb)
H_6 <- -p*d + sum(log(G_vb + c_vb))

# do the computation for c
L_7 <- -(1/2)*sum(log_c_vb) - sum(c_vb) - (p*d)*lgamma(1/2)
H_7 <- sum(-1 + log(1 + b_vb))



l <- L_1 + (L_2 - H_2) + L_H_3 + (L_4 - H_4) + (L_5 - H_5) +
     (L_6 - H_6) + (L_7 - H_7)


l2 <- (1/2)*(5*p*d+2) -(n*d/2) * log(2*pi) + ((n+p)/2)*log_tau_vb -
  sum(tau_vb * (nu_vb - colSums(m2_beta * b_vb) * sig2_inv_vb / 2 - nu)) + (1/2)*(p*d-2*eta_vb+1)*log_sig2_inv_vb -
  (1/2) * sig2_inv_vb * sum(tau_vb * colSums(b_vb * m2_beta)) +
  (1/2)*sum(log(sig2_beta_vb)) +
  sum((lambda-lambda_vb) * log_tau_vb + lambda_vb * (1 - (nu/nu_vb)) - lgamma(lambda) + lgamma(lambda_vb) +
        lambda*log(nu) - lambda_vb*log(nu_vb)) -
  a_inv_vb*(A^{-2}+sig2_inv_vb) - (2+2*p*d)*lgamma(1/2) +
  eta_vb + lgamma(eta_vb) - eta_vb*log(kappa_vb) -
   log(A) - log((A^{-2}) + sig2_inv_vb) -
     sum(c_vb*(b_vb+1))  -
     sum(log(G_vb + c_vb)) -
     sum(log(1 + b_vb))


  return(l2)

}

