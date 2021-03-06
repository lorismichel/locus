#' @title Horseshoe sparse regression with Cauchy priors on variances
#' @export
#'
horseshoe_core <- function(Y, X, d, n, p, list_hyper, alpha_vb, mu_beta_vb,
                           sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec,
                           full_output = FALSE) {

  # Y must have been centered, and X, standardized.

  with(list_hyper, { # list_init not used with the with() function to avoid
    # copy-on-write for large objects

    # second moment of the \beta's
    m2_beta  <- (mu_beta_vb ^ 2) + sig2_beta_vb
    # regression matrix
    mat_x_m1 <-  X %*% mu_beta_vb

    a_inv_vb <- 1 / (2 * A^{-2})
    b_inv_vb <- 1 / (2 * B^{-2})

    converged <- FALSE
    lb_old <- -Inf
    it <- 1

    # storage of ELBO along the algorithm
    ELBO <- c()

    while ((!converged) & (it <= maxit)) {

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))


      # % # update of sigma^{-2}
      lambda_vb  <- (p*d + 1) / 2
      nu_vb      <- update_nu_vb_horseshoe(b_inv_vb, m2_beta, alpha_vb,
                                           tau_vb, shared_prec)
      sig2_inv_vb <- lambda_vb / nu_vb

      # % # update of b^{-1}
      b_inv_vb <- 1 / (sig2_inv_vb + B^{-2})

      # % # # update of tau_{t}
      if(!shared_prec) {
        eta_vb <- (n+1)/2
      } else {
        eta_vb <- (n+p+1)/2
      }

      kappa_vb <- update_kappa_vb_horseshoe(Y, X,  mat_x_m1, alpha_vb, d, n,
                                            p, mu_beta_vb, m2_beta,
                                            sig2_inv_vb, a_inv_vb, shared_prec)

      tau_vb <- eta_vb / kappa_vb

      a_inv_vb <- 1/(A^{-2} + tau_vb)

      # % # update of the variance of the \beta's (inefficient replication of the tau value)
      if(!shared_prec) {
        sig2_beta_vb <- 1 / sweep(sig2_inv_vb * alpha_vb, 2,
                                  (n-1)*tau_vb,`+`)
      } else {
        sig2_beta_vb <- 1 / sweep((n-1) + (sig2_inv_vb * alpha_vb),
                                  2, tau_vb,`*`)
      }

      coreHorseShoeLoop(X, Y, mat_x_m1, mu_beta_vb, mu_beta_vb, sig2_beta_vb, tau_vb)

      # % # update of the G values
      if(!shared_prec) {
        G_vb <- (1/2) * sig2_inv_vb * m2_beta
      } else {
        G_vb <- (1/2) * sig2_inv_vb * sweep(m2_beta, 2, tau_vb,`*`)
      }

      # % # update of the b values
      alpha_vb <- (G_vb * Q_approx(G_vb))^{-1} - 1

      m2_beta <- (mu_beta_vb ^ 2)  +  sig2_beta_vb





      # % # computation of the lower bound
      lb_new <- lower_bound_horseshoe(Y, X, d, n, p, sig2_beta_vb, sig2_inv_vb, tau_vb,
                                      eta, kappa, lambda, nu, alpha_vb, mat_x_m1, mu_beta_vb,
                                      m2_beta, G_vb, eta_vb, kappa_vb, lambda_vb, nu_vb, a_inv_vb,
                                      b_inv_vb, A, B, shared_prec)

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
      create_named_list_(mu_beta_vb, sig2_beta_vb, m2_beta, sig2_inv_vb, tau_vb, alpha_vb,
                         eta_vb, kappa_vb, lambda_vb, nu_vb,
                         a_inv_vb, b_inv_vb, A, B, ELBO)
    } else {
      names_x <- colnames(X)
      names_y <- colnames(Y)

      rownames(mu_beta_vb) <- names_x
      colnames(mu_beta_vb) <- names_y
      rownames(alpha_vb) <- names_x
      colnames(alpha_vb) <- names_y


      create_named_list_(lb_opt, mu_beta_vb, sig2_beta_vb, tau_vb)
    }
  })

}


update_nu_vb_horseshoe <- function(b_inv_vb, m2_beta, alpha_vb, tau_vb, shared_prec) {

  # % # updating the \kappa
  if(!shared_prec) {
    as.numeric(b_inv_vb + (1/2)*sum(m2_beta * alpha_vb))
  } else {
    as.numeric(b_inv_vb + (1/2)*sum(tau_vb * colSums(m2_beta * alpha_vb)))
  }

}

update_kappa_vb_horseshoe <- function(Y_mat, X_mat, mat_x_m1, alpha_vb, d, n, p, m1_beta,
                                      m2_beta, sig2_inv_vb, a_inv_vb, shared_prec) {
  # put X_mat and Y_mat instead of X and Y to avoid conflicts with the function sapply,
  # which has also an "X" argument with different meaning...

  # X must be standardized as we use (X \hadamard X)^T \one_n = (n-1)\one_p
  if(!shared_prec) {
    kappa_vb <- a_inv_vb + (colSums(Y_mat ^ 2) - 2*colSums(Y_mat * (mat_x_m1)) +
                              +  colSums((n-1) * m2_beta) +
                              +  colSums(mat_x_m1^2) - (n - 1) * colSums(m1_beta^2))/ 2
  } else {

    kappa_vb <- a_inv_vb + (colSums(Y_mat ^ 2) - 2*colSums(Y_mat * (mat_x_m1)) +
                              +  colSums(((n-1) + sig2_inv_vb * alpha_vb) * m2_beta) +
                              +  colSums(mat_x_m1^2) - (n - 1) * colSums(m1_beta^2))/ 2

  }


}


# this function should be changed adequately
lower_bound_horseshoe <- function(Y, X, d, n, p, sig2_beta_vb, sig2_inv_vb, tau_vb,
                                  eta, kappa, lambda, nu,  alpha_vb, mat_x_m1, m1_beta,
                                  m2_beta, G_vb, eta_vb, kappa_vb, lambda_vb, nu_vb, a_inv_vb, b_inv_vb,
                                  A, B, shared_prec) {


  log_tau_vb <- digamma(eta_vb) - log(kappa_vb)
  log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
  log_a_inv <- digamma(1) - log(A^{-2} + tau_vb)
  log_b_inv <- digamma(1) - log(B^{-2} + sig2_inv_vb)


  if(!shared_prec) {
    L_1 <- -(n*d/2) * log(2*pi) + (n/2)*sum(log_tau_vb) -
      sum(tau_vb * (kappa_vb  - a_inv_vb))
  } else {
    L_1 <- -(n*d/2) * log(2*pi) + (n/2)*sum(log_tau_vb) -
      sum(tau_vb * (kappa_vb - colSums(m2_beta * alpha_vb) * sig2_inv_vb / 2 - a_inv_vb))
  }


  if(!shared_prec) {
    L_2 <- -(1/2)*(p*d)*log(2*pi) + (1/2)*(p*d)*log_sig2_inv_vb +
      - (1/2) * sig2_inv_vb * sum(alpha_vb * m2_beta)
  } else {
    L_2 <- -(1/2)*(p*d)*log(2*pi) + (1/2)*(p*d)*log_sig2_inv_vb +
      (1/2)*p*sum(log_tau_vb) -
      (1/2) * sig2_inv_vb * sum(tau_vb * colSums(alpha_vb * m2_beta))
  }

  H_2 <-  (1/2)*(-sum(log(sig2_beta_vb)) - p*d*(log(2*pi)+1))


  L_H_3 <- sum(((1/2)-eta_vb) * log_tau_vb) - sum(tau_vb * (a_inv_vb - kappa_vb)) - d*lgamma(1/2) + d*lgamma(eta_vb) +
    (1/2)*sum(log_a_inv) - sum(eta_vb*log(kappa_vb))

  L_H_4 <- ((1/2)-lambda_vb) * log_sig2_inv_vb - sig2_inv_vb * (b_inv_vb - nu_vb) -
    lgamma(1/2) + lgamma(lambda_vb) + (1/2) * log_b_inv -  lambda_vb*log(nu_vb)


  L_H_5 <- sum(-(1/2)*log_a_inv + a_inv_vb*tau_vb -lgamma(1/2) -
                 log(A^{-2} + tau_vb))

  L_H_6 <- -(1/2)*log_b_inv + b_inv_vb*sig2_inv_vb -lgamma(1/2) -
    log(B^{-2} + sig2_inv_vb)

  L_H_7 <-  -(p*d) * log(pi) + sum(log(Q_approx(G_vb)))  + sum(G_vb * alpha_vb)


  l <- L_1 + (L_2 - H_2) + L_H_3 + L_H_4 + L_H_5 + L_H_6 + L_H_7

  return(l)

}
