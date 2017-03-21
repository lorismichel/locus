locus_student_core_ <- function(Y, X, list_hyper, gam_vb, mu_beta_vb, sig2_beta_vb,
                        tau_vb, b_vb, c_vb, tol, maxit, batch, verbose, full_output = FALSE) {

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  # Y must have been centered, and X, standardized.

  with(list_hyper, { # list_init not used with the with() function to avoid
    # copy-on-write for large objects
    m1_beta <- mu_beta_vb * gam_vb
    m2_beta <- ((mu_beta_vb ^ 2) + sig2_beta_vb) * gam_vb

    mat_x_m1 <-  X %*% m1_beta

    rowsums_gam <- rowSums(gam_vb)
    sum_gam <- sum(rowsums_gam)

    lambda_vb <- nu_vb <- eta_vb <- kappa_vb <- NULL

    converged <- FALSE
    lb_old <- -Inf
    it <- 1

    while ((!converged) & (it <= maxit)) {

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % #
      lambda_vb <- update_lambda_vb_student(lambda, sum_gam)
      nu_vb <- update_nu_vb_student(nu, m2_beta, tau_vb, b_vb)

      sig2_inv_vb <- lambda_vb / nu_vb
      # % #

      # % #
      eta_vb <- update_eta_vb_student(n, eta, gam_vb)

      kappa_vb <- update_kappa_vb_student(Y, X, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb)


      tau_vb <- eta_vb / kappa_vb
      # % #

      sig2_beta_vb <- 1 / (sweep(sig2_inv_vb*b_vb + (n - 1),MARGIN = 2, tau_vb, `*`))
      # % #

      b_shape_vb <- 1 + gam_vb / 2
      b_scale_vb <- gam_vb*(c_vb + (1/2)*(sweep(sig2_inv_vb * m2_beta, 2, tau_vb,`*`)))

      b_vb <- (1/b_scale_vb)
      # % #

      c_shape_vb <- A + sum(gam_vb)
      c_scale_vb <- B + sum(gam_vb*b_vb)

      c_vb <- c_shape_vb / c_scale_vb
      # % #

      log_tau_vb <- digamma(eta_vb) - log(kappa_vb)
      log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
      log_b_vb <- digamma(b_shape_vb) - log(b_scale_vb)
      log_c_vb <- digamma(c_shape_vb) - log(c_scale_vb)

      vec_part_digam <- digamma(alpha + psi + d)

      if (batch) { # some updates are made batch-wise

        log_om_vb <- digamma(alpha + rowsums_gam) - vec_part_digam
        log_1_min_om_vb <- digamma(psi - rowsums_gam + d) - vec_part_digam

        for (j in 1:p) {
          mat_x_m1 <- mat_x_m1 - tcrossprod(X[, j], m1_beta[j, ])

          mu_beta_vb[j, ] <-  sig2_beta_vb[j,] * (tau_vb *
                                                crossprod(Y - mat_x_m1, X[, j]))

          log_gam_part_vb <-  log(sig2_beta_vb[j,]) / 2 +
                    (log_sig2_inv_vb  + log_tau_vb + log_b_vb[j,]) / 2 +
                   (mu_beta_vb[j, ] ^ 2) / (2 * sig2_beta_vb[j,])  +
                    log_om_vb[j] - log_1_min_om_vb[j] +
                   log_c_vb - b_vb[j,]*c_vb

          gam_vb[j, ] <- 1 / (1+exp(-log_gam_part_vb))


          m1_beta[j, ] <- mu_beta_vb[j, ] * gam_vb[j, ]

          mat_x_m1 <- mat_x_m1 + tcrossprod(X[, j], m1_beta[j, ])
        }

        rowsums_gam <- rowSums(gam_vb)

      }

      m2_beta <- ((mu_beta_vb ^ 2) + sig2_beta_vb) * gam_vb


      alpha_vb <- alpha + rowsums_gam
      psi_vb <- psi - rowsums_gam + d
      om_vb <- alpha_vb / (alpha_vb + psi_vb)

      sum_gam <- sum(rowsums_gam)


    #  lb_new <- lower_bound_student(Y, X, alpha, alpha_vb, psi, psi_vb, eta, gam_vb, kappa, lambda,
     #                        nu, sig2_beta_vb, sig2_inv_vb, tau_vb, b_vb, c_vb, m1_beta,
      #                       m2_beta, mat_x_m1, sum_gam, A, B)
#
      lb_new <- 0


      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Lower bound = ", format(lb_new), "\n\n", sep = ""))

    #  converged <- (abs(lb_new - lb_old) < tol)
       converged <- F
      lb_old <- lb_new
      it <- it + 1
    }



    if (verbose) {
      if (converged) {
        cat(paste("Convergence obtained after ", format(it),
                  " iterations with variational lower bound = ",
                  format(lb_new), ". \n\n", sep = ""))
      } else {
        cat("Maximal number of iterations reached before convergence. Exit.")
      }
    }

    lb_opt <- lb_new

    if (full_output) { # for internal use only

      create_named_list_(a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda,
                         nu, sig2_beta_vb, sig2_inv_vb, tau_vb, m1_beta,
                         m2_beta, mat_x_m1, sum_gam)

    } else {
      names_x <- colnames(X)
      names_y <- colnames(Y)

      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      names(om_vb) <- names_x

      create_named_list_(lb_opt, gam_vb, om_vb)
    }
  })

}


update_lambda_vb_student <- function(lambda, sum_gam) {

  lambda + sum_gam / 2

}

update_nu_vb_student <- function(nu, m2_beta, tau_vb, b_vb) {

  as.numeric(nu + crossprod(tau_vb, colSums(b_vb * m2_beta)) / 2)

}

update_eta_vb_student <- function(n, eta, gam_vb) {

  eta + n / 2 + colSums(gam_vb) / 2

}

update_kappa_vb_student <- function(Y, X, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb) {
  n <- nrow(Y)

  kappa + (colSums(Y^2) - 2 * colSums(Y * mat_x_m1)  +
             (n - 1) * colSums(m2_beta) + sig2_inv_vb * colSums(b_vb * m2_beta) +
             colSums(mat_x_m1^2) - (n - 1) * colSums(m1_beta^2))/ 2

}


lower_bound_student <- function(Y, X, alpha, alpha_vb, psi, psi_vb, eta, gam_vb, kappa, lambda, nu,
                         sig2_beta_vb, sig2_inv_vb, tau_vb, b_vb, c_vb, m1_beta, m2_beta,
                         mat_x_m1, sum_gam, A, B) {


  n <- nrow(Y)

  eta_vb <- update_eta_vb_student(n, eta, gam_vb)
  kappa_vb <- update_kappa_vb_student(Y, X, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb)

  lambda_vb <- update_lambda_vb_student(lambda, sum_gam)
  nu_vb <- update_nu_vb_student(nu, m2_beta, tau_vb, b_vb)

  b_shape_vb <- 1 + gam_vb/2
  b_scale_vb <- gam_vb*(c_vb + (1/2)*(sweep(sig2_inv_vb * m2_beta, 2, tau_vb,`*`)))

  b_vb <- b_shape_vb / b_scale_vb
  # % #

  c_shape_vb <- A + sum_gam
  c_scale_vb <- B + sum(b_vb)

  c_vb <- c_shape_vb / c_scale_vb
  # % #

  log_b_vb <- digamma(b_shape_vb) - log(b_scale_vb)
  log_c_vb <- digamma(c_shape_vb) - log(c_scale_vb)

  alpha_vb <- alpha + rowsums_gam
  psi_vb <- psi - rowsums_gam + d
  om_vb <- alpha_vb / (alpha_vb + psi_vb)

  log_tau_vb <- digamma(eta_vb) - log(kappa_vb)
  log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
  log_om_vb <- digamma(alpha_vb) - digamma(alpha + psi + d)
  log_1_min_om_vb <- digamma(psi_vb) - digamma(alpha + psi + d)

  L.1 <- sum(-n / 2 * log(2 * pi) + n / 2 * log_tau_vb -
             tau_vb * (kappa_vb - colSums(m2_beta) * sig2_inv_vb / 2 - kappa))
  H.1 <- sum((1/2)*gam_vb*(-log(2*pi)-log(sig2_beta_vb))-(1/2)*gam_vb)

  eps <- .Machine$double.eps # to control the argument of the log when gamma is very small
  L.2 <- sum(  -(p*d/2)*log(2*pi) +
             log_sig2_inv_vb * gam_vb / 2 +
             sweep(gam_vb, 2, log_tau_vb, `*`) / 2 +
             log_b_vb / 2 -
             sweep(b_vb * m2_beta, 2, tau_vb, `*`) * sig2_inv_vb / 2 )
            #+ sweep(gam_vb, 1, log_om_vb, `*`) +
            # sweep(1 - gam_vb, 1, log_1_min_om_vb, `*`) +
            # 1 / 2 * sweep(gam_vb, 2, log(sig2_beta_vb) + 1, `*`) -
            # gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))
  L.3 <- sum(gam_vb*(log_c_vb - c_vb*b_vb))

  H.3 <- sum( b_shape_vb*log_b_vb + b_shape_vb -lgamma(b_shape_vb) +
                b_shape_vb*log(b_scale_vb))

  L.4_H.2 <- sum((eta - eta_vb) * log_tau_vb -
                   (kappa - kappa_vb) * tau_vb + eta * log(kappa) -
                   eta_vb * log(kappa_vb) - lgamma(eta) + lgamma(eta_vb))

  L.5_H.4 <- sum((lambda - lambda_vb) * log_sig2_inv_vb -
             (nu - nu_vb) * sig2_inv_vb + lambda * log(nu) -
             lambda_vb * log(nu_vb) - lgamma(lambda) + lgamma(lambda_vb))

  H.5 <- sum( c_shape_vb*log_c_vb + c_shape_vb -lgamma(c_shape_vb) +
                c_shape_vb*log(c_scale_vb))

  L.6 <- sum(sweep(gam_vb, 1, log_om_vb, `*`) +
         sweep(1 - gam_vb, 1, log_1_min_om_vb, `*`))

  H.6 <- sum(gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

  L.7_H.7 <- sum((alpha - alpha_vb) * log_om_vb + (psi - psi_vb) * log_1_min_om_vb - lbeta(alpha, psi) +
             lbeta(alpha_vb, psi_vb))


  L.1 + L.2 + L.3 + L.6 + L.4_H.2 + L.5_H.4 +
  L.7_H.7 - H.1 - H.3 - H.5 -H.6
}
