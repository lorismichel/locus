require(expint)
locus_core_horseshoe <- function(Y, X, d, n, p, list_hyper, b_vb, sigma2_bv, mu_beta_vb,
                        sig2_beta_vb, tau_vb, tol, maxit, batch, verbose,
                        full_output = FALSE) {

  # Y must have been centered, and X, standardized.

  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects

    # second moment of the \beta's
    m2_beta <- sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`)

    # c_vb <- d_vb <- lambda_vb <- nu_vb <- a_vb <- b_vb <- NULL
    a_inv_vb <- A^{2} / 2
    converged <- FALSE
    lb_old <- -Inf
    it <- 1

    while ((!converged) & (it <= maxit)) {

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % # update of sigma^{-2}
      eta_vb <- (p*d+1)/2
      kappa_vb <- update_kappa_vb_horseshoe(a_inv_vb, m2_beta, b_vb)
      sig2_inv_vb <- eta_vb / kappa_vb

      # % # update of a^{-1}
      a_inv_vb <- 1 / (sig2_inv_vb + A^{-2})

      # % # # update of tau_{t}
      lambda_vb <- lambda + (n/2)
      nu_vb <- update_nu_vb_horseshoe(Y, X, d, n, p, mu_beta_vb, m2_beta, nu)
      tau_vb <- lambda_vb / nu_vb

      # % # update of the variance of the \beta's (inefficient replication of the tau value)
      sig2_beta_vb <- 1 / ((n-1)*matrix(tau_vb,ncol=d, nrow=p, byrow = T) + (sig2_inv_vb * b_vb))


      if (batch) { # some updates are made batch-wise

        mat_x_m1_j <-  X %*% mu_beta_vb

        for (j in 1:p) {

          mat_x_m1_j <- mat_x_m1_j - tcrossprod(X[, j], mu_beta_vb[j, ])

          # % # update of the \mu_beta
          mu_beta_vb[j, ] <- sig2_beta_vb[j,] * (tau_vb *
                                               crossprod(Y - mat_x_m1_j, X[, j]))

          mat_x_m1_j <- mat_x_m1_j + tcrossprod(X[, j], mu_beta_vb[j, ])

          m2_beta <- sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`)

          # % # update of the G values
          G_vb[j,] <- (1/2)*sig2_inv_vb * m2_beta[j,]

          # % # update of the b values
          b_vb[j,] <- (1 / (G_vb[j,] * exp(G_vb[j,]) * expint_E1(G_vb[j,]))) - 1
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
    lb_new <- lower_bound_horseshoe(Y, X, d, n, p, sig2_beta_vb, sig2_inv_vb, tau_vb,
                           eta, kappa, lambda, nu, b_vb,
                           m2_beta, G_vb)

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
      create_named_list_(mu_beta_vb, sig2_beta_vb, sig2_inv_vb, tau_vb, b_vb,
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


update_kappa_vb_horseshoe <- function(a_inv, m2_beta, b_vb) {

  # % # updating the \kappa
  as.numeric(a_inv + (1/2)*sum(m2_beta * b_vb))

}

update_nu_vb_ <- function(Y_mat, X_mat, d, n, p, m1_beta,
                             m2_beta, nu) {
  # put X_mat and Y_mat instead of X and Y to avoid conflicts with the function sapply,
  # which has also an "X" argument with different meaning...

  # X must be standardized as we use (X \hadamard X)^T \one_n = (n-1)\one_p
  nu_vb <- nu + colSums(Y_mat ^ 2) / 2 - colSums(Y_mat * (X_mat %*% m1_beta)) + ((n-1)/2) * colSums(m2_beta)

  # in the case that we observe more than one variable
  if (p > 1) {

    mat_x_list <- lapply(p:1, function(j) {
      X_mat[, j, drop = FALSE] %*% m1_beta[j,, drop = FALSE]
    }
    )

    cum_mat_x_list <- list(0)
    for (j in 1:(p-1)) {
      cum_mat_x_list[[j+1]] <- cum_mat_x_list[[j]] + mat_x_list[[j]]
    }
    cum_mat_x_list[[1]] <- NULL

    mix_x_sum <- sapply(p:2, function(j) {
      colSums(mat_x_list[[j]] * cum_mat_x_list[[j-1]])
    })

    if(d == 1) mix_x_sum <- t(mix_x_sum)


    nu_vb <- nu_vb + rowSums(mix_x_sum)

  }

  nu_vb

}



# this function should be changed adequately
lower_bound_horseshoe <- function(Y, X, d, n, p, sig2_beta_vb, sig2_inv_vb, tau_vb,
                         eta, kappa, lambda, nu,  b_vb,
                         m2_beta, G_vb) {

  # update for \tau_{t}
  lambda_vb <- lambda + (n/2)
  nu_vb <- update_nu_vb_(Y, X, d, n, p, mu_beta_vb, m2_beta, nu)

  # updates for \sigma^{-2}
  eta_vb <- (p*d+1)/2
  kappa_vb <- update_kappa_vb_(a_inv_vb, m2_beta, b_vb)

  # mean of a^{-1}
  a_inv_vb <- 1 / (sig2_inv_vb + A^{-2})

  # log values for |tau_{t} and |sig^{-2}
  log_tau_vb <- digamma(lambda_vb) - log(nu_vb)
  log_sig2_inv_vb <- digamma(eta_vb) - log(kappa_vb)

  # lower bound as in the paper
  l <- (-1)*(n*d/2)* log(2*pi) + p*d/2 + (n/2)*sum(log_tau_vb) - sum(log(sig2_beta_vb)) -
       crossprod(tau_vb,(nu_vb - nu)) + (p*d/2)*log_sig2_inv_vb - (1/2)* sum(b_vb * sig2_inv_vb * m2_beta) +
    sum((lambda_vb - lambda) * log_tau_vb + lambda_vb * (1 - (nu/nu_vb)) - lgamma(lambda) + lgamma(lambda_vb) +
          lambda*log(nu) - lambda_vb*log(nu_vb)) +
    ((1/2)-eta_vb)*log_sig2_inv_vb - a_inv_vb*(sig2_inv_vb + A^{-2}) -2*lgamma(1/2) + (1/2)*log(A^{-2}) +
    1 - log(sig2_inv_vb + A^{-2}) + eta_vb + lgamma(eta_vb) - eta_vb*log(kappa_vb) - (p*d*log(2*pi)) +
    sum(log(expint_E1(G_vb)) + G_vb) + p*d - sum(G_vb*exp(G_vb)*expint_E1(G_vb))

  return(l)

}

