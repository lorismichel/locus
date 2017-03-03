# Loris Michel for chair of stats EPFL
# implementation of the reduced variational schemes



TwoCompSpike <- function(X, Y, omega, var0, var1, tau, iter.max, error, batch=T) {

  # keep trace of the dimensions
  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(Y)

  # put the noise term in a vector format
  if(length(noise)!=1) {
    tau <- rep(tau,d)
  }

  # definition of the sigma who are fixed there
  Sigma0 <- 1/(tau*(n-1 + var0^{-1}))
  Sigma1 <- 1/(tau*(n-1 + var0^{-1}))



  while((nb.iter <= iter.max) && (err >= error)) {

    # update of the moments
    m1_beta  <- Gamma*Mu1 + (1-Gamma)*Mu0
    m2_beta <- Gamma*sweep(Mu1 ^ 2, 2, Sigma1, `+`) + (1-Gamma)*sweep(Mu0 ^ 2, 2, Sigma0, `+`)
    mat_x_m1 <-  X %*% m1_beta

    if (batch) { # some updates are made batch-wise

      for (j in 1:p) {

        mat_x_m1 <- mat_x_m1 - tcrossprod(X[, j], m1_beta[j, ])

        Mu1[j, ] <- Sigma1 * (tau_vb * crossprod(Y - mat_x_m1, X[, j]))
        Mu0[j, ] <- Sigma0 * (tau_vb * crossprod(Y - mat_x_m1, X[, j]))

        Gamma[j,] <- 1/(1+exp((1/2)*(Mu0[j,]^2)*Sigma0^{-2}+(1/2)*(Mu0[j,]^2)*Sigma0^{-2}+log(1-omega)+log(omega)))

        m1_beta[j, ] <- Mu1[j, ] * Gamma[j, ] + Mu0[j, ] * (1-Gamma[j, ])

        mat_x_m1 <- mat_x_m1 + tcrossprod(X[, j], m1_beta[j, ])
      }
    }
  }



}








    #update the gamma and the mu beta
    for(t in 1:d) {
      for(s in 1:p) {



        # uodate of
        Gamma[s,t] <- 1/(1+exp((1/2)*Mu0[s,t]^{2}Sigma0[s,t]^{-2}-(1/2)*Mu1[s,t]^{2}Sigma1[s,t]^{-2}+log(1-omega)-log(omega)))

        signal <- tau[t]*crossprod(X[s,], tcrossprod(Y[,t], -sum(mu_beta[-s,]*X[-s,])))
        Mu0[s,t] <- Sigma0*signal
        Mu1[s,t] <- Sigma1*signal

      }
    }
  }





}



L.twoCompSpike <- function(Y, X, a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda, nu,
                         sig2_beta_vb, sig2_inv_vb, tau_vb, m1_beta, m2_beta,
                         mat_x_m1, sum_gam) {


  A <- p*d - (n*d/2)*log(2*pi)+((n+pd)/2)*sum(log(tau))-sum(tau*colSums(y))/2
  B <- -sum(tau*tcrossprod(y,mat_x_m1))
  C <- sum(tau*(n-1)*m2_beta)
  D <- (1/2)*sum(Gamma*log(var1^{-1})+(1-Gamma)*log(var0^{-1}))
  E <- -(1/2)*sum(tau*colSums(var1^{-1}*Gamma*(Mu1^{2}+Sigma1)+var0^{-1}*(1-Gamma)*(Mu0^{2}+Sigma0)))
  G <- -(1/2)*sum(Gamma*log(Sigma1^{-1})+(1-Gamma)*log(Sigma0^{-1}))
  H <- sum(Gamma*log(omega)+(1-Gamma)*log(1-omega))
  J <- -sum(Gamma*log(Gamma+eps) + (1-Gamma)*log(1-Gamma+eps))

  A + B + C + D + E + G + H + J

}






TwoCompPrec <-
