#' Gather model hyperparameters provided by the user.
#'
#' This function must be used to provide hyperparameter values for the model
#' used in \code{\link{locus}}.
#'
#' The \code{\link{locus}} function can also be used with default
#' hyperparameter choices (without using \code{\link{set_hyper}}) by
#' setting its argument \code{list_hyper} to \code{NULL}.
#'
#' @param d Number of responses.
#' @param p Number of candidate predictors.
#' @param lambda Vector of length 1 or d providing the values of hyperparameter
#'   \eqn{\eta} for the prior distributions of the response residual precisions,
#'   \eqn{\tau} (vector of size d). If of length 1, the provided value is
#'   repeated d times.
#' @param nu Vector of length 1 or d providing the values of hyperparameter
#'   \eqn{\kappa} for the prior distributions of the response residual
#'   precisions, \eqn{\tau} (vector of size d). If of length 1, the provided
#'   value is repeated d times.
#' @param eta Value of hyperparameter \eqn{\lambda} for the prior
#'   distribution of \eqn{\sigma^{-2}}. \eqn{\sigma^2} represents the typical
#'   size of nonzero effects.
#' @param kappa Value of hyperparameter \eqn{\nu} for the prior distribution of
#'   \eqn{\sigma^{-2}}. \eqn{\sigma^2} represents the typical size of nonzero
#'   effects.
#' @param q Number of covariates. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#' @param phi Vector of length 1 or q providing the values of hyperparameter
#'   \eqn{\phi} for the prior distributions for the sizes of the nonzero
#'   covariate effects, \eqn{\zeta} (vector of size q). If of length 1, the
#'   provided value is repeated q times. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#' @param xi Vector of length 1 or q providing the values of hyperparameter
#'   \eqn{\xi} for the prior distributions for the sizes of the nonzero
#'   covariate effects, \eqn{\zeta} (vector of size q). If of length 1, the
#'   provided value is repeated q times. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#'
#' @return An object of class "\code{hyper}" preparing user hyperparameter in a
#'   form that can be passed to the \code{\link{locus}} function.
#'
#' @examples
#' user_seed <- 123
#' n <- 200; p <- 400; p0 <- 100; d <- 25; d0 <- 20
#' list_X <- generate_snps(n = n, p = p, user_seed = user_seed)
#' list_Y <- generate_phenos(n = n, d = d, var_err = 0.25, user_seed = user_seed)
#'
#' dat <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
#'                            ind_d0 = sample(1:d, d0), ind_p0 = sample(1:p, p0),
#'                            vec_prob_sh = 0.1, max_tot_pve = 0.9,
#'                            user_seed = user_seed)
#'
#' # a and b chosen so that each candidate predictor has a prior probability to
#' # be included in the model of 1/4.
#' list_hyper <- set_hyper(d, p, eta = 1, kappa = apply(dat$phenos, 2, var),
#'                               lambda = 1, nu = 1, a = 1, b = 4*d-1)
#'
#' vb <- locus(Y = dat$phenos, X = dat$snps, p0_av = p0, list_hyper = list_hyper,
#'             user_seed = user_seed)
#'
#' @seealso  \code{\link{set_init}}, \code{\link{locus}}
#'
#' @export
#'
set_hyper <- function(d, p, lambda, nu, eta, kappa,
                      q = NULL, phi = NULL, xi = NULL) {

  check_structure_(lambda, "vector", "double", c(1, d))
  check_positive_(lambda)
  if (length(eta) == 1) eta <- rep(eta, d)

  check_structure_(nu, "vector", "double", c(1,d))
  check_positive_(nu)
  if (length(kappa) == 1) kappa <- rep(kappa, d)

  check_structure_(eta, "vector", "double", 1)
  check_positive_(eta)

  check_structure_(kappa, "vector", "double", 1)
  check_positive_(kappa)

  check_structure_(a, "vector", "double", c(1, p))
  check_positive_(a)
  if (length(a) == 1) a <- rep(a, p)

  check_structure_(b, "vector", "double", c(1, p))
  check_positive_(b)
  if (length(b) == 1) b <- rep(b, p)

  if (!is.null(q)) {

    check_structure_(phi, "vector", "double", c(1, q))
    check_positive_(phi)
    if (length(phi) == 1) phi <- rep(phi, q)

    check_structure_(xi, "vector", "double", c(1, q))
    check_positive_(xi)
    if (length(xi) == 1) xi <- rep(xi, q)

  } else if (!is.null(phi) | !is.null(xi)) {
    stop("Provided q = NULL, not consitent with phi or xi being non-null.")
  }

  d_hyper <- d
  p_hyper <- p
  q_hyper <- q

  list_hyper <- create_named_list_(d_hyper, p_hyper, q_hyper, lambda, nu, eta,
                                   kappa, phi, xi)

  class(list_hyper) <- "hyper"

  list_hyper

}


auto_set_hyper_ <- function(Y, p, p_star, q = NULL) {

  d <- ncol(Y)

  # hyperparameter set using the data Y
  lambda <- 1 / median(apply(Y, 2, var)) #median to be consistent when doing permutations
  if (!is.finite(eta)) lambda <- 1e3
  lambda <- rep(lambda, d)
  nu <- rep(1, d)

  eta <- 1e-2
  kappa <- 1

  if (!is.null(q)) {
    phi <- rep(1, q)
    xi <- rep(1, q)
  } else {
    phi <- NULL
    xi <- NULL
  }

  d_hyper <- d
  p_hyper <- p
  q_hyper <- q

  list_hyper <- create_named_list_(d_hyper, p_hyper, q_hyper, eta, kappa, lambda,
                                   nu, a, b, phi, xi)

  class(list_hyper) <- "out_hyper"

  list_hyper

}

#' Gather initial variational parameters provided by the user.
#'
#' This function must be used to provide initial values for the variational
#' parameters used in \code{\link{locus}}.
#'
#' The \code{\link{locus}} function can also be used with default initial
#' parameter choices (without using \code{\link{set_init}}) by setting
#' its argument \code{list_init} to \code{NULL}.
#'
#' @param d Number of responses.
#' @param p Number of candidate predictors.
#' @param gam_vb Matrix of size p x d with initial values for the variational
#'   parameter yielding posterior probabilities of inclusion.
#' @param mu_beta_vb Matrix of size p x d with initial values for the variational
#'   parameter yielding regression coefficient estimates for predictor-response
#'   pairs included in the model.
#' @param sig2_beta_vb Vector of size d with initial values for the variational
#'   parameter yielding estimates of effect variances for predictor-response
#'   pairs included in the model. These values are the same for all predictors
#'   (as a result of the predictor variables being standardized before the
#'   variational algorithm).
#' @param tau_vb  Vector of size d with initial values for the variational
#'   parameter yielding estimates for the response residual precisions.
#' @param q Number of covariates. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#' @param mu_alpha_vb Matrix of size p x q with initial values for the
#'   variational parameter yielding regression coefficient estimates for
#'   covariate-response pairs. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#' @param sig2_alpha_vb Matrix of size p x q with initial values for the
#'   variational parameter yielding estimates of effect variances for
#'   covariate-response pairs. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#'
#' @return An object of class "\code{init}" preparing user initial values for
#'   the variational parameters in a form that can be passed to the
#'   \code{\link{locus}} function.
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 200; p <- 400; p0 <- 100; d <- 25; d0 <- 20
#' list_X <- generate_snps(n = n, p = p)
#' list_Y <- generate_phenos(n = n, d = d, var_err = 0.25)
#'
#' dat <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
#'                            ind_d0 = sample(1:d, d0), ind_p0 = sample(1:p, p0),
#'                            vec_prob_sh = 0.1, max_tot_pve = 0.9)
#'
#' # gam_vb chosen so that each candidate predictor has a prior probability to
#' # be included in the model of 1/4.
#' gam_vb <- matrix(rbeta(p * d, shape1 = 1, shape2 = 4*d-1), nrow = p)
#' mu_beta_vb <- matrix(rnorm(p * d), nrow = p)
#' tau_vb <- 1 / apply(dat$phenos, 2, var)
#' sig2_beta_vb <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb)
#'
#' list_init <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb)
#'
#' vb <- locus(Y = dat$phenos, X = dat$snps, p0_av = p0, list_init = list_init)
#'
#' @seealso  \code{\link{set_hyper}}, \code{\link{locus}}
#'
#' @export
#'
set_init <- function(d, p, b_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
                     q = NULL, mu_alpha_vb = NULL, sig2_alpha_vb = NULL) {

  check_structure_(b_vb, "matrix", "double", c(p, d))
  check_positive_(b_vb)

  check_structure_(mu_beta_vb, "matrix", "double", c(p, d))

  check_structure_(sig2_beta_vb, "vector", "double", d)
  check_positive_(sig2_beta_vb)

  check_structure_(tau_vb, "vector", "double", d)
  check_positive_(tau_vb)


  if (!is.null(q)) {

    check_structure_(mu_alpha_vb, "matrix", "double", c(q, d))

    check_structure_(sig2_alpha_vb, "matrix", "double", c(q, d))
    check_positive_(sig2_alpha_vb)

  } else if (!is.null(mu_alpha_vb) | !is.null(sig2_alpha_vb)) {
    stop(paste("Provided q = NULL, not consistent with mu_alpha_vb or ",
               "sig2_alpha_vb being non-null.", sep = ""))
  }

  d_init <- d
  p_init <- p
  q_init <- q

  list_init <- create_named_list_(d_init, p_init, q_init, gam_vb, mu_beta_vb,
                                  sig2_beta_vb, tau_vb, mu_alpha_vb, sig2_alpha_vb)

  class(list_init) <- "init"

  list_init
}


auto_set_init_ <- function(Y, p, p_star, user_seed, q = NULL) {

  d <- ncol(Y)

  if (!is.null(user_seed)) set.seed(user_seed)




  b_vb <- matrix(rexp(p * d), nrow = p)
  mu_beta_vb <- matrix(rnorm(p * d), nrow = p)


  tau_vb <- 1 / median(apply(Y, 2, var))
  if (!is.finite(tau_vb)) tau_vb <- 1e3
  tau_vb <- rep(tau_vb, d)

  sig2_inv_vb <- 1e-2
  sig2_beta_vb <- 1 / rgamma(d, shape = 2, rate = 1 / (sig2_inv_vb * b_vb + tau_vb))


  if (!is.null(q)) {
    mu_alpha_vb <- matrix(rnorm(q * d), nrow = q)
    zeta2_inv_vb <- rgamma(q, shape = 1, rate = 1)
    sig2_alpha_vb <- 1 / sapply(tau_vb,
                                function(tau_vb_t) {
                                  rgamma(q, shape = 2,
                                         rate = 1 / (zeta2_inv_vb * tau_vb_t))
                                } )
  } else {
    mu_alpha_vb <- NULL
    sig2_alpha_vb <- NULL
  }


  d_init <- d
  p_init <- p
  q_init <- q

  list_init <- create_named_list_(d_init, p_init, q_init, gam_vb, mu_beta_vb,
                                  sig2_beta_vb, tau_vb, mu_alpha_vb, sig2_alpha_vb)

  class(list_init) <- "out_init"

  list_init
}
