% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_hyper_init.R
\name{set_hyper}
\alias{set_hyper}
\title{Gather model hyperparameters provided by the user.}
\usage{
set_hyper(d, p, lambda, nu, a, b, eta, kappa, link = "identity",
  ind_bin = NULL, q = NULL, phi = NULL, xi = NULL, r = NULL,
  m0 = NULL)
}
\arguments{
\item{d}{Number of responses.}

\item{p}{Number of candidate predictors.}

\item{lambda}{Vector of length 1 providing the values of hyperparameter
\eqn{\lambda} for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma}
represents the typical size of nonzero effects.}

\item{nu}{Vector of length 1 providing the values of hyperparameter \eqn{\nu}
for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma} represents
the typical size of nonzero effects.}

\item{a}{Vector of length 1 or p providing the values of hyperparameter
\eqn{a} for the prior distributions for the proportion of responses
associated with each candidate predictor, \eqn{\omega} (vector of length p).
If of length 1, the provided value is repeated p times.}

\item{b}{Vector of length 1 or p providing the values of hyperparameter
\eqn{b} for the prior distributions for the proportion of responses
associated with each candidate predictor, \eqn{\omega} (vector of length p).
If of length 1, the provided value is repeated p times.}

\item{eta}{Vector of length 1 or d for \code{link = "identity"}, and of
length 1 or d_cont = d - length(ind_bin) (the number of continuous response
variables) for \code{link = "mix"}. Provides the values of
hyperparameter \eqn{\eta} for the prior distributions of the continuous
response residual precisions, \eqn{\tau}. If of length 1, the provided
value is repeated d, resp. d_cont, times. Must be \code{NULL} for
\code{link = "logit"} and \code{link = "probit"}.}

\item{kappa}{Vector of length 1 or d for \code{link = "identity"}, and of
length 1 or d_cont = d - length(ind_bin) (the number of continuous response
variables) for \code{link = "mix"}. Provides the values of hyperparameter
\eqn{\kappa} for the prior distributions of the response residual
precisions, \eqn{\tau}. If of length 1, the provided value is repeated d,
resp. d_cont, times. Must be \code{NULL} for \code{link = "logit"} and
\code{link = "probit"}.}

\item{link}{Response link. Must be "\code{identity}" for linear regression,
"\code{logit}" for logistic regression, "\code{probit}"
for probit regression, or "\code{mix}" for a mix of identity and probit
link functions (in this case, the indices of the binary responses must be
gathered in argument \code{ind_bin}, see below).}

\item{ind_bin}{If \code{link = "mix"}, vector of indices corresponding to
the binary variables in \code{Y}. Must be \code{NULL} if
\code{link != "mix"}.}

\item{q}{Number of covariates. Default is \code{NULL}, for \code{Z}
\code{NULL}.}

\item{phi}{Vector of length 1 or q providing the values of hyperparameter
\eqn{\phi} for the prior distributions for the sizes of the nonzero
covariate effects, \eqn{\zeta}. If of length 1, the provided value is
repeated q times. Default is \code{NULL}, for \code{Z} \code{NULL}.}

\item{xi}{Vector of length 1 or q providing the values of hyperparameter
\eqn{\xi} for the prior distributions for the sizes of the nonzero
covariate effects, \eqn{\zeta}. If of length 1, the provided value is
repeated q times. Default is \code{NULL}, for \code{Z} \code{NULL}.}

\item{r}{Number of variables representing external information on the
candidate predictors. Default is \code{NULL}, for \code{V} \code{NULL}.}

\item{m0}{Vector of length 1 or p providing the values of hyperparameter
\eqn{m0} for the prior distribution of \eqn{c0} linked to the proportion of
responses associated with each candidate predictor. If of length 1, the
provided value is repeated p times. Default is \code{NULL}, for \code{V}
\code{NULL}.}
}
\value{
An object of class "\code{hyper}" preparing user hyperparameter in a
  form that can be passed to the \code{\link{locus}} function.
}
\description{
This function must be used to provide hyperparameter values for the model
used in \code{\link{locus}}.
}
\details{
The \code{\link{locus}} function can also be used with default hyperparameter
choices (without using \code{\link{set_hyper}}) by setting its argument
\code{list_hyper} to \code{NULL}.
}
\examples{
user_seed <- 123
n <- 200; p <- 250; p0 <- 50; d <- 25; d0 <- 15
list_X <- generate_snps(n = n, p = p, user_seed = user_seed)
list_Y <- generate_phenos(n = n, d = d, var_err = 1, user_seed = user_seed)

# Continuous outcomes
#
dat_g <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
                             ind_d0 = sample(1:d, d0),
                             ind_p0 = sample(1:p, p0), vec_prob_sh = 0.1,
                             family = "gaussian", max_tot_pve = 0.9,
                             user_seed = user_seed)

# a and b chosen so that the prior mean number of responses associated with
# each candidate predictor is 1/4.
list_hyper_g <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
                          eta = 1, kappa = apply(dat_g$phenos, 2, var),
                          link = "identity")

# we take p0_av = p0 (known here); this choice may result in variable
# selections that are (too) conservative in some cases. In practice, often
# p0_av as a slightly overestimated guess of p0.
vb_g <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0,
              link = "identity", list_hyper = list_hyper_g,
              user_seed = user_seed)

# Continuous outcomes with covariates
#
q <- 4
Z <- matrix(rnorm(n * q), nrow = n)

list_hyper_g_z <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
                            eta = 1, kappa = apply(dat_g$phenos, 2, var),
                            link = "identity", q = q, phi = 1, xi = 1)

vb_g_z <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0, Z = Z,
                link = "identity", list_hyper = list_hyper_g_z,
                user_seed = user_seed)

# Continuous outcomes with external annotation
#
r <- 4
V <- matrix(rnorm(p * r), nrow = p)
bool_p0 <- rowSums(dat_g$pat) > 0
V[bool_p0, ] <- rnorm(sum(bool_p0) * r, mean = 2) # informative annotations

list_hyper_g_v <- set_hyper(d, p, lambda = 1, nu = 1, a = NULL, b = NULL,
                            eta = 1, kappa = apply(dat_g$phenos, 2, var),
                            link = "identity", r = r, m0 = 0)

vb_g_v <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0,  V = V,
                link = "identity", list_hyper = list_hyper_g_v,
                user_seed = user_seed)

# Binary outcomes
#
dat_b <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
                             ind_d0 = sample(1:d, d0),
                             ind_p0 = sample(1:p, p0), vec_prob_sh = 0.1,
                             family = "binomial", max_tot_pve = 0.9,
                             user_seed = user_seed)

list_hyper_logit <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
                              eta = NULL, kappa = NULL,
                              link = "logit")

vb_logit <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0,
                  link = "logit", list_hyper = list_hyper_logit,
                  user_seed = user_seed)

list_hyper_probit <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
                               eta = NULL, kappa = NULL,
                               link = "probit")

vb_probit <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0,
                   link = "probit", list_hyper = list_hyper_probit,
                   user_seed = user_seed)

# Binary outcomes with covariates
#
list_hyper_logit_z <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
                                eta = NULL, kappa = NULL,
                                link = "logit", q = q, phi = 1,
                                xi = 1)

vb_logit_z <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0, Z = Z,
                    link = "logit",
                    list_hyper = list_hyper_logit_z, user_seed = user_seed)

list_hyper_probit_z <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
                                 eta = NULL, kappa = NULL,
                                 link = "probit", q = q, phi = 1,
                                 xi = 1)

vb_probit_z <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0, Z = Z,
                     link = "probit",
                     list_hyper = list_hyper_probit_z, user_seed = user_seed)

# Mix of continuous and binary outcomes
#
Y_mix <- cbind(dat_g$phenos, dat_b$phenos)
ind_bin <- (d+1):(2*d)
p0_mix <- sum(rowSums(cbind(dat_g$pat, dat_b$pat)) > 0)

list_hyper_mix <- set_hyper(2*d, p, lambda = 1, nu = 1, a = 1, b = 8*d-1,
                            eta = 1, kappa = apply(dat_g$phenos, 2, var),
                            link = "mix", ind_bin = ind_bin)

vb_mix <- locus(Y = Y_mix, X = dat_b$snps, p0_av = p0_mix, link = "mix",
                ind_bin = ind_bin, list_hyper = list_hyper_mix,
                user_seed = user_seed)

list_hyper_mix_z <- set_hyper(2*d, p, lambda = 1, nu = 1, a = 1, b = 8*d-1,
                              eta = 1, kappa = apply(dat_g$phenos, 2, var),
                              link = "mix", ind_bin = ind_bin, q = q,
                              phi = 1, xi = 1)

vb_mix_z <- locus(Y = Y_mix, X = dat_b$snps, p0_av = p0_mix, Z = Z,
                  link = "mix", ind_bin = ind_bin,
                  list_hyper = list_hyper_mix_z, user_seed = user_seed)

}
\seealso{
\code{\link{set_init}}, \code{\link{locus}}
}
