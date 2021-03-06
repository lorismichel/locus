% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/locus.R
\name{locus}
\alias{locus}
\title{Fit sparse multivariate regression models using variational inference.}
\usage{
locus(Y, X, p0_av, Z = NULL, V = NULL, link = "identity",
  ind_bin = NULL, list_hyper = NULL, list_init = NULL, list_cv = NULL,
  list_blocks = NULL, user_seed = NULL, tol = 1e-04, maxit = 1000,
  batch = TRUE, save_hyper = FALSE, save_init = FALSE, verbose = TRUE)
}
\arguments{
\item{Y}{Response data matrix of dimension n x d, where n is the number of
observations and d is the number of response variables.}

\item{X}{Input matrix of dimension n x p, where p is the number of candidate
predictors. \code{X} cannot contain NAs. No intercept must be supplied.}

\item{p0_av}{Prior average number of predictors expected to be included in
the model. Must be \code{NULL} if \code{list_init} and \code{list_hyper}
are both non-\code{NULL} or if \code{list_cv} is non-\code{NULL}. Can also
be a vector of length p with entry s corresponding to the prior probability
that candidate predictor s is associated with at least one response.}

\item{Z}{Covariate matrix of dimension n x q, where q is the number of
covariates. Variables in \code{Z} are not subject to selection. \code{NULL}
if no covariate. Factor covariates must be supplied after transformation to
dummy coding. No intercept must be supplied.}

\item{V}{Annotation matrix of dimension p x r, where r is the number of
variables representing external information on the candidate predictors
which may make their selection more or less likely. \code{NULL} if no such
information.}

\item{link}{Response link. Must be "\code{identity}" for linear regression,
"\code{logit}" for logistic regression, "\code{probit}" for probit
regression, or "\code{mix}" for a mix of identity and probit link functions
(in this case, the indices of the binary responses must be gathered in
argument \code{ind_bin}, see below).}

\item{ind_bin}{If \code{link = "mix"}, vector of indices corresponding to
the binary variables in \code{Y}. Must be \code{NULL} if
\code{link != "mix"}.}

\item{list_hyper}{An object of class "\code{hyper}" containing the model
hyperparameters. Must be filled using the \code{\link{set_hyper}}
function or must be \code{NULL} for default hyperparameters.}

\item{list_init}{An object of class "\code{init}" containing the initial
variational parameters. Must be filled using the \code{\link{set_init}}
function or be \code{NULL} for a default initialization.}

\item{list_cv}{An object of class "\code{cv}" containing settings for
choosing the prior average number of predictors expected to be included in
the model, \code{p0_av}, by cross-validation. Must be filled using the
\code{\link{set_cv}} function or must be \code{NULL} for no
cross-validation. If non-\code{NULL}, \code{p0_av}, \code{list_init} and
\code{list_hyper} must all be \code{NULL}. Cross-validation only available
for \code{link = "identity"}.}

\item{list_blocks}{An object of class "\code{blocks}" containing settings for
parallel inference on a partitioned predictor space. Must be filled using
the \code{\link{set_blocks}} function or must be \code{NULL} for no
partitioning.}

\item{user_seed}{Seed set for reproducible default choices of hyperparameters
(if \code{list_hyper} is \code{NULL}) and initial variational parameters
(if \code{list_init} is \code{NULL}). Also used at the cross-validation
stage (if \code{list_cv} is non-\code{NULL}). Default is \code{NULL}, no
seed set.}

\item{tol}{Tolerance for the stopping criterion.}

\item{maxit}{Maximum number of iterations allowed.}

\item{batch}{If \code{TRUE} a fast batch updating scheme is used
(recommended).}

\item{save_hyper}{If \code{TRUE}, the hyperparameters used for the model are
saved as output.}

\item{save_init}{If \code{TRUE}, the initial variational parameters used for
the inference are saved as output.}

\item{verbose}{If \code{TRUE}, messages are displayed during execution.}
}
\value{
An object of class "\code{vb}" containing the following variational
  estimates and settings:
 \item{lb_opt}{Optimized variational lower bound for the marginal
               log-likelihood.}
 \item{gam_vb}{Posterior inclusion probability matrix of dimension p x d.
               Entry (s, t) corresponds to the posterior probability of
               association between candidate predictor s and response t.}
 \item{mu_alpha_vb}{Matrix of dimension q x d whose entries are the posterior
                    mean regression coefficients for the covariates provided
                    in \code{Z} (if \code{link = "logit"},
                    \code{link = "logit"} or
                    \code{link = "mix"} also for the intercept).
                    \code{NULL} if \code{Z} is \code{NULL}.}
 \item{om_vb}{Vector of length p containing the posterior mean of omega.
              Entry s controls the proportion of responses associated with
              candidate predictor s.}
 \item{rmvd_cst_x, rmvd_cst_z}{Vectors containing the indices of constant
                               variables in \code{X} (resp. \code{Z}) removed
                               prior to the analysis.}
 \item{rmvd_coll_x, rmvd_coll_z}{Vectors containing the indices of variables
                                 in \code{X} (resp. \code{Z}) removed prior
                                 to the analysis because collinear to other
                                 variables. The entry name indicates the
                                 corresponding variable kept in the analysis
                                 (i.e., that causing the collinearity for the
                                 entry in question).}
 \item{list_hyper, list_init}{If \code{save_hyper}, resp. \code{save_init},
                              \code{TRUE}, hyperparameters, resp. initial
                              variational parameters, used for inference are
                              saved as output.}
}
\description{
Variational approximation procedure fitting sparse multivariate regression
models for combined selection of predictors and associated responses in
high-dimensional set-ups. Dependence across responses linked to the same
predictors is captured through the model hierarchical structure.
The responses can be purely continuous, purely binary (logit or probit link
fits), or a mix of continuous and binary variables.
}
\details{
The continuous response variables in \code{Y} (if any) will be centered
before application of the variational algorithm, and the candidate predictors
and covariates resp. in \code{X} and \code{Z} will be standardized. An
intercept will be added if \code{link} is \code{"logit"}, \code{"probit"} or
\code{"mix"} (do not supply it in \code{X} or \code{Z}).
}
\examples{
user_seed <- 123; set.seed(user_seed)
n <- 200; p <- 300; p0 <- 50; d <- 40; d0 <- 30
list_X <- generate_snps(n = n, p = p)
list_Y <- generate_phenos(n = n, d = d, var_err = 1)

# Continuous outcomes
#
dat_g <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
                             ind_d0 = sample(1:d, d0),
                             ind_p0 = sample(1:p, p0),
                             vec_prob_sh = 0.1, family = "gaussian",
                             max_tot_pve = 0.9)

# we take p0_av = p0 (known here); this choice may result in variable
# selections that are (too) conservative in some cases. In practice, often
# p0_av as a slightly overestimated guess of p0.
vb_g <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0,
              link = "identity", user_seed = user_seed)

# Continuous outcomes with covariates
#
q <- 4
Z <- matrix(rnorm(n * q), nrow = n)
vb_g_z <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0,  Z = Z,
                link = "identity", user_seed = user_seed)

# Continuous outcomes with external annotation
#
r <- 4
V <- matrix(rnorm(p * r), nrow = p)
bool_p0 <- rowSums(dat_g$pat) > 0
V[bool_p0, ] <- rnorm(sum(bool_p0) * r, mean = 2) # informative annotations
vb_g_v <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0,  V = V,
                link = "identity", user_seed = user_seed)

# Binary outcomes
#
dat_b <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
                             ind_d0 = sample(1:d, d0),
                             ind_p0 = sample(1:p, p0),
                             vec_prob_sh = 0.1, family = "binomial",
                             max_tot_pve = 0.9)

vb_logit <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0,
                  link = "logit", user_seed = user_seed)

vb_probit <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0,
                   link = "probit", user_seed = user_seed)

# Binary outcomes with covariates
#
vb_logit_z <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0,  Z = Z,
                    link = "logit", user_seed = user_seed)

vb_probit_z <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0,  Z = Z,
                     link = "probit", user_seed = user_seed)

# Mix of continuous and binary outcomes
#
Y_mix <- cbind(dat_g$phenos, dat_b$phenos)
ind_bin <- (d+1):(2*d)
p0_mix <- sum(rowSums(cbind(dat_g$pat, dat_b$pat)) > 0)

vb_mix <- locus(Y = Y_mix, X = dat_b$snps, p0_av = p0, link = "mix",
                ind_bin = ind_bin, user_seed = user_seed)

# Mix of continuous and binary outcomes with covariates
#
vb_mix_z <- locus(Y = Y_mix, X = dat_b$snps, p0_av = p0,  Z = Z,
                  link = "mix", ind_bin = ind_bin, user_seed = user_seed)

}
\seealso{
\code{\link{set_hyper}}, \code{\link{set_init}},
  \code{\link{set_cv}}, \code{\link{set_blocks}}
}
