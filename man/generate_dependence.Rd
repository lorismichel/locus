% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_data.R
\name{generate_dependence}
\alias{generate_dependence}
\title{Generate pleiotropic associations between SNPs and phenotypes.}
\usage{
generate_dependence(list_snps, list_phenos, ind_d0, ind_p0, vec_prob_sh,
  family = "gaussian", pve_per_snp = NULL, max_tot_pve = NULL,
  user_seed = NULL)
}
\arguments{
\item{list_snps}{An object of class "sim_snps" containing simulated SNP data
and their corresponding sample minor allele frequencies. It must be
obtained from the function \code{\link{generate_snps}} or
\code{\link{replicate_real_snps}}.}

\item{list_phenos}{An object of class "sim_pheno" containing simulated
phenotypic data, their sample variance and block structure information.
It must be obtained from the function \code{\link{generate_phenos}} or
\code{\link{replicate_real_phenos}}.}

\item{ind_d0}{A vector of indices specifying the position of the "active"
phenotypes (i.e., which will be associated with at least one SNP). Must
range between 1 and \code{ncol(list_phenos$phenos)}.}

\item{ind_p0}{A vector of indices specifying the position of the "active"
SNPs (i.e., which will be associated with at least one phenotype). Must
range between 1 and \code{ncol(list_snps$snps)}.}

\item{vec_prob_sh}{Vector providing a set of probabilities with which an
active SNP is associated with an additional active phenotype in a given
phenotypic block. See Details section.}

\item{family}{Distribution used to generate the phenotypes. Must be either
"\code{gaussian}" or "\code{binomial}" for binary phenotypes.}

\item{pve_per_snp}{Average proportion of phenotypic variance explained by
each active SNP (for an active phenotype). Must be \code{NULL} if
\code{max_tot_pve} is provided. See Details section.}

\item{max_tot_pve}{Maximum proportion of phenotypic variance explained by the
active SNPs across all phenotypes. Must be \code{NULL} if
\code{pve_per_snp} is provided. See Details section.}

\item{user_seed}{Seed set for reproducibility. Default is \code{NULL}, no
seed set.}
}
\value{
An object of class "\code{sim_data}".
 \item{phenos}{Matrix containing the updated phenotypic data (whose variance
               is now partly explained by genetic effects).}
 \item{snps}{Matrix containing the original SNPs data.}
 \item{beta}{Matrix containing the generated effect sizes between the SNPs
            (rows) and phenotypes (columns).}
 \item{pat}{Matrix of booleans specifying the generated association pattern
            between the SNPs (rows) and phenotypes (columns).}
 \item{pve_per_snp}{Average proportion of phenotypic variance explained by
                    each active SNP (for an active phenotype).}
}
\description{
This function sets the association pattern and the effect sizes between SNP
and phenotype objects previously obtained from the functions
\code{\link{generate_snps}} or \code{\link{replicate_real_snps}}, and
\code{\link{generate_phenos}} or \code{\link{replicate_real_phenos}}. It
therefore adds a genetic contribution to the phenotypic data.
}
\details{
The user can provide using the argument \code{vec_prob_sh} a selection of
probabilities describing the propensity with which a given active SNP (i.e.,
associated with at least one phenotype) will be associated with active
phenotypes (i.e., associated with at least one SNP) of given phenotypic
blocks. More precisely, for each active SNP and each phenotypic block, a
value from this vector is selected uniformly at random; for instance a large
probability implies that the SNPs is highly likely to be associated with
each active phenotype in the block. If a single value is provided, all active
SNPs will have the same probability to be associated with active phenotypes
of all blocks.

The user can provide either argument \code{pve_per_snp}, specifying the
average proportion of phenotypic variance explained per active SNP for a
given active phenotype, or \code{max_tot_pve}, specifying the maximum value
for an active phenotype of its proportion of variance explained by the
cummulated genetic effects. If both \code{pve_per_snp} and \code{max_tot_pve}
are \code{NULL}, the proportion of phenotypic variance explained per SNP is
set to its maximum value so that the total proportion of variance explained
for the phenotypes are all below 1. Individual proportions of variance
explained are drawn from a Beta distribution with shape parameters 2 and 5,
putting more weights on smaller effects.

If family is "\code{binomial}", the phenotypes are generated from a probit
model, and the phenotypic variance explained by the SNPs is with respect to
the latent Gaussian variables involved in the probit model.
}
\examples{
user_seed <- 123; set.seed(user_seed)
n <- 500; p <- 5000; p0 <- 200; d <- 500; d0 <- 400

list_snps <- generate_snps(n = n, p = p)

cor_type <- "equicorrelated"; vec_rho <- runif(100, min = 0.25, max = 0.95)

list_phenos <- generate_phenos(n, d, var_err = 1, cor_type = cor_type,
                               vec_rho = vec_rho, n_cpus = 2)

# Gaussian phenotypes
dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = sample(1:d, d0),
                           ind_p0 = sample(1:p, p0), vec_prob_sh = 0.05,
                           family = "gaussian", max_tot_pve = 0.5)

# Binary phenotypes
dat_b <- generate_dependence(list_snps, list_phenos, ind_d0 = sample(1:d, d0),
                           ind_p0 = sample(1:p, p0), vec_prob_sh = 0.05,
                           family = "binomial", max_tot_pve = 0.5)

}
\seealso{
\code{\link{generate_snps}}, \code{\link{replicate_real_snps}},
  \code{\link{generate_phenos}}, \code{\link{replicate_real_phenos}}
}
