# Loris Michel, EPFL
# this is a first script to validate on a small sample both the fact that
# the ELBO is increasing, converging and that the results make sense
# with the horseshoe
require(locus)
require(gsl)
require(scales)
require(reshape2)

# repro
user_seed <- 121
set.seed(user_seed)
# define the dimensions of the problems
n <- 250; p <- 5000; p0 <- 100; d <- 50; d0 <- 40

score.HS_ga <- list()
score.HS_ga_prec <- list()

score.HS_exp_ga <- list()
score.HS_exp_ga_prec <- list()

score.HS <- list()
score.HS_prec <- list()

score.HS_exp <- list()
score.HS_exp_prec <- list()

score.HS_plus <- list()
score.HS_plus_prec <- list()

score.loc <- list()

labels <- list()

for (i in 1:100) {

list_X <- generate_snps(n = n, p = p,
                        cor_type = "equicorrelated",
                        user_seed = NULL,vec_rho = 0.75)

# params for phenos
vec_rho <- c(0.8,0.3,0.2,0.5)
var_err = 1
cor_type <- "equicorrelated"

list_Y <- generate_phenos(n = n, d = d, var_err = var_err,
                          cor_type = cor_type,
                          vec_rho = vec_rho,
                          user_seed = NULL)

vec_prob_sh <- 0.15
dat <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
                           ind_d0 = sample(1:d, d0), ind_p0 = sample(1:p, p0),
                           vec_prob_sh =  vec_prob_sh,
                           max_tot_pve = 0.3,user_seed = NULL)
X <- scale(dat$snps)
if(sum(is.na(X))>0) next
  #X <- X[,apply(X,2,function(x) any(!is.na(x)))]

Y <- scale(dat$phenos,center = T,scale = F)


alpha_vb <- matrix(1,ncol=d,nrow=ncol(X))
typical_size <- 0.01
mu_beta_vb <- matrix(0,nrow=ncol(X),ncol=d)
sig2_beta_vb <- matrix(1/(n-1 + typical_size^{-2}),nrow=ncol(X),ncol=d)
sigma2_vb <- typical_size^{2}
tol <- 10
maxit <- 200
batch = TRUE
verbose = TRUE
full_output = TRUE


c_vb <- matrix(1, nrow=ncol(X), ncol=d)
list_hyper <- list(lambda = 1,nu = typical_size,eta =  1 ,kappa = median(apply(Y, 2, var)),
                   A = rep(sqrt(median(apply(Y, 2, var))),d),B = typical_size)
tau_vb <- list_hyper$A
G_vb <- (1/2)*sweep(sig2_beta_vb,2,1/list_hyper$kappa,`*`)

# algorithms

# exponential HS with gamma prior on variances
HS_exp_ga <- horseshoe_core_exp_gamma(Y, X, d, n, p, list_hyper, alpha_vb, c_vb, mu_beta_vb,
                                           sig2_beta_vb, tau_vb, tol, maxit,  verbose, shared_prec = F,
                                           full_output = T)
HS_exp_ga_prec <- horseshoe_core_exp_gamma(Y, X, d, n, p, list_hyper, alpha_vb, c_vb, mu_beta_vb,
                                   sig2_beta_vb, tau_vb, tol, maxit,  verbose, shared_prec = T,
                                   full_output = T)

# HS with gamma prior on variances
HS_ga <-  horseshoe_core_gamma(Y, X, d, n, p, list_hyper = list_hyper, alpha_vb,  mu_beta_vb,
                            sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec =  F,
                            full_output = TRUE)
HS_ga_prec <-  horseshoe_core_gamma(Y, X, d, n, p, list_hyper = list_hyper, alpha_vb,  mu_beta_vb,
                               sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec =  T,
                               full_output = TRUE)

# HS with Cauchy priors
HS <- horseshoe_core(Y, X, d, n, p, list_hyper, alpha_vb, mu_beta_vb,
                     sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec = F,
                     full_output = TRUE)
HS_prec <- horseshoe_core(Y, X, d, n, p, list_hyper, alpha_vb, mu_beta_vb,
                          sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec = T,
                          full_output = TRUE)

# exponential HS with Cauchy prior on variances
HS_exp <- horseshoe_core_exp(Y, X, d, n, p, list_hyper, alpha_vb, c_vb, mu_beta_vb,
                                sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec = F,
                                full_output = T)
HS_exp_prec <- horseshoe_core_exp(Y, X, d, n, p, list_hyper, alpha_vb, c_vb, mu_beta_vb,
                             sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec = T,
                             full_output = T)

# plus HS with cauchy priors
HS_plus <- horseshoe_core_plus(Y, X, d, n, p, list_hyper=list_hyper, alpha_vb, c_vb, d_vb = matrix(1,ncol=d,nrow=ncol(X)), e_vb = matrix(1,ncol=d,nrow=ncol(X)), mu_beta_vb,
                              sig2_beta_vb, tau_vb, tol, maxit, verbose=T, shared_prec = F,
                              full_output = T)
HS_plus_prec <- horseshoe_core_plus(Y, X, d, n, p, list_hyper=list_hyper, alpha_vb, c_vb, d_vb = matrix(1,ncol=d,nrow=ncol(X)), e_vb = matrix(1,ncol=d,nrow=ncol(X)), mu_beta_vb,
                               sig2_beta_vb, tau_vb, tol, maxit, verbose=T, shared_prec = T,
                               full_output = T)

# locus
loc <- locus(Y = Y, X = X, p0_av =  p0*d0,
              link = "identity",tol = 10^{-2}, user_seed = NULL)

# compute the scores
score.HS_ga_prec[[i]] <- melt(1 / (n-1 + HS_ga_prec$sig2_inv_vb*HS_ga_prec$alpha_vb))$value
score.HS_ga[[i]] <- melt(1 / (n-1 + HS_ga$sig2_inv_vb*HS_ga$alpha_vb))$value

score.HS_exp_ga_prec[[i]] <- melt(1 / (n-1 + HS_exp_ga_prec$sig2_inv_vb*HS_exp_ga_prec$alpha_vb))$value
score.HS_exp_ga[[i]] <- melt(1 / (n-1 + HS_exp_ga$sig2_inv_vb*HS_exp_ga$alpha_vb))$value

score.HS_prec[[i]] <- melt(1 / (n-1 + HS_prec$sig2_inv_vb*HS_prec$alpha_vb))$value
score.HS[[i]] <- melt(1 / (n-1 + HS$sig2_inv_vb*HS$alpha_vb))$value

score.HS_exp_prec[[i]] <- melt(1 / (n-1 + HS_exp_prec$sig2_inv_vb*HS_exp_prec$alpha_vb))$value
score.HS_exp[[i]] <- melt(1 / (n-1 + HS_exp$sig2_inv_vb*HS_exp$alpha_vb))$value


score.HS_plus_prec[[i]] <- melt(1 / (n-1 + HS_plus_prec$sig2_inv_vb*HS_plus_prec$alpha_vb))$value
score.HS_plus[[i]] <- melt(1 / (n-1 + HS_plus$sig2_inv_vb*HS_plus$alpha_vb))$value

score.loc[[i]] <- melt(loc$gam)$value
labels[[i]] <- as.numeric(melt(dat$beta!=0)$value)

print(paste("finished iteration: ",i))

local.preds.HS_ga <- prediction(score.HS_ga[[i]],labels[[i]])
local.perf.HS_ga <- performance(local.preds.HS_ga,"tpr","fpr")

local.preds.HS_ga_prec <- prediction(score.HS_ga_prec[[i]],labels[[i]])
local.perf.HS_ga_prec <- performance(local.preds.HS_ga_prec,"tpr","fpr")

local.preds.HS_exp_ga <- prediction(score.HS_exp_ga[[i]],labels[[i]])
local.perf.HS_exp_ga <- performance(local.preds.HS_exp_ga,"tpr","fpr")

local.preds.HS_exp_ga_prec <- prediction(score.HS_exp_ga_prec[[i]],labels[[i]])
local.perf.HS_exp_ga_prec <- performance(local.preds.HS_exp_ga_prec,"tpr","fpr")

local.preds.HS <- prediction(score.HS[[i]],labels[[i]])
local.perf.HS <- performance(local.preds.HS,"tpr","fpr")

local.preds.HS_prec <- prediction(score.HS_prec[[i]],labels[[i]])
local.perf.HS_prec <- performance(local.preds.HS_prec,"tpr","fpr")

local.preds.HS_exp <- prediction(score.HS_exp[[i]],labels[[i]])
local.perf.HS_exp <- performance(local.preds.HS_exp,"tpr","fpr")

local.preds.HS_exp_prec <- prediction(score.HS_exp_prec[[i]],labels[[i]])
local.perf.HS_exp_prec <- performance(local.preds.HS_exp_prec,"tpr","fpr")

local.preds.HS_plus <- prediction(score.HS_plus[[i]],labels[[i]])
local.perf.HS_plus <- performance(local.preds.HS_plus,"tpr","fpr")

local.preds.HS_plus_prec <- prediction(score.HS_plus_prec[[i]],labels[[i]])
local.perf.HS_plus_prec <- performance(local.preds.HS_plus_prec,"tpr","fpr")

local.preds.locus <- prediction(score.loc[[i]],labels[[i]])
local.perf.locus <- performance(local.preds.locus,"tpr","fpr")

if(i==1) {
  plot(local.perf.HS_ga,lty=1,lwd=1,col=alpha("purple", 0.4))
  pl.HS_ga <- recordPlot()
  plot(local.perf.HS_ga_prec,lty=1,lwd=1,col=alpha("purple", 0.4))
  pl.HS_ga_prec <- recordPlot()

  plot(local.perf.HS_exp_ga,lty=1,lwd=1,col=alpha("brown", 0.4))
  pl.HS_exp_ga <- recordPlot()
  plot(local.perf.HS_exp_ga_prec,lty=1,lwd=1,col=alpha("brown", 0.4))
  pl.HS_exp_ga_prec <- recordPlot()

  plot(local.perf.HS,lty=1,lwd=1,col=alpha("blue", 0.4))
  pl.HS <- recordPlot()
  plot(local.perf.HS_prec,lty=1,lwd=1,col=alpha("blue", 0.4))
  pl.HS_prec <- recordPlot()

  plot(local.perf.HS_exp,lty=1,lwd=1,col=alpha("red", 0.4))
  pl.HS_exp <- recordPlot()
  plot(local.perf.HS_exp_prec,lty=1,lwd=1,col=alpha("red", 0.4))
  pl.HS_exp_prec <- recordPlot()

  plot(local.perf.HS_plus,lty=1,lwd=1,col=alpha("orange", 0.4))
  pl.HS_plus <- recordPlot()
  plot(local.perf.HS_plus_prec,lty=1,lwd=1,col=alpha("orange", 0.4))
  pl.HS_plus_prec <- recordPlot()

  plot(local.perf.locus,lty=1,lwd=1,col=alpha("orange", 0.4))
  pl.locus <- recordPlot()


} else {
  replayPlot(pl.HS_ga)
  plot(local.perf.HS_ga,lty=1,lwd=1,add=T,col=alpha("purple", 0.4))
  pl.HS_ga <- recordPlot()
  replayPlot(pl.HS_ga_prec)
  plot(local.perf.HS_ga_prec,lty=1,lwd=1,add=T,col=alpha("purple", 0.4))
  pl.HS_ga_prec <- recordPlot()
  replayPlot(pl.HS_exp_ga)
  plot(local.perf.HS_exp_ga,lty=1,lwd=1,add=T,col=alpha("brown", 0.4))
  pl.HS_exp_ga <- recordPlot()
  replayPlot(pl.HS_exp_ga_prec)
  plot(local.perf.HS_exp_ga_prec,lty=1,lwd=1,add=T,col=alpha("brown", 0.4))
  pl.HS_exp_ga_prec <- recordPlot()
  replayPlot(pl.HS)
  plot(local.perf.HS,lty=1,lwd=1,add=T,col=alpha("blue", 0.4))
  pl.HS <- recordPlot()
  replayPlot(pl.HS_prec)
  plot(local.perf.HS_prec,lty=1,lwd=1,add=T,col=alpha("blue", 0.4))
  pl.HS_prec <- recordPlot()
  replayPlot(pl.HS_exp)
  plot(local.perf.HS_exp,lty=1,lwd=1,add=T,col=alpha("red", 0.4))
  pl.HS_exp <- recordPlot()
  replayPlot(pl.HS_exp_prec)
  plot(local.perf.HS_exp_prec,lty=1,lwd=1,add=T,col=alpha("red", 0.4))
  pl.HS_exp_prec <- recordPlot()
  replayPlot(pl.HS_plus)
  plot(local.perf.HS_plus,lty=1,lwd=1,add=T,col=alpha("orange", 0.4))
  pl.HS_plus <- recordPlot()
  replayPlot(pl.HS_plus_prec)
  plot(local.perf.HS_plus_prec,lty=1,lwd=1,add=T,col=alpha("orange", 0.4))
  pl.HS_plus_prec <- recordPlot()
  replayPlot(pl.locus)
  plot(local.perf.locus,lty=1,lwd=1,add=T,col=alpha("orange", 0.4))
  pl.locus <- recordPlot()
}
#skip.streams(3)
}



score.HS_ga <- Filter(Negate(function(x) is.null(unlist(x))), score.HS_ga)
score.HS_ga_prec <- Filter(Negate(function(x) is.null(unlist(x))), score.HS_ga_prec)
score.HS_exp_ga <- Filter(Negate(function(x) is.null(unlist(x))), score.HS_exp_ga)
score.HS_exp_ga_prec <- Filter(Negate(function(x) is.null(unlist(x))), score.HS_exp_ga_prec)
score.HS <- Filter(Negate(function(x) is.null(unlist(x))), score.HS)
score.HS_prec <- Filter(Negate(function(x) is.null(unlist(x))), score.HS_prec)
score.HS_exp <- Filter(Negate(function(x) is.null(unlist(x))), score.HS_exp)
score.HS_exp_prec <- Filter(Negate(function(x) is.null(unlist(x))), score.HS_exp_prec)
score.HS_plus <- Filter(Negate(function(x) is.null(unlist(x))), score.HS_plus)
score.HS_plus_prec <- Filter(Negate(function(x) is.null(unlist(x))), score.HS_plus_prec)
score.loc <- Filter(Negate(function(x) is.null(unlist(x))), score.loc)

labels <- Filter(Negate(function(x) is.null(unlist(x))), labels)



preds.HS_ga <- prediction(score.HS_ga,labels)
preds.HS_ga_prec <- prediction(score.HS_ga_prec,labels)
preds.HS_exp_ga <- prediction(score.HS_exp_ga,labels)
preds.HS_exp_ga_prec <- prediction(score.HS_exp_ga_prec,labels)
preds.HS <- prediction(score.HS,labels)
preds.HS_prec <- prediction(score.HS_prec,labels)
preds.HS_exp <- prediction(score.HS_exp,labels)
preds.HS_exp_prec <- prediction(score.HS_exp_prec,labels)
preds.HS_plus <- prediction(score.HS_plus,labels)
preds.HS_plus_prec <- prediction(score.HS_plus_prec,labels)
preds.locus <- prediction(score.loc,labels)



perf.HS_ga <- performance(preds.HS_ga,"tpr","fpr")
perf.HS_ga_auc <- performance(preds.HS_ga,"auc")
perf.HS_ga_prec <- performance(preds.HS_ga_prec,"tpr","fpr")
perf.HS_ga_prec_auc <- performance(preds.HS_ga_prec,"auc")
perf.HS_exp_ga <- performance(preds.HS_exp_ga,"tpr","fpr")
perf.HS_exp_ga_auc <- performance(preds.HS_exp_ga,"auc")
perf.HS_exp_ga_prec <- performance(preds.HS_exp_ga_prec,"tpr","fpr")
perf.HS_exp_ga_prec_auc <- performance(preds.HS_exp_ga_prec,"auc")
perf.HS <- performance(preds.HS,"tpr","fpr")
perf.HS_auc <- performance(preds.HS,"auc")
perf.HS_prec <- performance(preds.HS_prec,"tpr","fpr")
perf.HS_prec_auc <- performance(preds.HS_prec,"auc")
perf.HS_exp <- performance(preds.HS_exp,"tpr","fpr")
perf.HS_exp_auc <- performance(preds.HS_exp,"auc")
perf.HS_exp_prec <- performance(preds.HS_exp_prec,"tpr","fpr")
perf.HS_exp_prec_auc <- performance(preds.HS_exp_prec,"auc")
perf.HS_plus <- performance(preds.HS_plus,"tpr","fpr")
perf.HS_plus_auc <- performance(preds.HS_plus,"auc")
perf.HS_plus_prec <- performance(preds.HS_plus_prec,"tpr","fpr")
perf.HS_plus_prec_auc <- performance(preds.HS_plus_prec,"auc")
perf.locus <- performance(preds.locus,"tpr","fpr")
perf.locus_auc <- performance(preds.locus,"auc")

save(perf.HS_ga_auc, perf.HS_ga_prec_auc,
     perf.HS_exp_ga_auc, perf.HS_exp_ga_prec_auc,
     perf.HS_auc, perf.HS_prec_auc,
     perf.HS_exp_auc, perf.HS_exp_prec_auc,
     perf.HS_plus_auc, perf.HS_plus_prec_auc,
     perf.locus_auc, file="perf_p5000_d50.Rdata")

replayPlot(pl.HS_ga)
plot(perf.HS_ga,avg="vertical",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T,col=alpha("purple", 0.4))
pl.HS_ga <- recordPlot()
replayPlot(pl.HS_ga_prec)
plot(perf.HS_ga_prec,avg="vertical",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T,col=alpha("purple", 0.4))
pl.HS_ga_prec <- recordPlot()
replayPlot(pl.HS_exp_ga)
plot(perf.HS_exp_ga,avg="vertical",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T,col=alpha("brown", 0.4))
pl.HS_exp_ga <- recordPlot()
replayPlot(pl.HS_exp_ga_prec)
plot(perf.HS_exp_ga_prec,avg="vertical",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T,col=alpha("brown", 0.4))
pl.HS_exp_ga_prec <- recordPlot()
replayPlot(pl.HS)
plot(perf.HS,avg="vertical",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T,col=alpha("blue", 0.4))
pl.HS <- recordPlot()
replayPlot(pl.HS_prec)
plot(perf.HS_prec,avg="vertical",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T,col=alpha("blue", 0.4))
pl.HS_prec <- recordPlot()
replayPlot(pl.HS_exp)
plot(perf.HS_exp,avg="vertical",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T,col=alpha("red", 0.4))
pl.HS_exp <- recordPlot()
replayPlot(pl.HS_exp_prec)
plot(perf.HS_exp_prec,avg="vertical",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T,col=alpha("red", 0.4))
pl.HS_exp_prec <- recordPlot()
replayPlot(pl.HS_plus)
plot(perf.HS_plus,avg="vertical",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T,col=alpha("orange", 0.4))
pl.HS_plus <- recordPlot()
replayPlot(pl.HS_plus_prec)
plot(perf.HS_plus_prec,avg="vertical",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T,col=alpha("orange", 0.4))
pl.HS_plus_prec <- recordPlot()
replayPlot(pl.locus)
plot(perf.locus,avg="vertical",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T,col=alpha("orange", 0.4))
pl.locus <- recordPlot()


png("pl.HS_ga.png")
replayPlot(pl.HS_ga)
dev.off()

png("pl.HS_ga_prec.png")
replayPlot(pl.HS_ga_prec)
dev.off()

png("pl.HS_exp_ga.png")
replayPlot(pl.HS_exp_ga)
dev.off()

png("pl.HS_exp_ga_prec.png")
replayPlot(pl.HS_exp_ga_prec)
dev.off()

png("pl.HS.png")
replayPlot(pl.HS)
dev.off()

png("pl.HS_prec.png")
replayPlot(pl.HS_prec)
dev.off()

png("pl.HS_exp.png")
replayPlot(pl.HS_exp)
dev.off()

png("pl.HS_exp_prec.png")
replayPlot(pl.HS_exp_prec)
dev.off()

png("pl.HS_plus.png")
replayPlot(pl.HS_plus)
dev.off()

png("pl.HS_plus_prec.png")
replayPlot(pl.HS_plus_prec)
dev.off()

png("pl.locus.png")
replayPlot(pl.locus)
dev.off()




