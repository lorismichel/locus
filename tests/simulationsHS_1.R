# Loris Michel, EPFL
library(scales)
cols <- cut(z, 6, labels = c("pink", "red", "yellow", "blue", "green", "purple"))
plot(c(1,2), c(1,2), main= "Fragment recruitment plot - FR-HIT",
     ylab = "Percent identity", xlab = "Base pair position",
     col = alpha("purple", 0.4), pch=16)
# first simu for assessing performances of HS vs locus
median.impute <- function(mat) {
  for(i in 1:ncol(mat)) {
    mat[is.na(mat[,i]),i] <- median(mat[,i],na.rm=T)
  }

  return(mat)
}

require(locus)
require(reshape2)
# define a seed for reproducibility

#RNGkind("L'Ecuyer-CMRG")
#skip.streams <- function(n) {
#  x <- .Random.seed
#  for (i in seq_len(n))
 #   x <- nextRNGStream(x)
 # assign('.Random.seed', x, pos=.GlobalEnv)
#}
user_seed <- 121
set.seed(user_seed)
# define the dimensions of the problems
n <- 50; p <- 1000; p0 <- 100; d <- 30; d0 <- 15

score.HS_ga <- list()
score.HS_exp_ga <- list()
score.HS <- list()
score.HS_exp <- list()
score.HS_plus <- list()
score.loc <- list()
labels <- list()

for (i in 1:3) {
# generate the data
list_X <- generate_snps(n = n, p = p,user_seed = NULL)
vec_rho <- runif(d, min = 0.25, max = 0.95)
vec_prob_sh <- 0.3
cor_type <- "equicorrelated"
list_Y <- generate_phenos(n = n, d = d, var_err = 1,cor_type = cor_type,vec_rho = vec_rho,
                          user_seed = NULL)

dat <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
                           ind_d0 = sample(1:d, d0), ind_p0 = sample(1:p, p0),
                           vec_prob_sh =  vec_prob_sh, max_tot_pve = 0.5,user_seed = NULL)
X <- scale(dat$snps)
if(sum(is.na(X))>0) next
  #X <- X[,apply(X,2,function(x) any(!is.na(x)))]

Y <- scale(dat$phenos,center = T,scale = F)


alpha_vb <- matrix(1,ncol=d,nrow=ncol(X))
typical_size <- 0.2
mu_beta_vb <- matrix(0,nrow=ncol(X),ncol=d)
sig2_beta_vb <- matrix(1/(n-1 + typical_size^{-2}),nrow=ncol(X),ncol=d)
sigma2_vb <- typical_size^{2}
tol <- 1
maxit <- 300
batch = TRUE
verbose = TRUE
full_output = TRUE

#if(sum(is.na(X))>0) {
#  X <- median.impute(X)
#}

c_vb <- matrix(1, nrow=ncol(X), ncol=d)

list_hyper <- list(lambda = 1,nu = typical_size,eta =  1 ,kappa = median(apply(Y, 2, var)),
                   A = rep(sqrt(median(apply(Y, 2, var))),d),B = typical_size)
tau_vb <- list_hyper$A
G_vb <- (1/2)*sweep(sig2_beta_vb,2,1/list_hyper$kappa,`*`)
# algorithms
HS_exp_ga <- horseshoe_core_exp_gamma(Y, X, d, n, p, list_hyper, alpha_vb, c_vb, mu_beta_vb,
                                   sig2_beta_vb, tau_vb, tol, maxit,  verbose, shared_prec = F,
                                   full_output = T)
HS_ga <-  horseshoe_core_gamma(Y, X, d, n, p, list_hyper = list_hyper, alpha_vb,  mu_beta_vb,
                            sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec =  F,
                            full_output = TRUE)
HS <- horseshoe_core(Y, X, d, n, p, list_hyper, alpha_vb, mu_beta_vb,
                     sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec = F,
                     full_output = TRUE)
HS_exp <- horseshoe_core_exp(Y, X, d, n, p, list_hyper, alpha_vb, c_vb, mu_beta_vb,
                                sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec = F,
                                full_output = T)
HS_plus <- horseshoe_core_plus(Y, X, d, n, p, list_hyper=list_hyper, alpha_vb, c_vb, d_vb = matrix(1,ncol=d,nrow=ncol(X)), e_vb = matrix(1,ncol=d,nrow=ncol(X)), mu_beta_vb,
                              sig2_beta_vb, tau_vb, tol, maxit, verbose=T, shared_prec = F,
                              full_output = T)
loc <- locus(Y = Y, X = X, p0_av =  vec_prob_sh*p0*d0,
              link = "identity", user_seed = NULL)

score.HS_ga[[i]] <- melt(1 / (1 + HS_ga$sig2_inv_vb*HS_ga$alpha_vb))$value

score.HS_exp_ga[[i]] <- melt(1 / (1 + HS_exp_ga$sig2_inv_vb*HS_exp_ga$alpha_vb))$value
score.HS[[i]] <- melt(1 / (1 + score.HS$sig2_inv_vb*score.HS$alpha_vb))$value
score.HS_exp[[i]] <- melt(1 / (1 + score.HS_exp$sig2_inv_vb*score.HS_exp$alpha_vb))$value
score.HS_plus[[i]] <- melt(1 / (1 + score.HS_plus$sig2_inv_vb*score.HS_plus$alpha_vb))$value
score.loc[[i]] <- melt(loc$gam)$value
labels[[i]] <- as.numeric(melt(dat$beta!=0)$value)
print(i)
local.preds.HS_ga <- prediction(score.HS_ga[[i]],labels[[i]])
local.perf.HS_ga <- performance(local.preds.HS_ga,"tpr","fpr")
local.preds.HS_exp_ga <- prediction(score.HS_exp_ga[[i]],labels[[i]])
local.perf.HS_exp_ga <- performance(local.preds.HS_exp_ga,"tpr","fpr")
if(i==1) {
  plot(local.perf.HS_ga,lty=1,lwd=1,col=alpha("purple", 0.4))
  pl.HS_ga <- recordPlot()
  plot(local.perf.HS_exp_ga,lty=1,lwd=1,col=alpha("brown", 0.4))
  pl.HS_exp_ga <- recordPlot()
} else {
  replayPlot(pl.HS_ga)
  plot(local.perf.HS_ga,lty=1,add=T,lwd=1,col=alpha("purple", 0.4))
  pl.HS_ga <- recordPlot()
  Sys.sleep(1)
  replayPlot(pl.HS_exp_ga)
  plot(local.perf.HS_exp_ga,lty=1,add=T,lwd=1,col=alpha("brown", 0.4))
  pl.HS_exp_ga <- recordPlot()
}
#skip.streams(3)
}

score.HS_ga <- score.HS_ga[-4]
score.HS_exp_ga <- score.HS_exp_ga[-4]
labels <- labels[-4]
preds.HS_ga <- prediction(score.HS_ga,labels)
preds.HS_exp_ga <- prediction(score.HS_exp_ga,labels)
#AUC.HS_ga <- performance(perf.HS_ga, "auc")
perf.HS_ga <- performance(preds.HS_ga,"tpr","fpr")
perf.HS_ga_auc <- performance(preds.HS_ga,"auc")

perf.HS_exp_ga <- performance(preds.HS_exp_ga,"tpr","fpr")
perf.HS_exp_ga_auc <- performance(preds.HS_exp_ga,"auc")
replayPlot(pl.HS_ga)
plot(perf.HS_ga,avg="vertical",col="purple",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T)
pl.HS_ga <- recordPlot()
replayPlot(pl.HS_exp_ga)
plot(perf.HS_exp_ga,avg="vertical",col="brown",lty=1,lwd=2,spread.estimate="stderror",plotCI.lwd=2,add=T)
pl.HS_exp_ga <- recordPlot()


png("~/Desktop/plot1.png")
replayPlot(pl.HS_ga)
dev.off()
png("~/Desktop/plot2.png")
replayPlot(pl.HS_exp_ga)
dev.off()
