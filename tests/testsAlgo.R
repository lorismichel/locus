user_seed <- 121; set.seed(user_seed)
 n <- 50; p <- 500; p0 <- 20; d <- 10; d0 <- 5
 list_X <- generate_snps(n = n, p = p)
 list_Y <- generate_phenos(n = n, d = d, var_err = 4)

 vec_prob_sh <- 0.2
 dat <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
                            ind_d0 = sample(1:d, d0), ind_p0 = sample(1:p, p0),
                            vec_prob_sh =  vec_prob_sh, max_tot_pve = 0.3)



 #beta <- matrix(rexp(n = p*d,rate = 1000),ncol=d,nrow=p)
 #noise <-  matrix(rnorm(n = n*d,mean = 0,sd = 1),ncol=d,nrow=n)
 #Y <- X %*% beta + noise
 alpha_vb <- matrix(1,ncol=d,nrow=p)
 #b_vb <- alpha_vb
 #lambda <- rep(1*10^{-1},d)
 #nu <- rep(1*10^{-1},d)
 #A <- rep(1, d)
 #B <- 1
 #list_hyper <- list(lambda = lambda, nu = nu, A = A, B = B)
 typical_size <- 0.2
 mu_beta_vb <- matrix(0,nrow=p,ncol=d)
 sig2_beta_vb <- matrix(1/(n-1 + typical_size^{-2}),nrow=p,ncol=d)
sigma2_vb <- typical_size^{2}
 tol <- 10^{-5}
 maxit <- 200
 batch = TRUE
 verbose = TRUE
 full_output = TRUE
 X <- scale(dat$snps)
 Y <- scale(dat$phenos,center = T,scale = F)
 G_vb <- (1/2)*sweep(sig2_beta_vb,2,1/list_hyper$kappa,`*`)
 c_vb <- matrix(1, nrow=p, ncol=d)

 list_hyper <- list(lambda = 1,nu = typical_size,eta =  1 ,kappa = median(apply(Y, 2, var)),
                    A = rep(sqrt(median(apply(Y, 2, var))),d),B = typical_size)
 tau_vb <- list_hyper$A
 exp_ga <- horseshoe_core_exp_gamma(Y, X, d, n, p, list_hyper, alpha_vb, c_vb, mu_beta_vb,
                                 sig2_beta_vb, tau_vb, tol, maxit,  verbose, shared_prec = F,
                                 full_output = T)

G_vb <- (1/2)*sweep(sig2_beta_vb,2,1/list_hyper$kappa,`*`)
ga <-  horseshoe_core_gamma(Y, X, d, n, p, list_hyper = list_hyper, alpha_vb,  mu_beta_vb,
                                       sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec =  F,
                                       full_output = TRUE)


G_vb <- matrix(1,ncol=d,nrow=p)
HS <- horseshoe_core(Y, X, d, n, p, list_hyper, alpha_vb, mu_beta_vb,
                                 sig2_beta_vb, tau_vb, tol, maxit, verbose, ,
                                 full_output = TRUE)


# cauchy <- locus_core_horseshoeCauchy(Y, X, d, n, p, list_hyper, b_vb, sigma2_bv, mu_beta_vb,
#                                      sig2_beta_vb, tau_vb, tol, maxit, batch, verbose, scheme = "noPrec",loop="c++",
#                                      full_output = T)
ExpCauchy <- horseshoe_core_exp(Y, X, d, n, p, list_hyper, alpha_vb, c_vb, mu_beta_vb,
                               sig2_beta_vb, tau_vb, tol, maxit, verbose, shared_prec = F,
                               full_output = T)


 HSplus <- horseshoe_core_plus(Y, X, d, n, p, list_hyper=list_hyper, alpha_vb, c_vb, d_vb = matrix(1,ncol=d,nrow=p), e_vb = matrix(1,ncol=d,nrow=p), mu_beta_vb,
                                            sig2_beta_vb, tau_vb, tol, maxit, verbose=T, shared_prec = F,
                                            full_output = T)


 #mod <- locus_core_horseshoeExp(Y, X, d = d, n = n, p = p, list_hyper, b_vb, c_vb, sigma2_vb, mu_beta_vb,
 #                                    sig2_beta_vb, tau_vb, tol, maxit, batch, verbose, scheme = "Prec",
 #                                    loop="c++",full_output = TRUE)

vb_g <- locus(Y = Y, X = X, p0_av =  vec_prob_sh*p0*d0,
               link = "identity", user_seed = user_seed)

#r#es<- locus_core_horseshoe(Y, X, d = d, n = n, p = p, list_hyper, b_vb, sigma2_vb, mu_beta_vb,
 #                                  sig2_beta_vb, tau_vb, tol, maxit = maxit, batch, verbose,
#                                   scheme="noPrec",loop = "c++",full_output = TRUE)

plot(exp_ga$ELBO,type="l")
plot(ga$ELBO,type="l")
plot(HS$ELBO,type="l")
plot(ExpCauchy$ELBO,type="l")
plot(HSplus$ELBO,type="l")
# roc curves
require(ROCR)
#score.HSEXP <- (1 / (n-1 + mod$sig2_inv_vb*mod$))
score.ga <- (1 / (1 + ga$sig2_inv_vb*ga$alpha_vb))
score.exp_ga <- (1 / (1 + exp_ga$sig2_inv_vb*exp_ga$alpha_vb))
score.HScauchy <- (1 / (1 + HS$sig2_inv_vb*HS$alpha_vb))
score.HSEXPcauchy <- (1 / (1 + ExpCauchy$sig2_inv_vb*ExpCauchy$alpha_vb))
score.HSplus <- (1 / (1 + HSplus$sig2_inv_vb*HSplus$alpha_vb))


preds.LOCUS <- prediction(as.numeric(vb_g$gam_vb),as.numeric(dat$beta != 0))
preds.HScauchy <- prediction(as.numeric(score.HScauchy),as.numeric(dat$beta != 0))
preds.HSEXPcauchy <- prediction(as.numeric(score.HSEXPcauchy),as.numeric(dat$beta != 0))
preds.ga <- prediction(as.numeric(score.ga),as.numeric(dat$beta != 0))
preds.exp_ga <- prediction(as.numeric(score.exp_ga),as.numeric(dat$beta != 0))
preds.HSplus <- prediction(as.numeric(score.HSplus),as.numeric(dat$beta != 0))

plot(performance(preds.LOCUS,"tpr","fpr"),col="blue")
plot(performance(preds.HScauchy ,"tpr","fpr"),col="pink",add=T)
plot(performance(preds.HSEXPcauchy ,"tpr","fpr"),col="black",add=T)
plot(performance(preds.ga ,"tpr","fpr"),col="orange",add=T)
plot(performance(preds.exp_ga ,"tpr","fpr"),col="purple",add=T)
plot(performance(preds.HSplus ,"tpr","fpr"),col="green",add=T)
#score.ga <- 1/ga$b_vb
#score.exp_ga <- 1/exp_ga$alpha_vb
#score.HSEXP <- 1 / mod$b_vb
#score.HSEXPcauchy <- (1 / (n-1 + ExpCauchy$sig2_inv_vb*ExpCauchy$))
#score.HS <- (1 / (n-1 + res$sig2_inv_vb*res$b_vb))

