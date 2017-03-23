user_seed <- 121; set.seed(user_seed)
 n <- 80; p <- 500; p0 <- 50; d <- 50; d0 <- 20
 list_X <- generate_snps(n = n, p = p)
 list_Y <- generate_phenos(n = n, d = d, var_err = 0.8)

 dat <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
                            ind_d0 = sample(1:d, d0), ind_p0 = sample(1:p, p0),
                            vec_prob_sh = 0.5, max_tot_pve = 0.3)


 #beta <- matrix(rexp(n = p*d,rate = 1000),ncol=d,nrow=p)
 #noise <-  matrix(rnorm(n = n*d,mean = 0,sd = 1),ncol=d,nrow=n)
 #Y <- X %*% beta + noise
 b_vb <- matrix(1,ncol=d,nrow=p)
 lambda <- rep(2*10^{-1},d)
 nu <- rep(1*10^{-1},d)
 A <- 10^{-2}
 list_hyper <- list(lambda = lambda, nu = nu, A = A)
 sigma2_vb <- 1000
 mu_beta_vb <- matrix(0,nrow=p,ncol=d)
 sig2_beta_vb <- matrix(10^{-5},nrow=p,ncol=d)
 tau_vb <- lambda / nu
 tol <- 1
 maxit <- 120
 batch = TRUE
 verbose = TRUE
 full_output = TRUE
 X <- scale(dat$snps)
 Y <- scale(dat$phenos,center = T,scale = F)
 G_vb <- b_vb
 c_vb <- matrix(1/2, nrow=p, ncol=d)

 mod <- locus_core_horseshoeExp(Y, X, d, n, p, list_hyper, b_vb, c_vb, sigma2_vb, mu_beta_vb,
                                     sig2_beta_vb, tau_vb, tol, maxit, batch, verbose, scheme = "Prec",
                                     loop="c++",full_output = TRUE)

vb_g <- locus(Y = dat$phenos, X = dat$snps, p0_av = 100,
               link = "identity", user_seed = user_seed)

res<- locus_core_horseshoe(Y, X, d = d, n = n, p = p, list_hyper, b_vb, sigma2_vb, mu_beta_vb,
                                   sig2_beta_vb, tau_vb, tol, maxit = maxit, batch, verbose,
                                   scheme="Prec",loop = "c++",full_output = TRUE)

# roc curves
require(ROCR)
score.HSEXP <- (1 / (n-1 + mod$sig2_inv_vb*mod$b_vb))
#score.HSEXP <- 1 / mod$b_vb
score.HS <- (1 / (n-1 + res$sig2_inv_vb*res$b_vb))
preds.HSEXP <- prediction(as.numeric(score.HSEXP),as.numeric(dat$beta != 0))
preds.LOCUS <- prediction(as.numeric(vb_g$gam_vb),as.numeric(dat$beta != 0))
preds.HS <- prediction(as.numeric(score.HS),as.numeric(dat$beta != 0))
plot(performance(preds.LOCUS,"tpr","fpr"),col="blue")
plot(performance(preds.HSEXP ,"tpr","fpr"),col="green",add=T)
plot(performance(preds.HS ,"tpr","fpr"),col="brown",add=T)


score <- 1 / (1 + mod$b_vb)
w <- cbind(dat$beta[,1],mod$mu_beta_vb[,1])
w <- cbind(dat$beta[,1],score.HSEXP[,1])
w <- cbind(dat$beta[,1],mod$b_vb[,1])
w <- cbind(dat$beta[,1],mod$sig2_beta_vb[,1])
w <- cbind(dat$beta[,1],G_vb[,1])
w <- cbind(dat$beta[,1],m2_beta[,1])
m <- res$mu_beta_vb / sqrt(res$sig2_beta_vb)
w <- cbind(dat$beta[,1],m[,1])

vb <- locus(Y = dat$phenos, X = dat$snps, p0_av = p0, user_seed = user_seed,batch = T,verbose = T)






 user_seed <- 123; set.seed(user_seed)
 n <- 200; p <- 300; p0 <- 50; d <- 40; d0 <- 30
 list_X <- generate_snps(n = n, p = p)
 list_Y <- generate_phenos(n = n, d = d, var_err = 0.25)

dat <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
                            ind_d0 = sample(1:d, d0), ind_p0 = sample(1:p, p0),
                            vec_prob_sh = 0.1, max_tot_pve = 0.9)

 vb <- locus(Y = dat$phenos, X = dat$snps, p0_av = p0, user_seed = user_seed)






 user_seed <- 123; set.seed(user_seed)
 n <- 200; p <- 300; p0 <- 50; d <- 40; d0 <- 30
 list_X <- generate_snps(n = n, p = p)
 list_Y <- generate_phenos(n = n, d = d, var_err = 0.25)

 dat <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
                            ind_d0 = sample(1:d, d0), ind_p0 = sample(1:p, p0),
                            vec_prob_sh = 0.1, max_tot_pve = 0.9)







 # test of new method
 b_vb <- matrix(1,ncol=40,nrow=300)
 gam_vb <- matrix(0.9,ncol=40,nrow=300)
 lambda <- 1
 nu <- 1
 eta <- rep(1,40)
 kappa <- rep(1,40)

 d <- 40
 A <- 100
 B <- 100
 list_hyper <- list(lambda = lambda,nu = nu, eta = eta, kappa = kappa,alpha = 1,psi = 1,A = A,B=B)
 sigma2_vb <- 1
 mu_beta_vb <- matrix(0,nrow=300,ncol=40)
 sig2_beta_vb <- matrix(rnorm(300*40,0,sqrt(sigma2_vb))^{2},nrow=300,ncol=40)
 tau_vb <- rep(1,40)
 tol <- 1
 maxit <- 6
 batch = TRUE
 verbose = TRUE
 full_output = TRUE
 X <- scale(dat$snps)
 Y <- scale(dat$phenos, center = T,scale = F)

 vb_g <- locus(Y = Y, X = X, p0_av = p0,
               link = "identity", user_seed = user_seed)

 a<-locus_student_core_(Y = Y, X = X, list_hyper = list_hyper, gam_vb = gam_vb, mu_beta_vb = mu_beta_vb, sig2_beta_vb = sig2_beta_vb,
                                 tau_vb = tau_vb, b_vb = b_vb, c_vb = A/B, tol = tol, maxit = maxit, batch = batch, verbose = verbose, full_output = FALSE)
 boxplot(a$gam_vb)

