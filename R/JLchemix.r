#' Latent Class Models for Joint Analysis of Disease Prevalence and High-dimensional Semicontinuous Chimical Biomarker Data
#' @import MCMCpack mvtnorm stats
#' @description This function fit a couple-based joint latent class model with an interaction between a
#' couple(e.g., female and male partners) and High-dimensional semicontinuous chemical biomarker
#' for each partner of the couple. This formulation introduces a dependence structure between the chemical
#' patterns within a couple and between the chemical patterns and the risk of desease.
#' A Bayesian framework examines the chemical biomarker profile from each member of the couple
#' and the risk of disease. The complex chemical mixtures on each couple link to disease risk through unobserved
#' latent classes. we posit that two sets of latent classes, each characterizing the chemical mixture patterns
#' of one partner of the couple, are linked to the risk of disease through a logistic model with main and
#' interaction  effects between latent classes. The semicontinuous chimical biomarker viarables (1/4 zeros and
#' right-skewed non-zero values) are processed through Tobit modeling framework. Markov chain Monte Carlo
#' algorithms was used to obtain posterior estimates of model parameters.The user supplies data and priors,
#' and a list of posterior estimates of model parameters is returned.
#'
#' @param nsim Number of simulations
#' @param nburn Burn in number
#' @param Yvariable Binary indicating dependent variable for the couple disease status, 1 for disease
#' @param X.f_mat chemical exposure variables for female individual
#' @param X.m_mat chemical exposure variables for male individual
#' @param covariate_f_mat subject-specific covariates such as age or smoking status for female individual
#' @param covariate_m_mat subject-specific covariates such as age or smoking status for male individual
# #' @param U.f_mat Binary nonzero measurement indicator for female individual
# #' @param U.m_mat Binary nonzero measurement indicator for male individual
# #' @param V.f_mat Nonzero measurement for female individual
# #' @param V.m_mat Nonzero measurement for male individual
#' @param seed_num The seed for the random number generator. If NA, the default seed 678443068 is used
#' @param K latent class number. Default is 3
#' @param alpha.0f Initial value for fixed effect coefficient of the latent class variable for female individual
#' @param alpha.1f Initial value for fixed effect coefficient of the latent class variable for female individual
#' @param alpha.0m Initial value for fixed effect coefficient of the latent class variable for male individual
#' @param alpha.1m Initial value for fixed effect coefficient of the latent class variable for male individual
#' @param b.0f Initial value for random effect coefficient of the latent class variable for female individual
#' @param b.1f Initial value for random effect coefficient of the latent class variable for female individual
#' @param b.0m Initial value for random effect coefficient of the latent class variable for male individual
#' @param b.1m Initial value for random effect coefficient of the latent class variable for male individual
#' @param beta.0 Initial value for Regression coefficients representing the association between the risk of Binary indicating variable for the couple and the latent class variables
#' @param beta.1f Initial value for Regression coefficients representing the association between the risk of Binary indicating variable for the couple and the latent class variables
#' @param beta.1m Initial value for Regression coefficients representing the association between the risk of Binary indicating variable for the couple and the latent class variables
#' @param beta.2 Initial value for Regression coefficients representing the association between the risk of Binary indicating variable for the couple and the latent class variables
#' @param eta.0f Initial value for coefficient for the association between the risk of Binary nonzero measurement indicator for female individual and Nonzero measurement on the log scale for female individual
#' @param eta.1f Initial value for coefficient for the association between the risk of Binary nonzero measurement indicator for female individual and Nonzero measurement on the log scale for female individual
#' @param eta.0m Initial value for coefficient for the association between the risk of Binary nonzero measurement indicator for male individual and Nonzero measurement on the log scale for male individual
#' @param eta.1m coefficient for the association between the risk of Binary nonzero measurement indicator for male individual and Nonzero measurement on the log scale for male individual
#' @param Sigma.b Initial value for Variance covariance matrix for Nonzero measurement specific shared random effects vector
#' @param tau2.f Initial value for variace of V.f_mat Nonzero measurement on the log scale for female individual
#' @param tau2.m variace of V.m_mat Nonzero measurement on the log scale for male individual
#' @param lambda.1f Initial value for Parameter vector for subject-specific covariates such as age, BMI or smoking status for female individual
#' @param lambda.2f Initial value for Parameter vector for subject-specific covariates such as age, BMI or smoking status for female individual
#' @param lambda.3f Initial value for Parameter vector for subject-specific covariates such as age, BMI or smoking status for female individual
#' @param lambda.4f Initial value for Parameter vector for subject-specific covariates such as age, BMI or smoking status for female individual
#' @param lambda.5f Initial value for Parameter vector for subject-specific covariates such as age, BMI or smoking status for female individual
#' @param lambda.1m Initial value for Parameter vector for subject-specific covariates such as age, BMI or smoking status for male individual
#' @param lambda.2m Initial value for Parameter vector for subject-specific covariates such as age, BMI or smoking status for male individual
#' @param lambda.3m Initial value for Parameter vector for subject-specific covariates such as age, BMI or smoking status for male individual
#' @param lambda.4m Initial value for Parameter vector for subject-specific covariates such as age, BMI or smoking status for male individual
#' @param lambda.5m Initial value for Parameter vector for subject-specific covariates such as age, BMI or smoking status for male individual
#' @param mu.b0 Initial value for prior distribution of beta0
#' @param mu.b1f Initial value for prior distribution of beta0 for female individual
#' @param mu.b1m Initial value for prior distribution of beta0 for male individual
#' @param mu.b2 Initial value for prior distribution of beta0
#' @param sig2.b0 Initial value for prior distribution of beta0
#' @param sig2.b1f Initial value for  prior distribution of beta0
#' @param sig2.b1m Initial value for  prior distribution of beta0
#' @param sig2.b2 Initial value for  prior distribution of beta0
#' @param mu.a0f Initial value for noninformative normal prior distribution of alpha for female individual
#' @param mu.a1f Initial value for noninformative normal prior distribution of alpha female individual
#' @param mu.a0m Initial value for noninformative normal prior distribution for male individual alpha
#' @param mu.a1m Initial value for noninformative normal prior distribution for male individual alpha
#' @param sig2.a0f Initial value for noninformative normal prior distribution for female individual alpha
#' @param sig2.a1f Initial value for noninformative normal prior distribution for female individual alpha
#' @param sig2.a0m Initial value for noninformative normal prior distribution for male individual alpha
#' @param sig2.a1m Initial value for noninformative normal prior distribution for male individual alpha
#' @param mu.e0f Initial value for noninformative normal prior distribution for female individual eta0
#' @param mu.e1f Initial value for noninformative normal prior distribution for female individual eta1
#' @param mu.e0m Initial value for noninformative normal prior distribution for male individual eta0
#' @param mu.e1m Initial value for noninformative normal prior distribution for male individual eta1
#' @param sig2.e0f Initial value for noninformative normal prior distribution for female individual eta0
#' @param sig2.e1f Initial value for noninformative normal prior distribution for female individual eta1
#' @param sig2.e0m Initial value for noninformative normal prior distribution for male individual eta0
#' @param sig2.e1m Initial value for noninformative normal prior distribution for male individual eta0
#' @param a.tau Initial value for prior distribution for inverse-Gamma distribution tau square
#' @param b.tau Initial value for prior distribution for inverse-Gamma distribution tau square
#' @param mu.lf Initial value for prior distribution for female individual lambda
#' @param mu.lm Initial value for prior distribution for male individual lambda
#' @param Sigma.lf Initial value for prior distribution for female individual lambda
#' @param Sigma.lm Initial value for prior distribution for male individual lambda
#' @param s.b0 Initial value for conditional posterior distribution of beta0 for adaptive Metropolis algorithm
#' @param s.b1f Initial value for conditional posterior distribution of female beta1 for adaptive Metropolis algorithm
#' @param s.b1m Initial value for conditional posterior distribution of male beta1 for adaptive Metropolis algorithm
#' @param s.b2 Initial value for conditional posterior distribution of beta2 for adaptive Metropolis algorithm
#' @param var.b0 Initial value for conditional posterior distribution of beta0 for adaptive Metropolis algorithm
#' @param var.b1f Initial value for conditional posterior distribution of female beta1 for adaptive Metropolis algorithm
#' @param var.b1m Initial value for conditional posterior distribution of male beta1 for adaptive Metropolis algorithm
#' @param var.b2 Initial value for conditional posterior distribution of beta2 for adaptive Metropolis algorithm
#' @param s.e0f Initial value for conditional posterior distribution of female eta0 for adaptive Metropolis algorithm
#' @param s.e1f Initial value for conditional posterior distribution of female eta1 for adaptive Metropolis algorithm
#' @param s.e0m Initial value for conditional posterior distribution of male eta0 for adaptive Metropolis algorithm
#' @param s.e1m Initial value for conditional posterior distribution of male eta1 for adaptive Metropolis algorithm
#' @param var.e0f Initial value for conditional posterior distribution of female eta0 for adaptive Metropolis algorithm
#' @param var.e1f Initial value for conditional posterior distribution of female eta1 for adaptive Metropolis algorithm
#' @param var.e0m Initial value for conditional posterior distribution of male eta0 for adaptive Metropolis algorithm
#' @param var.e1m Initial value for conditional posterior distribution of male eta1 for adaptive Metropolis algorithm
#' @param s.r Initial value for conditional posterior distribution of random effect b for adaptive Metropolis algorithm
#' @param cov.r Initial value for conditional posterior distribution of random effect b for adaptive Metropolis algorithm
#' @param s.lf Initial value for conditional posterior distribution of female lambda for adaptive Metropolis algorithm
#' @param s.lm Initial value for conditional posterior distribution of male lambda for adaptive Metropolis algorithm
#' @param eps Initial value for beta.0
#' @param Sigma.0 Initial value for Sigma.b
#' @param nu Initial value for Sigma.b
# #' @keywords
# #' @details
#'
# #' @export

#' @return A list of posterior estimates, latent class probability estimates and DIC
#' @references
#' Beom Seuk Hwang, Zhen Chen,  Germaine M. Buck Louis, and Paul S. Albert. (2018)
#' \emph{A Bayesian multi-dimensional couple-based latent risk model with an application to infertility}.
#' Biometrics, 75, 315--325. \url{https://doi.org/10.1111/biom.12972}
#' @examples
#' library(MCMCpack)
#' library(mvtnorm)
#' data(sampledata)
#' try1 <-lchemix:::JLchemix(Yvariable= sampledata[,1], X.f_mat = sampledata[,2:37],
#' X.m_mat = sampledata[,38:73], covariate_f_mat = sampledata[,74:78],
#' covariate_m_mat = sampledata[,79:83])


##########################################
################ 1st MCMC ################
##########################################
#  X.f_mat Data matrix for measurements of female individual, represented by U.f and V.f. #Only used for J
# pi.f Class membership probability for female individual
# pi.m Class membership probability for male individual
# alpha.f vector for alpha.0f and alpha.1f, fixed effects coefficients of the latent class variablesfor female individual

JLchemix<-function(nsim=700,
                nburn=500,
                Yvariable ,
                X.f_mat  ,
                X.m_mat  ,
                covariate_f_mat  ,
                covariate_m_mat  ,
#                U.f_mat = U.f,
#                U.m_mat = U.m,
#                V.f_mat = V.f,
#                V.m_mat = V.m,
                seed_num = 678443068,


                # maximum number of classes
                K=3,

                # 2-class model
                #pi.f=rep(1,K)/K,
                # 2-class model
                #pi.m=rep(1,K)/K ,

                alpha.0f=-6,
                alpha.1f=1,
                alpha.0m=-6,
                alpha.1m=1,


                beta.0=-4.5,
                beta.1f=0.5,
                beta.1m=0.5,
                beta.2=1,

                b.0f = c( -7.254321, -7.010908, -3.880426, -8.118554, -6.759431, -7.230387, -4.878213, -5.007503,
                          -6.856297, -4.489543, -8.068346, -5.106003, -6.452870, -4.495532, -5.517385, -6.268557,
                          -6.053185, -5.464252, -6.452734, -5.906192, -6.260919, -5.691590, -5.614410, -6.143081,
                          -5.047088, -5.301305, -5.913374, -5.756906, -6.054626, -6.873340, -6.849659, -7.004676,
                          -5.332403, -7.307853, -5.183989, -6.383119),
                b.0m = c(-5.079509, -7.272471, -6.521358, -6.128509, -5.696533, -6.374789, -5.777319, -4.826086,
                          -8.989660, -6.708051, -4.058636, -5.976921, -5.306304, -4.786394, -6.425503, -4.867691,
                          -6.091003, -7.162082, -5.652558, -6.368305, -4.130053, -6.062718, -6.272025, -6.096108,
                          -7.405642, -4.415025, -7.636100, -5.613704, -6.741973, -5.124894, -7.442003, -6.065161,
                          -6.127977, -6.397780, -5.039697, -5.912360),
                b.1f = c(0.06081535,0.27933553,  1.09567058,  1.30525030, -1.71363430,  3.56305321,  1.44386159,
                          1.11045081,  0.91064552, 2.33861095, 1.89400060,  1.73022487,  0.12834589,  0.98972580,
                          0.70672882,  1.82998947, 0.67922938, 0.91639977,  1.92834454,  1.01663356,  0.13329777,
                          1.94157555,  0.85399484,  1.57496884,  1.45862290,  1.78167318, -0.07774665,2.66127628,
                          -0.58454029,0.77320667, 0.97268288, 0.85628677,0.75621566,1.06990346, 1.14924001,
                          0.95091695),
                b.1m = c(1.68599796,  0.53185998,  1.14397009,  2.86619860,  4.08110658,  0.37164953,  1.54640353,
                          1.27457107,  3.80397148,  0.82511839,  2.99448637,  2.41212410,  1.47281944,  0.82053754,
                          2.01676277, -0.99610755, -0.03900866,  1.47895046,  1.15950092,  0.91964213,  1.17409250,
                          1.33819320,  1.36568913,  0.55896975,  1.41116326,  3.30927990,  1.51011954,  0.90399847,
                          0.40348049,  1.27435217,  1.53464794,  2.73497557,  1.91075938,  0.93572173,  2.02485773,
                          3.82461971),

                eta.0f=2,
                eta.1f=0.5,
                eta.0m=2,
                eta.1m=0.5,

                Sigma.b=diag(4) ,

                tau2.f=0.5,
                tau2.m=0.5,

#                L.f=sample(0:(K-1),I,replace=TRUE,prob=pi.f) ,
#                L.m=sample(0:(K-1),I,replace=TRUE,prob=pi.m) ,


                lambda.1f=0.5,
                lambda.2f=0,
                lambda.3f=0,
                lambda.4f=0,
                lambda.5f=0,
#                lambda.f=c(lambda.1f,lambda.2f,lambda.3f,lambda.4f,lambda.5f),
#                xlambda.f=covariate_f_mat%*%lambda.f,

                lambda.1m=0.5,
                lambda.2m=0,
                lambda.3m=0,
                lambda.4m=0,
                lambda.5m=0,

#                lambda.m=c(lambda.1m,lambda.2m,lambda.3m,lambda.4m,lambda.5m),
#                xlambda.m=covariate_m_mat%*%lambda.m,


#                mu.ij.f=matrix(1,I,J),
#                mu.ij.m=matrix(1,I,J),
                ### priors
                mu.b0=0,
                mu.b1f=0,
                mu.b1m=0,
                mu.b2=0,

                sig2.b0=100,
                sig2.b1f=1,
                sig2.b1m=1,
                sig2.b2=100,

                mu.a0f=0,
                mu.a1f=0,
                mu.a0m=0,
                mu.a1m=0,
#                mu.a=c(mu.a0f,mu.a1f,mu.a0m,mu.a1m),

                sig2.a0f=100,
                sig2.a1f=100,
                sig2.a0m=100,
                sig2.a1m=100,


                mu.e0f=0,
                mu.e1f=0,
                mu.e0m=0,
                mu.e1m=0,

                sig2.e0f=100,
                sig2.e1f=100,
                sig2.e0m=100,
                sig2.e1m=100,

                a.tau=1,
                b.tau=1,

                mu.lf=rep(0,covariate_num),
                mu.lm=rep(0,covariate_num),
                Sigma.lf=10*diag(covariate_num),
                Sigma.lm=10*diag(covariate_num),


#                iset=1:I,

                s.b0=2.4^2,
                s.b1f=2.4^2,
                s.b1m=2.4^2,
                s.b2=2.4^2,

                var.b0=1,
                var.b1f=1,
                var.b1m=1,
                var.b2=1,

                s.e0f=2.4^2,
                s.e1f=2.4^2,
                s.e0m=2.4^2,
                s.e1m=2.4^2,

                var.e0f=1,
                var.e1f=1,
                var.e0m=1,
                var.e1m=1,

                s.r=2.4^2,
                cov.r=diag(4),

                s.lf=2.4^2,
                s.lm=2.4^2,

                eps=0.01,
                nu=4,
                Sigma.0=diag(4)

){
  #ptm=proc.time()
  force(seed_num)
  set.seed(seed_num)
  #Create U matrix: if x !=0 then 1 else 0
  U.f = X.f_mat
  U.f[is.na(U.f)]=0
  U.f[U.f !=0]= 1
  U.m = X.m_mat
  U.m[is.na(U.m)]=0
  U.m[U.m !=0]= 1

    #Create V matrix: if x !=0 then x
  V.f = X.f_mat
  V.f[V.f ==0]= NA
  V.m = X.m_mat
  V.m[V.m ==0]= NA

#  force(I)
#  force(J)
  I=length(Yvariable)
  J=dim(X.f_mat)[2]
  covariate_num <- dim(covariate_m_mat)[2]
  #pi.f=rep(1,K)/K    # 2-class model
  #force(pi.f)
  #pi.m=rep(1,K)/K    # 2-class model
  #force(pi.m)
  # 2-class model
  pi.f=rep(1,K)/K
  # 2-class model
  pi.m=rep(1,K)/K

  alpha.f <- c(alpha.0f,alpha.1f)
  alpha.m <- c(alpha.0m,alpha.1m)
  alpha <- c(alpha.f,alpha.m)

  alpha.f=c(alpha.0f,alpha.1f)
  alpha.m=c(alpha.0m,alpha.1m)
  alpha=c(alpha.f,alpha.m)
  #This may problematic
  Sigma.b=diag(4)

  L.f=sample(0:(K-1),I,replace=TRUE,prob=pi.f)
  L.m=sample(0:(K-1),I,replace=TRUE,prob=pi.m)

  lambda.f=c(lambda.1f,lambda.2f,lambda.3f,lambda.4f,lambda.5f)
  xlambda.f=covariate_f_mat%*%lambda.f

  lambda.m=c(lambda.1m,lambda.2m,lambda.3m,lambda.4m,lambda.5m)
  xlambda.m=covariate_m_mat%*%lambda.m

  random.b=cbind(b.0f,b.1f,b.0m,b.1m)
  cov.lf=diag(covariate_num)
  cov.lm=diag(covariate_num)

#  force(mu.ij.f)     #=matrix(1,I,J)
#  force(mu.ij.m)     #=matrix(1,I,J)
  mu.ij.f=matrix(1,I,J)
  mu.ij.m=matrix(1,I,J)
  for (j in 1:J){
    mu.ij.f[,j]=b.0f[j]+b.1f[j]*L.f+xlambda.f
    mu.ij.m[,j]=b.0m[j]+b.1m[j]*L.m+xlambda.m
  }

  mu.ij.f1=mu.ij.f-mean(mu.ij.f)
  mu.ij.m1=mu.ij.m-mean(mu.ij.m)

  #force(mu.a)       #=c(mu.a0f,mu.a1f,mu.a0m,mu.a1m)
  mu.a=c(mu.a0f, mu.a1f, mu.a0m, mu.a1m)
  #force(Sigma.a)    #=diag(c(sig2.a0f,sig2.a1f,sig2.a0m,sig2.a1m),4)

  Sigma.a <- diag(c(sig2.a0f,sig2.a1f,sig2.a0m,sig2.a1m),4)

  ##############This seems questionable
  ##############This seems questionable
  iset=1:I
  #force(iset)  #=1:I
  #force(e.f)   #=rep(1,K)
  #force(e.m)   #=rep(1,K)

  #force(dff)     #=rep(0,K)
  #force(dmm)     #=rep(0,K)

  e.f=rep(1,K)
  e.m=rep(1,K)

  dff=rep(0,K)
  dmm=rep(0,K)

  ### output matrix
  pi.f_out=matrix(0.5,nsim,K)
  pi.m_out=matrix(0.5,nsim,K)

  beta.0_out=rep(1,nsim)
  beta.1f_out=rep(1,nsim)
  beta.1m_out=rep(1,nsim)
  beta.2_out=rep(1,nsim)

  alpha.0f_out=rep(1,nsim)
  alpha.1f_out=rep(1,nsim)
  alpha.0m_out=rep(1,nsim)
  alpha.1m_out=rep(1,nsim)

  eta.0f_out=rep(1,nsim)
  eta.1f_out=rep(1,nsim)
  eta.0m_out=rep(1,nsim)
  eta.1m_out=rep(1,nsim)

  b.0f_out=matrix(1,nsim,J)
  b.1f_out=matrix(1,nsim,J)
  b.0m_out=matrix(1,nsim,J)
  b.1m_out=matrix(1,nsim,J)

  Sigma.b_out=array(1,dim=c(4,4,nsim))

  sig.11_out=rep(1,nsim)
  sig.12_out=rep(1,nsim)
  sig.13_out=rep(1,nsim)
  sig.14_out=rep(1,nsim)
  sig.22_out=rep(1,nsim)
  sig.23_out=rep(1,nsim)
  sig.24_out=rep(1,nsim)
  sig.33_out=rep(1,nsim)
  sig.34_out=rep(1,nsim)
  sig.44_out=rep(1,nsim)

  tau2.f_out=rep(1,nsim)
  tau2.m_out=rep(1,nsim)

  lambda.f_out=matrix(1,nsim,covariate_num)
  lambda.m_out=matrix(1,nsim,covariate_num)

  L.f_out=matrix(0,nsim,K)
  L.m_out=matrix(0,nsim,K)

  loglike=rep(1,nsim)
  loglike.est=rep(1,nsim)


  ### traceplot matrix
  pi.f_T=matrix(0.5,(nsim-nburn)/10,K)
  pi.m_T=matrix(0.5,(nsim-nburn)/10,K)

  beta.0_T=rep(1,(nsim-nburn)/10)
  beta.1f_T=rep(1,(nsim-nburn)/10)
  beta.1m_T=rep(1,(nsim-nburn)/10)
  beta.2_T=rep(1,(nsim-nburn)/10)

  alpha.0f_T=rep(1,(nsim-nburn)/10)
  alpha.1f_T=rep(1,(nsim-nburn)/10)
  alpha.0m_T=rep(1,(nsim-nburn)/10)
  alpha.1m_T=rep(1,(nsim-nburn)/10)

  eta.0f_T=rep(1,(nsim-nburn)/10)
  eta.1f_T=rep(1,(nsim-nburn)/10)
  eta.0m_T=rep(1,(nsim-nburn)/10)
  eta.1m_T=rep(1,(nsim-nburn)/10)

  b.0f_T=matrix(1,(nsim-nburn)/10,J)
  b.1f_T=matrix(1,(nsim-nburn)/10,J)
  b.0m_T=matrix(1,(nsim-nburn)/10,J)
  b.1m_T=matrix(1,(nsim-nburn)/10,J)

  tau2.f_T=rep(1,(nsim-nburn)/10)
  tau2.m_T=rep(1,(nsim-nburn)/10)

  sig.11_T=rep(1,(nsim-nburn)/10)
  sig.12_T=rep(1,(nsim-nburn)/10)
  sig.13_T=rep(1,(nsim-nburn)/10)
  sig.14_T=rep(1,(nsim-nburn)/10)
  sig.22_T=rep(1,(nsim-nburn)/10)
  sig.23_T=rep(1,(nsim-nburn)/10)
  sig.24_T=rep(1,(nsim-nburn)/10)
  sig.33_T=rep(1,(nsim-nburn)/10)
  sig.34_T=rep(1,(nsim-nburn)/10)
  sig.44_T=rep(1,(nsim-nburn)/10)

  L.f_T=matrix(1,(nsim-nburn)/10,K)
  L.m_T=matrix(1,(nsim-nburn)/10,K)

  lambda.f_T=matrix(1,(nsim-nburn)/10,covariate_num)
  lambda.m_T=matrix(1,(nsim-nburn)/10,covariate_num)
  ## make Metropolis-Hastings function for beta.0
  MH_beta.0=function(beta.0.star,beta.0){

    Lambda.star=beta.0.star+beta.1f*L.f+beta.1m*L.m+beta.2*L.f*L.m
    Lambda=beta.0+beta.1f*L.f+beta.1m*L.m+beta.2*L.f*L.m

    part1=dnorm(beta.0.star,mu.b0,sqrt(sig2.b0))/dnorm(beta.0,mu.b0,sqrt(sig2.b0))
    part2=((exp(Lambda.star)/(1+exp(Lambda.star)))/(exp(Lambda)/(1+exp(Lambda))))^Yvariable
    part3=((1+exp(Lambda))/(1+exp(Lambda.star)))^(1-Yvariable)

    fn=part1*prod(part2*part3)

    return(fn)
  }

  ## make Metropolis-Hastings function for beta.1f
  MH_beta.1f=function(beta.1f.star,beta.1f){

    Lambda.star=beta.0+beta.1f.star*L.f+beta.1m*L.m+beta.2*L.f*L.m
    Lambda=beta.0+beta.1f*L.f+beta.1m*L.m+beta.2*L.f*L.m

    part1=dlnorm(beta.1f.star,mu.b1f,sqrt(sig2.b1f))/dlnorm(beta.1f,mu.b1f,sqrt(sig2.b1f))
    part2=((exp(Lambda.star)/(1+exp(Lambda.star)))/(exp(Lambda)/(1+exp(Lambda))))^Yvariable
    part3=((1+exp(Lambda))/(1+exp(Lambda.star)))^(1-Yvariable)

    fn=part1*prod(part2*part3)

    return(fn)
  }

  ## make Metropolis-Hastings function for beta.1m
  MH_beta.1m=function(beta.1m.star,beta.1m){

    Lambda.star=beta.0+beta.1f*L.f+beta.1m.star*L.m+beta.2*L.f*L.m
    Lambda=beta.0+beta.1f*L.f+beta.1m*L.m+beta.2*L.f*L.m

    part1=dlnorm(beta.1m.star,mu.b1m,sqrt(sig2.b1m))/dlnorm(beta.1m,mu.b1m,sqrt(sig2.b1m))
    part2=((exp(Lambda.star)/(1+exp(Lambda.star)))/(exp(Lambda)/(1+exp(Lambda))))^Yvariable
    part3=((1+exp(Lambda))/(1+exp(Lambda.star)))^(1-Yvariable)

    fn=part1*prod(part2*part3)

    return(fn)
  }

  ## make Metropolis-Hastings function for beta.2
  MH_beta.2=function(beta.2.star,beta.2){

    Lambda.star=beta.0+beta.1f*L.f+beta.1m*L.m+beta.2.star*L.f*L.m
    Lambda=beta.0+beta.1f*L.f+beta.1m*L.m+beta.2*L.f*L.m

    part1=dnorm(beta.2.star,mu.b2,sqrt(sig2.b2))/dnorm(beta.2,mu.b2,sqrt(sig2.b2))
    part2=((exp(Lambda.star)/(1+exp(Lambda.star)))/(exp(Lambda)/(1+exp(Lambda))))^Yvariable
    part3=((1+exp(Lambda))/(1+exp(Lambda.star)))^(1-Yvariable)

    fn=part1*prod(part2*part3)

    return(fn)
  }


  ## make Metropolis-Hastings function for b
  MH_b=function(U.f,U.m,V.f,V.m,b.0f.star,b.1f.star,b.0m.star,b.1m.star,b.0f,b.1f,b.0m,b.1m){   ## for specific j

    mu.j.f.star=b.0f.star+b.1f.star*L.f+xlambda.f
    mu.ij.f.star=mu.ij.f
    mu.ij.f.star[,j]=mu.j.f.star
    mu.j.f.star1=mu.j.f.star-mean(mu.ij.f.star)
    mu.j.f=b.0f+b.1f*L.f+xlambda.f
    mu.j.f1=mu.j.f-mean(mu.ij.f)

    mu.j.m.star=b.0m.star+b.1m.star*L.m+xlambda.m
    mu.ij.m.star=mu.ij.m
    mu.ij.m.star[,j]=mu.j.m.star
    mu.j.m.star1=mu.j.m.star-mean(mu.ij.m.star)
    mu.j.m=b.0m+b.1m*L.m+xlambda.m
    mu.j.m1=mu.j.m-mean(mu.ij.m)

    part1=dmvnorm(c(b.0f.star,b.1f.star,b.0m.star,b.1m.star),alpha,Sigma.b)/dmvnorm(c(b.0f,b.1f,b.0m,b.1m),alpha,Sigma.b)
    part2=((exp(eta.0f+eta.1f*mu.j.f.star1)/(1+exp(eta.0f+eta.1f*mu.j.f.star1)))/(exp(eta.0f+eta.1f*mu.j.f1)/(1+exp(eta.0f+eta.1f*mu.j.f1))))^U.f
    part3=((1+exp(eta.0f+eta.1f*mu.j.f1))/(1+exp(eta.0f+eta.1f*mu.j.f.star1)))^(1-U.f)
    part4=((exp(eta.0m+eta.1m*mu.j.m.star1)/(1+exp(eta.0m+eta.1m*mu.j.m.star1)))/(exp(eta.0m+eta.1m*mu.j.m1)/(1+exp(eta.0m+eta.1m*mu.j.m1))))^U.m
    part5=((1+exp(eta.0m+eta.1m*mu.j.m1))/(1+exp(eta.0m+eta.1m*mu.j.m.star1)))^(1-U.m)
    part6=(dlnorm(V.f,meanlog=mu.j.f.star,sdlog=sqrt(tau2.f))/dlnorm(V.f,meanlog=mu.j.f,sdlog=sqrt(tau2.f)))^U.f
    part7=(dlnorm(V.m,meanlog=mu.j.m.star,sdlog=sqrt(tau2.m))/dlnorm(V.m,meanlog=mu.j.m,sdlog=sqrt(tau2.m)))^U.m

    fn=part1*prod(part2*part3*part4*part5*part6*part7,na.rm=TRUE)

    return(fn)
  }

  ## make Metropolis-Hastings function for b for DIC
  MH_b.est=function(U.f,U.m,V.f,V.m,b.0f.star,b.1f.star,b.0m.star,b.1m.star,b.0f,b.1f,b.0m,b.1m){   ## for specific j

    mu.j.f.star=b.0f.star+b.1f.star*L.f+xlambda.f.est
    mu.ij.f.star=mu.ij.f
    mu.ij.f.star[,j]=mu.j.f.star
    mu.j.f.star1=mu.j.f.star-mean(mu.ij.f.star)
    mu.j.f=b.0f+b.1f*L.f+xlambda.f.est
    mu.j.f1=mu.j.f-mean(mu.ij.f)

    mu.j.m.star=b.0m.star+b.1m.star*L.m+xlambda.m.est
    mu.ij.m.star=mu.ij.m
    mu.ij.m.star[,j]=mu.j.m.star
    mu.j.m.star1=mu.j.m.star-mean(mu.ij.m.star)
    mu.j.m=b.0m+b.1m*L.m+xlambda.m.est
    mu.j.m1=mu.j.m-mean(mu.ij.m)

    part1=dmvnorm(c(b.0f.star,b.1f.star,b.0m.star,b.1m.star),alpha.est,Sigma.b.est)/dmvnorm(c(b.0f,b.1f,b.0m,b.1m),alpha.est,Sigma.b.est)
    part2=((exp(eta.0f.est+eta.1f.est*mu.j.f.star1)/(1+exp(eta.0f.est+eta.1f.est*mu.j.f.star1)))/(exp(eta.0f.est+eta.1f.est*mu.j.f1)/(1+exp(eta.0f.est+eta.1f.est*mu.j.f1))))^U.f
    part3=((1+exp(eta.0f.est+eta.1f.est*mu.j.f1))/(1+exp(eta.0f.est+eta.1f.est*mu.j.f.star1)))^(1-U.f)
    part4=((exp(eta.0m.est+eta.1m.est*mu.j.m.star1)/(1+exp(eta.0m.est+eta.1m.est*mu.j.m.star1)))/(exp(eta.0m.est+eta.1m.est*mu.j.m1)/(1+exp(eta.0m.est+eta.1m.est*mu.j.m1))))^U.m
    part5=((1+exp(eta.0m.est+eta.1m.est*mu.j.m1))/(1+exp(eta.0m.est+eta.1m.est*mu.j.m.star1)))^(1-U.m)
    part6=(dlnorm(V.f,meanlog=mu.j.f.star,sdlog=sqrt(tau2.f.est))/dlnorm(V.f,meanlog=mu.j.f,sdlog=sqrt(tau2.f.est)))^U.f
    part7=(dlnorm(V.m,meanlog=mu.j.m.star,sdlog=sqrt(tau2.m.est))/dlnorm(V.m,meanlog=mu.j.m,sdlog=sqrt(tau2.m.est)))^U.m

    fn=part1*prod(part2*part3*part4*part5*part6*part7,na.rm=TRUE)

    return(fn)
  }

  ## make Metropolis-Hastings function for eta.0f
  MH_eta.0f=function(eta.0f.star,eta.0f){

    part1=dnorm(eta.0f.star,mu.e0f,sqrt(sig2.e0f))/dnorm(eta.0f,mu.e0f,sqrt(sig2.e0f))
    part2=((exp(eta.0f.star+eta.1f*mu.ij.f1)/(1+exp(eta.0f.star+eta.1f*mu.ij.f1)))/(exp(eta.0f+eta.1f*mu.ij.f1)/(1+exp(eta.0f+eta.1f*mu.ij.f1))))^U.f
    part3=((1+exp(eta.0f+eta.1f*mu.ij.f1))/(1+exp(eta.0f.star+eta.1f*mu.ij.f1)))^(1-U.f)

    fn=part1*prod(part2*part3,na.rm=TRUE)

    return(fn)
  }

  ## make Metropolis-Hastings function for eta.1f
  MH_eta.1f=function(eta.1f.star,eta.1f){

    part1=dnorm(eta.1f.star,mu.e1f,sqrt(sig2.e1f))/dnorm(eta.1f,mu.e1f,sqrt(sig2.e1f))
    part2=((exp(eta.0f+eta.1f.star*mu.ij.f1)/(1+exp(eta.0f+eta.1f.star*mu.ij.f1)))/(exp(eta.0f+eta.1f*mu.ij.f1)/(1+exp(eta.0f+eta.1f*mu.ij.f1))))^U.f
    part3=((1+exp(eta.0f+eta.1f*mu.ij.f1))/(1+exp(eta.0f+eta.1f.star*mu.ij.f1)))^(1-U.f)

    fn=part1*prod(part2*part3,na.rm=TRUE)

    return(fn)
  }

  ## make Metropolis-Hastings function for eta.0m
  MH_eta.0m=function(eta.0m.star,eta.0m){

    part1=dnorm(eta.0m.star,mu.e0m,sqrt(sig2.e0m))/dnorm(eta.0m,mu.e0m,sqrt(sig2.e0m))
    part2=((exp(eta.0m.star+eta.1m*mu.ij.m1)/(1+exp(eta.0m.star+eta.1m*mu.ij.m1)))/(exp(eta.0m+eta.1m*mu.ij.m1)/(1+exp(eta.0m+eta.1m*mu.ij.m1))))^U.m
    part3=((1+exp(eta.0m+eta.1m*mu.ij.m1))/(1+exp(eta.0m.star+eta.1m*mu.ij.m1)))^(1-U.m)

    fn=part1*prod(part2*part3,na.rm=TRUE)

    return(fn)
  }

  ## make Metropolis-Hastings function for eta.1m
  MH_eta.1m=function(eta.1m.star,eta.1m){

    part1=dnorm(eta.1m.star,mu.e1m,sqrt(sig2.e1m))/dnorm(eta.1m,mu.e1m,sqrt(sig2.e1m))
    part2=((exp(eta.0m+eta.1m.star*mu.ij.m1)/(1+exp(eta.0m+eta.1m.star*mu.ij.m1)))/(exp(eta.0m+eta.1m*mu.ij.m1)/(1+exp(eta.0m+eta.1m*mu.ij.m1))))^U.m
    part3=((1+exp(eta.0m+eta.1m*mu.ij.m1))/(1+exp(eta.0m+eta.1m.star*mu.ij.m1)))^(1-U.m)

    fn=part1*prod(part2*part3,na.rm=TRUE)

    return(fn)
  }

  ## make Metropolis-Hastings function for lambda.f
  MH_lambda.f=function(lambda.f.star,lambda.f){

    xlambda.f.star=covariate_f_mat%*%lambda.f.star

    mu.ij.f.star=matrix(1,I,J)
    for (j in 1:J){
      mu.ij.f.star[,j]=b.0f[j]+b.1f[j]*L.f+xlambda.f.star
    }
    mu.ij.f.star1=mu.ij.f.star-mean(mu.ij.f.star)

    part1=dmvnorm(lambda.f.star,mu.lf,Sigma.lf)/dmvnorm(lambda.f,mu.lf,Sigma.lf)
    part2=((exp(eta.0f+eta.1f*mu.ij.f.star1)/(1+exp(eta.0f+eta.1f*mu.ij.f.star1)))/(exp(eta.0f+eta.1f*mu.ij.f1)/(1+exp(eta.0f+eta.1f*mu.ij.f1))))^U.f
    part3=((1+exp(eta.0f+eta.1f*mu.ij.f1))/(1+exp(eta.0f+eta.1f*mu.ij.f.star1)))^(1-U.f)
    part4=(dlnorm(V.f,meanlog=mu.ij.f.star,sdlog=sqrt(tau2.f))/dlnorm(V.f,meanlog=mu.ij.f,sdlog=sqrt(tau2.f)))^U.f

    fn=part1*prod(part2*part3*part4,na.rm=TRUE)

    return(fn)
  }

  ## make Metropolis-Hastings function for lambda.m
  MH_lambda.m=function(lambda.m.star,lambda.m){

    xlambda.m.star=covariate_m_mat%*%lambda.m.star

    mu.ij.m.star=matrix(1,I,J)
    for (j in 1:J){
      mu.ij.m.star[,j]=b.0m[j]+b.1m[j]*L.m+xlambda.m.star
    }
    mu.ij.m.star1=mu.ij.m.star-mean(mu.ij.m.star)

    part1=dmvnorm(lambda.m.star,mu.lm,Sigma.lm)/dmvnorm(lambda.m,mu.lm,Sigma.lm)
    part2=((exp(eta.0m+eta.1m*mu.ij.m.star1)/(1+exp(eta.0m+eta.1m*mu.ij.m.star1)))/(exp(eta.0m+eta.1m*mu.ij.m1)/(1+exp(eta.0m+eta.1m*mu.ij.m1))))^U.m
    part3=((1+exp(eta.0m+eta.1m*mu.ij.m1))/(1+exp(eta.0m+eta.1m*mu.ij.m.star1)))^(1-U.m)
    part4=(dlnorm(V.m,meanlog=mu.ij.m.star,sdlog=sqrt(tau2.m))/dlnorm(V.m,meanlog=mu.ij.m,sdlog=sqrt(tau2.m)))^U.m

    fn=part1*prod(part2*part3*part4,na.rm=TRUE)

    return(fn)
  }


  ##########################################
  ################ 1st MCMC ################
  ##########################################

  ptm=proc.time()


  for (gt in 1:nsim){

    ##################
    #### update beta
    ##################
    ## update beta.0
    if (gt>100){
      var.b0=var(beta.0_out[1:(gt-1)])
    }
    if (gt>1 & ((gt-1)%%100==0)){
      ar.b0=sum(diff(beta.0_out[(gt-100):(gt-1)])!=0)/99
      if (ar.b0<0.44){s.b0=s.b0/sqrt(2)}
      else if (ar.b0>0.44){s.b0=s.b0*sqrt(2)}
    }

    u.b0=runif(1)
    beta.0.star=rnorm(1,beta.0,sqrt(s.b0*var.b0+s.b0*eps))
    R.b0=MH_beta.0(beta.0.star,beta.0)

    if (u.b0<=R.b0){beta.0=beta.0.star}

    ## update beta.1f
    if (gt>100){
      var.b1f=var(beta.1f_out[1:(gt-1)])
    }
    if (gt>1 & ((gt-1)%%100==0)){
      ar.b1f=sum(diff(beta.1f_out[(gt-100):(gt-1)])!=0)/99
      if (ar.b1f<0.44){s.b1f=s.b1f/sqrt(2)}
      else if (ar.b1f>0.44){s.b1f=s.b1f*sqrt(2)}
    }

    repeat{
      beta.1f.star=rnorm(1,beta.1f,sqrt(s.b1f*var.b1f+s.b1f*eps))
      if (beta.1f.star>0)
        break
    }
    u.b1f=runif(1)
    R.b1f=MH_beta.1f(beta.1f.star,beta.1f)*pnorm(beta.1f,0,sqrt(s.b1f*var.b1f+s.b1f*eps))/pnorm(beta.1f.star,0,sqrt(s.b1f*var.b1f+s.b1f*eps))
    if (u.b1f<=R.b1f){beta.1f=beta.1f.star}

    ## update beta.1m
    if (gt>100){
      var.b1m=var(beta.1m_out[1:(gt-1)])
    }
    if (gt>1 & ((gt-1)%%100==0)){
      ar.b1m=sum(diff(beta.1m_out[(gt-100):(gt-1)])!=0)/99
      if (ar.b1m<0.44){s.b1m=s.b1m/sqrt(2)}
      else if (ar.b1m>0.44){s.b1m=s.b1m*sqrt(2)}
    }

    repeat{
      beta.1m.star=rnorm(1,beta.1m,sqrt(s.b1m*var.b1m+s.b1m*eps))
      if (beta.1m.star>0)
        break
    }
    u.b1m=runif(1)
    R.b1m=MH_beta.1m(beta.1m.star,beta.1m)*pnorm(beta.1m,0,sqrt(s.b1m*var.b1m+s.b1m*eps))/pnorm(beta.1m.star,0,sqrt(s.b1m*var.b1m+s.b1m*eps))
    if (u.b1m<=R.b1m){beta.1m=beta.1m.star}

    ## update beta.2
    if (gt>100){
      var.b2=var(beta.2_out[1:(gt-1)])
    }
    if (gt>1 & ((gt-1)%%100==0)){
      ar.b2=sum(diff(beta.2_out[(gt-100):(gt-1)])!=0)/99
      if (ar.b2<0.44){s.b2=s.b2/sqrt(2)}
      else if (ar.b2>0.44){s.b2=s.b2*sqrt(2)}
    }

    u.b2=runif(1)
    beta.2.star=rnorm(1,beta.2,sqrt(s.b2*var.b2+s.b2*eps))
    R.b2=MH_beta.2(beta.2.star,beta.2)

    if (u.b2<=R.b2){beta.2=beta.2.star}


    ##################
    #### update alpha
    ##################
    ## update alpha.f
    v2=solve(solve(Sigma.a)+J*solve(Sigma.b))
    v1=v2%*%(solve(Sigma.a)%*%mu.a+solve(Sigma.b)%*%colSums(random.b))

    alpha=rmvnorm(1,v1,v2)
    alpha.0f=alpha[,1]
    alpha.1f=alpha[,2]
    alpha.0m=alpha[,3]
    alpha.1m=alpha[,4]


    ####################
    #### update Sigma.b
    ####################

    A=random.b-rep(1,J)%*%alpha
    Sigma.b=riwish(J+nu,t(A)%*%A+Sigma.0)

    sig.11=Sigma.b[1,1]
    sig.12=Sigma.b[1,2]
    sig.13=Sigma.b[1,3]
    sig.14=Sigma.b[1,4]
    sig.22=Sigma.b[2,2]
    sig.23=Sigma.b[2,3]
    sig.24=Sigma.b[2,4]
    sig.33=Sigma.b[3,3]
    sig.34=Sigma.b[3,4]
    sig.44=Sigma.b[4,4]

    ####################
    #### update tau2
    ####################
    g1.f=a.tau+0.5*sum(U.f,na.rm=TRUE)
    g2.f=b.tau+0.5*sum((log(V.f[U.f==1])-mu.ij.f[U.f==1])^2,na.rm=TRUE)

    g1.m=a.tau+0.5*sum(U.m,na.rm=TRUE)
    g2.m=b.tau+0.5*sum((log(V.m[U.m==1])-mu.ij.m[U.m==1])^2,na.rm=TRUE)

    tau2.f=rinvgamma(1,shape=g1.f,scale=g2.f)
    tau2.m=rinvgamma(1,shape=g1.m,scale=g2.m)


    ##################
    #### update eta
    ##################
    ## update eta.0f
    if (gt>100){
      var.e0f=var(eta.0f_out[1:(gt-1)])
    }
    if (gt>1 & ((gt-1)%%100==0)){
      ar.e0f=sum(diff(eta.0f_out[(gt-100):(gt-1)])!=0)/99
      if (ar.e0f<0.44){s.e0f=s.e0f/sqrt(2)}
      else if (ar.e0f>0.44){s.e0f=s.e0f*sqrt(2)}
    }

    u.e0f=runif(1)
    eta.0f.star=rnorm(1,eta.0f,sqrt(s.e0f*var.e0f+s.e0f*eps))
    R.e0f=MH_eta.0f(eta.0f.star,eta.0f)

    if (u.e0f<=R.e0f){eta.0f=eta.0f.star}

    ## update eta.1f
    if (gt>100){
      var.e1f=var(eta.1f_out[1:(gt-1)])
    }
    if (gt>1 & ((gt-1)%%100==0)){
      ar.e1f=sum(diff(eta.1f_out[(gt-100):(gt-1)])!=0)/99
      if (ar.e1f<0.44){s.e1f=s.e1f/sqrt(2)}
      else if (ar.e1f>0.44){s.e1f=s.e1f*sqrt(2)}
    }

    u.e1f=runif(1)
    eta.1f.star=rnorm(1,eta.1f,sqrt(s.e1f*var.e1f+s.e1f*eps))
    R.e1f=MH_eta.1f(eta.1f.star,eta.1f)

    if (u.e1f<=R.e1f){eta.1f=eta.1f.star}

    ## update eta.0m
    if (gt>100){
      var.e0m=var(eta.0m_out[1:(gt-1)])
    }
    if (gt>1 & ((gt-1)%%100==0)){
      ar.e0m=sum(diff(eta.0m_out[(gt-100):(gt-1)])!=0)/99
      if (ar.e0m<0.44){s.e0m=s.e0m/sqrt(2)}
      else if (ar.e0m>0.44){s.e0m=s.e0m*sqrt(2)}
    }

    u.e0m=runif(1)
    eta.0m.star=rnorm(1,eta.0m,sqrt(s.e0m*var.e0m+s.e0m*eps))
    R.e0m=MH_eta.0m(eta.0m.star,eta.0m)

    if (u.e0m<=R.e0m){eta.0m=eta.0m.star}

    ## update eta.1m
    if (gt>100){
      var.e1m=var(eta.1m_out[1:(gt-1)])
    }
    if (gt>1 & ((gt-1)%%100==0)){
      ar.e1m=sum(diff(eta.1m_out[(gt-100):(gt-1)])!=0)/99
      if (ar.e1m<0.44){s.e1m=s.e1m/sqrt(2)}
      else if (ar.e1m>0.44){s.e1m=s.e1m*sqrt(2)}
    }

    u.e1m=runif(1)
    eta.1m.star=rnorm(1,eta.1m,sqrt(s.e1m*var.e1m+s.e1m*eps))
    R.e1m=MH_eta.1m(eta.1m.star,eta.1m)

    if (u.e1m<=R.e1m){eta.1m=eta.1m.star}


    ####################
    #### update lambda
    ####################
    ######## Female
    if (gt>100){
      cov.lf=cov(lambda.f_out[1:(gt-1),])
    }
    if (gt>1 & ((gt-1)%%100==0)){
      ar.lf=sum(diff(lambda.f_out[(gt-100):(gt-1),1])!=0)/99
      if (ar.lf<0.23){s.lf=s.lf/sqrt(2)}
      else if (ar.lf>0.23){s.lf=s.lf*sqrt(2)}
    }

    u.lf=runif(1)
    lambda.f.star=as.vector(rmvnorm(1,lambda.f,s.lf*cov.lf+s.lf*eps*diag(covariate_num)))
    R.lf=MH_lambda.f(lambda.f.star,lambda.f)

    if (u.lf<=R.lf){lambda.f=lambda.f.star}

    lambda.f[2:covariate_num]=0
    xlambda.f=covariate_f_mat%*%lambda.f

    ######## Male
    if (gt>100){
      cov.lm=cov(lambda.m_out[1:(gt-1),])
    }
    if (gt>1 & ((gt-1)%%100==0)){
      ar.lm=sum(diff(lambda.m_out[(gt-100):(gt-1),1])!=0)/99
      if (ar.lm<0.23){s.lm=s.lm/sqrt(2)}
      else if (ar.lm>0.23){s.lm=s.lm*sqrt(2)}
    }

    u.lm=runif(1)
    lambda.m.star=as.vector(rmvnorm(1,lambda.m,s.lm*cov.lm+s.lm*eps*diag(covariate_num)))
    R.lm=MH_lambda.m(lambda.m.star,lambda.m)

    if (u.lm<=R.lm){lambda.m=lambda.m.star}

    lambda.m[2:covariate_num]=0
    xlambda.m=covariate_m_mat%*%lambda.m


    #########################
    #### update latent class
    #########################
    ## update pi for female and male
    for (k in 1:K){
      dff[k]=sum(L.f==(k-1))+e.f[k]
      dmm[k]=sum(L.m==(k-1))+e.m[k]
    }
    pi.f=rdirichlet(1,dff)
    pi.m=rdirichlet(1,dmm)

    ## update latent class
    mu.ij.fL=matrix(1,I,J)
    mu.ij.mL=matrix(1,I,J)
    pr.latent.f=matrix(0,I,K)
    pr.latent.m=matrix(0,I,K)

    for (k in 1:K){
      ### for female
      Lamb.fL=beta.0+beta.1f*(k-1)+beta.1m*L.m+beta.2*(k-1)*L.m
      for (j in 1:J){
        mu.ij.fL[,j]=b.0f[j]+b.1f[j]*(k-1)+xlambda.f
      }
      mu.ij.fL1=mu.ij.fL-mean(mu.ij.fL)

      part1.f=(exp(Lamb.fL)/(1+exp(Lamb.fL)))^Yvariable*(1/(1+exp(Lamb.fL)))^(1-Yvariable)

      temp.f=(exp(eta.0f+eta.1f*mu.ij.fL1)/(1+exp(eta.0f+eta.1f*mu.ij.fL1)))^U.f*(1/(1+exp(eta.0f+eta.1f*mu.ij.fL1)))^(1-U.f)*
        (dlnorm(V.f,meanlog=mu.ij.fL,sdlog=sqrt(tau2.f)))^U.f
      part2.f=apply(temp.f,1,prod,na.rm=TRUE)

      pr.latent.f[,k]=pi.f[k]*part1.f*part2.f

      ### for male
      Lamb.mL=beta.0+beta.1f*L.f+beta.1m*(k-1)+beta.2*L.f*(k-1)
      for (j in 1:J){
        mu.ij.mL[,j]=b.0m[j]+b.1m[j]*(k-1)+xlambda.m
      }
      mu.ij.mL1=mu.ij.mL-mean(mu.ij.mL)

      part1.m=(exp(Lamb.mL)/(1+exp(Lamb.mL)))^Yvariable*(1/(1+exp(Lamb.mL)))^(1-Yvariable)

      temp.m=(exp(eta.0m+eta.1m*mu.ij.mL1)/(1+exp(eta.0m+eta.1m*mu.ij.mL1)))^U.m*(1/(1+exp(eta.0m+eta.1m*mu.ij.mL1)))^(1-U.m)*
        (dlnorm(V.m,meanlog=mu.ij.mL,sdlog=sqrt(tau2.m)))^U.m
      part2.m=apply(temp.m,1,prod,na.rm=TRUE)

      pr.latent.m[,k]=pi.m[k]*part1.m*part2.m
    }

    pr.latent.f=pr.latent.f/rowSums(pr.latent.f)
    pr.latent.m=pr.latent.m/rowSums(pr.latent.m)
    for (i in 1:I){
      L.f[i]=sample(0:(K-1),1,replace=TRUE,prob=pr.latent.f[i,])
      L.m[i]=sample(0:(K-1),1,replace=TRUE,prob=pr.latent.m[i,])
    }


    ############################
    #### update random effect b
    ############################
    for (j in 1:J){
      if (gt>100){
        cov.r=cov(cbind(b.0f_out[1:(gt-1),j],b.1f_out[1:(gt-1),j],b.0m_out[1:(gt-1),j],b.1m_out[1:(gt-1),j]))
      }
      if (gt>1 & ((gt-1)%%100==0)){
        ar.r=sum(diff(b.0f_out[(gt-100):(gt-1),j])!=0)/99
        if (ar.r<0.44){s.r=s.r/sqrt(2)}
        else if (ar.r>0.44){s.r=s.r*sqrt(2)}
      }

      u.r=runif(1)
      b.star=rmvnorm(1,c(b.0f[j],b.1f[j],b.0m[j],b.1m[j]),s.r*cov.r+s.r*eps*diag(4))
      b.0f.star=b.star[1]
      b.1f.star=b.star[2]
      b.0m.star=b.star[3]
      b.1m.star=b.star[4]

      R.r=MH_b(U.f[,j],U.m[,j],V.f[,j],V.m[,j],b.0f.star,b.1f.star,b.0m.star,b.1m.star,b.0f[j],b.1f[j],b.0m[j],b.1m[j])

      if (u.r<=R.r){b.0f[j]=b.0f.star;b.1f[j]=b.1f.star;b.0m[j]=b.0m.star;b.1m[j]=b.1m.star}
    }
    random.b=cbind(b.0f,b.1f,b.0m,b.1m)


    ###################################
    #### update mu.ij.f and mu.ij.m
    ###################################
    for (j in 1:J){
      mu.ij.f[,j]=b.0f[j]+b.1f[j]*L.f+xlambda.f
      mu.ij.m[,j]=b.0m[j]+b.1m[j]*L.m+xlambda.m
    }
    mu.ij.f1=mu.ij.f-mean(mu.ij.f)
    mu.ij.m1=mu.ij.m-mean(mu.ij.m)


    ########################################
    #### calculate 1st term of DIC
    ########################################
    Lamb=beta.0+beta.1f*L.f+beta.1m*L.m+beta.2*L.f*L.m

    pt1=sum(Yvariable*log((exp(Lamb)/(1+exp(Lamb))))-(1-Yvariable)*log(1+exp(Lamb)))
    pt2=sum(log(pi.f[(L.f+1)])+log(pi.m[(L.m+1)]))
    pt3=sum(U.f*log(exp(eta.0f+eta.1f*mu.ij.f1)/(1+exp(eta.0f+eta.1f*mu.ij.f1)))-(1-U.f)*log(1+exp(eta.0f+eta.1f*mu.ij.f1))+
              U.m*log(exp(eta.0m+eta.1m*mu.ij.m1)/(1+exp(eta.0m+eta.1m*mu.ij.m1)))-(1-U.m)*log(1+exp(eta.0m+eta.1m*mu.ij.m1)),na.rm=TRUE)
    pt4=sum(log(dlnorm(V.f,meanlog=mu.ij.f,sdlog=sqrt(tau2.f))),na.rm=TRUE)+
      sum(log(dlnorm(V.m,meanlog=mu.ij.m,sdlog=sqrt(tau2.m))),na.rm=TRUE)

    ptt=rep(0,J)
    for (j in 1:J){
      ptt[j]=log(dmvnorm(c(b.0f[j],b.1f[j],b.0m[j],b.1m[j]),alpha,Sigma.b))
    }
    pt5=sum(ptt)

    loglike[gt]=pt1+pt2+pt3+pt4+pt5


    #####################
    ## make output files
    #####################
    pi.f_out[gt,]=pi.f
    pi.m_out[gt,]=pi.m

    beta.0_out[gt]=beta.0
    beta.1f_out[gt]=beta.1f
    beta.1m_out[gt]=beta.1m
    beta.2_out[gt]=beta.2

    alpha.0f_out[gt]=alpha.0f
    alpha.1f_out[gt]=alpha.1f
    alpha.0m_out[gt]=alpha.0m
    alpha.1m_out[gt]=alpha.1m

    eta.0f_out[gt]=eta.0f
    eta.1f_out[gt]=eta.1f
    eta.0m_out[gt]=eta.0m
    eta.1m_out[gt]=eta.1m

    b.0f_out[gt,]=b.0f
    b.1f_out[gt,]=b.1f
    b.0m_out[gt,]=b.0m
    b.1m_out[gt,]=b.1m

    Sigma.b_out[,,gt]=Sigma.b

    sig.11_out[gt]=sig.11
    sig.12_out[gt]=sig.12
    sig.13_out[gt]=sig.13
    sig.14_out[gt]=sig.14
    sig.22_out[gt]=sig.22
    sig.23_out[gt]=sig.23
    sig.24_out[gt]=sig.24
    sig.33_out[gt]=sig.33
    sig.34_out[gt]=sig.34
    sig.44_out[gt]=sig.44

    tau2.f_out[gt]=tau2.f
    tau2.m_out[gt]=tau2.m

    lambda.f_out[gt,]=lambda.f
    lambda.m_out[gt,]=lambda.m

    for (k in 1:K){
      L.f_out[gt,k]=sum(L.f==(k-1))
      L.m_out[gt,k]=sum(L.m==(k-1))
    }
  }  ###### 1st MCMC Done!



  ## save data at every 10th iteration for traceplot
  for (t in 1:((nsim-nburn)/10)){
    pi.f_T[t,]=pi.f_out[(nburn+10*t),]
    pi.m_T[t,]=pi.m_out[(nburn+10*t),]

    beta.0_T[t]=beta.0_out[(nburn+10*t)]
    beta.1f_T[t]=beta.1f_out[(nburn+10*t)]
    beta.1m_T[t]=beta.1m_out[(nburn+10*t)]
    beta.2_T[t]=beta.2_out[(nburn+10*t)]

    alpha.0f_T[t]=alpha.0f_out[(nburn+10*t)]
    alpha.1f_T[t]=alpha.1f_out[(nburn+10*t)]
    alpha.0m_T[t]=alpha.0m_out[(nburn+10*t)]
    alpha.1m_T[t]=alpha.1m_out[(nburn+10*t)]

    eta.0f_T[t]=eta.0f_out[(nburn+10*t)]
    eta.1f_T[t]=eta.1f_out[(nburn+10*t)]
    eta.0m_T[t]=eta.0m_out[(nburn+10*t)]
    eta.1m_T[t]=eta.1m_out[(nburn+10*t)]

    b.0f_T[t,]=b.0f_out[(nburn+10*t),]
    b.1f_T[t,]=b.1f_out[(nburn+10*t),]
    b.0m_T[t,]=b.0m_out[(nburn+10*t),]
    b.1m_T[t,]=b.1m_out[(nburn+10*t),]

    tau2.f_T[t]=tau2.f_out[(nburn+10*t)]
    tau2.m_T[t]=tau2.m_out[(nburn+10*t)]

    lambda.f_T[t,]=lambda.f_out[(nburn+10*t),]
    lambda.m_T[t,]=lambda.m_out[(nburn+10*t),]

    sig.11_T[t]=sig.11_out[(nburn+10*t)]
    sig.12_T[t]=sig.12_out[(nburn+10*t)]
    sig.13_T[t]=sig.13_out[(nburn+10*t)]
    sig.14_T[t]=sig.14_out[(nburn+10*t)]
    sig.22_T[t]=sig.22_out[(nburn+10*t)]
    sig.23_T[t]=sig.23_out[(nburn+10*t)]
    sig.24_T[t]=sig.24_out[(nburn+10*t)]
    sig.33_T[t]=sig.33_out[(nburn+10*t)]
    sig.34_T[t]=sig.34_out[(nburn+10*t)]
    sig.44_T[t]=sig.44_out[(nburn+10*t)]

    L.f_T[t,]=L.f_out[(nburn+10*t),]
    L.m_T[t,]=L.m_out[(nburn+10*t),]
  }


  ## calculate posterior estimates
  beta.0.est=mean(beta.0_out[(nburn+1):nsim])
  beta.1f.est=mean(beta.1f_out[(nburn+1):nsim])
  beta.1m.est=mean(beta.1m_out[(nburn+1):nsim])
  beta.2.est=mean(beta.2_out[(nburn+1):nsim])

  beta.0.ci=quantile(beta.0_out[(nburn+1):nsim],probs=c(0.025,0.975))
  beta.1f.ci=quantile(beta.1f_out[(nburn+1):nsim],probs=c(0.025,0.975))
  beta.1m.ci=quantile(beta.1m_out[(nburn+1):nsim],probs=c(0.025,0.975))
  beta.2.ci=quantile(beta.2_out[(nburn+1):nsim],probs=c(0.025,0.975))

  alpha.0f.est=mean(alpha.0f_out[(nburn+1):nsim])
  alpha.1f.est=mean(alpha.1f_out[(nburn+1):nsim])
  alpha.0m.est=mean(alpha.0m_out[(nburn+1):nsim])
  alpha.1m.est=mean(alpha.1m_out[(nburn+1):nsim])
  alpha.est=c(alpha.0f.est,alpha.1f.est,alpha.0m.est,alpha.1m.est)

  alpha.0f.ci=quantile(alpha.0f_out[(nburn+1):nsim],probs=c(0.025,0.975))
  alpha.1f.ci=quantile(alpha.1f_out[(nburn+1):nsim],probs=c(0.025,0.975))
  alpha.0m.ci=quantile(alpha.0m_out[(nburn+1):nsim],probs=c(0.025,0.975))
  alpha.1m.ci=quantile(alpha.1m_out[(nburn+1):nsim],probs=c(0.025,0.975))

  eta.0f.est=mean(eta.0f_out[(nburn+1):nsim])
  eta.1f.est=mean(eta.1f_out[(nburn+1):nsim])
  eta.0m.est=mean(eta.0m_out[(nburn+1):nsim])
  eta.1m.est=mean(eta.1m_out[(nburn+1):nsim])

  eta.0f.ci=quantile(eta.0f_out[(nburn+1):nsim],probs=c(0.025,0.975))
  eta.1f.ci=quantile(eta.1f_out[(nburn+1):nsim],probs=c(0.025,0.975))
  eta.0m.ci=quantile(eta.0m_out[(nburn+1):nsim],probs=c(0.025,0.975))
  eta.1m.ci=quantile(eta.1m_out[(nburn+1):nsim],probs=c(0.025,0.975))

  b.0f.est=colMeans(b.0f_out[(nburn+1):nsim,])
  b.1f.est=colMeans(b.1f_out[(nburn+1):nsim,])
  b.0m.est=colMeans(b.0m_out[(nburn+1):nsim,])
  b.1m.est=colMeans(b.1m_out[(nburn+1):nsim,])

  Sigma.b.est=rowMeans(Sigma.b_out[,,(nburn+1):nsim],dims=2)

  pi.f.est=colMeans(pi.f_out[(nburn+1):nsim,])
  pi.m.est=colMeans(pi.m_out[(nburn+1):nsim,])

  tau2.f.est=mean(tau2.f_out[(nburn+1):nsim])
  tau2.m.est=mean(tau2.m_out[(nburn+1):nsim])

  tau2.f.ci=quantile(tau2.f_out[(nburn+1):nsim],probs=c(0.025,0.975))
  tau2.m.ci=quantile(tau2.m_out[(nburn+1):nsim],probs=c(0.025,0.975))

  lambda.f.est=colMeans(lambda.f_out[(nburn+1):nsim,])
  lambda.m.est=colMeans(lambda.m_out[(nburn+1):nsim,])
  lambda.f.ci=apply(lambda.f_out[(nburn+1):nsim,],2,quantile,probs=c(0.025,0.975))
  lambda.m.ci=apply(lambda.m_out[(nburn+1):nsim,],2,quantile,probs=c(0.025,0.975))

  xlambda.f.est=covariate_f_mat%*%lambda.f.est
  xlambda.m.est=covariate_m_mat%*%lambda.m.est

  sig.11.est=mean(sig.11_out[(nburn+1):nsim])
  sig.12.est=mean(sig.12_out[(nburn+1):nsim])
  sig.13.est=mean(sig.13_out[(nburn+1):nsim])
  sig.14.est=mean(sig.14_out[(nburn+1):nsim])
  sig.22.est=mean(sig.22_out[(nburn+1):nsim])
  sig.23.est=mean(sig.23_out[(nburn+1):nsim])
  sig.24.est=mean(sig.24_out[(nburn+1):nsim])
  sig.33.est=mean(sig.33_out[(nburn+1):nsim])
  sig.34.est=mean(sig.34_out[(nburn+1):nsim])
  sig.44.est=mean(sig.44_out[(nburn+1):nsim])

  sig.11.ci=quantile(sig.11_out[(nburn+1):nsim],probs=c(0.025,0.975))
  sig.12.ci=quantile(sig.12_out[(nburn+1):nsim],probs=c(0.025,0.975))
  sig.13.ci=quantile(sig.13_out[(nburn+1):nsim],probs=c(0.025,0.975))
  sig.14.ci=quantile(sig.14_out[(nburn+1):nsim],probs=c(0.025,0.975))
  sig.22.ci=quantile(sig.22_out[(nburn+1):nsim],probs=c(0.025,0.975))
  sig.23.ci=quantile(sig.23_out[(nburn+1):nsim],probs=c(0.025,0.975))
  sig.24.ci=quantile(sig.24_out[(nburn+1):nsim],probs=c(0.025,0.975))
  sig.33.ci=quantile(sig.33_out[(nburn+1):nsim],probs=c(0.025,0.975))
  sig.34.ci=quantile(sig.34_out[(nburn+1):nsim],probs=c(0.025,0.975))
  sig.44.ci=quantile(sig.44_out[(nburn+1):nsim],probs=c(0.025,0.975))
  for (gt in 1:nsim){

    #########################
    #### update latent class
    #########################
    ## update latent class
    mu.ij.fL.est=matrix(1,I,J)
    mu.ij.mL.est=matrix(1,I,J)
    pr.latent.f=matrix(0,I,K)
    pr.latent.m=matrix(0,I,K)

    for (k in 1:K){
      ### for female
      Lamb.fL.est=beta.0.est+beta.1f.est*(k-1)+beta.1m.est*L.m+beta.2.est*(k-1)*L.m
      for (j in 1:J){
        mu.ij.fL.est[,j]=b.0f[j]+b.1f[j]*(k-1)+xlambda.f.est
      }
      mu.ij.fL1.est=mu.ij.fL.est-mean(mu.ij.fL.est)

      part1.f=(exp(Lamb.fL.est)/(1+exp(Lamb.fL.est)))^Yvariable*(1/(1+exp(Lamb.fL.est)))^(1-Yvariable)

      temp.f=(exp(eta.0f.est+eta.1f.est*mu.ij.fL1.est)/(1+exp(eta.0f.est+eta.1f.est*mu.ij.fL1.est)))^U.f*(1/(1+exp(eta.0f.est+eta.1f.est*mu.ij.fL1.est)))^(1-U.f)*
        (dlnorm(V.f,meanlog=mu.ij.fL.est,sdlog=sqrt(tau2.f.est)))^U.f
      part2.f=apply(temp.f,1,prod,na.rm=TRUE)

      pr.latent.f[,k]=pi.f.est[k]*part1.f*part2.f

      ### for male
      Lamb.mL.est=beta.0.est+beta.1f.est*L.f+beta.1m.est*(k-1)+beta.2.est*L.f*(k-1)
      for (j in 1:J){
        mu.ij.mL.est[,j]=b.0m[j]+b.1m[j]*(k-1)+xlambda.m.est
      }
      mu.ij.mL1.est=mu.ij.mL.est-mean(mu.ij.mL.est)

      part1.m=(exp(Lamb.mL.est)/(1+exp(Lamb.mL.est)))^Yvariable*(1/(1+exp(Lamb.mL.est)))^(1-Yvariable)

      temp.m=(exp(eta.0m.est+eta.1m.est*mu.ij.mL1.est)/(1+exp(eta.0m.est+eta.1m.est*mu.ij.mL1.est)))^U.m*(1/(1+exp(eta.0m.est+eta.1m.est*mu.ij.mL1.est)))^(1-U.m)*
        (dlnorm(V.m,meanlog=mu.ij.mL.est,sdlog=sqrt(tau2.m.est)))^U.m
      part2.m=apply(temp.m,1,prod,na.rm=TRUE)

      pr.latent.m[,k]=pi.m.est[k]*part1.m*part2.m
    }

    pr.latent.f=pr.latent.f/rowSums(pr.latent.f)
    pr.latent.m=pr.latent.m/rowSums(pr.latent.m)
    for (i in 1:I){
      L.f[i]=sample(0:(K-1),1,replace=TRUE,prob=pr.latent.f[i,])
      L.m[i]=sample(0:(K-1),1,replace=TRUE,prob=pr.latent.m[i,])
    }


    ############################
    #### update random effect b
    ############################
    for (j in 1:J){
      if (gt>100){
        cov.r=cov(cbind(b.0f_out[1:(gt-1),j],b.1f_out[1:(gt-1),j],b.0m_out[1:(gt-1),j],b.1m_out[1:(gt-1),j]))
      }
      if (gt>1 & ((gt-1)%%100==0)){
        ar.r=sum(diff(b.0f_out[(gt-100):(gt-1),j])!=0)/99
        if (ar.r<0.44){s.r=s.r/sqrt(2)}
        else if (ar.r>0.44){s.r=s.r*sqrt(2)}
      }

      u.r=runif(1)
      b.star=rmvnorm(1,c(b.0f[j],b.1f[j],b.0m[j],b.1m[j]),s.r*cov.r+s.r*eps*diag(4))
      b.0f.star=b.star[1]
      b.1f.star=b.star[2]
      b.0m.star=b.star[3]
      b.1m.star=b.star[4]

      R.r=MH_b.est(U.f[,j],U.m[,j],V.f[,j],V.m[,j],b.0f.star,b.1f.star,b.0m.star,b.1m.star,b.0f[j],b.1f[j],b.0m[j],b.1m[j])

      if (u.r<=R.r){b.0f[j]=b.0f.star;b.1f[j]=b.1f.star;b.0m[j]=b.0m.star;b.1m[j]=b.1m.star}
    }
    random.b=cbind(b.0f,b.1f,b.0m,b.1m)


    ###################################
    #### update mu.ij.f and mu.ij.m
    ###################################
    for (j in 1:J){
      mu.ij.f[,j]=b.0f[j]+b.1f[j]*L.f+xlambda.f.est
      mu.ij.m[,j]=b.0m[j]+b.1m[j]*L.m+xlambda.m.est
    }
    mu.ij.f1=mu.ij.f-mean(mu.ij.f)
    mu.ij.m1=mu.ij.m-mean(mu.ij.m)


    ########################################
    #### calculate 2nd term of DIC
    ########################################
    Lamb.est=beta.0.est+beta.1f.est*L.f+beta.1m.est*L.m+beta.2.est*L.f*L.m

    pt1=sum(Yvariable*log((exp(Lamb.est)/(1+exp(Lamb.est))))-(1-Yvariable)*log(1+exp(Lamb.est)))
    pt2=sum(log(pi.f.est[(L.f+1)])+log(pi.m.est[(L.m+1)]))
    pt3=sum(U.f*log(exp(eta.0f.est+eta.1f.est*mu.ij.f1)/(1+exp(eta.0f.est+eta.1f.est*mu.ij.f1)))-(1-U.f)*log(1+exp(eta.0f.est+eta.1f.est*mu.ij.f1))+
              U.m*log(exp(eta.0m.est+eta.1m.est*mu.ij.m1)/(1+exp(eta.0m.est+eta.1m.est*mu.ij.m1)))-(1-U.m)*log(1+exp(eta.0m.est+eta.1m.est*mu.ij.m1)),na.rm=TRUE)
    pt4=sum(log(dlnorm(V.f,meanlog=mu.ij.f,sdlog=sqrt(tau2.f.est))),na.rm=TRUE)+
      sum(log(dlnorm(V.m,meanlog=mu.ij.m,sdlog=sqrt(tau2.m.est))),na.rm=TRUE)

    ptt=rep(0,J)
    for (j in 1:J){
      ptt[j]=log(dmvnorm(c(b.0f[j],b.1f[j],b.0m[j],b.1m[j]),alpha.est,Sigma.b.est))
    }
    pt5=sum(ptt)

    loglike.est[gt]=pt1+pt2+pt3+pt4+pt5

  }  ###### 2nd MCMC Done!


  ## calculate DIC
  dev1=-2*mean(loglike[(nburn+1):nsim])
  dev2=-2*mean(loglike.est[(nburn+1):nsim])

  pD=dev1-dev2
  DIC=dev1+pD
  ## return posterior estimates
 mcmc2<- list(
   beta.0.est,    beta.1f.est,    beta.1m.est,    beta.2.est,
   beta.0.ci,    beta.1f.ci,   beta.1m.ci,   beta.2.ci,
   alpha.est,
   alpha.0f.ci,   alpha.1f.ci,    alpha.0m.ci,    alpha.1m.ci,
   eta.0f.est,   eta.1f.est, eta.0m.est,   eta.1m.est,
   eta.0f.ci,   eta.1f.ci,   eta.0m.ci,   eta.1m.ci,
   b.0f.est,   b.1f.est,   b.0m.est,   b.1m.est,
   Sigma.b.est,
   pi.f.est,   pi.m.est,
   tau2.f.est,   tau2.m.est,   tau2.f.ci,   tau2.m.ci,
   lambda.f.est,   lambda.m.est,   lambda.f.ci,   lambda.m.ci,
   xlambda.f.est,   xlambda.m.est,
   sig.11.est,   sig.12.est,   sig.13.est,   sig.14.est,   sig.22.est,
   sig.23.est,   sig.24.est,   sig.33.est,   sig.34.est,   sig.44.est,
   sig.11.ci,   sig.12.ci,   sig.13.ci,   sig.14.ci,   sig.22.ci,
   sig.23.ci,   sig.24.ci,   sig.33.ci,   sig.34.ci,   sig.44.ci,

   random.b,  mu.ij.fL.est, mu.ij.mL.est, pr.latent.f, pr.latent.m ,     mu.ij.f1,
                mu.ij.m1,    DIC, loglike.est)

}
