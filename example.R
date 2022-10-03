#library(SynIntTestExact)

setwd("/Volumes/GoogleDrive/My Drive/Research_Duke/Shape Constrained Inference/Code/Updated Codes/final-package-management/SynInt-main-test-exact")

Rcpp::sourceCpp("src/cpp_MCMCSampler_adaptive_optimized.cpp")
source("/Volumes/GoogleDrive/My Drive/Research_Duke/Shape Constrained Inference/Code/Updated Codes/final-package-management/SynInt-main-test-exact/R/cpp_optimized_SIDsampler.R", echo=TRUE)

r = 1
set.seed(2001 + r)

n = 500
ntot = 1000
p = 2
p_cov = 5

############ Generate the data #################

gamma0 = 1
sigma0sq = 100

Xwhole = matrix(runif(ntot*p), nrow = ntot, ncol = p)
Zwhole = matrix(rnorm(ntot*p_cov), nrow = ntot, ncol = p)

#true_int_whole = gamma0 * (Xwhole[,1]^2) * (Xwhole[,2]^2)
true_int_whole = gamma0 * ((Xwhole[,1] * Xwhole[,2]) - (2*(Xwhole[,1]*Xwhole[,2])^2))

true_sur_whole = rep(0.5, ntot) +
  (Xwhole[,1]^2) +
  (Xwhole[,2]^2) +
  true_int_whole

y_whole = true_sur_whole + rnorm(n = ntot, mean = 0, sd = sqrt(sigma0sq))

Xtrain = Xwhole[1:n,]
ytrain = y_whole[1:n]
Ztrain = Zwhole[1:n,]

Xtest = Xwhole[-c(1:n),]
ytest = y_whole[-c(1:n)]
Ztest = Zwhole[-c(1:n),]

true_sur_test = true_sur_whole[-c(1:n)]  
true_int_test = true_int_whole[-c(1:n)]

MC = 10000

SIMmodel = SIMsampler(y = ytrain,
                      X = Xtrain,
                      Z = NULL,
                      MC = MC,
                      cutoff = MC,
                      me_integral_constraint = TRUE,
                      accept_low = 0.5,
                      accept_high = 0.7,
                      accept_scale = 0.8,
                      K_IE = 3,
                      eps_MALA = 0.001)
