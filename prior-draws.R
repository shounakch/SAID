library(splines2)
library(ggplot2)
library(MASS)

setwd("/Volumes/GoogleDrive/My Drive/Research_Duke/Shape Constrained Inference/Code/Updated Codes/final-package-management/SynInt-main-test-exact")

Rcpp::sourceCpp("src/cpp_MCMCSampler_adaptive_optimized.cpp")
source("/Volumes/GoogleDrive/My Drive/Research_Duke/Shape Constrained Inference/Code/Updated Codes/final-package-management/SynInt-main-test-exact/R/cpp_optimized_SIDsampler.R", echo=TRUE)

penmatt<-function(M)
{
  
  Amat = matrix(0, nrow = M, ncol = M)
  
  Amat[1,] = c(1, rep(0, M-1))
  
  for(i in 2:M)
  {
    
    Amat[i,] = c(rep(0, i-2), -1, 1, rep(0, M-i))
    
  }
  
  return(t(Amat) %*% Amat)
  
}

qform<-function(x, A)
{
  
  b = t(x) %*% A %*% x
  
  return(as.numeric(b))
  
}

##### Generate X's

r = 1
set.seed(2001 + r)

n = 1000
p = 2

X = matrix(runif(n*p), nrow = n, ncol = p)

##### Construct B-spline matrices

#K_IE = 3

nspl = 6

M1 = bSpline(x = X[,1], df = nspl)
M2 = bSpline(x = X[,2], df = nspl)

A1 = (t(M1) %*% M1) / n
A2 = (t(M2) %*% M2) / n

##### Other hyperparameters

pen_param = 1
SigmaInv = penmatt(nspl)

MC = 20000

criterion_stor = rep(0, MC)

#### Initialize parameters and start MCMC

theta1 = rnorm(nspl, mean = 0, sd = 0.1)
theta2 = rnorm(nspl, mean = 0, sd = 0.1)
phi1 = rnorm(nspl, mean = 0, sd = 0.1)
phi2 = rnorm(nspl, mean = 0, sd = 0.1)

for(m in 1:MC)
{
  
  ## Generate theta1
  
  C11 = qform(phi1, A2) * qform(theta2, A1) * qform(phi2, A2)
  
  theta1_var = chol2inv(chol(SigmaInv + 
                       (2 * pen_param * C11 * A1)))
  
  theta1 = mvrnorm(n = 1, mu = rep(0, nspl), Sigma = theta1_var)
  
  ## Generate theta2
  
  C21 = qform(theta1, A1) * qform(phi1, A2) * qform(phi2, A2)
  
  theta2_var = chol2inv(chol(SigmaInv + 
                               (2 * pen_param * C21 * A1)))
  
  theta2 = mvrnorm(n = 1, mu = rep(0, nspl), Sigma = theta2_var)
  
  ## Generate phi1
  
  C12 = qform(theta1, A1) * qform(theta2, A1) * qform(phi2, A2)
  
  phi1_var = chol2inv(chol(SigmaInv + 
                               (2 * pen_param * C12 * A2)))
  
  phi1 = mvrnorm(n = 1, mu = rep(0, nspl), Sigma = phi1_var)
  
  ## Generate phi2
  
  C22 = qform(theta1, A1) * qform(phi1, A2) * qform(theta1, A1)
  
  phi2_var = chol2inv(chol(SigmaInv + 
                             (2 * pen_param * C22 * A2)))
  
  phi2 = mvrnorm(n = 1, mu = rep(0, nspl), Sigma = phi2_var)
  
  ## Compute h(\beta) at training points
  
  P1 = (M1 %*% theta1)^2
  P2 = (M2 %*% phi1)^2
  N1 = (M1 %*% theta2)^2
  N2 = (M2 %*% phi2)^2
  
  h = (P1 * P2) - (N1 * N2)
  
  hplus = sapply(h, max, 0)
  hminus = sapply(-h, max, 0)
  
  criterion_stor[m] = mean(hplus) * mean(hminus)
  
  ## Print iterate number
  
  if(m %% 1000 == 0)
  {
    
    print(paste("Iterate: ", m, sep = ""))
    
  }
  
}

## Look at histograms

df = as.data.frame(criterion_stor)

# ggpl = ggplot(df, aes(x=criterion_stor)) + geom_histogram(binwidth = 0.1)
# print(ggpl)

hist(criterion_stor, probability = TRUE, nclass = 50)

max(criterion_stor[10000:20000])

mean(criterion_stor[10000:20000] <= 0.001)

## Plot criterion vs kappa

df_plot = data.frame(matrix(nrow = 6, ncol = 3))
colnames(df_plot) = c("kappa", "max_crit", "prob_crit")
df_plot$kappa = c(0, 0.1, 1, 10, 100, 1000)
df_plot$max_crit = c(939.03, 26.05, 7.06, 1.59, 0.35, 0.05)
df_plot$prob_crit = c(0.128, 0.144, 0.179, 0.290, 0.537, 0.865)

pl1 = ggplot(df_plot, aes(x = kappa, y = prob_crit)) + geom_line() + geom_point()
pl1 = pl1 + labs(x = "Penalty Parameter", y = "Proportion") + ylim(c(0,1))
pl1 = pl1 + theme(text = element_text(size = 20))

print(pl1)

write.csv(df_plot, file = paste(getwd(), "/prior_plot_data.csv", sep = ""))
