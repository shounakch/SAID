\name{SAIDsampler}
\alias{SAIDsampler}
\title{Synergistic Antagonistic Interaction Detection posterior sampler}

\description{
  This function implements the Synergistic Antagonistic Interaction Detection (SAID) framework in Chattopadhyay et al. (2022+). The dose response surface is decomposed into main effects and pairwise interactions. The pairwise interactions are shrunk towards the space of synergistic or antagonistic interactions. 
}

\usage{
SAIDsampler(y, 
           X, 
           Z = NULL, 
           K_ME = 5, 
           K_IE = 2,
           a_lamb = 0.001,
           b_lamb = 0.001,
           eps_MALA = rep(0.01, choose(dim(X)[2], 2)), 
           c_HMC = 1.0, 
           L_HMC = 5, 
           MC = 10000, 
           zero_ind = rep(1, choose(dim(X)[2],2)), 
           me_integral_constraint = TRUE, 
           cutoff = 0.5*MC,
           accept_low = 0.65,
           accept_high = 0.9,
           accept_scale = 0.8)
}

\arguments{
  \item{y}{n by 1 vector of continuous responses.}
  \item{X}{n by p matrix of values taken by exposure variables.}
  \item{Z}{n by p_cov matrix of values taken by covariate variables. By default Z = NULL meaning there are no covariates.}
  \item{K_ME}{Number of knots taken to fit main effects. Default value is 5.}
  \item{K_IE}{Number of knots taken to fit interacton effects. Default value is 2.}
  \item{a_lamb, b_lamb}{Hyperparameters of the gamma prior for main effect scales.}
  \item{eps_MALA}{p*(p-1)/2 by 1 vector of step-sizes while performing HMC. Default sets all entries to 0.01.}
  \item{c_HMC}{Mass value while performing HMC. Default value is 1.}
  \item{L_HMC}{Number of HMC steps per MCMC iteration. Default is 5.}
  \item{MC}{Number of Markov chain Monte Carlo samples to be drawn from the posterior distribution. Default is 10000.}
  \item{zero_ind}{p*(p-1)/2 by 1 vector with entries = 0 or 1 according to whether the corresponding pairwise interaction is eliminated from the model or included, respectively. Default value rep(1, choose(p,2)) includes all pairwise interactions.}
  \item{me_integral_constraint}{If TRUE, then fits main effects using integral constraint; if FALSE, fits main effects using origin-start constraint. TRUE by default.}
  \item{cutoff}{Adaptive choosing of HMC step-size done till first cutoff * MC samples.}
  \item{accept_low}{Lower acceptance limit.}
  \item{accept_high}{Upper acceptance limit.}
  \item{accept_scale}{Perturbation factor in adaptive HMC phase.}
}
  

\details{
  A list is returned with the following components:
  1. SIM_model = List with values described below.
  2. ME_list
  3. ME_knots
  4. IE_knots
  5. data
  6. init_values
  7. zero_ind
}

\value{
  \item{intercept}{MC by 1 vector of sampled intercept values.}
  \item{error_var}{MC by 1 vector of sampled error variances.}
  \item{ME_coeff}{MCMC samples of main effect coefficients.}
  \item{ME_scale}{MCMC samples of main effect scales.}
  \item{IE_pos_coeff_part1}{MCMC samples of first part of positive function.}
  \item{IE_pos_coeff_part2}{MCMC samples of second part of positive function.}
  \item{IE_neg_coeff_part1}{MCMC samples of first part of negative function.}
  \item{IE_neg_coeff_part2}{MCMC samples of second part of negative function.}
  \item{ME_mat}{Returns the ME matrix of evaluated spline basis functions.}
  \item{IE_list}{Returns the IE list of spline basis matrices.}
  \item{cov_effects}{MCMC samples of covariate effects.}
  \item{Accept_Prop}{HMC acceptance indicators.}
  \item{HMC_epsilon}{HMC vector of step-sizes after sampling.}
  \item{ME_list}{Contains ME matrices.}
  \item{ME_knots}{Contains ME knots.}
  \item{IE_knots}{Contains IE knots.}
  \item{data}{List returning original data (y, X, Z) along with n, p, and MC.}
}

\references{
  [to be added]
}

\author{Shounak Chattopadhyay}

\examples{
[to be added]
}
