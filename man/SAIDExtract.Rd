\name{SAIDExtract}
\alias{SAIDExtract}
\title{Extract Main Effects and Interactions from SAID Object.}

\description{
  This function extracts estimates of main effects, interaction effects, and the overall dose response surface of the exposures provided a fitted SAID object.
}

\usage{
SAIDExtract(whole_model_SAID, Xtest, good_index, me_integral_constraint,
                      zero_ind = rep(1, dim(Xtest)[2]))
}

\arguments{
  \item{whole_model_SAID}{Fitted SAID object.}
  \item{Xtest}{ntest by p matrix of exposure test values.}
  \item{good_index}{Indices of MCMC samples to be used.}
  \item{me_integral_constraint}{TRUE if integral constraint used, FALSE if origin constraint used.}
  \item{zero_ind}{Vector of 0s and 1s with 1 indicating an interaction is selected in the model and 0 indicating it is not. By default, all interactions are assumed to have been selected.}
}
  

\details{
  A list is returned with the following components:
}

\value{
\item{whole_sur_samples}{Returns MCMC samples of the dose response surface evaluated at test exposure points.}
\item{int_sur_samples}{Returns MCMC samples of each individual interaction evaluated at test exposure points.}
\item{ME_spline_list}{Returns the spline matrix used to fit main effects.}
\item{IE_spline_list}{Returns the spline matrix used to fit interaction effects.}
\item{int_sur_train_samples}{Returns MCMC samples of the dose response surface evaluated at training exposure points.}
}

\references{
  [1. Chattopadhyay, Engel, and Dunson (2022+)]
}

\author{Shounak Chattopadhyay}

\examples{
[to be added]
}
