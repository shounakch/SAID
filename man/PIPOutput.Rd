\name{PIPOutput}
\alias{PIPOutput}
\title{SAID Posterior Inclusion Probability}

\description{
  This function implements the S2M (Sequential 2-Means) framework in Li and Pati (2017) to carry out variable selection based on posterior samples obtained from SAID. 
}

\usage{
PIPOutput(int_max, int_min, length.grid = 1000, lower.thres = 0.00001)
}

\arguments{
  \item{int_max}{MC by K matrix of integrals for positive part of interactions. MC = number of MC samples and K = number of interactions.}
  \item{int_min}{MC by K matrix of integrals for negative part of interactions.}
  \item{length.grid}{Length of grid of cutoff values for S2M.}
  \item{lower.thres}{Controls lower range of cutoff values.}
}
  

\details{
  A list is returned with the following components:
}

\value{
\item{PIPInt}{Returns the posterior inclusion probabilities (PIPs) for each interaction.}
\item{PSPInt}{Returns the posterior synergistic probabilities (PSPs) for each interaction.}
\item{PAPInt}{Returns the posterior antagonistic probabilities (PAPs) for each interaction.}
}

\references{
  [1. Li and Pati (2017)
  2. Chattopadhyay, Engel, and Dunson (2022+)]
}

\author{Shounak Chattopadhyay}

\examples{
[to be added]
}
