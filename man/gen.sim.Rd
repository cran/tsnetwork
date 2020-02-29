\name{gen.sim}
\alias{gen.sim}

\title{Generate ordinal multivariate time series data}

\description{
Generating data from time series chain graphical model with vector autoregressive model of order 1.
}

\usage{
gen.sim(t = NULL, n = NULL, p = NULL, k = NULL, network = NULL, prob0 = NULL)
}

\arguments{
	\item{t}{The number of time points.}
	\item{n}{The number of sample size (observations). The default value is 200.}
	\item{p}{The number of variables. The default value is 90.}
	\item{k}{The number of categories in the data. The default value is 5.}
	\item{network}{Specifies the type of network structure. This could be random, scale-free, hub
or user defined structures. Details on simultions from the various network structures can be found in the R package flare.}
	\item{prob0}{Initial sparsity level.}
}

\value{
An object with S3 class "gen.sim" is returned: \cr \cr
 \item{data}{A simulated time-series data.}
 \item{true_theta}{Sparse precision matrix.}
 \item{true_gamma}{Sparse autoregressive coefficients matrix.}
 \item{true_sigma}{Covariance matrix.}
} 

\author{ 
  Pariya Behrouzi and Fentaw Abegaz \cr
  Maintainer: Pariya Behrouzi <pariya.behrouzi@gmail.com>
}

\references{
1. Behrouzi, P., and Wit, E. C. (2017c). netgwas: An R Package for Network-Based Genome-Wide Association Studies. arXiv preprint, arXiv:1710.01236.\cr
}

\seealso{
\code{\link{tsnetwork}}
}

\examples{
dat <- gen.sim(t = 2, n = 10, p = 5, k = 2, network = "hub")
sim <- dat$data1
prec.true <- dat$theta
autoR.true<- dat$gamma
}