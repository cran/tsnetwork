\name{tsnetwork}
\alias{tsnetwork}
\title{
Time Series Chain Graphical Models for mixed discrete-continuous time-series data. 
}
\description{
Reconstructing instantaneous (undirected) and dynamic (directed) networks from 
 repeated multivariate mixed discrete-continuous or ordinal time series data. 
 This function computes sparse an autoregressive coefficient and a precision matrices for time
 series chain graphical models.
 }
 
\usage{ tsnetwork(dat, lower, upper, penalty= "scad", n_lam1= NULL, 
      lam1_ratio= NULL, n_lam2= NULL, lam2_ratio= NULL, em_tol= NULL,
      em_iter= NULL, iter_Mstep= NULL, pen_diag_gamma= FALSE,
      ncores = 1)
}

\arguments{
	\item{dat}{Longitudinal data format}
	\item{lower}{Lower boundry of the data. Can be cimouted internally. Deafult is NULL. }
	\item{upper}{Upper boundry of the data. Can be cimouted internally. Deafult is NULL.}
	\item{penalty}{This specifies the type of penalty function to be used. SCAD penalty function is applied if penalty = "scad" and GLASSO is applied if penalty = "lasso"}
	\item{n_lam1}{The number of regularization parameters for the instantaneous interactions.}
	\item{lam1_ratio}{Determines the sequence of lam1.}
	\item{n_lam2}{The number of regularization parameters for the dynamics interactions.}
	\item{lam2_ratio}{Determines the sequence of lam2.}
	\item{em_tol}{A value to meet the convergence criteria of the EM algorithm. Default value is 0.01}
	\item{em_iter}{The number of EM iterations. The default value is 10.}
	\item{iter_Mstep}{The number of iterations in the M-step to guarantee the convergence. The default value is 5.}
	\item{pen_diag_gamma}{Penalazing the diagonal elements of the autoregressive matrix.}
	\item{ncores}{The number of cores to use for the calculations.}
	}
\author{
Pariya Behrouzi \cr
Maintainer: Pariya Behrouzi <pariya.behrouzi@gmail.com>
}
\references{
\cite{Pariya Behrouzi, Fentaw Abegaz and Ernst Wit (2018). Dynamic Chain Graph Models for Ordinal Time Series Data. 
\emph{Arxiv}. \bold{14}, 3: 586-599.}
  
\cite{Fentaw Abegaz and Ernst Wit  (2013). Sparse time series chain graphical models for reconstructing
genetic networks. \emph{Biostatistics}. \bold{14}, 3: 586-599.}
}
\keyword{package}

\examples{
simulate <- gen.sim(t = 3, n = 2, p = 3, k = 3, network = "scale-free") 
sim.dat <- simulate$dat
out <- tsnetwork(dat =sim.dat, lower= NULL, upper= NULL, penalty= "lasso",
                 n_lam1= 1, lam1_ratio= NULL, n_lam2= 1, lam2_ratio= NULL, em_tol = NULL,
                 em_iter= 1, iter_Mstep = 1, pen_diag_gamma= FALSE, ncores = 1)

# Estimated sparse precision (undirected) and autoregression (directed) matrices
undirected <- out$theta[ , , 1, 1]
directed <- out$gamma[ , , 1, 1]

oldpar <- par(no.readonly =TRUE)
par(mfrow=c(1,2))
plotG(undirected, mod="undirected", main= "Estimated precision matrix", label=TRUE)     
plotG(directed, mod="directed", main ="Estimated autoregression coef. matrix", label=TRUE) 
par(oldpar)
}
