\name{estim_herit}
\alias{estim_herit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of heritability in linear mixed models
}
\description{
The function computes heritability when the random effects follow either a Gaussian distribution or a mixture of a Dirac and a Gaussian distribution.
 It also provides a confidence interval of heritability when the random effects are Gaussian.
 }
\usage{
estim_herit(Y,W)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
vector of observations of size n
}
  \item{W}{
matrix with n rows and N columns
}
}

\value{

\item{heritability}{Heritability}
\item{CI_low}{Lower bound of the confidence interval}
\item{CI_up}{Upper bound of the confidence interval}
\item{standard_deviation }{Standard deviation}
}
\references{
The method is developped in the paper "Heritability estimation in high dimensional linear mixed models" (A.Bonnet,E.Gassiat,C.L??vy-Leduc,2014)
}
\author{
Anna Bonnet
}

\examples{
library(HiLMM)
data(Y)
data(W)
estim_herit(Y,W)$heritability 
estim_herit(Y,W)$CI_low
estim_herit(Y,W)$CI_up
estim_herit(Y,W)$standard_deviation
}
\keyword{heritability}
