\name{estim_herit}
\alias{estim_herit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the heritability in high dimensional linear mixed models
}
\description{
The function estimates the heritability in a linear mixed model in which
the random effect part is a mixture of a Gaussian and a Dirac distribution.
It also provides a confidence interval for the heritability. For further
details on the considered model and on the definition of the heritability, we refer
the reader to the paper quoted in the references section.}
\usage{
estim_herit(Y, W)
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
\item{IC_inf}{Lower bound of the confidence interval}
\item{IC_sup}{Upper bound of the confidence interval}
\item{standard_deviation }{Standard deviation}
}
\references{
This function implements the method described in the paper "Heritability estimation in high dimensional linear mixed models" by A. Bonnet, E. Gassiat, C. Levy-Leduc, 2014 (submitted and availavable on arXiv:1404.3397).
}
\author{
Anna Bonnet
}

\examples{
library(HiLMM)
data(vec_Y)
data(mat_W)
Y=data.matrix(vec_Y,rownames.force=NA)
W=data.matrix(mat_W,rownames.force=NA)
resultat=estim_herit(Y,W)
resultat$heritability 
resultat$IC_inf
resultat$IC_sup
resultat$standard_deviation
}
\keyword{ heritability }
\keyword{ linear mixed model }

