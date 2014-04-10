\name{HiLMM-package}
\alias{HiLMM-package}
\alias{HiLMM}
\docType{package}
\title{
Estimation of heritability in high dimensional Linear Mixed Models
}

\description{
The package computes the heritability when the random effects follow either a Gaussian distribution or a mixture of a Dirac and a Gaussian distribution. It also provides a confidence interval of the heritability when the random effects are Gaussian. For further details, see the paper quoted in the references section.
}
\details{
\tabular{ll}{
Package: \tab HiLMM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-10\cr
License: \tab GPL-2\cr
}
}
\author{
Anna Bonnet

Maintainer: Anna Bonnet <anna.bonnet@agroparistech.fr>

}
\references{
This function implements the method described in the paper "Heritability estimation in high dimensional linear mixed models" by A. Bonnet, E. Gassiat, C. Levy-Leduc, 2014 (submitted and availavable on arXiv:1404.3397).
}
\keyword{ heritability }
\keyword{ linear mixed model }
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
