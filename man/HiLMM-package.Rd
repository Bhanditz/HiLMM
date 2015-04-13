\name{HiLMM-package}
\alias{HiLMM-package}
\docType{package}
\title{
Estimation of heritability in linear mixed models
}

\description{
The package contains two functions: estim_herit computes heritability when the random effects follow either a Gaussian distribution or a mixture of a Dirac and a Gaussian distribution. It also provides a confidence interval of heritability when the random effects are Gaussian. The function data_simu allows the user to generate dataset which are compatible wit the arguments of HiLMM.
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
 Anna Bonnet <anna.bonnet@agroparistech.fr>

}
\references{
The method is developped in the paper "Heritability estimation in sparse high dimensional linear mixed models" (A.Bonnet,E.Gassiat,C.Levy-Leduc,2014)
}
\keyword{heritability }
\examples{
library(HiLMM)
data_HiLMM=data_simu(100,1000,0.7,1)
Y=data_HiLMM$Y
W=data_HiLMM$W
estim_herit(Y,W)$heritability 
estim_herit(Y,W)$CI_low
estim_herit(Y,W)$CI_up
estim_herit(Y,W)$standard_deviation
}
