\name{data_simu}
\alias{data_simu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to generate a dataset which is compatible with the arguments of the function HiLMM
}
\description{
This function generates a dataset which contains a matrix corresponding 
to the genetic information and a vector of observations. This dataset corresponds 
to the type of arguments which is expected in the funtion HiLMM.
}
\usage{
data_simu(n, N, eta_star, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
 size of the vector of observations that the user wants to generate.
}
  \item{N}{
   number of columns of the genetic information matrix 
  that the user wants to generate.
}
  \item{eta_star}{
 value of the heritability.
}
  \item{q}{
 proportion of non zero components in the random effects associated 
to the genetic information matrix.
}
}

\value{
  \item{Y}{vector of observations of size n}
  \item{W}{Matrix of size n x N which contains entries with 0,1 and 2}
}


\examples{
library(HiLMM)
data_sim=data_simu(100,1000,0.5,0.5)
 Y=data_sim$Y
 W=data_sim$W
## The function is currently defined as
function (n, N, eta_star, q) 
{
    sigma_u = 1
    P = runif(N, 0.1, 0.5)
    W = matrix(0, n, N)
    for (j in 1:N) {
        W[, j] = rbinom(n, 2, P[j])
    }
    nb_comp_non_zero = q * N
    sigma_e = sqrt(q * N * sigma_u^2 * (1 - eta_star)/eta_star)
    b = sample(1:N, nb_comp_non_zero)
    a1 = sort(b)
    u = rnorm(nb_comp_non_zero, 0, sigma_u)
    e = rnorm(n, 0, sigma_e)
    U = matrix(0, N)
    U[a1] = u
    Z = scale(W, center = TRUE, scale = TRUE)
    Y = Z \%*\% U + e
    list(Z = Z, Y = Y)
  }
}
