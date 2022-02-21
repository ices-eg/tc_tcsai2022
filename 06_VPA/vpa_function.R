#' Virtual Population Analysis
#'
#' Run simple VPA model, with no tuning index.
#'
#' @param C catch-at-age matrix.
#' @param M natural mortality rate, a scalar.
#' @param Fterm terminal F, either a scalar or a vector of same length as the
#'        number of ages.
#' @param Fages number of ages to calculate F for oldest age. For example, if
#'        the catch-at-age matrix contains ages up to 15 and \code{Fages=5},
#'        then F for the oldest age will be set as the average of ages 10-14.
#'
#' @return
#' List containing the model results (\code{N}, \code{F}, \code{Z}), as well as
#' the input objects passed by the user (\code{C}, \code{M}, \code{Fterm},
#' \code{Fages}).
#'
#' @export

vpa <- function(C, M, Fterm, Fages)
{
  ## Prepare matrices
  C <- as.matrix(C)
  T <- nrow(C)
  A <- ncol(C)
  N <- F <- Z <- matrix(NA_real_, nrow=T, ncol=A, dimnames=dimnames(C))

  ## Set F in terminal year
  F[T,] <- Fterm
  Z <- F + M

  ## Calculate N in terminal year
  N <- C*Z / (F*(1-exp(-Z)))

  ## Calculate N and F up to terminal year,
  ## assuming F[oldest] = avg(preceding ages)
  for(t in (T-1):1)
  {
    for(a in 1:(A-1))
    {
      N[t,a] <- N[t+1,a+1] * exp(M) + C[t,a] * exp(M/2)
      F[t,a] <- log(N[t,a] / N[t+1,a+1]) - M
    }
    F[t,A] <- mean(F[t,A-(1:Fages)])
    Z[t,] <- F[t,] + M
    N[t,A] <- C[t,A]*Z[t,A] / (F[t,A]*(1-exp(-Z[t,A])))
  }

  list(N=N, F=F, Z=Z, C=C, M=M, Fterm=Fterm, Fages=Fages)
}
