## Read data, prepare matrices
catch <-
  as.matrix(
    read.csv("06_VPA/cod_catch.csv",
      header = TRUE,
      check.names = FALSE,
      row.names = 1
    )
  ) / 1000

nyears <- nrow(catch)
nages <- ncol(catch)

N <- F <- Z <- matrix(NA, nrow = nyears, ncol = nages, dimnames = dimnames(catch))

## Assume F in terminal year = 0.1 and M = 0.2
F[nyears,] <- 0.1
M <- matrix(c(rep(0.1, 2), rep(0.3, 8)), nrow = nyears, ncol = nages, byrow = TRUE, dimnames = dimnames(catch))
Z <- F + M

## Calculate N in terminal year
N <- catch * Z / (F * (1-exp(-Z)))

verbose <- TRUE

## Calculate N and F up to terminal year,
## assuming F[oldest] = avg(3 preceding ages)
for(y in (nyears-1):1)
{
  for(a in 1:(nages-1))
  {
    N[y, a] <- N[y + 1, a + 1] * exp(M[y, a]) + catch[y, a] * exp(M[y, a] / 2)
    F[y, a] <- log(N[y, a] / N[y + 1, a + 1]) - M[y, a]
  }
  F[y, nages] <- mean(F[y, nages-(1:3)])
  Z[y, ] <- F[y, ] + M[y, ]
  N[y, nages] <-
    catch[y, nages] * Z[y, nages] / (F[y, nages] * (1 - exp(-Z[y, nages])))

  # little but of debuging type code
  if (verbose) {
    print(N)
    ans <- readline("press key to continue (Q to run to end)")
    if (ans == "Q") verbose <- FALSE
  }
}
