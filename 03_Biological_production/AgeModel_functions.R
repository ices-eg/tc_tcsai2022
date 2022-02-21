################################################################################
## Functions

AgeModel <- function(Ninit, M, Fmort, mat, w, Amax, Tmax, alpha, beta, v)
{
  N <- matrix(nrow=Amax, ncol=Tmax)
  SSB <- numeric(Tmax)

  ## Year 1
  N[1,1] <- Ninit
  for (a in 1:(Amax - 1)) {
    N[a + 1, 1] <- N[a, 1] * exp(-M[a] - Fmort[a])
  }
  SSB[1] <- sum(N[,1] * mat * w)

  ## Later years
  for(t in 1:(Tmax-1))
  {
    N[1,t+1] <- (alpha*SSB[t]) / (beta+SSB[t]) + v*runif(1,-0.5,0.5)
    for (a in 1:(Amax - 1)) {
      N[a + 1, t + 1] <- N[a, t] * exp(-M[a] - Fmort[a])
    }
    SSB[t+1] <- sum(N[,t+1] * mat * w)
  }

  list(alpha=alpha, beta=beta, N=N, SSB=SSB, Rec=N[1,])
}

RecPlot <- function(SR, Slim=NULL, Rlim=NULL)
{
  ## SR contains elements SSB, Rec, alpha, beta

  if(is.null(Slim))
    Slim <- c(0, 1.1*max(SR$SSB))
  if(is.null(Rlim))
    Rlim <- c(0, 1.1*max(SR$Rec))

  SSB.curve <- seq(Slim[1], Slim[2], length=100)
  Rec.curve <- (alpha*SSB.curve) / (beta+SSB.curve)

  opar <- par(mfrow=c(2,1), oma=c(0,0,0,3))
  plot(
    SR$SSB,
    ylim = Slim, type = "l", lwd = 2, col = "blue",
    xlab = "Year", ylab = "SSB (t)"
  )
  par(new=TRUE)
  plot(
    SR$Rec,
    ylim = Rlim, type = "h", lwd = 2, col = "red",
    ann = FALSE, axes = FALSE
  )
  axis(4)
  mtext("Recruitment at age 1 (1000s)", side=4, line=3)

  plot(
    SSB.curve, Rec.curve,
    xlim = Slim, ylim = Rlim, type = "l",
    xlab = "SSB (t)", ylab = "Recruitment at age 1 (1000s)", lwd = 2, col = "red"
  )
  points(SR$SSB, SR$Rec, type="o", pch=16, col="blue")
  par(opar)
}

################################################################################
## Simulation

Amax <- 15
Tmax <- 25
Ninit <- 1000
M <- rep(0.2, Amax)
Fmort <- c(0, 0, rep(0.3,Amax-2))
mat <- c(0, 0, rep(1,Amax-2))
w <- 10 * (1-exp(-0.5*((1:Amax)+0.1)))^3
alpha <- 1000
beta <- 1000

pop <- AgeModel(Ninit, M, Fmort, mat, w, Amax, Tmax, alpha, beta, v=1000)
round(pop$N)

RecPlot(
  AgeModel(Ninit, M, Fmort, mat, w, Amax, Tmax, alpha, beta, v = 1000),
  Slim = c(0, 20000), Rlim = c(0, 2000)
)
