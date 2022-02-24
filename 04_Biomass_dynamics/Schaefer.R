Schaefer <- function(par, data, verbose = FALSE) {
  # extract parameters
  r <- exp(par[["logr"]])
  K <- exp(par[["logK"]])
  Binit <- exp(par[["logBinit"]])
  q <- exp(par[["logq"]])

  # extract data
  year <- data$Year
  C <- data$Catch
  I <- data$Index

  # useful variables
  n <- length(year)
  B <- numeric(n)

  # create population from catch data and model parameters
  B[1] <- Binit
  for (i in 1:(n - 1))
  {
    B[i + 1] <- max(B[i] + r * B[i] * (1 - B[i] / K) - C[i], 1e-4)
  }
  Ifit <- q * B

  # get residuals and sums of squares
  res <- log(I) - log(Ifit) # log(I / Ifit)
  RSS <- sum(res^2)

  if (verbose) {
    # make useful summaries
    pars <- c(r = r, K = K, Binit = Binit, q = q)
    refpts <- c(HRmsy = 0.5 * r, Bmsy = 0.5 * K, MSY = 0.25 * r * K)

    list(B = B, HR = C / B, Ifit = Ifit, res = res, pars = pars, refpts = refpts, RSS = RSS)
  } else {
    RSS
  }
}

plot_shaefer <- function(fit, data, main) {
  oldpar <- par(mfrow = c(2, 2))

  plot(data$Year, fit$Ifit,
    ylim = c(0, max(fit$Ifit)), type = "l", lwd = 4,
    col = "gray", xlab = "Year", ylab = "Biomass index",
    main = paste0(main, ": Fit to data")
  )
  points(Index ~ Year, data)

  plot(data$Year, fit$B,
    type = "l", ylim = c(0, max(fit$B)), lwd = 2,
    xlab = "Year", ylab = "Biomass and catch", main = paste0(main, ": Biomass and catch")
  )
  points(Catch ~ Year, data, type = "h", lwd = 6)

  plot(data$Year, fit$HR,
    ylim = c(0, max(fit$HR)), type = "l",
    lwd = 2, xlab = "Year", ylab = "Harvest rate", main = paste0(main, ": Harvest rate")
  )

  plot(data$Year, fit$res,
    xlab = "Year", ylab = "log residuals",
    main = paste0(main, ": Residuals")
  )
  abline(h = 0)

  par(oldpar)
}

################################################################################
## Norther Shelf Haddock


haddock <- read.csv("04_Biomass_dynamics/haddock.csv", header=TRUE)

init <-
  c(
    logr = log(1.000),
    logK = log(8 * mean(haddock$Catch)),
    logBinit = log(4 * mean(haddock$Catch)),
    logq = log(haddock$Index[1] / (4 * mean(haddock$Catch)))
  )

Schaefer(par = init, haddock)
opt <- optim(init, Schaefer, data = haddock)
opt
fit <- Schaefer(opt$par, haddock, verbose = TRUE)


plot_shaefer(fit, haddock, main = "Haddock")

fit$pars
fit$refpts



# different initial values?
init2 <-
  c(
    logr = log(1),
    logK = log(950),
    logBinit = log(500),
    logq = log(3)
  )

opt2 <- nlminb(init2, Schaefer, data = haddock, control = list(eval.max = 1e6, iter.max = 1e6))
opt2

fit2 <- Schaefer(opt2$par, haddock, verbose = TRUE)
plot_shaefer(fit2, haddock, main = "Haddock")




################################################################################
## South Atlantic albacore

albacore <- read.table("04_Biomass_dynamics/albacore.dat", header=TRUE)
init <- c(logr=log(0.5), logK=log(200), logBinit=log(100), logq=log(0.5))

Schaefer(par=init, albacore)
optim(init, Schaefer, data=albacore)
est <- optim(init, Schaefer, data=albacore)$par
fit <- Schaefer(est, albacore, verbose=TRUE)

plot_shaefer(fit, albacore, main = "Albacore: Fit to data")

par(mfrow=c(2,2))

plot(albacore$Year, fit$Ifit, ylim=c(0,90), yaxs="i", type="l", lwd=4,
     col="gray", xlab="Year", ylab="Biomass index",
     main="Albacore: Fit to data")
points(Index~Year, albacore)

plot(albacore$Year, fit$B, type="l", ylim=c(0,300), yaxs="i", lwd=2,
     xlab="Year", ylab="Biomass and catch", main="Albacore: Biomass and catch")
points(Catch~Year, albacore, type="h", lwd=6)

plot(albacore$Year, fit$HR, ylim=c(0,0.35), yaxs="i", type="l",
     lwd=2, xlab="Year", ylab="Harvest rate", main="Albacore: Harvest rate")

fit$pars
fit$refpts

################################################################################
## Georges Bank winter flounder

flounder <- read.table("flounder.dat", header=TRUE)

K.init <- 8 * mean(flounder$Catch)
B.init <- 0.5 * K.init
q.init <- flounder$Index[1] / B.init
init <- c(logr=log(0.5), logK=log(K.init),
          logBinit=log(B.init), logq=log(q.init))

Schaefer(par=init, flounder)
optim(init, Schaefer, data=flounder)
optim(init, Schaefer, data=flounder, method="Nelder-Mead",
      control=list(maxit=1e5, reltol=1e-10))
nlminb(init, Schaefer, data=flounder, control=list(eval.max=1e4, iter.max=1e4))
est <- nlminb(init, Schaefer, data=flounder,
              control=list(eval.max=1e4, iter.max=1e4))$par
fit <- Schaefer(est, flounder, verbose=TRUE)

par(mfrow=c(2,2))

plot(flounder$Year, fit$Ifit, ylim=c(0,8), yaxs="i", lwd=4, col="gray",
     type="l", xlab="Year", ylab="Biomass index", main="Flounder: Fit to data")
points(Index~Year, flounder)

plot(flounder$Year, fit$B, type="l", ylim=c(0,15), yaxs="i", lwd=2,
     xlab="Year", ylab="Biomass and catch", main="Flounder: Biomass and catch")
points(Catch~Year, flounder, type="h", lwd=6)

plot(flounder$Year, fit$HR, ylim=c(0,0.6), yaxs="i", type="l",
     lwd=2, xlab="Year", ylab="Harvest rate", main="Flounder: Harvest rate")

t(fit$pars)
t(fit$refpts)


################################################################################
## Cod

cod <- read.csv("04_Biomass_dynamics/cod.csv", header=TRUE)
init <- c(logr     = log(1.000),
          logK     = log(8*mean(cod$Catch)),
          logBinit = log(4*mean(cod$Catch)),
          logq     = log(cod$Index[1] / (4*mean(cod$Catch)))
          )

Schaefer(par=init, cod)
opt <- optim(init, Schaefer, data=cod)
opt
fit <- Schaefer(opt$par, cod, verbose=TRUE)
exp(opt$par)

plot_shaefer(fit, cod, main = "Cod")

fit$pars
fit$refpts

# use contstraints ?
lower <- c(logr=log(.8), logK=log(500), logBinit=log(300), logq=log(.2))
upper <- c(logr=log(1.5), logK=log(1000), logBinit=log(700), logq=log(.9))

opt <- nlminb(init, Schaefer, data=cod, lower = lower, upper = upper, control=list(eval.max=1e4, iter.max=1e4))
opt

exp(opt$par)

fit <- Schaefer(opt$par, cod, verbose=TRUE)

plot_shaefer(fit, cod, main = "Cod")

fit$pars
fit$refpts


# initialise at nlminb solution?
opt <- optim(init, Schaefer, data=cod, method = "BFGS")
opt

opt <- optim(opt$par, Schaefer, data=cod, method = "Nelder-Mead")
opt

fit <- Schaefer(opt$par, cod, verbose=TRUE)

plot_shaefer(fit, cod, main = "Cod")

fit$pars
fit$refpts



# initialise at solver solution ?
init <- c(logr=log(1.02798251784361), logK=log(920.324078858206), logBinit=log(460.961477373726), logq=log(0.729428157796693))

opt <- optim(init, Schaefer, data=cod)
opt
fit <- Schaefer(opt$par, cod, verbose=TRUE)
exp(opt$par)

plot_shaefer(fit, cod, main = "Cod")

fit$pars
fit$refpts


# try a different optimiser - this is the best fit... but is it sensible?
set.seed(15) # this optimiser is simulation based and the fits can depend on the random seed
             # if the model does not fit that well
opt <- optim(init, Schaefer, data=cod, method = "SANN")
opt
fit <- Schaefer(opt$par, cod, verbose=TRUE)
exp(opt$par)

plot_shaefer(fit, cod, main = "Cod")

fit$pars
fit$refpts
