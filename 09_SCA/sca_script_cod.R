source("09_SCA/sca_function.R")

## Read data
C <-
  as.matrix(read.table("09_SCA/nscod_catage.dat",
    header = TRUE,
    check.names = FALSE, row.names = 1
  ))
I <- as.matrix(read.table("09_SCA/nscod_survey.dat", header = TRUE,
                          check.names = FALSE, row.names = 1))
M <- as.matrix(read.table("09_SCA/nscod_natmort.dat", header = TRUE,
                          check.names = FALSE, row.names = 1))
data <- list(C = C, I = I, M = M)

## Set initial parameter values

par <- c(logNa = rep(8, ncol(C)),
         logNt = rep(8, nrow(C)),
         logFa = rep(0, ncol(C) - 1),
         logFt = rep(0, nrow(C)),
         logQ = rep(-5, ncol(I)))

################################################################################

## Fit model

sca(par, data, full = TRUE)
sca(par, data)


opt1 <- optim(par, sca, data = data)
opt1

optim(par, sca,
  data = data,
  control = list(maxit = 1000)
)


opt2 <- optim(par, sca, data = data, method = "BFGS")
opt2

# or we can use nlminb
# ?nlminb
opt4 <- nlminb(par, sca, data = data)
opt4

opt5 <- nlminb(par, sca, data = data,
               control = list(eval.max = 1000, iter.max = 1000))
opt5

# summarise fits to check objective value and convergence
opt4$value <- opt4$objective
opt5$value <- opt5$objective
opts <- list(opt1 = opt1, opt2 = opt2, opt3 = opt3, opt4 = opt4, opt5 = opt5)
sapply(opts,
       function(x) {
         c(value = x$value,
           convergence = x$convergence)
       })

# lets go with BFGS and maximum iterations of 1000

## final run

run <- optim(par = par, fn = sca, data = data, method = "BFGS", hessian = TRUE)

run

# evaluate the model at the optimised parameter values
predictions <- sca(run$par, data, full = TRUE)

## View results

par(mfrow = c(2, 2))

## 1 Population
round(predictions$N)
Year <- as.integer(rownames(predictions$N))
plot(apply(predictions$N, 2, median), ylim = c(0, 400),
     yaxs = "i", type = "l", lty = 3,
     main = "Population in 2016 (bars) vs.\n median population (line)",
     xlab = "Age", ylab = "Individuals")
points(c(tail(predictions$N, 1)), type = "h", lwd = 6)

## 2 Recruitment
barplot(predictions$N[, 1], ylab = "Individuals at age 1", main = "Recruitment")

## 3 Selectivity
round(predictions$F, 2)
plot(colMeans(predictions$F) / max(colMeans(predictions$F)),
     ylim = c(0, 1.05), yaxs = "i", type = "l",
     main = "Selectivity", xlab = "Age", ylab = "Average F at age")

## 4 Fbar
Fbar2.4 <- rowMeans(predictions$F[, 2:4])
plot(Year[-length(Year)], Fbar2.4,
     ylim = c(0, 1.2), yaxs = "i", type  =  "l",
     main = "Fbar (2-4)", ylab = "Average F at ages 2-4")


# get some errors out (WARNING, probably not quite right, but here as an example)
library(MASS)
Sigma <- solve(2 * run$hessian)
par_sim <- mvrnorm(1000, run$par, Sigma)

sims <- lapply(1:1000, function(i) sca(par_sim[i,], data, full = TRUE))

Fbar2.4_sim <-
  sapply(
    sims,
    function(x) {
      rowMeans(x$F[, 2:4])
    }
  )

# plot a few curves to see the uncertainty in the relationship
matplot(Year[-length(Year)], Fbar2.4_sim[, sample(1:ncol(Fbar2.4_sim), 100)],
  type = "l", lty = 1, col = grey(0.5, alpha = 0.5),
  ylim = c(0, max(Fbar2.4_sim)), # set y limits
  main = "Fbar (2-4)", ylab = "Average F at ages 2-4", xlab = "Year"
)

# overlay fit
lines(Year[-length(Year)], Fbar2.4, lty = 1, lwd = 2)

# add confidence intervals
lines(Year[-length(Year)], apply(Fbar2.4_sim, 1, quantile, 0.025), col = "red", lty = 2)
lines(Year[-length(Year)], apply(Fbar2.4_sim, 1, quantile, 0.975), col = "red", lty = 2)
