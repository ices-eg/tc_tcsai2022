source("sca_function.R")

## Read data

library(dplyr)

C <-
  read.table("nscod_catage.dat",
    header = TRUE,
    check.names = FALSE, row.names = 1
  ) %>%
  tibble()



I <- as.matrix(read.table("nscod_survey.dat", header = TRUE,
                          check.names = FALSE, row.names = 1))
M <- as.matrix(read.table("nscod_natmort.dat", header = TRUE,
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

opt2 <- optim(par, sca, data = data, method = "BFGS")
opt2

opt3 <- optim(par, sca, data = data,
              method = "BFGS", control = list(maxit = 1000))
opt3

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

run <- optim(par = par, fn = sca, data = data, method = "BFGS",
             control = list(maxit = 1000), hessian = TRUE)


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
