source("06_VPA/vpa_function.R")

## Read catch at age
catch <- read.csv("06_VPA/cod_catch.csv", header = TRUE, check.names = FALSE, row.names = 1)
Year <- as.numeric(row.names(catch))

catch <- cbind(catch[, 1:5], rowSums(catch[, -(1:5)]))

## Run model
model <- vpa(catch, M=0.2, Fterm=0.1, Fages=3)

## View results
par(mfrow=c(2,2))

## Fishing mortality
round(model$F, 3)

## Recruitment
barplot(model$N[,1], ylab="Recruitment at age 1", main="Recruitment (N age 1)")

## Selectivity
round(model$F, 3)
plot(colMeans(model$F)/max(colMeans(model$F)), ylim=c(0,1.05),
     type="l", main="Selectivity", xlab="Age", ylab="Mean F at age")

## Fbar
Fbar <- rowMeans(model$F)
plot(Year, Fbar, ylim=c(0, max(Fbar)), main="Mean F ages 1-10",
     ylab="Mean F (ages 1-10)", type="l")

## SSB

## Read weigths and maturity at age
wt <- read.csv("06_VPA/cod_weights.csv", header = TRUE, check.names = FALSE, row.names = 1)
mat <- read.csv("06_VPA/cod_maturity.csv", header = TRUE, check.names = FALSE, row.names = 1)

wt <- wt[, 1:6]
mat <- mat[, 1:6]

ssb <- rowSums(model$N * wt * mat) / 1000
plot(Year, ssb, ylim = c(0, max(ssb)), main="Spawning stock biomass (SSB)",
     ylab="SSB (kt)", type="l")
