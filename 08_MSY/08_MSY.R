# Calculating Maximum sustainable yeild reference points in R

# Colin Millar modified from Iago Mosqueira

# This script shows the steps followed to Calculate
# yeild per recruit reference points from the a VPA fit to
# the data in 'bluefin_catage.dat'

#==============================================================================
#
# Develop an age-based projection to estimate MSY for western Atlantic
# Bluefin tuna, using the stock-recruit relationship from
# assignment 1 and YPR-SPR from assignment 7.
#
#==============================================================================


#------------------------------------------------------------------------------
# first lets rerun the SR model from Monday:
#------------------------------------------------------------------------------

# load data from comma-separated file to data.frame
herring <- read.csv(file = '../02_Model_fitting/north_sea_herring_SR.csv', header = TRUE)
names(herring) <- c("yc", "ssb", "rec")

# define beverton and holt stock recruit model
bevholt <- function(b, S) {
  b[1] * S / (b[2] + S)
}

# define the ssq function to minimise
ssq <- function(b) {
  pred <- bevholt(b, herring$ssb)
  sum(
    (log(herring$rec) - log(pred))^2
  )
}

# run the optimiser
bhfit <- optim(c(1000000, 300000), ssq)


#------------------------------------------------------------------------------
# then lets rerun the assessment from yesterday:
#------------------------------------------------------------------------------

# get vpa function
source("../06_VPA/vpa_function.R")

# Read catch at age, change units
catage <- read.table("../06_VPA/bluefin.dat", header=TRUE, check.names=FALSE, row.names=1)
catage <- catage / 1000

# Run model
vpafit <- vpa(catage, 0.14, 0.1, 5)


#------------------------------------------------------------------------------
# (1) Start with a yearclass of recruits
#   (344687 at age-1, b1 asymptotic recruitment, Rmax).
#------------------------------------------------------------------------------

R <- bhfit$par[1]

#------------------------------------------------------------------------------
# (2) Assume natural mortality (M=0.14)
#------------------------------------------------------------------------------

M <- 0.14

#------------------------------------------------------------------------------
# (3) Fully-recruited F=0.5 (Fmax, for now)
#------------------------------------------------------------------------------

F = 0.5

#------------------------------------------------------------------------------
# (4) Assume a schedule of fishing mortality (table in previous slide)
#------------------------------------------------------------------------------

Fsteps <- seq(0, 1, by = 0.1)

#------------------------------------------------------------------------------
# (5) For year 1 of the projection assume that each yearclass in the first year was 344,687 at age-1) to calculate abundance at age for ages 2 to 20 (same as your YPR, SPR calculations, but with a different age-1 recruitment):
#------------------------------------------------------------------------------

# set up N and Pa vector
a <- 1:20
Pa <- numeric(20)
Na <- numeric(20)

# fill in ages 1 to 15
Pa[1:15] <- vpafit$F["2007",]
Pa[16:20] <- vpafit$F["2007",15]
# scale
Pa <- Pa / max(Pa)

# assign age 1 and propagate down the ages
Na[1] <- R
for(i in seq_along(Na)[-1]) {
  Na[i] <- Na[i-1] * exp(-(Pa[i-1] * F + M))
}

#------------------------------------------------------------------------------
# (6) Calculate SSB (/1000 to convert kg to tons!):
#------------------------------------------------------------------------------

# set up growth parameters
Linf <- 315; k <- 0.089; t0 <- -1.13
# calculate mean length at age
La <- Linf * (1 - exp(-k * (a - t0)))
# calculate mean weight at age
wa <- 0.00002861 * La^2.929

# set up maturity
ma <- numeric(20)
ma[1:8] <- 0
ma[9:20] <- 1

# calculate ssb (in tonnes)
SSB <- sum(Na * ma * wa) / 1000

#------------------------------------------------------------------------------
# (7) Calculate catch-at-age:
#------------------------------------------------------------------------------

Ca <- Na *
      Pa * F / (Pa * F + M) *
      (1 - exp(-M - Pa*F))

#------------------------------------------------------------------------------
# (8) Calculate yield (/1000 to convert kg to tons):
#------------------------------------------------------------------------------

Y <- sum(Ca * wa) / 1000

#------------------------------------------------------------------------------
# (9)   Calculate recruitment in year 2:
#   9.1. calculate a random deviation from predicted recruitment for each year
#        using the root mean square residual from assignment 1: =(rand()-0.5)*58721
#   9.2. copy those random values as values so that the values do not change.
#   9.3. calculate recruitment:
#        = max(0.001,'b1'*S/('b2'+S))+'dev_value‘
#            where b1=344687 and b2=42478
#------------------------------------------------------------------------------

# predict recruitement
Rnew <- bevholt(bhfit$par, SSB)

# simulate some noise / deviation
# one option is from our boostrap, but lets use
# a parametric version this time
# the model fit says that the residuals have a variance of:
#    bhfit$value = 4.200643
dev <- rnorm(1, mean = 0, sd = sqrt(bhfit$value))

# add the noise into the Recruitment estimate
Rnew <- Rnew * exp(dev)

#------------------------------------------------------------------------------
# (10) Calculate abundance out to age-20 as a function of
#      natural mortality, fishing mortality,
#      abundance in the previous age and the previous year:
#------------------------------------------------------------------------------

# set up next years N vector
N2a <- numeric(20)

# assign age 1 and propagate down the ages
N2a[1] <- Rnew
for(i in seq_along(Na)[-1]) {
  N2a[i] <- Na[i-1] * exp(-(Pa[i-1] * F + M))
}

#------------------------------------------------------------------------------
# (11) Repeat for 50 years:
#------------------------------------------------------------------------------

# a function would be usefull here:
project <- function(Na) {
  # set up next years N vector
  N2a <- numeric(20)

  # calculate ssb (in tonnes)
  SSB <- sum(Na * ma * wa) / 1000

  # simulate recruitment
  Rnew <- bevholt(bhfit$par, SSB) *
          exp(rnorm(1, mean = 0, sd = sqrt(bhfit$value)))

  # assign age 1 and propagate down the ages
  N2a[1] <- Rnew
  for(i in seq_along(Na)[-1]) {
    N2a[i] <- Na[i-1] * exp(-(Pa[i-1] * F + M))
  }

  N2a
}

# check it works
project(Na)

# Set up a matrix to hold N
Nay <- matrix(NA, 20, 50)

# fill in the first year:
Nay[,1] <- Na

# loop over and fill in later years
for (y in 2:50) {
  Nay[,y] <- project(Nay[,y-1])
}


#------------------------------------------------------------------------------
# At F=0.5:
#  What is the ‘equilibrium’ recruitment (average R for years 26-50)?
#  What is the ‘equilibrium’ spawning biomass (average SSB for years 26-50)?
#  What is the ‘equilibrium’ yield (average Y for years 26-50)?
#------------------------------------------------------------------------------

# equilibrium recruiment
mean(Nay[1,26:50])

# equilibrium SSB
mean( colSums(Nay * ma * wa) / 1000 )

# equilibrium yeild
Cay <- Nay *
       Pa * F / (Pa * F + M) *
       (1 - exp(-M - Pa*F))

mean( colSums(Cay * wa) / 1000 )


#------------------------------------------------------------------------------
# What fishing mortality produces the greatest equilibrium yield (Fmsy)?
# What is the estimate of MSY?
# What is the estimate of SSBmsy?
#------------------------------------------------------------------------------

# Now we need to build another function!
# In the simplest case we will need:
#   F
#
#   and we will 'hard code' the other options

msy <- function(F) {

  # first lets set up all the things we need
  #-----------------------------------------------------
  # set up Pa vector
  a <- 1:20
  Pa <- numeric(20)

  # fill in ages 1 to 15
  Pa[1:15] <- vpafit$F["2007",]
  Pa[16:20] <- vpafit$F["2007",15]
  # scale
  Pa <- Pa / max(Pa)

  # set up growth parameters
  Linf <- 315; k <- 0.089; t0 <- -1.13
  # calculate mean length at age
  La <- Linf * (1 - exp(-k * (a - t0)))
  # calculate mean weight at age
  wa <- 0.00002861 * La^2.929

  # set up maturity
  ma <- numeric(20)
  ma[1:8] <- 0
  ma[9:20] <- 1

  # Set up a matrix to hold N
  Nay <- matrix(NA, 20, 50)

  # now the modelling:
  #-----------------------------------------------------

  # (1) Start with a yearclass of recruits
  R <- bhfit$par[1]

  # (2) Assume natural mortality (M=0.14)
  M <- 0.14

  # (11) predict out for 50 years:
  # a (local) function would be usefull here:
  project <- function(Na) {
    # set up next years N vector
    N2a <- numeric(20)

    # calculate ssb (in tonnes)
    SSB <- sum(Na * ma * wa) / 1000

    # simulate recruitment
    Rnew <- bevholt(bhfit$par, SSB) *
            exp(rnorm(1, mean = 0, sd = sqrt(bhfit$value)))

    # assign age 1 and propagate down the ages
    N2a[1] <- Rnew
    for(i in seq_along(Na)[-1]) {
      N2a[i] <- Na[i-1] * exp(-(Pa[i-1] * F + M))
    }

    N2a
  }

  # fill in the first year:
  Nay[1,1] <- R
  for(i in 2:20) {
    Nay[i,1] <- Nay[i-1,1] * exp(-(Pa[i-1] * F + M))
  }

  # loop over and fill in later years
  for (y in 2:50) {
    Nay[,y] <- project(Nay[,y-1])
  }

  # summarise the equilibrium state of the stock
  # equilibrium recruiment
  R <- mean(Nay[1,26:50])

  # equilibrium SSB
  SSB <- mean( colSums(Nay * ma * wa) / 1000 )

  # equilibrium yeild
  Cay <- Nay *
         Pa * F / (Pa * F + M) *
         (1 - exp(-M - Pa*F))

  Y <- mean( colSums(Cay * wa) / 1000 )

  # format the results and return
  out <- c(R, SSB, Y)
  names(out) <- c("R", "SSB", "Y")
  out
}

# check it!
msy(0.5)


#------------------------------------------------------------------------------
# Reminder of questions:
#   What fishing mortality produces the greatest equilibrium yield (Fmsy)?
#   What is the estimate of MSY?
#   What is the estimate of SSBmsy?
#------------------------------------------------------------------------------

# create a quick function, suitable for an optimiser
msy_optim <- function(F) {
  msy(F)["Y"]
}

# test it works!
msy_optim(1)

# now find the F that maximises equilibrium yeild
opt <- optimize(msy_optim, interval = c(0,1), maximum = TRUE)
Fmsy <- opt$maximum


# what is MSY and SSBmsy?
msy(Fmsy)


