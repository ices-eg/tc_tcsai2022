# Calculating Yeild per recruit reference points in R

# Colin Millar

# This script shows the steps followed to Calculate
# yeild per recruit reference points from the a VPA fit to
# the data in 'bluefin_catage.dat'

#==============================================================================
#
# Use results from VPA assignment to estimate Fmax and F40% for
# western Atlantic bluefin tuna.
#
#   M=0.14
#   Assume 2007 selectivity pattern from VPA (less influenced by 2013 F assumption).
#   Maturity at age 9 (0% mature at ages 1-8, 100% mature at ages 9+; ICCAT 2014)
#   Growth k=0.089 Linf = 315, t0 = -1.13
#   Weight = 0.00002861 x Length^2.929
#
#==============================================================================

#------------------------------------------------------------------------------
# first lets rerun the assessment from yesterday:
#------------------------------------------------------------------------------

# get vpa function
source("../06_VPA/vpa_function.R")

# Read catch at age, change units
catage <- read.table("../06_VPA/bluefin.dat", header=TRUE, check.names=FALSE, row.names=1)
catage <- catage / 1000

# Run model
vpafit <- vpa(catage, 0.14, 0.1, 5)

# inspect contents
str(vpafit)


#------------------------------------------------------------------------------
# (1) Calculate partial recruitment for ages 1 to 20: Pa=Fa/Ffull, assuming P16+=P15
#------------------------------------------------------------------------------

# create a container for partial F
a <- 1:20
Pa <- numeric(20)

# fill in ages 1 to 15
Pa[1:15] <- vpafit$F["2007",]
Pa[16:20] <- vpafit$F["2007",15]
# scale
Pa <- Pa / max(Pa)

#------------------------------------------------------------------------------
# (2) Calculate mean lengh at age: La= L∞[1-e-k(a-to)]
#   Growth k=0.089 Linf = 315, t0 = -1.13
#------------------------------------------------------------------------------

# set up growth parameters
Linf <- 315; k <- 0.089; t0 <- -1.13

# calculate mean length at age
La <- Linf * (1 - exp(-k * (a - t0)))

#------------------------------------------------------------------------------
# (3) Calculate mean weight at age: wa = alpha * L^beta
#   Weight = 0.00002861 x Length^2.929
#------------------------------------------------------------------------------

wa <- 0.00002861 * La^2.929

#------------------------------------------------------------------------------
# (4) Assume an arbitrary number of recruits (N1=R=1000)
#------------------------------------------------------------------------------

R <- 1000

#------------------------------------------------------------------------------
# (5) Assume F=0 (for now)
#------------------------------------------------------------------------------

F <- 0

#------------------------------------------------------------------------------
# (6) Calculate abundance at age: Na+1=Na exp(-(Pa*F + M))
#     M = 0.14
#------------------------------------------------------------------------------

# set up N vector
Na <- numeric(20)

# define M
M <- 0.14

# assign age 1 and propagate down the ages
Na[1] <- R
for(i in seq_along(Na)[-1]) {
  Na[i] <- Na[i-1] * exp(-(Pa[i-1] * F + M))
}

# have a quick peek
plot(a, Na, type = "b", ylim=c(0,1000),
     main = "Population when F=0 (no fishing)")

#------------------------------------------------------------------------------
# (7) Calculate catch at age:
#------------------------------------------------------------------------------

Ca <- Na *
      Pa * F / (Pa * F + M) *
      (1 - exp(-M - Pa*F))

#------------------------------------------------------------------------------
# (8) Calculate yield
#------------------------------------------------------------------------------

Y <- sum(Ca * wa)

#------------------------------------------------------------------------------
# (9) Calculate yield per recruit:
#------------------------------------------------------------------------------

YPR <- Y / R

#------------------------------------------------------------------------------
# (10) Calculate spawning biomass
#    Maturity at age 9
#    (0% mature at ages 1-8,
#     100% mature at ages 9+; ICCAT 2014)
#------------------------------------------------------------------------------

# set up maturity
ma <- numeric(20)
ma[1:8] <- 0
ma[9:20] <- 1

# calculate ssb
SSB <- sum(Na * ma * wa)

#------------------------------------------------------------------------------
# (11) Calculate spawning biomass per recruit
#------------------------------------------------------------------------------

SPR <- SSB / R

#------------------------------------------------------------------------------
# (12) Save YPR and SPR results for the assumed F.
#------------------------------------------------------------------------------

# okay, so now we need to put this together into a recipe (function)
# In the simplest case we will need:
#   F
#
#   and we will 'hard code' the other options

ypr <- function(F) {

  # first lets set up all the things we need
  #-----------------------------------------------------
  a <- 1:20
  Pa <- numeric(20)
  Na <- numeric(20)

  # define M
  M <- 0.14

  # set up growth parameters
  Linf <- 315; k <- 0.089; t0 <- -1.13
  # weight params
  alpha <- 0.00002861; beta <- 2.929

  # set up maturity
  ma <- numeric(20)
  ma[1:8] <- 0
  ma[9:20] <- 1

  # now the modelling:
  #-----------------------------------------------------

  # (1) Calculate partial recruitment for ages 1 to 20: Pa=Fa/Ffull, assuming P16+=P15
  Pa[1:15] <- vpafit$F["2007",]
  Pa[16:20] <- vpafit$F["2007",15]
  # scale
  Pa <- Pa / max(Pa)

  # (2) Calculate mean lengh at age: La= L∞[1-e-k(a-to)]
  La <- Linf * (1 - exp(-k * (a - t0)))

  # (3) Calculate mean weight at age: wa = alpha * L^beta
  wa <- alpha * La^beta

  # (4) Assume an arbitrary number of recruits (N1=R=1000)
  R <- 1000

  # (6) Calculate abundance at age: Na+1=Na exp(-(Pa*F + M))
  Na[1] <- R
  for(i in seq_along(Na)[-1]) {
    Na[i] <- Na[i-1] * exp(-(Pa[i-1] * F + M))
  }

  # (7) Calculate catch at age:
  Ca <- Na *
        Pa * F / (Pa * F + M) *
        (1 - exp(-M - Pa*F))

  # (8) Calculate yield
  Y <- sum(Ca * wa)

  # (9) Calculate yield per recruit:
  YPR <- Y / R

  # (10) Calculate spawning biomass
  SSB <- sum(Na * ma * wa)

  # (11) Calculate spawning biomass per recruit
  SPR <- SSB / R

  # format the results and return
  out <- c(YPR, SPR)
  names(out) <- c("YPR", "SPR")
  out
}


# test the function
ypr(F=0)

# should be the same as:
c(YPR, SPR)

#------------------------------------------------------------------------------
# (13) Increase F by 0.1
#------------------------------------------------------------------------------

ypr(0.1)

#------------------------------------------------------------------------------
# (14) Repeat steps 6-12 up to F=1.
#------------------------------------------------------------------------------

# set up all the Fs we want to try
?seq
Fsteps <- seq(0, 1, by = 0.1)

# calculate YPR and SPR for each F
results <- sapply(Fsteps, ypr)
results <- data.frame(Fsteps, t(results))

#------------------------------------------------------------------------------
# (15) Plot YPR and SPR as a function of F
#------------------------------------------------------------------------------

# two plots on one page
par(mfrow = c(2,1))

# ypr plot
plot(Fsteps, results$YPR, type = "l",
     ylab = "Yeild per recruit", las = 1)

# spr plot
plot(Fsteps, results$SPR, type = "l",
     ylab = "Spawners per recruit", las = 1)

#------------------------------------------------------------------------------
# (16) Find Fmax
#------------------------------------------------------------------------------

# create a quick function, suitable for an optimiser
ypr_optim <- function(F) {
  ypr(F)["YPR"]
}

# test it works!
ypr_optim(1)

# now find the F that maximises YPR
opt <- optimize(ypr_optim, interval = c(0,1), maximum = TRUE)
Fmax <- opt$maximum

#------------------------------------------------------------------------------
# (17) According to the VPA results and the YPR results, what was the
# status of the tuna fishery in 2007?
#------------------------------------------------------------------------------

plot(as.integer(rownames(vpafit$F)), vpafit$F[,3], type = "l",
     main = "Fmax", xlab="Year", ylab="Fishing mortality rate")
abline(h = Fmax, col = "blue")


#------------------------------------------------------------------------------
# (18) Extra Credit: What is F40%?
#------------------------------------------------------------------------------

# use linear interpolation to get the value of F when SPR is 0.4 SPR when F = 0
results
SPR0 <- results$SPR[results$Fsteps==0.0]
apprx <- approx(results$SPR, Fsteps, xout = 0.4 * SPR0)
F40 <- apprx$y

# spr plot
plot(Fsteps, results$SPR, type = "b", main = "F40%",
     xlab = "Fishing mortality rate", ylab = "Spawners per recruit", las = 1)
# a line showing F40
lines(c(F40, F40), c(0, 0.4 * SPR0), col = "red")

