################################################################################
## 1  Mortality

Amax <- 15
age <- 1:Amax
Ninit <- 1000
Nvec <- numeric(Amax)
M <- rep(0.2, Amax)
Fmort <- c(0, 0, rep(0.3, Amax-2))

Nvec[1] <- Ninit
for (a in 1:(Amax-1)) {
  Nvec[a+1] <- Nvec[a] * exp(-M[a]-Fmort[a])
}


data.frame(age, Nvec, M, Fmort)

################################################################################
## 2  Growth

Winf <- 10
k <- 0.5
t0 <- -0.1
b <- 3

w <- Winf * (1-exp(-k*(age-t0)))^b
B <- Nvec * w

data.frame(age, Nvec, w, B)
sum(B)

################################################################################
## 3  Recruitment

mat <- c(0, 0, rep(1,Amax-2))
SSBvec <- B * mat
sum(SSBvec)

alpha <- 1000
beta <- 1000
Tmax <- 25
year <- 1:Tmax

N <- matrix(nrow=Amax, ncol=Tmax)
SSB <- matrix(nrow=Amax, ncol=Tmax)

N[,1] <- Nvec
SSB[,1] <- N[,1] * mat * w

N[1,2] <- (alpha*sum(SSB[,1])) / (beta+sum(SSB[,1]))
for(a in 1:(Amax-1))
  N[a+1,2] <- N[a,1] * exp(-M[a]-Fmort[a])
SSB[,2] <- N[,2] * mat * w

for(t in 1:(Tmax-1))
{
  N[1,t+1] <- (alpha*sum(SSB[,t])) / (beta+sum(SSB[,t]))
  for(a in 1:(Amax-1))
    N[a+1,t+1] <- N[a,t] * exp(-M[a]-Fmort[a])
  SSB[,t+1] <- N[,t+1] * mat * w
}

round(N)

rec_age <- 1
SR <-
  data.frame(
    year = year[1:(Tmax - rec_age)],
    SSB = colSums(SSB)[1:(Tmax - rec_age)],
    Rec = N[1, (rec_age+1):Tmax]
  )

SSB.curve <- seq(0, 12000, by=500)
Rec.curve <- (alpha*SSB.curve) / (beta+SSB.curve)

par(mfrow=c(2,1), oma=c(0,0,0,3))
plot(SR$year, SR$SSB, ylim=c(0,12000), type="l", lwd=2, col="blue", xlab="Year",
     ylab="SSB (t)")
lines(SR$year, SR$Rec * 10, lwd = 2, col = "red")
axis(4, at=seq(0,12000,by=2000), labels=seq(0,1200,by=200))
mtext("Recruitment at age 1 (1000s)", side=4, line=3)

plot(SSB.curve, Rec.curve, ylim=c(0,1000), type="l", xlab="SSB (t)",
     ylab="Recruitment at age 1 (1000s)", lwd=2, col="red")
points(SR$SSB, SR$Rec, type="o", pch=16, col="blue")

data.frame(SSB.curve, Rec.curve)

################################################################################
## 4  Recruitment variability

v <- 1000

## Year 1
N[,1] <- Nvec
SSB[,1] <- N[,1] * mat * w

## Later years
for(t in 1:(Tmax-1))
{
  N[1,t+1] <- (alpha*sum(SSB[,t])) / (beta+sum(SSB[,t])) + v*runif(1,-0.5,0.5)
  for(a in 1:(Amax-1))
    N[a+1,t+1] <- N[a,t] * exp(-M[a]-Fmort[a])
  SSB[,t+1] <- N[,t+1] * mat * w
}

round(N, 3)

rec_age <- 1
SR <-
  data.frame(
    year = year[1:(Tmax - rec_age)],
    SSB = colSums(SSB)[1:(Tmax - rec_age)],
    Rec = N[1, (rec_age+1):Tmax]
  )

SSB.curve <- seq(0, 15000, by=500)
Rec.curve <- (alpha*SSB.curve) / (beta+SSB.curve)

plot(year, SR$SSB, ylim=c(0,15000), type="l", lwd=2, col="blue", xlab="Year",
     ylab="SSB (t)")
lines(year, SR$Rec*10, lwd=2, col="red")
axis(4, at=seq(0,12000,by=2000), labels=seq(0,1200,by=200))
mtext("Recruitment at age 1 (1000s)", side=4, line=3)

plot(SSB.curve, Rec.curve, xlim=c(0,15000), ylim=c(0,1500), type="l",
     xlab="SSB (t)", ylab="Recruitment at age 1 (1000s)", lwd=2, col="red")
points(SR$SSB, SR$Rec, type="o", pch=16, col="blue")
