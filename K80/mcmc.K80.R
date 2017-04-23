# Clean up your R space
rm(list=ls())

#######################
# PART 1: Introduction
#######################

# The data are the 12S rRNA alignment of human and orangutang, with 948 base pairs and 
# 90 differences (84 transitions and 6 transversions).
# Example data from Table 1.3, p.7 of Yang (2014) Molecular Evolution: A Statistical 
# Approach. Oxford University Press.

n <- 948 # length of alignment in bp
ns <- 84 # total number of transitions (23+30+10+21)
nv <- 6  # total number of transversions (1+0+2+0+2+1+0+0)
 
# log-likelihood function, f(D|d,k), using Kimura's (1980) substitution model
# see p.8 in Yang (2014)

k80.lnL <- function(d, k, n=948, ns=84, nv=6) {
 
 p0 <- .25 + .25 * exp(-4*d/(k+2)) + .5 * exp(-2*d*(k+1)/(k+2))
 p1 <- .25 + .25 * exp(-4*d/(k+2)) - .5 * exp(-2*d*(k+1)/(k+2))
 p2 <- .25 - .25 * exp(-4*d/(k+2))
 
 return ((n - ns - nv) * log(p0/4) +
        ns * log(p1/4) + nv * log(p2/4))
}

dim <- 100  # dimension for the plot
d.v <- seq(from=0, to=0.3, len=dim) # vector of d values
k.v <- seq(from=0, to=100, len=dim) # vector of k values
dk <- expand.grid(d=d.v, k=k.v)
 
par(mfrow=c(1, 3))

# Prior surface, f(D)f(k)
Pri <- matrix(dgamma(dk$d, 2, 20) * dgamma(dk$k, 2, .1),
       ncol=dim)
 
image(d.v, k.v, -Pri, las=1, col=heat.colors(50),
      main="Prior", xlab="distance, d",
      ylab="kappa, k", cex.main=2.0,
      cex.lab=1.5, cex.axis=1.5)
 
contour(d.v, k.v, Pri, nlev=10, drawlab=FALSE, add=TRUE)
 
# Likelihood surface, f(D|d,k)
lnL <- matrix(k80.lnL(d=dk$d, k=dk$k), ncol=dim)
 
# for numerical reasons, we scale the likelihood to be 1
# at the maximum, i.e. we subtract max(lnL)
L <- exp(lnL - max(lnL))  
 
image(d.v, k.v, -L, las=1, col=heat.colors(50),
      main="Likelihood", xlab="distance, d",
      ylab="kappa, k", cex.main=2.0,
      cex.lab=1.5, cex.axis=1.5)
 
contour(d.v, k.v, L, nlev=10,
        drawlab=FALSE, add=TRUE) # creates a contour plot
# Unscaled posterior surface, f(d)f(k)f(D|d,k)
Pos <- Pri * L
 
image(d.v, k.v, -Pos, las=1, col=heat.colors(50),
      main="Posterior", xlab="distance, d",
      ylab="kappa, k", cex.main=2.0,
      cex.lab=1.5, cex.axis=1.5)
 
contour(d.v, k.v, Pos, nlev=10,
        drawlab=FALSE, add=TRUE)
        
        
##########################################
# PART 2: Markov Chain Monte Carlo (MCMC)
##########################################

# We now obtain the posterior distribution by MCMC sampling.
# In most practical problems, constant z cannot be calculated (either
# analytically or numerically), and so the MCMC algorithm becomes necessary.

# Function that returns the logarithm of the unscaled posterior:
# f(d) * f(k) * f(D|d,k)
# by default we set the priors as: 
# f(d) = Gamma(d | 2, 20) and f(k) = Gamma(k | 2, .1)
ulnPf <- function(d, k, n=948, ns=84, nv=6, a.d=2, b.d=20, a.k=2, b.k=.1) {
  # The normalizing constant in the prior densities can be ignored:
  lnpriord <- (a.d - 1)*log(d) - b.d * d
  lnpriork <- (a.k - 1)*log(k) - b.k * k
  
  # log-Likelihood (K80 model):
  expd1 <- exp(-4*d/(k+2));
  expd2 <- exp(-2*d*(k+1)/(k+2))
  p0 <- .25 + .25 * expd1 + .5 * expd2
  p1 <- .25 + .25 * expd1 - .5 * expd2
  p2 <- .25 - .25 * expd1
  lnL <- ((n - ns - nv) * log(p0/4) + ns * log(p1/4) + nv * log(p2/4))
  
  # Return unnormalised posterior:
  return (lnpriord + lnpriork + lnL)
}

# Draft MCMC algorithm:
# 1. Set initial states for d and k.
# 2. Propose a new state d* (from an appropriate proposal density).
# 3. Accept or reject the proposal with probability
#      min(1, p(d*)p(x|d*) / p(d)p(x|d))
#      If the proposal is accepted set d=d*, otherwise d=d.
# 4. Save d.
# 5. Repeat 2-4 for k.
# 6. Go to step 2.

mcmcf <- function(init.d, init.k, N, w.d, w.k) {
  # init.d and init.k are the initial states
  # N is the number of MCMC iterations.
  # w.d is the width of the sliding-window proposal for d.
  # w.k is the width of the sliding-window proposal for k.
  
  # We keep the visited states (d, k) in sample.d and sample.k 
  # for easy plotting. In practical MCMC applications, these 
  # are usually written into a file.
  sample.d <- sample.k <- numeric(N+1)
  
  # STEP 1: Initialise the MCMC chain
  d <- init.d;  sample.d[1] <- init.d
  k <- init.k;  sample.k[1] <- init.k
  ulnP <- ulnPf(d, k)
  acc.d <- 0;  acc.k <- 0 # number of acceptances
  
  for (i in 1:N) {
    # STEP 2: Propose new state d*
    # we use a uniform sliding window of width w with reflection
    # to propose new values d* and k* 
    # propose d* and accept/reject the proposal
    dprop <- d + runif(1, -w.d/2, w.d/2)
    # reflect if dprop is negative
    if (dprop < 0) dprop <- -dprop
    
    ulnPprop <- ulnPf(dprop, k)
    lnalpha <- ulnPprop - ulnP
    
    # STEP 3: Accept or reject the proposal
    # if ru < alpha accept proposed d*:
    if (lnalpha > 0 || runif(1) < exp(lnalpha)){
      d <- dprop;  ulnP <- ulnPprop; 
      acc.d  <- acc.d + 1
    }
    # else reject and stay where we are (so that nothing needs 
    # to be done).
    
    # STEP 4: Repeat steps 2-3 for k
    # propose k* and accept/reject the proposal
    kprop <- k + runif(1, -w.k/2, w.k/2)
    # reflect if kprop is negative
    if (kprop < 0) kprop <- -kprop
    
    ulnPprop <- ulnPf(d, kprop)
    lnalpha <- ulnPprop - ulnP
    # if ru < alpha accept proposed k*:
    if (lnalpha > 0 || runif(1) < exp(lnalpha)){
      k <- kprop;  ulnP <- ulnPprop
      acc.k  <- acc.k + 1
    }
    # else reject and stay where we are.
    
    # STEP 5: Save chain state
    sample.d[i+1] <- d;  sample.k[i+1] <- k
  }
  
  # print out the proportion of times
  # the proposals were accepted
  print("Acceptance proportions (d, k):")
  print(c(acc.d/N, acc.k/N))
  
  # return vector of d and k visited during MCMC
  
  return (list(d=sample.d, k=sample.k))
}

# Test run-time:
system.time(mcmcf(0.2, 20, 1e4, .12, 180)) # about 0.3s
# Run again and save MCMC output:
dk.mcmc <- mcmcf(0.2, 20, 1e4, .12, 180) 

par(mfrow=c(1,3))
# trace plot of d
plot(dk.mcmc$d, ty='l', xlab="Iteration", ylab="d", main="Trace of d")
# trace plot of k
plot(dk.mcmc$k, ty='l', xlab="Iteration", ylab="k", main="Trace of k")
 
# We can also plot the joint sample of d and k
# (points sampled from posterior surface)
plot(dk.mcmc$d, dk.mcmc$k, pch='.', xlab="d", ylab="k", main="Joint of d and k")

##########################################
# PART 3: Efficiency of the MCMC chain
##########################################

# Values sampled in an MCMC chain are autocorrelated because new states
# are either the previous state or a modification of it.
# The efficiency of an MCMC chain is closely related to the autocorrelation.
# Intuitively, if the autocorrelation is high, the chain will be inefficient,
# i.e. we will need to run the chain for a long time to obtain a good
# approximation to the posterior distribution.
# The efficiency of a chain is defined as:
# eff = 1 / (1 + 2(r1 + r2 + r3 + ...))
# where ri is the correlation for lag i.

# run a very long chain (1e6 generations take about
# 40s in my MacBook Air) to calculate efficiency
dk.mcmc2 <- mcmcf(0.2, 20, 1e6, .12, 180)
 
# R's acf function (for AutoCorrelation Function) 
par(mfrow=c(1,2))
acf(dk.mcmc2$d)
acf(dk.mcmc2$k)
 
# Define efficiency function
eff <- function(acf) 1 / (1 + 2 * sum(acf$acf[-1]))
 
# the efficiencies are roughly 22% and 20% for d and k respectively:
eff(acf(dk.mcmc2$d)) # [1] 0.2255753 # mcmcf(0.2, 20, 1e7, .12, 180)
eff(acf(dk.mcmc2$k)) # [1] 0.2015054 # mcmcf(0.2, 20, 1e7, .12, 180)


#To illustrate inefficient chains, we will run our MCMC again by using a proposal density
#with a too large step size for d, and another with a too small step size for k.

# The window width for the d proposal density is too large,
# while it is too small for k
dk.mcmc3 <- mcmcf(0.2, 20, 1e4, 3, 5)
 
par(mfrow=c(1,2))
# because proposal width for d is too large,
# chain gets stuck at same values of d:
plot(dk.mcmc3$d, ty='l', main="Trace of d", cex.main=2.0,
     cex.lab=1.5, cex.axis=1.5, ylab="d")
 
# whereas proposal width for k is too small,
# so chain moves slowly:
plot(dk.mcmc3$k, ty='l', main="Trace of k", cex.main=2.0,
     cex.lab=1.5, cex.axis=1.5, ylab="k")
     
dk.mcmc4 <- mcmcf(0.2, 20, 1e6, 3, 5)
 
# Efficiencies are roughly 1.5% for d, and 0.35% for k:
eff(acf(dk.mcmc4$d, lag.max=2e3)) # [1] 0.01530385  # mcmcf(0.2, 20, 1e7, 3, 5)
eff(acf(dk.mcmc4$k, lag.max=2e3)) # [1] 0.003493112 # mcmcf(0.2, 20, 1e7, 3, 5)

# plot the traces for efficient (part 2) and inefficient chains

par(mfrow=c(2,2))
 
plot(dk.mcmc$d, ty='l', las=1, ylim=c(.05,.2),
     main="Trace of d, efficient chain", xlab='',
     ylab="Distance, d", cex.main=2.0, cex.lab=1.5)
plot(dk.mcmc3$d, ty='l', las=1, ylim=c(.05,.2),
     main="Trace of d, inefficient chain", xlab='',
     ylab='', cex.main=2.0, cex.lab=1.5)
plot(dk.mcmc$k, ty='l', las=1, ylim=c(0,100),
     main="Trace of k, efficient chain",
     xlab='', ylab="ts/tv ratio, k",
     cex.main=2.0, cex.lab=1.5)
plot(dk.mcmc3$k, ty='l', las=1, ylim=c(0,100),
     main="Trace of k, inefficient chain",
     xlab='', ylab='', cex.main=2.0, cex.lab=1.5)
     
####################################
# PART 4: Checking for convergence
####################################

# We now illustrate the concept of burn-in
# We will run a chain with a high starting value,
# and another with a low starting value.
dk.mcmc.l <- mcmcf(0.01, 20, 1e4, .12, 180)
dk.mcmc.h <- mcmcf(0.4, 20, 1e4, .12, 180)

plot(dk.mcmc.l$d, xlim = c(1,200), ylim = c(0,0.4), ty = "l")
lines(dk.mcmc.h$d, col="red")

# We use the low chain to calculate the mean
# and sd of d. We could have used the high chain
# as well.
mean.d <- mean(dk.mcmc.l$d)
sd.d <- sd(dk.mcmc.l$d)
abline(h = mean.d + 2 * c(-sd.d, sd.d), lty = 2)
# The horizontal dashed lines indicate approximately
# the 95% CI. Notice how the chains move from either 
# the high or low starting values towards the
# stationary phase (the area within the dashed lines).
# The area before it reaches stationarity is the burn-in.


# Efficient chain (good proposal step sizes)
dk.mcmc.b <- mcmcf(0.05, 5, 1e4, .12, 180)
# Inefficient chain (bad proposal step sizes)
dk.mcmc3.b <- mcmcf(0.05, 5, 1e4, 3, 5)

# plot and compare histograms
par(mfrow=c(1,2))
bks <- seq(from=0, to=150, by=1)
 
hist(dk.mcmc.b$k, prob=TRUE, breaks=bks, border=NA,
     col=rgb(0, 0, 1, .5), las=1, xlab="kappa",
     xlim=c(0,100), ylim=c(0,.055))
 
hist(dk.mcmc$k, prob=TRUE, breaks=bks, border=NA,
     col=rgb(.5, .5, .5, .5), add=TRUE)
 
hist(dk.mcmc3.b$k, prob=TRUE, breaks=bks, border=NA,
     col=rgb(0, 0, 1, .5), las=1, xlab="kappa",
     xlim=c(0,100), ylim=c(0,.055))
 
hist(dk.mcmc3$k, prob=TRUE, breaks=bks, border=NA,
     col=rgb(.5, .5, .5, .5), add=TRUE)

# to calculate the posterior
# posterior means (similar for efficient chains):
mean(dk.mcmc$d); mean(dk.mcmc.b$d)
mean(dk.mcmc$k); mean(dk.mcmc.b$k)
 
# posterior means (not so similar for the inefficient chains):
mean(dk.mcmc3$d); mean(dk.mcmc3.b$d)
mean(dk.mcmc3$k); mean(dk.mcmc3.b$k)

# efficient chain, standard error of the means
sqrt(1/1e4 * var(dk.mcmc$d) / 0.23) # roughly 2.5e-4
sqrt(1/1e4 * var(dk.mcmc$k) / 0.20) # roughly 0.2
 
# inefficient chain, standard error of the means
sqrt(1/1e4 * var(dk.mcmc3$d) / 0.015) # roughly 9.7e-4
sqrt(1/1e4 * var(dk.mcmc3$k) / 0.003) # roughly 1.6

# plot densities (smoothed histograms) for the efficient and inefficient chains
par(mfrow=c(1,2)); adj <- 1.5
# Efficient chains:
plot(density(dk.mcmc.b$k, adj=adj), col="blue", las=1,
     xlim=c(0, 100), ylim=c(0, .05), xaxs="i", yaxs="i")
lines(density(dk.mcmc$k, adj=adj), col="black")
 
# Inefficient chains:
plot(density(dk.mcmc3.b$k, adj=adj), col="blue", las=1,
     xlim=c(0, 100), ylim=c(0, .05), xaxs="i", yaxs="i")
lines(density(dk.mcmc3$k, adj=adj), col="black")

