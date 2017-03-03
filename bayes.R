# Examples for Nature Ecology and Evolution manuscript
# October 2016

rm(list=ls())
# The data are the 12S rRNA alignment of human and orangutang, with 948 base pairs and 
# 90 differences (84 transitions and 6 transversions).
# Example data from Table 1.3, p.7 of Yang (2014) Molecular Evolution: A Statistical 
# Approach. Oxford University Press.

n <- 948 
ns <- 84 # 23+30+10+21
nv <- 6  # 1+0+2+0+2+1+0+0
  
# proportions of transversional and transitional differences
V <- nv/n; S <- ns/n

# log-likelihood function, f(D|d,k), using Kimura's (1980) substitution model
# see p.8 in Yang (2014)
k80.lnL <- function(d, k, n=948, ns=84, nv=6) {
  p0 <- .25 + .25 * exp(-4*d/(k+2)) + .5 * exp(-2*d*(k+1)/(k+2))
  p1 <- .25 + .25 * exp(-4*d/(k+2)) - .5 * exp(-2*d*(k+1)/(k+2))
  p2 <- .25 - .25 * exp(-4*d/(k+2))
  return ((n - ns - nv) * log(p0/4) + ns * log(p1/4) + nv * log(p2/4))
}

# function that returns the logarithm of the unscaled posterior, f(d) * f(k) * f(D|d,k)
# by default we set the priors as f(d) = Gamma(d | 2, 20) and f(k) = Gamma(k | 2, .1)
ulnP <- function(d, k, n=948, ns=84, nv=6, a.d=2, b.d=20, a.k=2, b.k=.1)
  return ( dgamma(d, a.d, b.d, log=TRUE) + dgamma(k, a.k, b.k, log=TRUE) + k80.lnL(d, k, n, ns, nv))
# prior, likelihood and posterior surfaces:
dim <- 100  # dimension for the plot
d.v <- seq(from=0, to=0.3, len=dim) # vector of d values
k.v <- seq(from=0, to=100, len=dim) # vector of k values
dk <- expand.grid(d=d.v, k=k.v)

par(mfrow=c(1, 3))
# prior surface, f(d) * f(k)
# d ~ Gamma(2, 20); k ~ Gamma(2, .1)
Pri <- matrix(dgamma(dk$d, 2, 20) * dgamma(dk$k, 2, .1), ncol=dim)
image(d.v, k.v, -Pri, las=1, col=heat.colors(50), main="Prior",xlab="distance, d", 
      ylab="kappa, k", cex.main=2.0, cex.lab=1.5, cex.axis=1.5)
contour(d.v, k.v, Pri, nlev=10, drawlab=FALSE, add=TRUE)
# Likelihood surface, f(D|d,k)
lnL <- matrix(k80.lnL(d=dk$d, k=dk$k), ncol=dim)
# for numerical reasons, we scale the likelihood to be 1 at the maximum,
# i.e. we subtract max(lnL)
L <- exp(lnL - max(lnL))  
image(d.v, k.v, -L, las=1, col=heat.colors(50), main="Likelihood", 
      xlab="distance, d", ylab="kappa, k", cex.main=2.0, cex.lab=1.5, cex.axis=1.5)
contour(d.v, k.v, L, nlev=10, drawlab=FALSE, add=TRUE)

# Posterior surface (unscaled), f(d) * f(k) * f(D|d,k)
Pos <- Pri * L
image(d.v, k.v, -Pos, las=1, col=heat.colors(50), main="Posterior", 
      xlab="distance, d", ylab="kappa, k", cex.main=2.0, cex.lab=1.5, cex.axis=1.5)
contour(d.v, k.v, Pos, nlev=10, drawlab=FALSE, add=TRUE)

# dev.print(pdf, file="~/Dropbox/work/BayesNatEcolEvol/densities.pdf")

# #############################################
# Markov chain Monte Carlo (MCMC)
# #############################################

# We now obtain the posterior distribution by MCMC sampling.
# In most practical problems, constant z cannot be calculated (either
# analytically or numerically), and so the MCMC algorithm becomes necessary.

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
  # N is the number of 'generations' the algorithm is run for.
  # w.d and w.k are the 'widths' of the proposal densities for d and k.
  d <- k <- numeric(N+1) # Here we will keep the visited values of d and k.
  d[1] <- init.d; k[1] <- init.k
  acc.d <- acc.k <- 0 # number of acceptances
  
  for (i in 1:N) {
    # here we use uniform densities with reflection to propose new d* and k*
    # propose and accept/reject new d
    d.prop <- abs(d[i] + runif(1, -w.d/2, w.d/2))
    p.ratio <- exp(ulnP(d.prop, k[i]) - ulnP(d[i], k[i]))
    alpha <- min(c(1, p.ratio))
    # if ru < alpha accept proposed d*:
    if (runif(1, 0, 1) < alpha) { d[i+1] <- d.prop; acc.d <- acc.d + 1 }
    # else reject it:
    else d[i+1] <- d[i]
    
    # propose and accept/reject new k
    k.prop <- abs(k[i] + runif(1, -w.k/2, w.k/2))
    p.ratio <- exp(ulnP(d[i], k.prop) - ulnP(d[i], k[i]))
    alpha <- min(c(1, p.ratio))
    # if ru < alpha accept proposed k*:
    if (runif(1, 0, 1) < alpha) { k[i+1] <- k.prop; acc.k  <- acc.k + 1 }
    # else reject it:
    else k[i+1] <- k[i]
  }
  
  # print out the proportion of times the proposals were accepted
  print("Acceptance proportions (d, k):")
  print(c(acc.d/N, acc.k/N))
  return (list(d=d, k=k)) # return vector of d and k visited during MCMC
}

# Test the algorithm
# Note the algorithm is random, so repetitions of it will give different results
dk.mcmc <- mcmcf(0.2, 20, 1e4, .12, 180) # system.time(mcmcf(0.2, 20, 1e4, .12, 180)) # 0.7s
dk.mcmc.b <- mcmcf(0.2, 20, 1e4, .12, 180) # run again to test for convergence
par(mfrow=c(1,3))
plot(dk.mcmc$d, ty='l', main="Trace of d", cex.main=2.0, cex.lab=1.5, cex.axis=1.5)  # trace plot of d
plot(dk.mcmc$k, ty='l', main="Trace of k", cex.main=2.0, cex.lab=1.5, cex.axis=1.5)  # trace plot of k
plot(dk.mcmc$d, dk.mcmc$k, pch='.', main="Joint sample of d and k", cex.main=2.0, cex.lab=1.5, cex.axis=1.5) # joint sample of d and k (points sampled from posterior surface)

# Efficiency of the chain:
# Values sampled in an MCMC chain are autocorrelated because new states
# are either the previous state or a modification of it.
# The efficiency of an MCMC chain is closely related to the autocorrelation.
# Intuitively, if the autocorrelation is high, the chain will be inefficient,
# i.e. we will need to run the chain for a long time to obtain a good
# approximation to the posterior distribution.
# The efficiency of a chain is defined as:
# eff = 1 / (1 + 2(r1 + r2 + r3 + ...))
# where ri is the correlation for lag i.

# run a very long chain (1e6 generations take about 1.2min in my MacBook Air)
# to calculate efficiency
dk.mcmc2 <- mcmcf(0.2, 20, 1e6, .12, 180)

# R's acf function (for AutoCorrelation Function) will calculate the autocorrelation
# of a sample, and produce the autocorrelation plot
par(mfrow=c(1,2))
acf(dk.mcmc2$d)
acf(dk.mcmc2$k)

# Define efficiency function
eff <- function(acf) 1 / (1 + 2 * sum(acf$acf[-1]))

# the efficiencies are roughly 22% and 20% for d and k respectively:
eff(acf(dk.mcmc2$d)) # [1] 0.2255753 # mcmcf(0.2, 20, 1e7, .12, 180)
eff(acf(dk.mcmc2$k)) # [1] 0.2015054 # mcmcf(0.2, 20, 1e7, .12, 180)

# Ineficient chain
# The window width for the d proposal density is too large,
# while it is too samll for k
dk.mcmc3 <- mcmcf(0.2, 20, 1e4, 3, 5)
dk.mcmc3.b <- mcmcf(0.2, 20, 1e4, 3, 5)

par(mfrow=c(1,2))
# because proposal width for d is too large, chain gets stuck at same values of d:
plot(dk.mcmc3$d, ty='l', main="Trace of d", cex.main=2.0, cex.lab=1.5, cex.axis=1.5)
# whereas proposal width for k is too small, so chain moves slowly:
plot(dk.mcmc3$k, ty='l', main="Trace of k", cex.main=2.0, cex.lab=1.5, cex.axis=1.5)

# run a very long inefficient chain (using 1e7 is better but slow)
dk.mcmc4 <- mcmcf(0.2, 20, 1e6, 3, 5)

# Efficiencies are roughly 1.5% for d, and 0.35% for k:
eff(acf(dk.mcmc4$d, lag.max=2e3)) # [1] 0.01530385  # mcmcf(0.2, 20, 1e7, 3, 5)
eff(acf(dk.mcmc4$k, lag.max=2e3)) # [1] 0.003493112 # mcmcf(0.2, 20, 1e7, 3, 5)

# All trace plots together, for efficient and inefficient chains:
par(mfrow=c(2,2))
plot(dk.mcmc$d, ty='l', las=1, ylim=c(.05,.2), 
     main="Trace of d, efficient chain", ylab="Distance, d", xlab='',
     cex.main=2.0, cex.lab=1.5)
plot(dk.mcmc3$d, ty='l', las=1, ylim=c(.05,.2),
     main="Trace of d, inefficient chain", ylab='', xlab='',
     cex.main=2.0, cex.lab=1.5)
plot(dk.mcmc$k, ty='l', las=1, ylim=c(0,100),
     main="Trace of k, efficient chain", ylab="ts/tv ratio, k", xlab='',
     cex.main=2.0, cex.lab=1.5)
plot(dk.mcmc3$k, ty='l', las=1, ylim=c(0,100),
     main="Trace of k, inefficient chain", ylab='', xlab='',
     cex.main=2.0, cex.lab=1.5)
# the efficient chains on the left have good mixing
# whereas the inefficient chains on the right have poor mixing
#dev.print(pdf, file="~/Dropbox/work/BayesNatEcolEvol/traces.pdf")

# Checking for convergence:
# histograms:
par(mfrow=c(1,2))
bks <- seq(from=0, to=105, by=1)
# Efficient chains converge quickly, and histograms overlap well:
hist(dk.mcmc.b$k, prob=TRUE, breaks=bks, border=NA, col=rgb(0, 0, 1, .5),
     las=1, xlab="kappa", xlim=c(0,100), ylim=c(0,.055), 
     main = "Histogram of efficient chains", cex.main=2.0, cex.lab=1.5)
hist(dk.mcmc$k, prob=TRUE, breaks=bks, border=NA, col=rgb(.5, .5, .5, .5), add=TRUE)
# Inefficient chains need to be run for a very long time to converge, here they
# were not run long enough, so that the histograms do not overlap well:
hist(dk.mcmc3.b$k, prob=TRUE, breaks=bks, border=NA, col=rgb(0, 0, 1, .5), las=1, 
     xlab="kappa", xlim=c(0,100), ylim=c(0,.055), 
     main = "Histogram of inefficient chains", cex.main=2.0, cex.lab=1.5)
hist(dk.mcmc3$k, prob=TRUE, breaks=bks, border=NA, col=rgb(.5, .5, .5, .5), add=TRUE)
#dev.print(pdf, file="~/Dropbox/work/BayesNatEcolEvol/hists.pdf")

# We can also calculate the posterior means
# posterior means (similar for efficient chains):
mean(dk.mcmc$d); mean(dk.mcmc.b$d)  
mean(dk.mcmc$k); mean(dk.mcmc.b$k)  
# posterior means (not so similar for the inefficient chains):
mean(dk.mcmc3$d); mean(dk.mcmc3.b$d) 
mean(dk.mcmc3$k); mean(dk.mcmc3.b$k)

# densities (smooth histograms)
par(mfrow=c(1,2)); adj <- 1.5
# Efficient chains:
plot(density(dk.mcmc.b$k, adj=adj), col="blue", las=1, 
     xlim=c(0, 100), ylim=c(0, .05), xaxs="i", yaxs="i", main = "Density for efficient chains", 
     cex.main=2.0, cex.lab=1.5)
lines(density(dk.mcmc$k, adj=adj), col="black")
# Inefficient chains:
plot(density(dk.mcmc3.b$k, adj=adj), col="blue", las=1, 
     xlim=c(0, 100), ylim=c(0, .05), xaxs="i", yaxs="i",
     main = "Density for inefficient chains", cex.main=2.0, cex.lab=1.5)
lines(density(dk.mcmc3$k, adj=adj), col="black")

#dev.print(pdf, file="~/Dropbox/work/BayesNatEcolEvol/kdensities.pdf")

