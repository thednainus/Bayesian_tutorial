# Solution to EXERCISE 7:
lnL <- function(d, x, n)  x*log(3/4 - 3*exp(-4*d/3)/4) + (n - x) *log(1/4 + 3*exp(-4*d/3)/4)
r.lnprior <- function(r, mu=0.10/22.45) -r/mu              # -log(mu) term deleted
t.lnprior <- function(t, mu=22.45, sd=5.6) -((t-mu)/sd)^2  # constant deleted
lnL.rt <- function(r, t, x, n) lnL(d=r*t, x=x, n=n)
rt.lnpost <- function(r, t, x=90, n=948)  r.lnprior(r) + t.lnprior(t) + lnL.rt(r, t, x, n)


rt.mcmcf <- function(init.r, init.t, N, w.r, w.t) {
  # init.r and init.t are the initial states
  # N is the number of 'generations' the algorithm is run for.
  # w.r and w.t are the step sizes of the proposal densities.
  r <- t <- numeric(N+1)
  r[1] <- init.r; t[1] <- init.t
  rnow <- init.r; tnow <- init.t
  ap.r <- ap.t <- 0 # acceptance proportions
  for (i in 1:N) {
    # Propose, and accept or reject new r:
    rnew <- rnow + runif(1, -w.r/2, w.r/2)
    if(rnew < 0) rnew <- -rnew;
    lnalpha <- rt.lnpost(rnew, tnow, x=90, n=948) - rt.lnpost(rnow, tnow, x=90, n=948)
    # if ru < alpha accept proposal:
    if (lnalpha>0 || runif(1, 0, 1) < exp(lnalpha)) { rnow <- rnew; ap.r  <- ap.r + 1 }
    # else reject it, so that rnow = rnow.
    
    # Propose, and accept or reject new t:
    tnew <- tnow + runif(1, -w.t/2, w.t/2)
    if(tnew < 0) tnew <- -tnew;
    lnalpha <- rt.lnpost(rnow, tnew, x=90, n=948) - rt.lnpost(rnow, tnow, x=90, n=948)
    if (lnalpha > 0 || runif(1, 0, 1) < exp(lnalpha)) { tnow <- tnew; ap.t  <- ap.t + 1 }
    # else reject it so that tnow = tnow.
    r[i+1] <- rnow;  t[i+1] <- tnow;   # take the sample
  }
  # print out the acceptance proportions
  print(c(ap.r/N, ap.t/N))
  return (list(r=r, t=t)) # return vector of d's visited during MCMC
}

# Test the chain:
rt.1 <- rt.mcmcf(0.01, 40, 1e2, 0.001, 10)

plot(rt.1$r, rt.1$t, ty='b', pch=19, cex=.5, xlab="rate", ylab="time")

# Do a longer run and finetune:
rt.2 <- rt.mcmcf(0.01, 40, 1e4, 0.001, 10)
# If the acceptance proportion is high, step is too small and needs to be increased
# On the other hand, if acceptance proportion is low, decrease step size
# What step sizes above lead to 30% and 30% acceptance proportions?

# Look at the trace files:
par(mfrow=c(2,1))
plot(rt.2$r, ty='l', main="rate")
plot(rt.2$t, ty='l', main="time")

require(MASS)  # This package has additional statistical functions
zz <- kde2d(rt.2$r, rt.2$t)
par(mfrow=c(1,1))
image(zz, xlab="rate", ylab="time")
contour(zz, add=TRUE)

# Add a first few points to see MCMC progress:
lines(rt.2$r[1:200], rt.2$t[1:200], ty='b', pch=19, cex=.5, col="blue")