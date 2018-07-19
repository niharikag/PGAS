library(smcUtils)
PGAS <- function(numMCMC, y, prior, N, qinit, rinit, q0, r0)
# Runs the PGAS algorithm,
# F. Lindsten and M. I. Jordan T. B. Schon, "Ancestor sampling for
# Particle Gibbs", Proceedings of the 2012 Conference on Neural
# Information Processing Systems (NIPS), Lake Taho, USA, 2012.
#
# The function returns the sample paths of (q, r, x_{1:T}).
{
  T = length(y)
  q = matrix(0, numMCMC, 1)
  r = matrix(0, numMCMC, 1)
  X = matrix(0, numMCMC, T)
  # Initialize the parameters
  q[1] = qinit
  r[1] = rinit
  # Initialize the state by running a PF
  param = c(q[1], r[1])
  res = conditionalParticleFilter(param, y, N, X)
  X[1, ] <- res$xHatFiltered

  # Run MCMC loop
  for(k in 2:numMCMC)
  {
    # Sample the parameters (inverse gamma posteriors)
    err_q = X[k-1,2:T] - stateTransFunc(X[k-1,1:T-1], 1:(T-1))
    q[k] = rinvgamma(1, prior.a + (T-1)/2, prior.b + err_q*t(err_q)/2)
    err_r = y - transferFunc(X[k-1,])
    r[k] = rinvgamma(1, prior.a + T/2, prior.b + err_r*t(err_r)/2)
    # Run CPF-AS
    param <- c(q[k], r[k])
    res = conditionalParticleFilter(param, y, N, X[k-1, ])
    X[k, ] <- res$xHatFiltered
  }
  return(list(q = q, r = r, x = X))
}
