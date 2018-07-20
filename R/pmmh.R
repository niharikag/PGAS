library(invgamma)

PMMH <- function(numMCMC, y, prior, prop, N, qinit, rinit, q0, r0)
# Runs the PMMH algorithm,
#
#   C. Andrieu, A. Doucet and R. Holenstein, "Particle Markov chain Monte
#   Carlo methods" Journal of the Royal Statistical Society: Series B,
#   2010, 72, 269-342.
#
# The function returns the sample paths of (q, r, x_{1:T}).
{
  T = length(y)
  q = matrix(0, numMCMC, 1)
  r = matrix(0, numMCMC,1)
  X = matrix(0, numMCMC,T)
  loglik = matrix(0, numMCMC,1)
  # Initialize the parameters
  prior.a <- prior[1]
  prior.b <- prior[2]
  prop.sigma_q <- prop[1]
  prop.sigma_r <- prop[2]
  q[1] = qinit
  r[1] = rinit
  # Initialize the state & likelihood by running a PF
  param = c(q[1], r[1])
  res = particleFilter(param, y, N)
  # Draw
  J <- which(runif(1) < cumsum(res$w[,T]))[1]
  X[1, ] = res$particles[J,]
  loglik[1] <- res$logLikelihood
  # Run MCMC loop

  for(k in 2:numMCMC)
  {
    # Propose a parameter
    q_prop = q[k-1] + prop.sigma_q*rnorm(1)
    r_prop = r[k-1] + prop.sigma_r*rnorm(1)

    if(q_prop <= 0 || r_prop <= 0) # Prior probability 0, reject
    {
        accept = FALSE
    }
    else
    {
        # Run a PF to evaluate the likelihood
        param <- c(q_prop, r_prop)
        res = particleFilter(param, y, N)
        J <- which(runif(1) < cumsum(res$w[,T]))[1]
        X[k, ] = res$particles[J,]
        #X[k, ] = res$xHatFiltered
        loglik_prop <- res$logLikelihood

        acceptprob = exp(loglik_prop - loglik[k-1]) # Likelihood contribution
        acceptprob = acceptprob *  # Prior contribution
          dinvgamma(q_prop,prior.a,prior.b)*dinvgamma(r_prop,prior.a,prior.b) /
          (dinvgamma(q[k-1],prior.a,prior.b)*dinvgamma(r[k-1],prior.a,prior.b))
        q[k] = q_prop
        r[k] = r_prop
        loglik[k] = loglik_prop

        accept = runif(1) < acceptprob
    }

    if(!accept)
    {
        q[k] = q[k-1]
        r[k] = r[k-1]
        X[k,] = X[k-1,]
        loglik[k] = loglik[k-1]
    }
  }
  return(list(q = q, r = r, x = X))
}
