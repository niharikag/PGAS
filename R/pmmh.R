# The PMMH algorithm as defined in,
#
#   C. Andrieu, A. Doucet and R. Holenstein, "Particle Markov chain Monte
#   Carlo methods" Journal of the Royal Statistical Society: Series B,
#   2010, 72, 269-342.
#
# Input:
#   param - state parameters
#   y - measurements
#   x0 - initial state
#   M - number of MCMC runs
#   N - number of particles
#   resamplingMethod - resampling methods:
#     multinomical and systematics resampling methods are supported
# Output:
#       The function returns the sample paths of (q, r, x_{1:T})
PMMH <- function(param, y, x0, prior, prop,
                 M = 1000, N = 100, resamplingMethod = "multi")
{
  # Stop, if input parameters are NULL
  if(is.null(param) || is.null(y) || is.null(x0))
  {
    stop("Error: the input parameters are NULL")
  }

  # Number of states
  T <- length(y)
  #Initialize the state parameters
  f <- param$f # state transition function
  g <- param$g # tranfer function
  QInit <- param$Q # process noise variance
  RInit <- param$R # measurement noise variance

  q = matrix(0, M, 1)
  r = matrix(0, M,1)
  X = matrix(0, M,T)
  loglik = matrix(0, M,1)
  # Initialize the parameters
  prior.a <- prior[1]
  prior.b <- prior[2]
  prop.sigma_q <- prop[1]
  prop.sigma_r <- prop[2]
  q[1] = QInit
  r[1] = RInit
  # Initialize the state & likelihood by running a PF
  param <- list(f = f, g = g, Q = q[1], R = r[1])
  res = particleFilter(param = param, y = y, x0 = x0, N = N)
  # Draw
  J <- which(runif(1) < cumsum(res$w[,T]))[1]
  X[1, ] = res$particles[J,]
  loglik[1] <- res$logLikelihood
  # Run MCMC loop

  for(k in 2:M)
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
        param <- list(f = f, g = g, Q = q_prop, R = r_prop)
        res = particleFilter(param = param, y = y, x0 = x0, N = 100)
        J <- which(runif(1) < cumsum(res$w[,T]))[1]
        X[k, ] = res$particles[J,]
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
