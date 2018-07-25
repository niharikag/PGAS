# Implementation of PGAS algorithm introduced in
#
#  Lindsten, Fredrik, Michael I. Jordan, and Thomas B. Schön.
#    "Particle Gibbs with ancestor sampling."
#    The Journal of Machine Learning Research 15.1 (2014): 2145-2184.
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
PGAS <- function(param, y, x0, prior,
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
  r = matrix(0, M, 1)
  X = matrix(0, M, T)
  # Initialize the parameters
  prior.a <- prior[1]
  prior.b <- prior[2]
  q[1] = QInit
  r[1] = RInit
  # Initialize the state by running a PF
  param <- list(f = f, g = g, Q = q[1], R = r[1])
  res = conditionalParticleFilter(param = param, y = y, x0 = x0, X = X[1,], N = N)
  J <- which(runif(1) < cumsum(res$w[,T]))[1]
  X[1, ] = res$particles[J,]
  # Run MCMC loop
  for(k in 2:M)
  {
    # Sample the parameters (inverse gamma posteriors)
    err_q = X[k-1,2:T] - f(X[k-1,1:T-1], 1:(T-1))
    err_q = sum(err_q^2)
    q[k] = rinvgamma(1, prior.a + (T-1)/2, prior.b + err_q/2)
    err_r = y - g(X[k-1,])
    err_r <- sum(err_r^2)
    r[k] = rinvgamma(1, prior.a + T/2, prior.b + err_r/2)
    # Run CPF-AS
    param <- list(f = f, g = g, Q = q[k], R = r[k])
    res = conditionalParticleFilter(param = param, y = y, x0 = x0, X = X[k-1,], N = 100)
    J <- which(runif(1) < cumsum(res$w[,T]))[1]
    X[k, ] = res$particles[J,]
  }
  return(list(q = q, r = r, x = X))
}
