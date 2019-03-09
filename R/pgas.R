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
PGAS <- function(param, y, prior, x0=0, M = 1000,
                 N = 100, resamplingMethod = "multi")
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

  q = rep(0, M)
  r = rep(0, M)
  X = matrix(0, M, T)
  # Initialize the parameters
  prior.a <- prior[1]
  prior.b <- prior[2]
  q[1] = QInit
  r[1] = RInit
  # Initialize the state by running a PF
  param <- list(f = f, g = g, Q = q[1], R = r[1])
  result = APF(param = param, y = y, x0 = x0, N = 100)
  X[1, ] = result$x

  # Run MCMC loop
  for(k in 2:M)
  {
    # Sample the parameters (inverse gamma posteriors)
    err_q = X[k-1,2:T] - f(X[k-1,1:(T-1)], 1:(T-1))
    err_q = sum(err_q^2)
    q[k] = 1/rgamma(1, prior.a + (T-1)/2, prior.b + err_q/2)
    err_r = y - g(X[k-1,])
    err_r <- sum(err_r^2)
    r[k] = 1/rgamma(1, prior.a + T/2, prior.b + err_r/2)
    # Run CPF-AS
    param <- list(f = f, g = g, Q = q[k], R = r[k])
    result = CPF_AS(param = param, y = y, x0 = x0, x_ref = X[k-1,], N = N)
    X[k, ] = result$x
  }
  return(list(q = q, r = r, x = X))
}
