######################################################################################
#
# Author : Niharika Gauraha
#          Uppsala University, Uppsala
#          Email : niharika.gauraha@farmbio.uu.se
#
#
#####################################################################################

# Conditional particle filter with Ancestor Sampling or
# Conditional SMC with AS (CSMC-AS)
# Input:
#   param   - state parameters
#   y       - measurements
#   x0      - initial state
#   x_ref   - reference trajecory
#   N       - number of particles
# Output:
#   x_star  - sample from target distribution


CPF_AS <- function(param, y, x0, x_ref, N = 100)
{
  # Stop, if input parameters are NULL
  if(is.null(param) || is.null(y) || is.null(x0) || is.null(x_ref))
  {
    stop("Error: the input parameters are NULL")
  }

  # Number of states
  T <- length(y)
  #Initialize the state parameters
  f <- param$f # state transition function
  g <- param$g # tranfer function
  Q <- param$Q # process noise variance
  R <- param$R # measurement noise variance

  # Initialize variables
  particles <- matrix(0, nrow = N, ncol = T)
  normalisedWeights <- matrix(0, nrow = N, ncol = T)

  # Init state, at t=0
  particles[, 1] = x0  # Deterministic initial condition
  particles[N, 1] = x_ref[1]  # Set the Nth particle to the reference particle

  # weighting step at t=0
  logweights = dnorm(y[1], mean = g(particles[,1]), sd = sqrt(R), log = TRUE)
  const = max(logweights)  # Subtract the maximum value for numerical stability
  weights = exp(logweights - const)
  normalisedWeights[,1] = weights/sum(weights)  # Save the normalized weights
  newAncestors = 1:N

  for(t in 2:T)
  {
    newAncestors <- multinomial.resample(normalisedWeights[, t - 1])
    xpred <- f(particles[, t-1], t-1)

    # propogation step
    particles[,t] = xpred[newAncestors] + sqrt(Q)*rnorm(N)
    particles[N,t] = x_ref[t] # set the nth particle deterministically

    #Ancestor sampling
    m = dnorm(x_ref[t], mean = xpred, sd = sqrt(Q), log = TRUE)
    const = max(m)  # Subtract the maximum value for numerical stability
    weights = exp(m - const)
    weights <- normalisedWeights[,t-1]*weights
    w_as = weights/sum(weights)  # Save the normalized weights
    newAncestors[N] <-  sample(1:N, prob=w_as, size= 1) # sample the nth ancestor
    # update the trajectory
    particles[, t-1] = particles[newAncestors, t-1]

    # weighting step
    logweights = dnorm(y[t], mean = g(particles[,t]), sd = sqrt(R), log = TRUE)
    const = max(logweights)  # Subtract the maximum value for numerical stability
    weights = exp(logweights - const)
    normalisedWeights[,t] = weights/sum(weights)  # Save the normalized weights

  }

  #return(list(particles=particles, normalisedWeights=normalisedWeights))
  J <- sample(1:N, prob=normalisedWeights[,T], size= 1)
  #J <- which(runif(1) < cumsum(normalisedWeights[,T]))[1]

  return(particles=particles[J,])
}

# Given theta estimate states using PGAS
PGAS <- function(param, y, x0, x_ref, N = 100, iter = 100)
{
  x_star = x_ref
  for (i in 1:iter) {
    x_star = CPF_AS(param = param, y = y, x_ref = x_ref, x0 = 0, N = 100)
    x_ref = x_star
  }
  return(x_star)
}
