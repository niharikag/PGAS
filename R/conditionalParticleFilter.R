# Conditional particle filter with ancestor sampling
# Input:
#   param - state parameters
#   y - measurements
#   x0 - initial state
#   X - conditioned particles
#   N - number of particles
#   resamplingMethod - resampling methods:
#     multinomical and systematics resampling methods are supported
# Output:
#   particles - particles
#   logLikelihood - log likelihood
#   normalisedWeights - normalised Weights
conditionalParticleFilter <- function(param, y, x0, X, N = 100, resamplingMethod = "multi")
{
  # Stop, if input parameters are NULL
  if(is.null(param) || is.null(y) || is.null(x0) || is.null(X))
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
  ancestorIndices <- matrix(0, nrow = N, ncol = T)
  normalisedWeights <- matrix(0, nrow = N, ncol = T)
  logLikelihood <- 0

  # Init state
  ancestorIndices[, 1] <- 1:N
  particles[ ,1] <- 0
  normalisedWeights[, 1] = 1 / N
  particles[, 1] = x0  # Deterministic initial condition
  particles[N, 1] = X[1]  # Set the Nth particle according to the conditioning

  for(t in 1:T)
  {
    if(t != 1)
    {
      # Resample
      if(resamplingMethod == "systematic")
      {
        newAncestors <- systematic.resample(normalisedWeights[, t - 1])
      }
      else
      {
        newAncestors <- multinomial.resample(normalisedWeights[, t - 1])
      }

      xpred = f(particles[, t-1], t-1)
      particles[,t] = xpred[newAncestors] + sqrt(Q)*rnorm(N)
      particles[N,t] = X[t]  # Set the Nth particle according to the conditioning
      # Ancestor sampling
      m = exp(-1/(2*Q)*(X[t]-xpred)^2)
      w_as = normalisedWeights[,t-1]*m
      w_as = w_as/sum(w_as)
      newAncestors[N] = which(runif(1) < cumsum(w_as))[1]
      # Store the ancestor indices
      ancestorIndices[, t] <- newAncestors
    }
    # Compute importance weights
    ypred = g(particles[,t])
    logweights = -1/(2*R)*(y[t] - ypred)^2  # (up to an additive constant)
    const = max(logweights)  # Subtract the maximum value for numerical stability
    weights = exp(logweights - const)
    normalisedWeights[,t] = weights/sum(weights)  # Save the normalized weights
  }

  # Generate the trajectories from ancestor indices
  ind = ancestorIndices[,T]
  for(t in (T-1):1)
  {
    particles[,t] = particles[ind,t];
    ind = ancestorIndices[ind,t]
  }

  return(list(particles = particles,
              weights = normalisedWeights,
              logLikelihood = logLikelihood))
}

