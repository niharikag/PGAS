# Utility file
# defines basic particle filter and
# conditional particle filter with ancestor sampling


### to be removed, should be passed as parameter
transferFunc <- function(xt)
{
  yt = xt^2/20
}
### to be removed, it should be passed as parameter
stateTransFunc <- function(xt,t)
{
  xt1 = 0.5*xt + 25*xt/(1+xt^2) + 8*cos(1.2*t)
  return(xt1)
}

particleFilter <- function(param, y, N = 100, resamplingMethod = "multi")
  # Particle filter
  # Input:
  #   y - measurements
  #   param - contains process noise variance and measurement noise variance
  #   N - number of particles
  #   resamplingMethod -  resampling method
  # Output:
  #   xHatFiltered - mean states
  #   logLikelihood - log likelihood
  #   normalisedWeights - normalised Weights
{
  T <- length(y)
  #Initialize the parameters
  #f <- param[1] , state transition function
  #h <- param[2], tranfer function
  Q <- param[1] # process noise variance
  R <- param[2] # measurement noise variance

  # Initialize variables
  particles <- matrix(0, nrow = N, ncol = T)
  ancestorIndices <- matrix(0, nrow = N, ncol = T)
  normalisedWeights <- matrix(0, nrow = N, ncol = T)
  xHatFiltered <-rep(0, T)
  logLikelihood <- 0

  ancestorIndices[, 1] <- 1:N
  particles[ ,1] <- 0 # Deterministic initial condition
  xHatFiltered[1] <- 0
  normalisedWeights[, 1] = 1 / N

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
      # Store the ancestor indices
      ancestorIndices[, t] <- newAncestors

      xpred = stateTransFunc(particles[, t-1],t-1)
      particles[,t] = xpred[newAncestors] + sqrt(Q)*rnorm(N)
      xHatFiltered[t] <- mean(particles[, t])
    }
    # Compute importance weights
    ypred = transferFunc(particles[,t])
    logweights = -1/2*log(2*pi*R) - 1/(2*R)*(y[t] - ypred)^2
    const = max(logweights) # Subtract the maximum value for numerical stability
    weights = exp(logweights-const)
    # Compute loglikelihood
    logLikelihood = logLikelihood + const + log(sum(weights)) - log(N)
    normalisedWeights[,t] = weights/sum(weights) # Save the normalized weights
  }

  return(list(xHatFiltered = xHatFiltered,
              logLikelihood = logLikelihood,
              particles = particles,
              weights = normalisedWeights))
}


#--------------------------------------------------------------------------
conditionalParticleFilter <- function(param, y, N, X, resamplingMethod = "multi")
{
  # Conditional particle filter with ancestor sampling
  # Input:
  #   y - measurements
  #   q - process noise variance
  #   r - measurement noise variance
  #   N - number of particles
  #   X - conditioned particles - if not provided, un unconditional PF is run
  conditioning = TRUE
  T <- length(y)
  #Initialize the parameters
  #f <- param[1] , state transition function
  #h <- param[2], tranfer function
  Q <- param[1] # process noise variance
  R <- param[2] # measurement noise variance

  # Initialize variables
  particles <- matrix(0, nrow = N, ncol = T)
  ancestorIndices <- matrix(0, nrow = N, ncol = T)
  normalisedWeights <- matrix(0, nrow = N, ncol = T)
  xHatFiltered <-rep(0, T)
  logLikelihood <- 0

  ancestorIndices[, 1] <- 1:N
  particles[ ,1] <- 0
  xHatFiltered[1] <- 0
  normalisedWeights[, 1] = 1 / N

  particles[, 1] = 0  # Deterministic initial condition
  particles[N, 1] = X[1]  # Set the 1st particle according to the conditioning

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
      # Store the ancestor indices
      ancestorIndices[, t] <- newAncestors

      xpred = stateTransFunc(particles[, t-1], t-1)
      particles[,t] = xpred[newAncestors] + sqrt(Q)*rnorm(N,1)
      if(conditioning)
      {
        particles[N,t] = X[t]  # Set the N:th particle according to the conditioning
        # Ancestor sampling
        m = exp(-1/(2*Q)*(X[t]-xpred)^2)
        w_as = normalisedWeights[,t-1]*m
        w_as = w_as/sum(w_as)
        newAncestors[N] = N
      }
      xHatFiltered[t] <- mean(particles[, t])
    }
    # Compute importance weights
    ypred = transferFunc(particles[,t])
    logweights = -1/(2*R)*(y[t] - ypred)^2  # (up to an additive constant)
    const = max(logweights)  # Subtract the maximum value for numerical stability
    weights = exp(logweights-const)
    normalisedWeights[,t] = weights/sum(weights)  # Save the normalized weights
  }

  return(list(xHatFiltered = xHatFiltered,
              logLikelihood = logLikelihood,
              particles = particles,
              weights = normalisedWeights))
}

