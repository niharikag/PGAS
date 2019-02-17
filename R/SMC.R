# SMC (APF): Auxiliary particle filter


BPF <- function(param, y, x0, N = 100)
{
  # Stop, if input parameters are NULL
  if(is.null(param) || is.null(y) || is.null(x0) )
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
  # weighting step
  logweights = dnorm(y[1], mean = g(particles[,1]), sd = sqrt(R), log = TRUE)
  const = max(logweights)  # Subtract the maximum value for numerical stability
  weights = exp(logweights - const)
  normalisedWeights[,1] = weights/sum(weights)  # Save the normalized weights

  for (t in 2:T) {
    # resampling step
    newAncestors <- multinomial.resample(normalisedWeights[, t - 1])
    particles[, t-1] = particles[newAncestors, t-1]

    # propogation step
    particles[,t] = f(particles[, t-1], t-1) + sqrt(Q)*rnorm(N)

    # weighting step
    logweights = dnorm(y[t], mean = g(particles[,t]), sd = sqrt(R), log = TRUE)
    const = max(logweights)  # Subtract the maximum value for numerical stability
    weights = exp(logweights - const)
    normalisedWeights[,t] = weights/sum(weights)  # Save the normalized weights
  }

  J <- sample(1:N, prob=normalisedWeights[,T], size= 1)
  #J <- which(runif(1) < cumsum(normalisedWeights[,T]))[1]

  return(particles=particles[J,])
}


APF <- function(param, y, x0, N = 100)
{
  # Stop, if input parameters are NULL
  if(is.null(param) || is.null(y) || is.null(x0) )
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

  # weighting step
  logweights = dnorm(y[1], mean = g(particles[,1]), sd = sqrt(R), log = TRUE)
  const = max(logweights)  # Subtract the maximum value for numerical stability
  weights = exp(logweights - const)
  normalisedWeights[,1] = weights/sum(weights)  # Save the normalized weights

  for (t in 2:T) {
    # resampling step
    newAncestors <- multinomial.resample(normalisedWeights[, t - 1])


    # propogation step
    particles[,t] = f(particles[newAncestors, t-1], t-1) + sqrt(Q)*rnorm(N)
    num = dnorm(y[t], mean = g(particles[,t]), sd = sqrt(R))
    deno = dnorm(y[t-1], mean = g(particles[newAncestors,t-1]), sd = sqrt(R))

    particles[, t-1] = particles[newAncestors, t-1]

    #particles[,t] = particles[,t] * num/deno
    # weighting step
    logweights = dnorm(y[t], mean = g(particles[,t]), sd = sqrt(R), log = TRUE)
    const = max(logweights)  # Subtract the maximum value for numerical stability
    weights = exp(logweights - const)
    normalisedWeights[,t] = weights/sum(weights)  # Save the normalized weights
  }

  J <- sample(1:N, prob=normalisedWeights[,T], size= 1)
  #J <- which(runif(1) < cumsum(normalisedWeights[,T]))[1]

  return(particles=particles[J,])
}
