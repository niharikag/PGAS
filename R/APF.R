######################################################################################
#
# Author : Niharika Gauraha
#          Uppsala University, Uppsala
#          Email : niharika.gauraha@farmbio.uu.se
#
#
#####################################################################################

# APF: Auxiliary particle filter
# Input:
#   param   - state parameters
#   y       - measurements
#   x0      - initial state
#   N       - number of particles
# Output:
#   x_star  - sample from target distribution

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

  for (t in 2:T) {

    # weighting step
    temp_p = f(particles[,t-1],t-1)
    logweights = dnorm(y[t], mean = g(temp_p), sd = sqrt(R), log = TRUE)
    const = max(logweights)
    # Subtract the maximum value for numerical stability
    unnorm_weights = exp(logweights - const)
    temp_wights = unnorm_weights/sum(unnorm_weights)  # Save the normalized weights

    # resampling step
    newAncestors <- multinomial.resample(temp_wights)

    # propogation step
    particles[, t-1] = particles[newAncestors, t-1]
    particles[,t] = f(particles[, t-1], t-1) + sqrt(Q)*rnorm(N)

    # weighting step
    logweights = dnorm(y[t], mean = g(particles[,t]), sd = sqrt(R), log = TRUE)
    const = max(logweights)
    # Subtract the maximum value for numerical stability
    weights = exp(logweights - const)/temp_wights[newAncestors]
    normalisedWeights[,t] = weights/sum(weights)  # Save the normalized weights

    newAncestors <- multinomial.resample(normalisedWeights[, t])
    particles[, t] = particles[newAncestors, t]
  }

  J <- sample(1:N, prob=normalisedWeights[,T], size= 1)
  x_star <- particles[J,]
  return(x_star)
}
