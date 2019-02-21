# Conditional particle filter or Conditional SMC (CSMC)
# Input:
#   param - state parameters
#   y - measurements
#   x0 - initial state
#   x_ref - reference trajecory
#   N - number of particles
# Output:
#   x_star - sample from target distribution

# Note: CPF with the simplest possible proposal density
CPF <- function(param, y, x0, x_ref, N = 100, resamplingMethod = "multi",
                plotGeneology = NULL, lengthGeneology = 10-1)
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
  B <- matrix(0, nrow = N, ncol = T) # for ancestral geneology

  # Init state, at t=0
  particles[, 1] = x0  # Deterministic initial condition
  particles[N, 1] = x_ref[1]  # Set the Nth particle to the reference particle

  for(t in 2:T)
  {
    # weighting step
    temp_p = f(particles[,t-1],t-1)
    logweights = dnorm(y[t], mean = g(temp_p), sd = sqrt(R), log = TRUE)
    const = max(logweights)
    # Subtract the maximum value for numerical stability
    w = exp(logweights - const)
    w = w/sum(w)  # Save the normalized weights

    # resampling step
    if(resamplingMethod == 'systematic')
    {
      newAncestors <- systematic.resample(w)
    }
    else if(resamplingMethod == 'stratified')
    {
      newAncestors <- stratified.resample(w)
    }
    else{
      newAncestors <- multinomial.resample(w)
    }

    # set the Nth ancestor
    newAncestors[N] = N
    B[, t-1] = newAncestors

    # propogation step
    #particles[, t-1] = particles[newAncestors, t-1]
    particles[,t] = f(particles[newAncestors, t-1], t-1) + sqrt(Q)*rnorm(N)
    #particles[N,t] = x_ref[t]

    # weighting step
    logweights = dnorm(y[t], mean = g(particles[,t]), sd = sqrt(R), log = TRUE)
    const = max(logweights)
    # Subtract the maximum value for numerical stability
    new_weights = exp(logweights - const)/w[newAncestors]
    normalisedWeights[,t] = new_weights/sum(new_weights)  # Save the normalized weights

    # resampling step
    if(resamplingMethod == 'systematic')
    {
      newAncestors <- systematic.resample(normalisedWeights[, t])
    }
    else if(resamplingMethod == 'stratified')
    {
      newAncestors <- stratified.resample(normalisedWeights[, t])
    }
    else{
      newAncestors <- multinomial.resample(normalisedWeights[, t])
    }

    particles[, t] = particles[newAncestors, t]
    particles[N,t] = x_ref[t]
  }

  B[,T] <- 1:N

  if(!is.null(plotGeneology)){
    if(plotGeneology == 'all'){
      particleGeneologyAll(particles, B, lengthGeneology = lengthGeneology)
    }
    else{
      particleGeneology(particles, B, lengthGeneology = lengthGeneology)
    }
  }

  x_star = rep(0,T)
  #J <- sample(1:N, prob=normalisedWeights[,T], size= 1)
  J <- which(runif(1) < cumsum(normalisedWeights[,T]))[1]
  for (t in 1:T) {
    x_star[t] = particles[B[J,t],t]
  }
  return(x_star)
}

# Given theta estimate states using PG
iteratedCPF <- function(param, y, x0, x_ref, N = 100, resamplingMethod = "multi",
                        plotGeneology = NULL, lengthGeneology = 10-1, iter = 100)
{
  #x_star = x_ref
  for (i in 1:(iter-1)) {
    x_star = CPF(param = param, y = y, x_ref = x_ref, x0 = 0, N = 100)
    x_ref = x_star
  }
  x_star = CPF(param = param, y = y, x_ref = x_ref, x0 = 0, N = 100,
               plotGeneology = plotGeneology, lengthGeneology = lengthGeneology)
  return(x_star)
}
