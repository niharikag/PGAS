# Implementation of blocked Particle Gibbs algorithm introduced in
#
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
#       The function returns the sample paths of (x_{1:T})

#library(foreach)
#library(doParallel)
source("R/utils.R")
source("R/auxiliaryParticleFilter.R")
require(plotly)
require(smcUtils)

blockedCPF <- function(param, y, x0, x_ref, N = 100, resamplingMethod = "multi",
                       AS = FALSE, x_last=NULL)
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

  # set resampling method
  if(resamplingMethod == 'systematic')
  {
    resampling <- systematic.resample
  }
  else if(resamplingMethod == 'stratified')
  {
    resampling <- stratified.resample
  }
  else{
    resampling <- multinomial.resample
  }

  # Initialize variables
  particles <- matrix(0, nrow = N, ncol = T)
  normalisedWeights <- matrix(0, nrow = N, ncol = T)
  B <- matrix(0, nrow = N, ncol = T) # for ancestral geneology
  logLikelihood <- 0

  # Init state, at t=0
  particles[, 1] = x0  # Deterministic initial condition
  particles[N, 1] = x_ref[1]  # Set the Nth particle to the reference particle
  normalisedWeights[, 1] = 1/N

  for(t in 2:T)
  {
    newAncestors <- resampling(normalisedWeights[, t-1])
    # weighting step
    x_pred = f(particles[,t-1],t-1)
    logweights = dnorm(y[t], mean = g(x_pred[newAncestors]), sd = sqrt(R), log = TRUE)
    const = max(logweights)
    # Subtract the maximum value for numerical stability
    w = exp(logweights - const)
    w = w/sum(w)  # Save the normalized weights


    ancestors <- resampling(w)
    newAncestors <- newAncestors[ancestors]
    newAncestors[N] = N

    # propogation step
    particles[,t] = x_pred[newAncestors] + sqrt(Q)*rnorm(N)
    particles[N,t] = x_ref[t] # set the nth particle deterministically

    if(AS){
      #Ancestor sampling
      m = dnorm(x_ref[t], mean = x_pred, sd = sqrt(Q), log = TRUE)
      const = max(m)  # Subtract the maximum value for numerical stability
      w_as = exp(m - const)
      #w_as <- normalisedWeights[,t-1]*w_as
      w_as <- w*w_as
      w_as = w_as/sum(w_as)  # Save the normalized weights
      newAncestors[N] <-  which(runif(1) < cumsum(w_as))[1]
    }

    B[, t-1] = newAncestors

    # weighting step
    logweights = dnorm(y[t], mean = g(particles[,t]), sd = sqrt(R), log = TRUE)
    max_weight = max(logweights)
    # Subtract the maximum value for numerical stability
    new_weights = exp(logweights - max_weight)/w[ancestors]
    normalisedWeights[,t] = new_weights/sum(new_weights)  # Save the normalized weights
  }

  # Reweighting step
  if(!is.null(x_last)){
    logweights = dnorm(x_last, mean = f(particles[,T],T), sd = sqrt(Q), log = TRUE)
    max_weight = max(logweights)
    weights_1 = exp(logweights - max_weight)*normalisedWeights[,T]
    normalisedWeights[,T] = weights_1/sum(weights_1)  # Save the normalized weights
    temp=f(particles[B[, T-1],T-1], T-1)
    logweights = dnorm(y[T], mean = g(temp), sd =sqrt(R), log = TRUE)
    max_weight = max(logweights)
    # Subtract the maximum value for numerical stability
    new_weights = exp(logweights - max_weight)*weights_1
    normalisedWeights[,T] = new_weights/sum(new_weights)  # Save the    normalized weights
  }

  B[,T] <- 1:N

  x_star = rep(0,T)

  J <- which(runif(1) < cumsum(normalisedWeights[,T]))[1]
  for (t in 1:T) {
    x_star[t] = particles[B[J,t],t]
  }

  return(x_star)
}

blockedPG <- function(param, y, x0=0, L = 4, P=2, N = 100, M = 1000,
                      resamplingMethod = "multi", AS=FALSE)
{
  # Stop, if input parameters are NULL
  if(is.null(param) || is.null(y) || is.null(x0))
  {
    stop("Error: the input parameters are NULL")
  }

  P = P+1 # initial state is considered to be fixed

  # Number of states
  T <- length(y)
  #Initialize the state parameters
  f <- param$f # state transition function
  g <- param$g # tranfer function
  Q <- param$Q # process noise variance
  R <- param$R # measurement noise variance
  #L Block Size
  #P overlap
  # Initialize the state by running an APF

  pf = APF(stateTransFunc, transferFunc, Q, R, x0)
  pf$generateWeightedParticles(y,N)
  x_ref = pf$sampleStateTrajectory()

  if(P >= L/2){
    stop("Error: Overlap cannot exceed block size")
  }

  internalStates = L-2*P
  startID = seq(1,T,L-P)
  numBlocks = length(startID);
  x_ref_new = rep(0, T)
  for(m in 2:M)
  {
    for(i in 1:numBlocks){
      s = startID[i]
      u = min(s+L-1,T)

      if(i == numBlocks){# last block
        x_init <- x_ref[s]
        x_last = NULL
      }else if(i > 1) { #Intermediate block
        x_init <- x_ref[s]
        x_last = x_ref[u+1]
      }else{ # First block
        x_init <- x0
        x_last = x_ref[u+1]
      }
      x_temp = blockedCPF(param, y[s:u], x_init, x_ref[s:u], N = 100, resamplingMethod = "multi",
              AS = AS, x_last=x_last)
      if(s==1){
        x_ref_new[1:u] = x_temp
      }else{
          x_ref_new[(s+1):u] = x_temp[2:length(x_temp)]
      }
    }
    x_ref <- x_ref_new
  }
  return(x_ref)
}

demo_blockedPG <- function()
{
  # Set up some parameters
  N = 100 # Number of particles
  T = 100 # Length of data record

  # define functions
  stateTransFunc = function(xt, t)  0.5*xt + 25*xt/(1+xt^2) + 8*cos(1.2*t)
  transferFunc = function(x) x^2/20

  # Generate data
  Q = 0.1  # True process noise variance
  R = 1 # True measurement noise variance
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = generateData(param = param, x0 = 0, T = T)
  x <- res$x
  y <- res$y
  x0 = 0
  x_ref = rep(0, N)
  L=50
  P = 0
  res_iPG = blockedPG(param = param, L=L, P = P, y = y, x0 = x0, N = 100, M=100,
                      AS=FALSE, resamplingMethod = 'stratified')

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res_iPG,
            name = 'CPF_AS Filtered States', type = 'scatter', mode = 'lines+markers')
}
