# Implementation of iPMCMC algorithm introduced in
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

library(foreach)
library(doParallel)

iPG <- function(param, y, x0=0, nNodes = 4, N = 100, M = 1000,
                resamplingMethod = "multi")
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
  Q <- param$Q # process noise variance
  R <- param$R # measurement noise variance

  nNode_SMC = nNodes/2 # no of nodes running SMC
  nNode_CSMC = nNodes - nNode_SMC # no of nodes running CSMC
  X_smc = matrix(0, nNode_SMC, T)
  X_csmc = matrix(0, nNode_CSMC, T)
  ll_smc = rep(0, nNode_SMC)
  ll_csmc = rep(0, nNode_CSMC)
  x_refs = matrix(0, nNode_CSMC, T)

  # Initialize the state by running an APF
  #param <- list(f = f, g = g, Q = Q, R = R)

  cl<-makeCluster(nNodes)
  clusterExport(cl, c( "APF", "CPF_AS", "param", "y", "x0"))
  registerDoParallel(cl)
  # initialize reference particles for all CSMS nodes
  res <- foreach(k = 1:nNode_CSMC, .combine = "cbind", .packages = c("smcUtils")) %dopar%
  {
    APF(param = param, y = y, x0 = x0,  N = N)
  }

  stopCluster(cl)
  x_refs[1, ] = unlist(res[1])
  index = 2
  for (i in 1:(nNode_CSMC-1)) {
    x_refs[index, ] = unlist(res[i*2+1])
    index = index+1
  }

  # Run MCMC loop
  for(m in 2:M)
  {
    cl<-makeCluster(nNodes)
    clusterExport(cl, c( "APF", "CPF_AS", "param", "y", "x0"))
    registerDoParallel(cl)

    res <- foreach(k = 1:nNode_SMC, .combine = "cbind", .packages = c("smcUtils")) %dopar%
    {
      APF(param = param, y = y, x0 = x0,  N = N)
    }

    stopCluster(cl)

    X_smc[1, ] = unlist(res[1])
    ll_smc[1] = exp(unlist(res[2]))
    index = 2
    for (i in 1:(nNode_SMC-1)) {
      X_smc[index, ] = unlist(res[i*2+1])
      ll_smc[index] <- exp(unlist(res[i*2+2]))
      index = index+1
    }

    cl<-makeCluster(nNodes)
    clusterExport(cl, c( "APF", "CPF_AS", "param", "y", "x0"))
    registerDoParallel(cl)

    res <- foreach(k = 1:nNode_CSMC, .combine = "cbind", .packages = c("smcUtils")) %dopar%
    {
      CPF_AS(param = param, y = y, x0 = x0,  x_ref = x_refs[k, ], N = N)
    }
    stopCluster(cl)

    X_csmc[1, ] = unlist(res[1])
    ll_csmc[1] = exp(unlist(res[2]))
    index = 2
    for (i in 1:(nNode_CSMC-1)) {
      X_csmc[index, ] = unlist(res[i*2+1])
      ll_csmc[index] <- exp(unlist(res[i*2+2]))
      index = index+1
    }


    # TO DO: weights
    weights = rep(0, nNode_SMC+1)
    weights[1:nNode_SMC] = ll_smc

    for (i in 1:nNode_CSMC) {
      weights[nNode_SMC+1] = ll_csmc[i]
      norm_weights  = weights/sum(weights)
      index = systematic.resample(norm_weights, num.samples=1)
      if(index > nNode_SMC){
        x_refs[i,] = X_csmc[i, ]
      }
      else{
        x_refs[i,] = X_smc[index, ]
        ll_csmc[i] = ll_smc[index]
      }
    }
  }

  norm_weights  = ll_csmc/sum(ll_csmc)
  index = systematic.resample(norm_weights, num.samples=1)
  return(x_refs[index,])
}



demo_iPG <- function()
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
  res_iPG = iPG(param = param, y = y, x0 = x0, N = N, M=10)

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res_iPG,
               name = 'CPF_AS Filtered States', type = 'scatter', mode = 'lines+markers')
}
