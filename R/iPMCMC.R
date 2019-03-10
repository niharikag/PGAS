# Implementation of iPMCMC algorithm introduced in:
#
# Rainforth, Tom, et al. "Interacting particle Markov chain Monte Carlo."
# International Conference on Machine Learning. 2016.
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


  #clusterExport(cl, c( "APF", "param", "y", "x0","stateTransFunc","transferFunc", "Q", "R"))
  list_pf = c()
  for (i in 1:nNode_CSMC) {
    pf = APF(stateTransFunc, transferFunc, Q, R, x0)
    list_pf = c(list_pf,pf)
  }

  #cl<-makeCluster(nNode_CSMC)
  #clusterExport(cl, c( "list_pf"))
  #registerDoParallel(cl)
  # initialize reference particles for all CSMS nodes
  #x_refs = foreach(k = 1:nNode_CSMC, .combine = "rbind", .packages = c("smcUtils")) %dopar%
  for (k in 1:nNode_CSMC)
  {
    #pf = unlist(list_pf[k])
    #pf$generateWeightedParticles(y)
    list_pf[[k]]$generateWeightedParticles(y)
    x_refs[k, ] = list_pf[[k]]$sampleStateTrajectory()
  }
  #stopCluster(cl)

  list_pf = c()
  for (i in 1:nNode_SMC) {
    pf = APF(stateTransFunc, transferFunc, Q, R, x0)
    list_pf = c(list_pf,pf)
  }

  list_cpf = c()
  for (i in 1:nNode_CSMC) {
    pf = CPF(stateTransFunc, transferFunc, Q, R, x0)
    list_cpf = c(list_cpf,pf)
  }

  # Run MCMC loop
  for(m in 2:M)
  {
    #cl<-makeCluster(nNode_SMC)
    #clusterExport(cl, c( "list_pf"))
    #registerDoParallel(cl)
    # simulate for all SMC nodes
    #res_pf = foreach(k = 1:nNode_SMC, .combine = "rbind", .packages = c("smcUtils")) %dopar%
    for (k in 1:nNode_SMC)
    {
      list_pf[[k]]$generateWeightedParticles(y)
      X_smc[k, ] = list_pf[[k]]$sampleStateTrajectory()
      ll_smc[k] = list_pf[[k]]$getLogLikelihood()
      #list(x_ref = xRef, ll=x_ll)
    }
    #stopCluster(cl)

    #for (i in 1:nNode_SMC) {
    #  X_smc[i, ] = unlist(res_pf[i])
    #  ll_smc[i] <- exp(unlist(res_pf[i+nNode_SMC]))
    #}

    #cl<-makeCluster(nNodes)
    #clusterExport(cl, c( "list_cpf"))
    #registerDoParallel(cl)

    #res_cpf <- foreach(k = 1:nNode_CSMC, .combine = "rbind", .packages = c("smcUtils")) %dopar%
    for (k in 1:nNode_CSMC)
    {
      list_cpf[[k]]$generateWeightedParticles(y, x_refs[k,])
      X_csmc[k,] = list_cpf[[k]]$sampleStateTrajectory()
      ll_csmc[k] = list_cpf[[k]]$getLogLikelihood()
      #list(x_ref = xRef, ll=x_ll)
    }
    #stopCluster(cl)

    #for (i in 1:nNode_CSMC) {
    #  X_csmc[i, ] = unlist(res_cpf[i])
    #  ll_csmc[i] <- exp(unlist(res_cpf[i+nNode_CSMC]))
    #}

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
