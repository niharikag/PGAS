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

  list_pf = c()
  for (i in 1:nNode_CSMC) {
    pf = APF(f, g, Q, R, x0)
    list_pf = c(list_pf,pf)
  }

  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores, type="FORK")
  registerDoParallel(cl)
  x_refs = foreach(k = 1:nNode_CSMC, .combine = "rbind", .packages = c("smcUtils")) %dopar%
  {
    list_pf[[k]]$generateWeightedParticles(y)
    list_pf[[k]]$sampleStateTrajectory()
  }
  stopCluster(cl)

  list_pf = c()
  for (i in 1:nNode_SMC) {
    pf = APF(f, g, Q, R, x0)
    list_pf = c(list_pf,pf)
  }

  list_cpf = c()
  for (i in 1:nNode_CSMC) {
    pf = CPF(f, g, Q, R, x0)
    list_cpf = c(list_cpf,pf)
  }

  # Run MCMC loop
  for(m in 2:M)
  {
    cl <- makeCluster(no_cores, type="FORK")
    registerDoParallel(cl)
    res_pf = foreach(k = 1:nNode_SMC, .combine = "rbind", .packages = c("smcUtils")) %dopar%
    {
      list_pf[[k]]$generateWeightedParticles(y)
      xRef = list_pf[[k]]$sampleStateTrajectory()
      x_ll = list_pf[[k]]$getLogLikelihood()
      list(x_ref = xRef, ll=x_ll)
    }
    stopCluster(cl)

    for (i in 1:nNode_SMC) {
      X_smc[i, ] = unlist(res_pf[i])
      ll_smc[i] <- exp(unlist(res_pf[i+nNode_SMC]))
    }

    cl <- makeCluster(no_cores, type="FORK")
    registerDoParallel(cl)
    res_cpf <- foreach(k = 1:nNode_CSMC, .combine = "rbind", .packages = c("smcUtils")) %dopar%
    {
      list_cpf[[k]]$generateWeightedParticles(y, x_refs[k,])
      xRef = list_cpf[[k]]$sampleStateTrajectory()
      x_ll = list_cpf[[k]]$getLogLikelihood()
      list(x_ref = xRef, ll=x_ll)
    }
    stopCluster(cl)

    for (i in 1:nNode_CSMC) {
      X_csmc[i, ] = unlist(res_cpf[i])
      ll_csmc[i] <- exp(unlist(res_cpf[i+nNode_CSMC]))
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
