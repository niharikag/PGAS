# Demo of the Particle Filter (PF), Conditional Particle Filter
# with ancenstor sampling (CPF-AS), Particle Gibbs with Ancestor
# Sampling (PGAS) and the interacting Particle MCMC  (iPMCMC)
# algorithms.
#
# The script generates a batch of data y_{1:T} from the the standard
# nonlinear time series model,
#
#   x_{t+1} = 0.5*x_t + 25*x_t/(1+x_t^2) + 8*cos(1.2*t) + v_t,
#   y_t = 0.05*x_t^2 + e_t,
#
# with v_t ~ N(0,q) and e_t ~ N(0,r). The process noise and measurement
# noise variances (q, r) are treated as known parameters for PF and CPF-AS
# and they are treated as unkown parameters with inverse Gamma
# priors for PGAS.
#
#
require(plotly)
require(smcUtils)
source("R/baseParticleFilter.R")
source("R/utils.R")
source("R/conditionalParticleFilter.R")
source("R/auxiliaryParticleFilter.R")
source("R/particleGibbs.R")
source("R/bootstrapParticleFilter.R") # for bootstrap PF
source("R/GPSSM.R")

demoBPF <- function()
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
  x0 = 0
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = generateData(param = param, x0 = 0, T = T)
  x <- res$x
  y <- res$y

  x_ref = rep(0, N)
  pf = BPF(stateTransFunc, transferFunc, Q, R, x0)
  pf$generateWeightedParticles(y)
  pf$plotGeneology()
  x_multi = pf$sampleStateTrajectory()

  print(pf$getLogLikelihood())
  pf$generateWeightedParticles(y, resamplingMethod = "systematic")
  x_sys = pf$sampleStateTrajectory()

  pf$generateWeightedParticles(y, resamplingMethod = "stratified")
  x_stra = pf$sampleStateTrajectory()

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  p<-add_lines(p, x = c(1:T), y = x_multi,
              name = 'SMC Filtered States', type = 'scatter', mode = 'lines+markers')
  p <- add_lines(p, x = c(1:T), y = x_sys,
            name = 'SMC+Systematic Filtered States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = x_stra,
            name = 'SMC+Stratified Filtered States', type = 'scatter', mode = 'lines+markers')
}

demoAPF <- function()
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
  x0 = 0
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = generateData(param = param, x0 = 0, T = T)
  x <- res$x
  y <- res$y

  x_ref = rep(0, N)
  pf = APF(stateTransFunc, transferFunc, Q, R, x0)
  pf$generateWeightedParticles(y)
  pf$plotGeneology()
  x_multi = pf$sampleStateTrajectory()
  print(pf$getLogLikelihood())

  pf$generateWeightedParticles(y, resamplingMethod = "systematic")
  x_sys = pf$sampleStateTrajectory()

  pf$generateWeightedParticles(y, resamplingMethod = "stratified")
  x_stra = pf$sampleStateTrajectory()

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  p<-add_lines(p, x = c(1:T), y = x_multi,
               name = 'APF Filtered States', type = 'scatter', mode = 'lines+markers')
  p <- add_lines(p, x = c(1:T), y = x_sys,
                 name = 'APF+Systematic Filtered States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = x_stra,
            name = 'APF+Stratified Filtered States', type = 'scatter', mode = 'lines+markers')
}

demoCPF <- function()
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
  x0 = 0
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = generateData(param = param, x0 = 0, T = T)
  x <- res$x
  y <- res$y

  x_ref = rep(0, N)
  pf = CPF(stateTransFunc, transferFunc, Q, R, x0)
  pf$generateWeightedParticles(y, x_ref)
  pf$plotGeneology()
  x_multi = pf$sampleStateTrajectory()
  print(pf$getLogLikelihood())

  pf$generateWeightedParticles(y, x_ref, resamplingMethod = "systematic")
  x_sys = pf$sampleStateTrajectory()

  pf$generateWeightedParticles(y, x_ref, resamplingMethod = "stratified")
  x_stra = pf$sampleStateTrajectory()

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  p<-add_lines(p, x = c(1:T), y = x_multi,
               name = 'CPF Filtered States', type = 'scatter', mode = 'lines+markers')
  p <- add_lines(p, x = c(1:T), y = x_sys,
                 name = 'CPF+Systematic Filtered States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = x_stra,
            name = 'CPF+Stratified Filtered States', type = 'scatter', mode = 'lines+markers')

  x_pg = pf$iteratedCPF(y, x_ref, resamplingMethod = "stratified")

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = x_pg,
            name = 'iterated CPF Filtered States', type = 'scatter', mode = 'lines+markers')

}


demoPG <- function()
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
  x0 = 0
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = generateData(param = param, x0 = 0, T = T)
  x <- res$x
  y <- res$y

  x_ref = rep(0, N)
  numMCMC = 1000
  pf = ParticleGibbs(stateTransFunc, transferFunc, X_init=0)
  pf$simulate(y, x_ref, resamplingMethod = 'systematic', M = numMCMC)
  pf$plotGeneology()
  x_multi = pf$sampleStateTrajectory()
  print(pf$getLogLikelihood())

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = x_multi,
               name = 'SMC Filtered States', type = 'scatter', mode = 'lines+markers')

  burnin = floor(numMCMC * .3)
  q = pf$sampleProcessNoise()
  r = pf$sampleMeasurementNoise()
  hist(q[burnin:numMCMC], main = "Distribution of the process noise variance",
       xlab = "", ylab = "", freq = FALSE, breaks = 100)
  hist(r[burnin:numMCMC], main = "Distribution of the measurement noise variance",
       xlab = "", ylab = "", freq = FALSE,  breaks = 100)

}

demoPGAS <- function()
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
  x0 = 0
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = generateData(param = param, x0 = 0, T = T)
  x <- res$x
  y <- res$y

  x_ref = rep(0, N)
  numMCMC = 1000
  pf = ParticleGibbs(stateTransFunc, transferFunc, X_init=0, ancestorSampling=TRUE)
  pf$simulate(y, x_ref, resamplingMethod = 'systematic', M = numMCMC)
  pf$plotGeneology()
  x_multi = pf$sampleStateTrajectory()
  print(pf$getLogLikelihood())

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = x_multi,
            name = 'SMC Filtered States', type = 'scatter', mode = 'lines+markers')

  burnin = floor(numMCMC * .3)
  q = pf$sampleProcessNoise()
  r = pf$sampleMeasurementNoise()
  hist(q[burnin:numMCMC], main = "Distribution of the process noise variance",
       xlab = "", ylab = "", freq = FALSE, breaks = 100)
  hist(r[burnin:numMCMC], main = "Distribution of the measurement noise variance",
       xlab = "", ylab = "", freq = FALSE,  breaks = 100)

}


demoGPSSM <- function(){
  T = 500
  # Define the true functions f and g
  f <- function(x, t){
    return(tanh(2 * x))
  }

  g <- function(x){
    return (x)
  }

  # Set true noise covariance matrix
  Q = 0.1
  R = 0.1

  # Simulate trajectory
  param <- list(f = f, g = g, Q = Q, R = R)
  res = generateData(param = param, x0 = 1, T = T)

  x = res$x
  y = res$y
  xv = seq(-3, 3, len=100)

  plot(xv, f(xv, 0), main = 'Nonlinear function (to be learned)', type="l")

  hist(x)

  plot(x, type = "l", main = 'Data as a time series')

  # Priors for hyperparameters
  ell_prior = function(ell){
    return(dnorm(log(ell), mean = 0, sd = 0.1, log = TRUE))
  }

  # Priors for hyperparameters
  Sf_prior = function(Sf){
    return(dnorm(log(Sf), mean = 0, sd = 10, log = TRUE))
  }

  #GP covariance function prior:
  S_SE = function(w, ell){
    return(sqrt(2 * pi * ell^2) * exp(-(w^2 / 2)* (ell^2)))
  }

  S_f = function(lx, Sfx, lambda_x){
    return((1 / Sfx) * diag(1 / S_SE(sqrt(lambda_x), lx)))
  }

  gpssm = GPSSM(stateTransFunc=NULL, transFunc=g, processNoise=0, observationNoise=R,
                X_init=1, X_ref=NULL, ancestorSampling=FALSE)

  x_gpssm = gpssm$fit(y, m=16, lQ = 10, LambdaQ = 1,
                      L = 4, N=10, K_MH = 5, K=100)
  xv = seq(-3, 3, len=100)
  phi_x_xv = gpssm$phi_x(m=16, L=4, xv)
  f_m = gpssm$f %*% phi_x_xv

  # plot functions first
  p <-plot_ly(x = xv, y = f(xv,1),
              name = 'True function', type = 'scatter', mode = 'lines')
  add_lines(p, x = xv, y = f_m[1,],
            name = 'GPSSM learned function', type = 'scatter', mode = 'lines')

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = x_gpssm,
            name = 'GPSSM learned States', type = 'scatter', mode = 'lines')

  hist(gpssm$Qr)
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
