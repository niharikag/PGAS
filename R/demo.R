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
