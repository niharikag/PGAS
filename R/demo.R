# Demo of the Particle Filter (PF), Conditional Particle Filter
# with ancenstor sampling (CPF-AS), Gibbs with Ancestor
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
# priors for PGAS and PMMH algorithms.
#
#
require(plotly)
require(smcUtils)
source("R/baseParticleFilter.R")
source("R/utils.R")
source("R/conditionalParticleFilter.R")
source("R/APF.R")
source("R/particleGibbs.R")
source("R/CPF_AS.R")
source("R/pgas.R")
source("R/SMC.R") # for bootstrap PF


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
  bpf = BPF(stateTransFunc, transferFunc, Q, R, x0)
  bpf$generateWeightedParticles(y)
  bpf$plotGeneology()
  x_multi = bpf$sampleStateTrajectory()

  bpf$generateWeightedParticles(y, resamplingMethod = "systematic")
  x_sys = bpf$sampleStateTrajectory()

  bpf$generateWeightedParticles(y, resamplingMethod = "stratified")
  x_stra = bpf$sampleStateTrajectory()

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  p<-add_lines(p, x = c(1:T), y = x_multi,
              name = 'SMC Filtered States', type = 'scatter', mode = 'lines+markers')
  p <- add_lines(p, x = c(1:T), y = x_sys,
            name = 'SMC+Systematic Filtered States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = x_stra,
            name = 'SMC+Stratified Filtered States', type = 'scatter', mode = 'lines+markers')
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
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = generateData(param = param, x0 = 0, T = T)
  x <- res$x
  y <- res$y

  x_ref = rep(0, N)
  res_cpf = iteratedCPF(param = param, y = y, x0 = 0, N = 10, x_ref = x_ref,
                plotGeneology='all', lengthGeneology = 19)
  res_cpf_multi = iteratedCPF(param = param, y = y, x0 = 0, N = 10, x_ref = x_ref)
  res_cpf_sys = iteratedCPF(param = param, y = y, x0 = 0, N = 10, x_ref = x_ref, resamplingMethod = "systematic")
  res_cpf_strat = iteratedCPF(param = param, y = y, x0 = 0, N = 10, x_ref = x_ref, resamplingMethod = "stratified")

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  p<-add_lines(p, x = c(1:T), y = res_cpf_multi,
               name = 'SMC Filtered States', type = 'scatter', mode = 'lines+markers')
  p <- add_lines(p, x = c(1:T), y = res_cpf_sys,
                 name = 'SMC+Systematic Filtered States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res_cpf_strat,
            name = 'SMC+Stratified Filtered States', type = 'scatter', mode = 'lines+markers')
}

demoCPF_AS <- function()
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
  res_cpfas = iteratedCPFAS(param = param, y = y, x0 = 0, N = N, x_ref = x_ref,
                plotGeneology='all', lengthGeneology = 19)
  res_cpfas_multi = iteratedCPFAS(param = param, y = y, x0 = 0, N = N, x_ref = x_ref)
  res_cpfas_sys = iteratedCPFAS(param = param, y = y, x0 = 0, N = N, x_ref = x_ref, resamplingMethod = "systematic")
  res_cpfas_strat = iteratedCPFAS(param = param, y = y, x0 = 0, N = N, x_ref = x_ref, resamplingMethod = "stratified")

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  p<-add_lines(p, x = c(1:T), y = res_cpfas_multi,
               name = 'CPF_AS Filtered States', type = 'scatter', mode = 'lines+markers')
  p <- add_lines(p, x = c(1:T), y = res_cpfas_sys,
                 name = 'CPF_AS+Systematic Filtered States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res_cpfas_strat,
            name = 'CPF_AS+Stratified Filtered States', type = 'scatter', mode = 'lines+markers')
}

demoAPF <- function()
{
  # Set up some parameters
  N = 100  # Number of particles
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

  #res_bpf = BPF(param = param, y = y, x0 = 0, N = 100)
  APF(param = param, y = y, x0 = 0, N = 100, resamplingMethod = "stratified", plotGeneology = 'all', lengthGeneology = 40)
  res_apf_multi = APF(param = param, y = y, x0 = 0, N = N, resamplingMethod = "multi")
  res_apf_sys = APF(param = param, y = y, x0 = 0, N = N, resamplingMethod = "systematic")
  res_apf_strat = APF(param = param, y = y, x0 = 0, N = N, resamplingMethod = "stratified")

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  p<-add_lines(p, x = c(1:T), y = res_apf_multi$x,
               name = 'APF Filtered States', type = 'scatter', mode = 'lines+markers')
  p <- add_lines(p, x = c(1:T), y = res_apf_sys$x,
                 name = 'APF+Systematic Filtered States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res_apf_strat$x,
            name = 'APF+Stratified Filtered States', type = 'scatter', mode = 'lines+markers')
}

demoPG <- function()
{
  # Set up some parameters
  N = 10               # Number of particles used in PGAS
  T = 100                # Length of data record
  numMCMC = 1000         # Number of iterations in the MCMC samplers
  burnin = 3000           # Number of interations to burn

  # define functions
  stateTransFunc = function(xt, time_t)  0.5*xt + 25*xt/(1+xt^2) + 8*cos(1.2*time_t)
  transferFunc = function(x) x^2/20

  # Generate data
  Q = .10  # True process noise variance
  R = 1 # True measurement noise variance
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = generateData(param = param, x0 = 0, T = T)
  x <- res$x
  y <- res$y
  #y = c(-0.82121,0.54702,0.92023,-0.87655,9.35938,14.26957,10.74371,
  #      2.74862,1.09014,8.70915)
  #x = c(0.00000,1.90620,4.43713,0.97950,13.68442,16.69604,13.85551,
  #        6.52509,   -1.87921,  -14.32580)
  # Hyperparameters for the inverse gamma priors (uninformative)
  prior = c(.01, 0.01)

  # Run the particle Gibbs filter
  cat("Running particle Gibbs ")
  param <- list(f = stateTransFunc, g = transferFunc, Q = .1, R = .1)
  x_ref = x #rep(0, T)
  numMCMC=1000
  burnin=300
  res_pg = PG(param = param, y = y, x0 = 0, prior = prior,
              N = N, M = numMCMC)

  p <-plot_ly(x = c(1:T), y = x, name = 'Real States', type = 'scatter',
              mode = 'lines+markers', line = list(color = "black"))
  add_lines(p, x = c(1:T), y = res_pg$x[numMCMC,],
            name = 'PG Filtered States', type = 'scatter',
            mode = 'markers', line = list(color = "red"))

  hist(res_pg$q[burnin:numMCMC], main = "Distribution of the process noise variance",
       xlab = "", ylab = "", freq = FALSE, breaks = 100)
  hist(res_pg$r[burnin:numMCMC], main = "Distribution of the measurement noise variance",
       xlab = "", ylab = "", freq = FALSE,  breaks = 100)

  cat("Running PGAS : ")
  param <- list(f = stateTransFunc, g = transferFunc, Q = .1, R = .1)
  #iPGAS <- iteratedPGAS(param, y, x0 = 0, x_ref= x_ref, N=10, iter = 1000)
  res_pgas = PGAS(param, y, x0 = 0, prior = prior, M = numMCMC, N = 10)
  p <-plot_ly(x = c(1:T), y = x, name = 'Real States', type = 'scatter',
              mode = 'lines+markers', line = list(color = "black"))
  add_lines(p, x = c(1:T), y = res_pgas$x[numMCMC,], name = 'PGAS States',
            type = 'scatter', mode = 'markers', line = list(color = "green"))
  burnin =300
  # plot histrograma of the process noise variance and the measurement variance
  hist(res_pgas$q[burnin:numMCMC], main = "Distribution of the process noise variance",
       breaks = 100, freq = FALSE)
  hist(res_pgas$r[burnin:numMCMC], main = "Distribution of the measurement noise variance", freq = FALSE)

  numMCMC=100
  param <- list(f = stateTransFunc, g = transferFunc, Q = 1, R = .1)
  res_pg = PG(param = param, y = y, x0 = 0, prior = prior, N = 5, M = numMCMC)
  res_pg$q

  numMCMC =1000
  param <- list(f = stateTransFunc, g = transferFunc, Q = .01, R = .1)
  res_pgas = PGAS(param=param, y=y, prior = prior, M = numMCMC, N = 5)
  round(res_pgas$q,2)
  round(res_pgas$r,2)
}
