# For matlab implementaion, see http://www.it.uu.se/katalog/freli660/software
# Demo of the Particle Filter (PF), Conditional Particle Filter
# with ancenstor sampling (CPF-AS), Gibbs with Ancestor
# Sampling (PGAS) and the Particle Marginal Metropolis-Hastings (PMMH)
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
#require(PGAS)
require(smcUtils)

generateData <- function(param, x0, T)
{
  #Initialize the state parameters
  f <- param$f # state transition function
  g <- param$g # tranfer function
  Q <- param$Q # process noise variance
  R <- param$R # measurement noise variance

  x = rep(0, T)
  y = rep(0, T)
  x[1] = x0  # Initial state
  y[1] = 0 #g(x[1]) + sqrt(R)*rnorm(1)

  for(t in 2:T)
  {
    x[t] = f(x[t-1],t-1) + sqrt(Q)*rnorm(1)
    y[t] = g(x[t]) + sqrt(R)*rnorm(1)
  }
  return(list(x = x, y = y))

}

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
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = generateData(param = param, x0 = 0, T = T)
  x <- res$x
  y <- res$y

  x_ref = rep(0, N)
  res_bpf = BPF(param = param, y = y, x0 = 0, N = 100)
  res_pg = PG(param = param, y = y, x_ref = x_ref, x0 = 0, N = 100, iter = 1000)
  res_pgas = PGAS(param = param, y = y, x_ref = x_ref, x0 = 0, N = 100, iter = 1000)

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  p<-add_lines(p, x = c(1:T), y = res_bpf,
              name = 'SMC Filtered States', type = 'scatter', mode = 'lines+markers')
  p <- add_lines(p, x = c(1:T), y = res_pg,
            name = 'PG Filtered States', type = 'scatter', mode = 'lines+markers')

  add_lines(p, x = c(1:T), y = res_pgas,
            name = 'PGAS Filtered States', type = 'scatter', mode = 'lines+markers')
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

  #x_ref = rep(0, T)
  x_ref = x
  res = APF(param = param, y = y, x0 = rep(0, T), N = 100)
  res_bpf = BPF(param = param, y = y, x0 = 0, N = 100)
  #J <- which(runif(1) < cumsum(res$normalisedWeights[,T]))[1]
  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  p<-add_lines(p, x = c(1:T), y = res_bpf,
               name = 'BPF Filtered States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res,
            name = 'APF Filtered States', type = 'scatter', mode = 'lines+markers')
}

demo <- function()
{
  # Set up some parameters
  N1 = 5                 # Number of particles used in PGAS
  N2 = 500               # Number of particles used in PMMH
  T = 100                # Length of data record
  numMCMC = 3000         # Number of iterations in the MCMC samplers
  burnin = 300           # Number of interations to burn

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

  # Hyperparameters for the inverse gamma priors (uninformative)
  prior = c(0.01, 0.01)

  #cat("First plot true states and observed states ")
  #p <-plot_ly(x = c(1:T), y = x,
  #            name = 'Real States', type = 'scatter', mode = 'lines+markers')
  #add_lines(p, x = c(1:T), y = y,
  #          name = 'Observed States', type = 'scatter', mode = 'lines+markers')

  # Run the particle filter
  cat("Running particle filter ")
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  x_ref = rep(0, T)
  res = conditionalParticleFilter(param = param, y = y, x0 = 0, x_ref = x_ref, N = 100)

  res = particleFilter(param = param, y = y, x0 = 0, N = 100)
  #p <-plot_ly(x = c(1:T), y = x,
  #            name = 'Real States', type = 'scatter', mode = 'lines+markers')
  #J <- which(runif(1) < cumsum(res$w[,T]))[1]
  #add_lines(p, x = c(1:T), y = res$particles[J,],
  #          name = 'Filtered States', type = 'scatter', mode = 'lines+markers')

  cat("Running conditional particle filter ")
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = conditionalParticleFilter(param = param, y = y, x0 = 0, X = x, N = 100)
  #J <- which(runif(1) < cumsum(res$w[,T]))[1]
  #p <-plot_ly(x = c(1:T), y = x,
  #            name = 'Real States', type = 'scatter', mode = 'lines+markers')
  #add_lines(p, x = c(1:T), y = res$particles[J,],
  #            name = 'CPF_AS Filtered States', type = 'scatter', mode = 'lines+markers')

  cat("Running PGAS : ")
  param <- list(f = stateTransFunc, g = transferFunc, Q = 1, R = 0.1)
  res = PGAS(param, y, x0 = 0, prior = prior, M = numMCMC, N = N1)
  #p <-plot_ly(x = c(1:T), y = x,
  #            name = 'Real States', type = 'scatter', mode = 'lines+markers')
  #add_lines(p, x = c(1:T), y = res$x[N1,], name = 'PGAS States',
  #          type = 'scatter', mode = 'lines+markers')
  # plot histrograma of the process noise variance and the measurement variance
  #hist(res$q[burnin:numMCMC], main = "Distribution of the process noise variance", freq = FALSE)
  #hist(res$r[burnin:numMCMC], main = "Distribution of the measurement noise variance", freq = FALSE)

  cat("Running PMMH : ")
  # Proposal for PMMH (Gaussian random walk)
  prop = c(.1, .1)
  param <- list(f = stateTransFunc, g = transferFunc, Q = .1, R = 1)
  res = PMMH(param, y, x0 = 0, prior, prop, N = N2, M = numMCMC)
  #p <-plot_ly(x = c(1:T), y = x,
  #            name = 'Real States', type = 'scatter', mode = 'lines+markers')
  #add_lines(p, x = c(1:T), y = res$x[N2,],
  #          name = 'PMMH States', type = 'scatter', mode = 'lines+markers')
  #plot histrograms of the process noise variance and the measurement variance
  #hist(res$q[burnin:numMCMC], main = "Distribution of the process noise variance", freq = FALSE)
  #hist(res$r[burnin:numMCMC], main = "Distribution of the measurement noise variance", freq = FALSE)
}
