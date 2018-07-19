# Demo of the Particle Gibbs with Ancestor
# Sampling (PGAS) and the Particle Marginal Metropolis-Hastings (PMMH)
# algorithms, presented in,
#
#   [1] F. Lindsten and M. I. Jordan T. B. Sch?n, "Ancestor sampling for
#   Particle Gibbs", Proceedings of the 2012 Conference on Neural
#   Information Processing Systems (NIPS), Lake Taho, USA, 2012.
#
# and
#
#   [2] C. Andrieu, A. Doucet and R. Holenstein, "Particle Markov chain Monte
#   Carlo methods" Journal of the Royal Statistical Society: Series B,
#   2010, 72, 269-342.
#
# respectively.
#
# The script generates a batch of data y_{1:T} from the the standard
# nonlinear time series model,
#
#   x_{t+1} = 0.5*x_t + 25*x_t/(1+x_t^2) + 8*cos(1.2*t) + v_t,
#   y_t = 0.05*x_t^2 + e_t,
#
# with v_t ~ N(0,q) and e_t ~ N(0,r). The process noise and measurement
# noise variances (q,r) are treated as unkown parameters with inverse Gamma
# priors. The PGAS and PMMH algorithms are then executed independently to
# find the posterior parameter distribution p(q, r | y_{1:T}).
#
#
library(plotly)

demo <- function()
{
  # Set up some parameters
  N1 = 5                 # Number of particles used in PGAS
  N2 = 500               # Number of particles used in PMMH
  T = 100                # Length of data record
  numMCMC = 3000         # Number of iterations in the MCMC samplers
  burnin = 300           # Number of interations to burn

  # Generate data
  q0 = 0.1  # True process noise variance
  r0 = 1 # True measurement noise variance
  res = generateData(T, q0, r0)
  x0 <- res$x
  y0 <- res$y

  # Hyperparameters for the inverse gamma priors (uninformative)
  prior.a = 0.01
  prior.b = 0.01

  # Parameter proposal for PMMH (Gaussian random walk)
  prop.sigma_q = 1
  prop.sigma_r = 1

  # Initialization for the parameters
  qinit = 1
  rinit = 0.1


  # Run the particle filter
  cat("Running particle filter ")
  param <- c(1, 0.1)
  res = particleFilter(param = param, y = y0, N = 100)
  plot_ly(x = c(1:T), y = y0,
          name = 'Simulated data', type = 'scatter', mode = 'lines+markers')
  p <-plot_ly(x = c(1:T), y = x0,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res$xHatFiltered,
            name = 'Filtered States', type = 'scatter', mode = 'lines+markers')

  cat("Running conditional particle filter ")
  param <- c(1, 0.1)
  res = conditionalParticleFilter(param = param, y = y0, N = 100, x0)
  p <-plot_ly(x = c(1:T), y = x0,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res$xHatFiltered,
            name = 'Filtered States', type = 'scatter', mode = 'lines+markers')

  # Run the algorithms
  cat("Running PGAS (N=#i). Progress: ")
  res = PGAS(numMCMC, y0, prior, N1, qinit, rinit, q0, r0)
  p <-plot_ly(x = c(1:T), y = x0,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res$xHatFiltered,
            name = 'Filtered States', type = 'scatter', mode = 'lines+markers')


  cat("Running PMMH (N=#i). Progress: ")
  res = PMMH(numMCMC, y0, prior, prop, N2, qinit, rinit, q0, r0)
  p <-plot_ly(x = c(1:T), y = x0,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res$xHatFiltered,
            name = 'Filtered States', type = 'scatter', mode = 'lines+markers')



}
