\name{PGAS}
\alias{PGAS}
\title{
An Implementation of Particle Gibbs with ancestor sampling
}
\description{
An Implementation of Particle Gibbs with ancestor sampling using inverse gamma priors on noise variables.
}
\usage{
PGAS(param, y, x0, prior, M = 1000, N = 100, resamplingMethod = "multi")
}
\arguments{
  \item{param}{
  state parameters
}
  \item{y}{
  measurements
}
  \item{x0}{
  initial state
}
  \item{prior}{
  Hyperparameters for the inverse gamma priors (uninformative)
}
  \item{M}{
  number of MCMC runs
}
  \item{N}{
  number of particles
}
  \item{resamplingMethod}{
  resampling methods: 'multi' multinomial and 'systematic' systematic resampling methods are currently supported.
}
}

\value{
  \item{q }{Sample path of the process noise Q}
  \item{r }{Sample path of the measurement noise R}
  \item{x }{Sample path of the states}
}
\references{
[1] Lindsten, Fredrik, Michael I. Jordan, and Thomas B. Schön. "Particle Gibbs with ancestor sampling." The Journal of Machine Learning Research 15.1 (2014): 2145-2184.
}
\author{
 Niharika Gauraha
}

\seealso{
 \code{\link{PMMH}},
 \code{\link{particleFilter}},
 \code{\link{conditionalParticleFilter}}
}
\examples{
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

  for(t in 1:T)
  {
    if(t < T)
    {
      x[t+1] = stateTransFunc(x[t],t) + sqrt(Q)*rnorm(1)
    }
    y[t] = transferFunc(x[t]) + sqrt(R)*rnorm(1)
  }
  return(list(x = x, y = y))
}

stateTransFunc = function(xt, t)  0.5*xt + 25*xt/(1+xt^2) + 8*cos(1.2*t)
  transferFunc = function(x) x^2/20

  # Set up some parameters
  N1 = 5                 # Number of particles used in PGAS
  N2 = 500               # Number of particles used in PMMH
  T = 100                # Length of data record
  numMCMC = 3000         # Number of iterations in the MCMC samplers
  burnin = 300           # Number of interations to burn

  # Generate data
  Q = 0.1  # True process noise variance
  R = 1 # True measurement noise variance
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = generateData(param = param, x0 = 0, T = T)
  x <- res$x
  y <- res$y

  # Hyperparameters for the inverse gamma priors (uninformative)
  prior = c(0.01, 0.01)

  cat("First plot true states and observed states ")
  #p <-plot_ly(x = c(1:T), y = x,
  #            name = 'Real States', type = 'scatter', mode = 'lines+markers')
  #add_lines(p, x = c(1:T), y = y,
  #          name = 'Observed States', type = 'scatter', mode = 'lines+markers')

  cat("Running PGAS : ")
  param <- list(f = stateTransFunc, g = transferFunc, Q = 1, R = 0.1)
  res = PGAS(param, y, x0 = 0, prior = prior, M = numMCMC, N = N1)
  #p <-plot_ly(x = c(1:T), y = x,
  #            name = 'Real States', type = 'scatter', mode = 'lines+markers')
  #add_lines(p, x = c(1:T), y = res$x[N1,], name = 'PGAS States',
  #          type = 'scatter', mode = 'lines+markers')
  # plot histrograma of the process noise variance and the measurement variance
  hist(res$q[burnin:numMCMC], main = "Distribution of the process noise variance", freq = FALSE)
  hist(res$r[burnin:numMCMC], main = "Distribution of the measurement noise variance", freq = FALSE)


}
