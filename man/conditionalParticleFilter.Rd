\name{conditionalParticleFilter}
\alias{conditionalParticleFilter}
\title{
An implementation of conditional particle filter with ancestor sampling
}
\description{
An implementation of conditional particle filter with ancestor sampling
}
\usage{
conditionalParticleFilter(param, y, x0, X, N = 100, resamplingMethod = "multi")
}
%- maybe also 'usage' for other objects documented here.
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
 \item{X}{
  conditioned particles
}
  \item{N}{
  number of particles
}

\item{resamplingMethod}{
  resampling methods: 'multi' multinomial and 'systematic' systematic resampling methods are currently supported.
}
}
\value{
  \item{particles }{particles}
  \item{logLikelihood }{log Likelihood}
  \item{normalisedWeights}{ normalised Weights}
}
\author{
Niharika Gauraha
}

\seealso{
 \code{\link{PGAS}},
 \code{\link{PMMH}},
 \code{\link{particleFilter}}
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
  T = 100                # Length of data record

  # Generate data
  Q = 0.1  # True process noise variance
  R = 1 # True measurement noise variance
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = generateData(param = param, x0 = 0, T = T)
  x <- res$x
  y <- res$y

  cat("First plot true states and observed states ")
  #p <-plot_ly(x = c(1:T), y = x,
  #            name = 'Real States', type = 'scatter', mode = 'lines+markers')
  #add_lines(p, x = c(1:T), y = y,
  #          name = 'Observed States', type = 'scatter', mode = 'lines+markers')

  cat("Running conditional particle filter ")
  param <- list(f = stateTransFunc, g = transferFunc, Q = Q, R = R)
  res = conditionalParticleFilter(param = param, y = y, x0 = 0, X = x, N = 100)
  J <- which(runif(1) < cumsum(res$w[,T]))[1]
  #p <-plot_ly(x = c(1:T), y = x,
  #            name = 'Real States', type = 'scatter', mode = 'lines     +markers')
  #add_lines(p, x = c(1:T), y = res$particles[J,],
  #          name = 'CPF_AS Filtered States', type = 'scatter', mode =    'lines+markers')

}
