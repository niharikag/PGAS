generateData <- function(T, q, r)
{
  x = rep(0, T)
  y = rep(0, T)
  x[1] = 0  # Initial state

  for(t in 1:T)
  {
      if(t < T)
      {
          x[t+1] = stateTransFunc(x[t],t) + sqrt(q)*rnorm(1)
      }
      y[t] = transferFunc(x[t]) + sqrt(r)*rnorm(1)
  }
  return(list(x = x, y = y))
}

