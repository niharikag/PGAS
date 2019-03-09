library(MCMCpack)
library(mvtnorm)
require(plotly)
require(smcUtils)
source("R/utils.R")
source("R/APF.R")
source("R/CPF_AS.R")



# Model: g and R are known, and to learn f as a GPR,
# Q and hyper-parameters of GPR, applied reduced rank GPSSM

# Priors for hyperparameters
ell_prior <- function(ell){
  return(dnorm(log(ell), mean = 0, sd = 0.1, log = TRUE))
}

# Priors for hyperparameters
Sf_prior <- function(Sf){
  return(dnorm(log(Sf), mean = 0, sd = 10, log = TRUE))
}

# eigen functions of Laplace operator
phi_x <- function(m, L, x)
{
  jv_x = c(1:m) # generate a sequence from 1 to m

  if(length(x) > 1)
    temp = matrix(0, m, length(x))
  else{
    temp = matrix(0, m, 1)
  }
  for(i in 1:m)
  {
    temp[i,] = L^(-1 / 2) * sin(pi * (jv_x[i] * (x + L)) / (2 * L))
  }
  return(temp)
}

#GP covariance function prior:
S_SE <- function(w, ell){
  return(sqrt(2 * pi * ell^2) * exp(-(w^2 / 2)* (ell^2)))
}

S_f <- function(lx, Sfx, lambda_x){
  return((1 / Sfx) * diag(1 / S_SE(sqrt(lambda_x), lx)))
}


GPSSM <- function(y, g, R, m=16, lQ = 10, LambdaQ = 1,
                  L = 4, N=10, K_MH = 5, K=100){
  # Priors for Q: lQ and LambdaQ

  jv_x = c(1:m) # generate a sequence from 1 to m
  # eigenvalues of Laplace operator
  lambda_x = (pi * jv_x / (2 * L))^2

  # initialization variables
  Qi = rep(0, K)
  Qi[1] = 0.1
  x_prim = rep(0, T)
  Ai = matrix(0, nrow = m, ncol = K)
  Ai[, 1] = 0.1 * rnorm(m)
  ell_f = rep(0,K)
  ell_f[1] = 0.1
  Sf_f = rep(0, K)
  Sf_f[1] = 20

  for (k in 1:(K-1))
  {
      f_i <- function(x, t){
          return(Ai[,k] %*% phi_x(m, L, x))
      }

      param <- list(f = f, g = g, Q = Qi[k], R = R)
      if (k == 1){
          x_prim = APF(param, y=y, x0=1, N=N)
      }
      else{
          x_prim = CPF_AS(param, y=y, x_ref=x_prim, x0=1, N=N)
      }
      print(k)
      #print('Sampling. k = ', str(k), '/', str(K))
      # Compute statistics
      xiX = x_prim[2:T]
      zX = matrix(0, m, T - 1)
      SigmaX = matrix(0, m, m)
      PsiX =matrix(0,1, m)
      for (t in 1:(T-1)){
        zX[, t] = phi_x(m, L, x_prim[t])
        SigmaX = SigmaX + zX[, t] %*% t(zX[, t])
        PsiX = PsiX + xiX[t] * t(zX[, t])
      }

      PhiX = sum(xiX * xiX)

      VX = S_f(ell_f[k], Sf_f[k], lambda_x)

      SigmaX_bar = SigmaX + VX #+ diag(.001 * ones(m))
      PsiX_bar = PsiX
      PhiX_bar = PhiX

      GammaX_star = PsiX_bar %*% ginv(SigmaX_bar)
      # print(GammaX_star)
      SigmaX_bar_inv =ginv(SigmaX_bar)
      PiX_star = PhiX_bar - sum(GammaX_star* PsiX_bar)
      #print(matmul(GammaX_star, PsiX_bar))
      #print(PiX_star)

      X = rnorm(m)
      Qi[k + 1] = riwish(S = LambdaQ + PiX_star, v=T - 1 + lQ)
      GamX = GammaX_star + Qi[k+1]* X %*% SigmaX_bar_inv
      Ai[, k + 1] = GamX

      p_S_f = dmvnorm(GamX, GammaX_star, SigmaX_bar/Qi[k+1], log=TRUE) +
      diwish(Qi[k + 1], S=LambdaQ + PiX_star, v=T - 1 + lQ)
      #print("Qi", Qi[k + 1])
      #print(p_S_f)
      ell_f[k + 1] = ell_f[k]
      Sf_f[k + 1] = Sf_f[k]

      for (k_mh in 1:K_MH){
        ell_prop = ell_f[k + 1] + 0.1 * rnorm(1)  # Random walk proposal
        Sf_prop = Sf_f[k + 1] + rnorm(1)  # Random walk proposal

        if ((ell_prop > 0) && (Sf_prop > 0)){
          V_prop = S_f(ell_prop, Sf_prop, lambda_x)

          SigmaX_bar_prop = SigmaX + V_prop
          PsiX_bar_prop = PsiX
          PhiX_bar_prop = PhiX

          GammaX_star_prop = PsiX_bar_prop %*% ginv(SigmaX_bar_prop)
          SigmaX_bar_inv_prop = ginv(SigmaX_bar_prop)
          PiX_star_prop = PhiX_bar_prop - GammaX_star_prop %*% t(PsiX_bar_prop)
          #print(GammaX_star_prop.shape)

          p_S_prop = dmvnorm(GamX, GammaX_star_prop, SigmaX_bar_inv_prop/Qi[k+1],log=TRUE) +
            diwish(Qi[k + 1], S=LambdaQ + PiX_star_prop, v =T - 1 + lQ)

          dv = runif(1)
          dl = min(1, exp(p_S_prop + ell_prior(ell_prop) + Sf_prior(Sf_prop) - p_S_f -
                           ell_prior(ell_f[k + 1]) - Sf_prior(Sf_f[k + 1])))


          if (dv < dl){
            ell_f[k + 1] = ell_prop
            Sf_f[k + 1] = Sf_prop
          }
        }
      }
  }


  # Remove burn - in
  burn_in = min(floor(K / 2), 100)
  Ar = Ai[, (burn_in + 1): K]
  Qr = Qi[ (burn_in + 1): K]

  covA = cov(Ar)
  meanA = rowMeans(Ar)

  #f_std =  sqrt(diag( matmul(matmul(phi_x_xv, covA) , phi_x(xv))))
  # Show true function
  #plot(xv, f(xv,1), main="true", ylab = c(-3,3), type="l", col="red")
  # Show learned function
  #lines(xv, f_m, col = "blue")
  # Compare the true states and the learned states
  #plot(x, main = "true states", type="l", col="red")
  #lines(x_prim, main="learned states", col ="blue")
  return(list(x=x_prim, f_m = meanA))
}


# for reproducible results, setting the seed
#random.seed(10)

# Number of measurements to simulate
T = 100



demoGPSSM <- function(){
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

  res_gpssm = GPSSM(y, g, R, m=16, L=4, N=10, K_MH=5, K=100)

  xv = seq(-3, 3, len=100)
  phi_x_xv = phi_x(m=16, L=4, xv)
  f_m = res_gpssm$f %*% phi_x_xv

  # plot functions first
  p <-plot_ly(x = xv, y = f(xv,1),
              name = 'True function', type = 'scatter', mode = 'lines')
  add_lines(p, x = xv, y = f_m[1,],
            name = 'GPSSM learned function', type = 'scatter', mode = 'lines')

  plot(xv, f(xv,1), main="true", ylab = c(-3,3), type="l", col="red")
  # Show learned function
  lines(xv, res_gpssm$f, col = "blue")

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res_gpssm$x,
            name = 'SMC Filtered States', type = 'scatter', mode = 'lines')
}

demo2GPSSM <- function(){
  # Set true noise covariance matrix
  N = 100 # Number of particles
  T = 100 # Length of data record

  # define functions
  f = function(xt, t)  0.5*xt + 25*xt/(1+xt^2) + 8*cos(1.2*t)
  g = function(x) x^2/20
  Q = 0.1  # True process noise variance
  R = 1 # True measurement noise variance
  param <- list(f = f, g = g, Q = Q, R = R)
  res = generateData(param = param, x0 = 1, T = T)

  x = res$x
  y = res$y
  xv = seq(-3, 3, len=100)

  plot(xv, f(xv, 0), main = 'Nonlinear function (to be learned)', type="l")

  hist(x)

  plot(x, type = "l", main = 'Data as a time series')

  res_gpssm = GPSSM(y, g, R, m=16, L=4, N=10, K_MH=5, K=1000)

  xv = seq(-3, 3, len=100)
  phi_x_xv = phi_x(m=16, L=4, xv)
  f_m = res_gpssm$f %*% phi_x_xv

  # plot functions first
  p <-plot_ly(x = xv, y = f(xv,1),
              name = 'True function', type = 'scatter', mode = 'lines')
  add_lines(p, x = xv, y = f_m[1,],
            name = 'GPSSM learned function', type = 'scatter', mode = 'lines')

  p <-plot_ly(x = c(1:T), y = x,
              name = 'Real States', type = 'scatter', mode = 'lines+markers')
  add_lines(p, x = c(1:T), y = res_gpssm$x,
            name = 'GPSSM learned States', type = 'scatter', mode = 'lines')
}
