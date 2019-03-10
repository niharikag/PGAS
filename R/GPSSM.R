# Implementation of iPMCMC algorithm introduced in:
#
# Svensson, Andreas, et al. "Computationally efficient Bayesian learning of
# Gaussian process state space models." Artificial Intelligence and Statistics. 2016.
#

library(MCMCpack)
library(mvtnorm)
require(plotly)
require(smcUtils)

# Model: Assuming g and R are known, and to learn f as a GPR,
# Q and hyper-parameters of GPR, using reduced rank GP and PG algorithms

GPSSM <- setRefClass(
  "GPSSM", contains = "CPF",
  fields = list(
    Qr = 'numeric',
    Ar = 'matrix',
    ell_f = 'numeric',
    Sf_f = 'numeric',
    ell_prior = 'ANY',
    Sf_prior = 'ANY',
    S_SE = 'ANY',
    S_f = 'ANY'
  ),
  methods = list(
    initialize = function(stateTransFunc, transFunc, processNoise=0, observationNoise=0,
                          X_init=0, X_ref=NULL, ancestorSampling=FALSE,
                          ell_p=NULL, sf_p=NULL, s_se=NULL, s_f=NULL)
    {
      "This method is called when you create an instance of the class."

      ell_prior <<- ell_p
      Sf_prior <<- sf_p
      S_SE <<- s_se
      S_f <<- s_f

      if(is.null(ell_p))
      {
        # Priors for hyperparameters
        ell_prior <<- function(ell){
          return(dnorm(log(ell), mean = 0, sd = 0.1, log = TRUE))
        }
      }

      if(is.null(sf_p) )
      {
        # Priors for hyperparameters
        Sf_prior <<- function(Sf){
          return(dnorm(log(Sf), mean = 0, sd = 10, log = TRUE))
        }
      }

      if(is.null(s_se) )
      {
        #GP covariance function prior:
        S_SE <<- function(w, ell){
          return(sqrt(2 * pi * ell^2) * exp(-(w^2 / 2)* (ell^2)))
        }
      }

      if(is.null(s_f) )
      {
        #GP hyper-param prior:
        S_f <<- function(lx, Sfx, lambda_x){
          return((1 / Sfx) * diag(1 / S_SE(sqrt(lambda_x), lx)))
        }
      }

      if(is.null(stateTransFunc) )
      {
        stateTransFunc=1 # dummy function
      }

      callSuper(stateTransFunc, transFunc, processNoise,
                observationNoise,  X_init, X_ref, ancestorSampling)
    },
    # eigen functions of Laplace operator
    phi_x = function(m, L, x)
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
    },

    # fit GPSSM
    fit = function(y, m=16, lQ = 10, LambdaQ = 1,
                    L = 4, N=10, K_MH = 5, K=100){
      T <<- length(y)
      # local variables
      Qi = rep(0, K)
      Ai = matrix(0, nrow = m, ncol = K)
      ell_f <<- rep(0,K)
      Sf_f <<- rep(0, K)

      # initialization variables
      Qi[1] = 0.1
      x_prim = rep(0, T)
      Ai[, 1] = 0.1 * rnorm(m)
      ell_f[1] <<- 0.1
      Sf_f[1] <<- 20

      jv_x = c(1:m) # generate a sequence from 1 to m
      # eigenvalues of Laplace operator
      lambda_x = (pi * jv_x / (2 * L))^2

      for (k in 1:(K-1))
      {

        f_i <- function(x, t){
          return(Ai[,k] %*% phi_x(m, L, x))
        }

        if (k == 1){
          pf = APF(f_i, g, Qi[k], R, x0)
          pf$generateWeightedParticles(y, N)
          x_prim = pf$sampleStateTrajectory()
        }
        else{
          setParameters(f_i, g, Qi[k], R)
          generateWeightedParticles(y, x_prim)
          x_prim = sampleStateTrajectory()
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
        ell_f[k + 1] <<- ell_f[k]
        Sf_f[k + 1] <<- Sf_f[k]

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
              ell_f[k + 1] <<- ell_prop
              Sf_f[k + 1] <<- Sf_prop
            }
          }
        }
      }

      # Remove burn - in
      burn_in = min(floor(K / 2), 100)
      Ar <<- Ai[, (burn_in + 1): K]
      Qr <<- Qi[ (burn_in + 1): K]

      covA = cov(Ar)
      meanA = rowMeans(Ar)

      f <<- meanA
      return(x_prim)
    },
    sampleProcessNoise = function()
    {
      return(Qr)
    }
  )
)
