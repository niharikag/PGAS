# Implementation of PG and PGAS algorithm introduced in:
#
#  Andrieu, Christophe, Arnaud Doucet, and Roman Holenstein.
#  "Particle markov chain monte carlo methods." Journal of the
#  Royal Statistical Society: Series B (Statistical Methodology) 72.3 (2010): 269-342.
# &
# Lindsten, Fredrik, Michael I. Jordan, and Thomas B. Schön. "Particle Gibbs with ancestor sampling."
# The Journal of Machine Learning Research 15.1 (2014): 2145-2184.
#
# Input:
#   param - state parameters
#   y - measurements
#   x0 - initial state
#   M - number of MCMC runs
#   N - number of particles
#   resamplingMethod - resampling methods:
#     multinomical and systematics resampling methods are supported
# Output:
#       The function returns the sample paths of (q, r, x_{1:T})

ParticleGibbs <- setRefClass(
  "ParticleGibbs", contains = "CPF",
  fields = list(
    q = 'matrix',
    r = 'matrix'
  ),
  methods = list(
    initialize = function(stateTransFunc, transFunc, processNoise=0,
                          observationNoise=0, X_init=0, X_ref=NULL, ancestorSampling=FALSE)
    {
      "This method is called when you create an instance of the class."
      if(is.null(stateTransFunc) || is.null(transFunc) )
      {
        stop("Error: the input parameters are NULL")
      }

      callSuper(stateTransFunc, transFunc, processNoise,
                observationNoise,  X_init, X_ref, ancestorSampling)
    },
    simulate = function(y, X_ref, nParticles = 100, resamplingMethod = 'multi',
                        QInit=0.1, RInit=0.1, prior_a = 0.01, prior_b = 0.01, M=100)
    {
      if(Q!=0 && R!=0){
        # noise varianced are known, sample only state trajectory
        return(iteratedCPF(y, X_ref, nParticles, resamplingMethod, M))
      }
      else{
        # Number of states
        T <<- length(y)
        #QInit <<- Q_init # process noise variance
        #RInit <<- R_init # measurement noise variance
        q <<- matrix(0, M, 1)
        r <<-matrix(0, M, 1)
        X = matrix(0, M, T)
        q[1] <<- QInit
        r[1] <<- RInit
        Q <<- QInit
        R <<- RInit
        # Initialize the state by running a PF
        generateWeightedParticles(y, X_ref, nParticles, resamplingMethod)
        X[1, ] = sampleStateTrajectory()

        # Run MCMC loop
        for(k in 2:M)
        {
          # Sample the parameters (inverse gamma posteriors)
          err_q = X[k-1,2:T] - f(X[k-1,1:(T-1)], 1:(T-1))
          err_q = sum(err_q^2)
          q[k] <<- 1/rgamma(1, shape= prior_a + (T-1)/2, rate = prior_b + err_q/2)
          err_r = y - g(X[k-1,])
          err_r <- sum(err_r^2)
          r[k]  <<-  1/rgamma(1, prior_a + T/2, prior_b + err_r/2)
          # Run CPF-AS
          Q <<- q[k]
          R <<- r[k]
          # Initialize the state by running a PF
          generateWeightedParticles(y, X[k-1,])
          X[k, ] = sampleStateTrajectory()
        }
        return(X[M, ])
      }
    },
    sampleProcessNoise = function()
    {
      return(q)
    },
    sampleMeasurementNoise = function()
    {
      return(r)
    }

  )
)


