# Conditional particle filter or Conditional SMC (CSMC)
# Input:
#   param - state parameters
#   y - measurements
#   x0 - initial state
#   x_ref - reference trajecory
#   N - number of particles
# Output:
#   x_star - sample from target distribution
CPF <- setRefClass(
  "CPF", contains = "BaseParticleFilter",
  fields = list(
    AS='logical'
  ),
  methods = list(
    initialize = function(stateTransFunc, transFunc, processNoise=1,
                          observationNoise=1, X_init=0, X_ref=NULL, ancestorSampling=FALSE)
    {
      "This method is called when you create an instance of the class."
      if(is.null(stateTransFunc) || is.null(transFunc) )
      {
        stop("Error: the input parameters are NULL")
      }
      AS <<- ancestorSampling
      callSuper(stateTransFunc, transFunc, processNoise,
                observationNoise,  X_init, X_ref)
    },
    generateWeightedParticles = function(y, X_ref, AS = FALSE, nParticles = 100, resamplingMethod = 'multi')
    {
      #print("generate weighted particles")
      if(is.null(y) || is.null(X_ref) )
      {
        stop("Error: the input parameters are NULL")
      }

      # set resampling method
      if(resamplingMethod == 'systematic')
      {
        .self$resampling <- systematic.resample
      }else if(resamplingMethod == 'stratified')
      {
        .self$resampling <- stratified.resample
      }else{
        .self$resampling <- multinomial.resample
      }

      # Number of states
      T <<- length(y)
      N <<- nParticles
      x_ref <<- X_ref

      # Initialize variables
      particles <<- matrix(0, nrow = N, ncol = T)
      normalisedWeights <<- matrix(0, nrow = N, ncol = T)
      B <<- matrix(0, nrow = N, ncol = T) # for ancestral geneology
      logLikelihood <<- 0

      # Init state, at t=0
      particles[, 1] <<- x0  # Deterministic initial condition
      particles[N, 1] <<- x_ref[1]  # Set the Nth particle to the reference particle
      normalisedWeights[, 1] <<- 1/N

      for(t in 2:T)
      {
        # weighting step
        newAncestors <- resampling(normalisedWeights[, t-1])
        xpred = f(particles[,t-1], t-1)
        logweights = dnorm(y[t], mean = g(xpred[newAncestors]), sd = sqrt(R), log = TRUE)
        max_weight = max(logweights)
        # Subtract the maximum value for numerical stability
        w = exp(logweights - max_weight)
        w = w/sum(w)  # Save the normalized weights

        # accumulate the log-likelihood
        logLikelihood <<- logLikelihood + max_weight +
          log(sum(w)) - log(N)

        ancestors = resampling(w)
        newAncestors = newAncestors[ancestors]
        # set the Nth ancestor
        newAncestors[N] = N

        # ancestor resampling
        if(AS){
          # ancestor sampling
          m = dnorm(x_ref[t], mean = xpred[newAncestors], sd = sqrt(Q), log = TRUE)
          const = max(m)  # Subtract the maximum value for numerical stability
          w_as = exp(m - const)
          w_as <- w*w_as
          w_as = w_as/sum(w_as)  # Save the normalized weights
          newAncestors[N] <-  which(runif(1) < cumsum(w_as))[1]
        }
        B[, t-1] <<- newAncestors

        # propogation step
        particles[,t] <<- xpred[newAncestors] + sqrt(Q)*rnorm(N)
        particles[N,t] <<- x_ref[t]

        # weighting step
        logweights = dnorm(y[t], mean = g(particles[,t]), sd = sqrt(R), log = TRUE)
        max_weight = max(logweights)
        # Subtract the maximum value for numerical stability
        new_weights = exp(logweights - max_weight)/w[ancestors]
        normalisedWeights[,t] <<- new_weights/sum(new_weights)  # Save the normalized weights
      }
      B[,T] <<- 1:N
    },
    # Given theta estimate states using PG
    iteratedCPF = function(y, X_ref, nParticles = 100, resamplingMethod = 'multi', iter = 100)
    {
      print("generate weighted particles")
      if(is.null(y) || is.null(X_ref) )
      {
        stop("Error: the input parameters are NULL")
      }
      x_ref <<- X_ref
      #x_star = x_ref
      for (i in 1:(iter-1)) {
        generateWeightedParticles(y, x_ref, nParticles = 100, resamplingMethod = 'multi')
        x_ref <<- sampleStateTrajectory()
      }
      return(x_ref)
    }
  )
)



