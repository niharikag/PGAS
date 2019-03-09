# SMC (APF): Auxiliary particle filter
BPF <- setRefClass(
  "BPF", contains = "BaseParticleFilter",
  #fields = list( ),
  methods = list(
    initialize = function(stateTransFunc, transFunc, processNoise=1,
                          observationNoise=1, X_init=0, X_ref=NULL)
    {
      "This method is called when you create an instance of the class."
      if(is.null(stateTransFunc) || is.null(transFunc) )
      {
        stop("Error: the input parameters are NULL")
      }

      callSuper(stateTransFunc, transFunc, processNoise,
                observationNoise,  X_init, X_ref)
    },
    generateWeightedParticles = function(y, nParticles = 100, resamplingMethod = 'multi')
    {
      print("generate weighted particles")
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

      # Initialize variables
      particles <<- matrix(0, nrow = N, ncol = T)
      normalisedWeights <<- matrix(0, nrow = N, ncol = T)
      B <<- matrix(0, nrow = N, ncol = T) # for ancestral geneology
      logLikelihood <<- 0

      # Init state, at t=0
      particles[, 1] <<- x0  # Deterministic initial condition
      normalisedWeights[, 1] <<- 1/N
      B[, 1] <<- 1:N

      for (t in 2:T) {
        newAncestors <- resampling(normalisedWeights[,t-1])
        #particles[, t-1] = particles[newAncestors, t-1]
        B[, t-1] <<- newAncestors

        # propogation step
        particles[,t] <<- f(particles[newAncestors, t-1], t-1) + sqrt(Q)*rnorm(N)

        # weighting step
        logweights = dnorm(y[t], mean = g(particles[,t]), sd = sqrt(R), log = TRUE)
        max_weight = max(logweights)  # Subtract the maximum value for numerical stability
        new_weights = exp(logweights - max_weight)
        normalisedWeights[,t] <<- new_weights/sum(new_weights)  # Save the normalized weights

        # accumulate the log-likelihood
        logLikelihood <<- logLikelihood + max_weight +
          log(sum(new_weights)) - log(N)
      }

      B[,T] <<- 1:N

    }
  )
)

