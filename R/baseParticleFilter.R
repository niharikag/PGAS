#############################################
### Base SVMPlus class
#############################################

.BaseParticleFilter <- setRefClass(
  "BaseParticleFilter",
  fields = list(
    f = "ANY",
    g = "ANY",
    Q = "numeric",
    R = "numeric",
    y = 'numeric',
    x0 = "numeric",
    x_ref = "ANY",
    resampling = "ANY",
    particles = "matrix",
    normalisedWeights = "matrix",
    B = "matrix",
    T = "numeric",
    N = "numeric",
    logLikelihood = "numeric"
  ),
  methods = list(
    initialize = function(stateTransFunc, transFunc, processNoise,
                          observationNoise, X_init, X_ref=NULL)
    {
      "This method is called when you create an instance of the class."
      .self$f <<- stateTransFunc
      .self$g <<- transFunc
      .self$Q <<- processNoise
      .self$R <<- observationNoise
      .self$x0 <<- X_init
      .self$x_ref <<- X_ref
      resampling <<- NULL
      T <<- 0
      N <<- 10
      logLikelihood <<- 0
    },
    generateWeightedParticles = function(samplingMethod)
    {
      print("generate weighted particles")
    },
    plotGeneology = function(type='all', lengthGeneology=10)
    {
      if(T==0){
        stop("Error: call generateWeightedParticles first")
      }
      print("plotGeneology")
      if(type == 'all'){
        particleGeneologyAll(particles, B, lengthGeneology = lengthGeneology)
      }
      else{
        particleGeneology(particles, B, lengthGeneology = lengthGeneology)
      }

    },
    sampleStateTrajectory = function()
    {
      if(T==0){
        stop("Error: call generateWeightedParticles first")
      }

      x_star = rep(0,T)

      J <- which(runif(1) < cumsum(normalisedWeights[,T]))[1]
      for (t in 1:T) {
        x_star[t] = particles[B[J,t],t]
      }

      return(x_star)

    },
    getLogLikelihood = function()
    {
      if(T==0){
        stop("Error: call generateWeightedParticles first")
      }

      return(logLikelihood)
    },
    setParameters = function(stateTransFunc, transFunc, processNoise, observationNoise)
    {
      "set parameters"
      f <<- stateTransFunc
      g <<- transFunc
      Q <<- processNoise
      R <<- observationNoise
    }
  )
)



