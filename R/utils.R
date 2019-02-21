# plot geneology of the survived particles
particleGeneology = function(particles, B, lengthGeneology = T-1)
{
  N = dim(particles)[1]
  T = dim(particles)[2]
  x_matrix = matrix(0,N,T)
  startIndex = T-lengthGeneology

  mIndex = 0
  for (t in 1:T) {
    x_matrix[, t] = t
  }

  #plot all the particles first
  matplot(x_matrix[,startIndex:T], particles[,startIndex:T], pch=19, cex = .4, col = "black",
          xlab = "time-stamp", ylab = "particles")

  # plot geneology
  x_star = rep(0, lengthGeneology+1)
  for (j in 1:N) {
    index=1
    for (t in startIndex:T) {
      x_star[index] = particles[B[j,t],t]
      index = index+1
    }
    lines(c(startIndex:T),x_star, type = "l", col = "grey")
  }
}

# plot geneology of all the particles generated (survived and died)
particleGeneologyAll = function(particles, B, lengthGeneology = T-1)
{
  N = dim(particles)[1]
  T = dim(particles)[2]
  x_matrix = matrix(0,N,T)
  startIndex = T-lengthGeneology

  for (t in startIndex:T) {
    x_matrix[, t] = t
  }


  #plot all the particles first
  matplot(x_matrix[,startIndex:T], particles[,startIndex:T], pch=19, cex = .4, col = "black",
          xlab = "time-stamp", ylab = "particles")


  # plot geneology
  x_star = rep(0, T)
  for (i in T:(startIndex+1)) {
    index=1
    for (j in 1:N) {
      x_star[i] = particles[j, i]
      for (t in startIndex:(i-1)) {
        x_star[t] = particles[B[j,t],t]
      }
      lines(c(startIndex:i),x_star[startIndex:i], type = "l", col = "grey")
    }
  }

  # plot geneology of survived
  x_star = rep(0, lengthGeneology+1)
  for (j in 1:N) {
    index=1
    for (t in startIndex:T) {
      x_star[index] = particles[B[j,t],t]
      index = index+1
    }
    lines(c(startIndex:T),x_star, type = "l", col = "red")
  }

}
