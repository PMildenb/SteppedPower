library(lme4)
library(Matrix)
library(microbenchmark)

Clusters <- rep(2,4)
construct.X <- function(Clusters){
  nCl <- length(Clusters)
  X   <- diag(nCl)
  X[upper.tri(X)] <- 1
  X   <- X[rep(1:nCl,Clusters),]
  return(X)
}
construct.X(rep(2,4))

triangle
construct.V <- function()

  
## power function for unequal variances in swd

## Covariance Matrix
V <- matrix(c(1,.5,0,.5,1,0,0,0,1),nrow = 3)

## Design matrix
X <- t(matrix(c(1,1,0,1,0,0,1,0,1),nrow = 3))

## covariance for beta_hat
microbenchmark({
  C <- solve(t(X) %*% solve(V) %*% X)
},{
  C2 <- chol2inv(chol(t(X) %*% chol2inv(chol(V)) %*% X))
})
  C2
C

t(X) %*%solve(V) %*% X




###########################################################################
View(swCRTdesign::swPwr)
View(lmer)
View(mkLmerDevfun)
View(lFormula)
View(optimizeLmer)
View(lme4:::optwrap)
