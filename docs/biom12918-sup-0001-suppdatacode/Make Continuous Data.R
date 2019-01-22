#############################################################
# Simulate correlated Gaussian responses
# in a balanced cohort stepped wedge CRT

# Feb 2018

# INPUT
# n: Number of clusters
# m: Cluster size
# t: Number of periods
# delta: Effect size
# s2: Dispersion parameter or total variance (assume to be 1
#     and suppressed in the inputs)
# beta: Vector of period effects
# alpha: Vector of correlations 
#        alpha_0=within period correlation
#        alpha_1=inter-period correlation
#        alpha_2=within-individual correlation
#############################################################
contGEN<-function(n,m,t,delta,beta,alpha){
  
  require(mvtnorm)
  
  ########################################################
  # Create block exchangeable correlation matrix.
  ########################################################
  
  bxch<-function(alpha){
    rho<-alpha[1]
    nu<-alpha[2]
    lambda<-alpha[3]
    bm1<-(1-lambda-rho+nu)*diag(t*m)
    bm2<-(lambda-nu)*kronecker(matrix(1,t,t),diag(m))
    bm3<-(rho-nu)*kronecker(diag(t),matrix(1,m,m))
    bm4<-nu*matrix(1,t*m,t*m)
    return(bm1+bm2+bm3+bm4)
  }
  
  ########################################################
  # returns variance matrix of Gaussian variables with 
  # dispersion parameter s2=1 and corr matrix r[,].
  ########################################################
  
  cor2var<-function(r,s2){
    return(s2*r)
  }
  
  # Create treatment sequences
  trtSeq<-matrix(0,t-1,t)
  trtSeq[upper.tri(trtSeq)]<-1
  g<-n/(t-1) # number of clusters per step
  
  # Simulate correlated Gaussian outcomes
  s2<-1
  y<-NULL
  r<-bxch(alpha)
  v<-cor2var(r,s2)   # v <- cov matrix
  for(i in 1:(t-1)){
    u_c<-c(beta+delta*trtSeq[i,])
    u<-rep(u_c,each=m)
    y<-cbind(y,t(rmvnorm(g,u,v))) # simulate data matrix
  }
  
  # Return simulated data matrix
  return(y)
}

set.seed(062017)
n<-20
m<-20
t<-5
delta<-0.4
beta<-cumsum(c(0,0.1,0.1/2,0.1/(2^2),0.1/(2^3)))
alpha<-c(0.03,0.015,0.2)
# Generate outcome
y<-contGEN(n,m,t,delta,beta,alpha)
y<-c(y)

# marginal mean design matrix including period and treatment indicators
X<-NULL
trtSeq<-matrix(0,t-1,t)
trtSeq[upper.tri(trtSeq)]<-1
g<-n/(t-1) # number of clusters per step
for(i in 1:(t-1)){
  for(j in 1:g){
    X<-rbind(X,kronecker(cbind(diag(t),trtSeq[i,]),rep(1,m)))}
}

# colume names of design matrix 
colnames(X)<-c("period1","period2","period3","period4","period5","treatment")
cluster<-rep(1:n,each=t*m)         # create cluster id
ind<-rep(rep(1:m,t),n)             # create individual id
period<-rep(rep(1:t,each=m),n)     # create period label

simdata_cont<-data.frame(cbind(y,ind,cluster,period,X))
setwd("D:/Research/CRT Methodology/SWDSampSize/Latex/Submission/R Code")
write.csv(simdata_cont, file = "simdata_cont.csv", row.names = FALSE)

