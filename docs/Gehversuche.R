require(SteppedPower)
require(swCRTdesign)
# rechts auf Help clicken hilft.

# Waves: du meinst sequence group
# waves sind in panel studien oder GHS Erhebungsaktivitäten die in derselben Periode stattfinden.
# bei swdesigns wären das in dem Bsp. dann 6 wellen = perioden.
# ich finde waves daher nicht so glücklich.

wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=c(12,8,10,9,14))
wlsMixedPower(mu0=0,mu1=1) # geht nicht

wlsMixedPower(mu0=0, mu1=1, Cl=rep(2,5), sigma=2, tau=0.33, N=c(12,8,10,9,14)) # geht nicht
wlsMixedPower(mu0=0, mu1=1, Cl=rep(2,5), sigma=2, tau=0.33, N=c(12,8,10,9,14,12,8,10,9,14)) #geht
wlsMixedPower(mu0=0, mu1=1, Cl=rep(2,5), sigma=2, tau=0.33, N=8) #geht
wlsMixedPower(mu0=0, mu1=1, Cl=rep(2,5), sigma=2, tau=0.33, N=rep(8,10)) #geht ist dasselbe
# random slopes besser auch random intervention effect (time constant, cluster-specific)


wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10)
##
##
## ... with auto-regressive cluster effect
sp<-wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=2, tauAR=0.5, N=10,verbose=TRUE)
wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, tauAR=1, N=10)# passt

# AR means (cov-matrix of cluster effects)_ij = tau^2*rho^|i-j|, rho=tauAR
# [rho laut SAS proc mixed repeated statement doc]
# tauAR defines the correlation of cluster-period effects a_it between adjacent time points

## ... with varying cluster size
wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=c(12,8,10,9,14))
wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33,
              N=matrix(c(12,8,10,9,14,
                         11,8,10,9,13,
                         11,7,11,8,12,
                         10,7,10,8,11,
                         9,7, 9,7,11,
                         9,6, 8,7,11),5,6))
## N gives the sample size per cluster per period. It may be specified as a number representing constant sample sizes that are constant
## across periods and clusters, or a vector of cluster specific sample sizes that are constant across periods or a matrix with cluster
## and period specific sample sizes with rows representing clusters and columns representing periods

##
##
## ... with random treatment effect (with sd=0.2), which is correlated with
## the cluster effect with rho=0.25
##
##  with time constant random treatment effect b_ik, (b_i0=0) cor(a_i, b_i1)=rho
## note that random cluster effect variance is representative for controls (k=0)
## variance of cluster effects in intervention group (k=1) is tau^2+eta^2+2*tau*eta*rho
## Hence rho=-eta/tau/2 leads to equal variance of cluster under both treatment conditions.
##
## [Die Erläuterung stimmt im Moment nur für time-constant effects]
##
sp<-wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, eta=.2, rho=.25, N=10,verbose=TRUE)
h<-sp$CovarianceMatrix[1:6,1:6]
hh<- swPwr(swDsn(rep(1,5)),distn="gaussian",n=10, mu0=0, mu1=1,sigma=2, tau=0.33, eta=.2, rho=.25, gamma=0)
hh

sp<-wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, eta=.2, rho=-.2/.33/2, N=10,verbose=TRUE)
h<-sp$CovarianceMatrix[1:6,1:6]
##
##
##
## ... with missing observations (a.k.a. incomplete stepped wedge design)
#  with empty cells
wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10, incomplete=3)
wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10,
              incomplete=matrix(c(1,1,1,0,0,
                                  1,1,1,1,0,
                                  1,1,1,1,1,
                                  1,1,1,1,1,
                                  0,1,1,1,1,
                                  0,0,1,1,1),5,6))

wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10, incomplete=3) # same as above?
wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10,
              incomplete=matrix(c(1,1,1,0,0,
                                  0,1,1,1,0,
                                  1,0,1,1,1,
                                  1,1,0,1,1,
                                  0,1,1,0,1,
                                  0,0,1,1,0),5,6))  # geht.


# incomplete	integer, either a vector (only for SWD) or a matrix. A vector defines the number of periods before and
# after the switch from control to intervention that are observed.
# A matrix consists of 1's for observed clusterperiods and 0's for unobserved clusterperiods.
#
# Für das Vektormärchen bitte ein Beispiel. Das versteh ich nicht.
#
# Kommentar: im Beispiel übergibts du numeric xx<-matrix(c(1L,1L,1L,0L),2,2) wäre integer (long integer) aber das spielt keine Rolle oder?

sp<- wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10,  incomplete=c(1,1,1,0,0)) # geht nicht
sp<- wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10,  incomplete=c(2,2,2,2,2)) # geht nicht
sp<- wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5), sigma=2, tau=0.33, N=10,  incomplete=2,verbose=TRUE) # geht nicht


##
##
## ... with two levels of clustering. This arises if the patients are
## observed over the whole  study period

# Das ist raffiniert. Du meinst das irgendwie, aber für mich ein interessantes Rätsel.
# Was ist dann der Fall wo man jeden Patienten in genau einer Period genau fünfmal beobachtet?
#   Auch two levels of clustering?
#   aus der Doku:
# psi	 random individuum effect. Leads to a closed cohort setting

#   random individuum effect kenne ich als random subject effect or random subject specific intercept
#
# CovMat numeric, a positive-semidefinite matrix of dimension (#Clusters \cdot timepoints) *#Cluster* \cdot *timepoints* rows/columns. If 'CovMat' is given, 'sigma', 'tau', 'eta' and 'rho' are ignored.
#
#   geht kürzer
#
#   CovMat numeric, a positive-semidefinite matrix with (#Clusters \cdot #timepoints) rows and columns. If 'CovMat' is given, 'sigma', 'tau', 'eta' and 'rho' are ignored.




## (often referred to as closed cohort design) or if subclusters exist
## (such as wards within clinics). For
mod_aggr  <- wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5),
                           sigma=2, tau=0.33, psi=.5,
                           N=10, incomplete=3, verbose=TRUE)
mod_indiv <- wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5),
                           sigma=2, tau=0.33, psi=.5,
                           N=10, incomplete=3, verbose=TRUE, INDIV_LVL=TRUE)
mod_aggr
mod_indiv
## Compare covariance matrices of first cluster
h<-mod_aggr$CovarianceMatrix[1:6,1:6] ; mod_indiv$CovarianceMatrix[1:60,1:60]
## indexed cluster1 subject 1 t=1,..., cluster 1 subject 2 t=1,.... etc
verstehe : diagonal: tau^2+sigma^2+psi^2.
within subject off diagona: tau^2+psi^2,
within cluster between subject (=off diagonal blocks) tau^2

6 x 6 sparse Matrix of class "dgCMatrix"

[1,] 0.5339 0.1339 0.1339 0.1339 0.1339 0.1339
[2,] 0.1339 0.5339 0.1339 0.1339 0.1339 0.1339
[3,] 0.1339 0.1339 0.5339 0.1339 0.1339 0.1339
[4,] 0.1339 0.1339 0.1339 0.5339 0.1339 0.1339
[5,] 0.1339 0.1339 0.1339 0.1339    Inf 0.1339
[6,] 0.1339 0.1339 0.1339 0.1339 0.1339    Inf

Stimmt das denn??????
 # auch schön: Now. two patients sharing a cluster and a period have more similarity than two patients sharing a cluster but
#  not a period.

  mod_i <- wlsMixedPower(mu0=0, mu1=1, Cl=rep(1,5),
                             sigma=2, tau=0.33, psi=.5,gamma=.4,
                             N=10, incomplete=3, verbose=TRUE, INDIV_LVL=TRUE)
mod_i$CovarianceMatrix[1:12,1:12]

##
##
##
## longitudinal parallel design, with 5 time periods, 3 clusters in treatment
## and control arm each.
wlsMixedPower(mu0=0, mu1=1, Cl=c(3,3), sigma=2, tau=0.33, N=10,
              dsntype="parallel", timepoints=5)
##
##
##
## ... with one baseline period and four parallel periods
wlsMixedPower(mu0=0, mu1=1, Cl=c(3,3), sigma=2, tau=0.33, N=10,
              dsntype="parallel_baseline", timepoints=c(1,4))
##
##
##
## cross-over design with two timepoints before and two after the switch
wlsMixedPower(mu0=0, mu1=1, Cl=c(3,3), sigma=2, tau=0.33, N=10,
              dsntype="crossover", timepoints=c(2,2))
##
##
##
## stepped wedge design with 32 Individuals in 8 waves, binomial outcome,
## 50% incidence under control, 25% incidence under interventional treatment.
## cluster effect sd = 0.5 (ICC of 1/3 under control),
## every individual is its own cluster.
## ... with incidences defined conditional on cluster effect=0

# du meinst: 32 individuals each observed once in each periods, were the intra individual correlation is modelled as a
subject specific random intercept.
"cluster effect =0"  meint:
  logit(p_it) =mu+ a_i + b_t + theta*X_it, a_i~N(0,tau^2), mu=logit(mu0)
  .5= 1/(1+exp(-mu)), .25=1/(1+exp(-mu-theta))

# was heißt jetzt cluster effect sd???

sp<-wlsMixedPower(mu0=0.5, mu1=0.25, Cl=rep(4,8), tau=0.5, N=1,
              family="binomial",verbose=TRUE)
##
##
## ... with  marginally defined incidences
wlsMixedPower(mu0=0.5, mu1=0.25, Cl=rep(4,8), tau=0.5, N=1,
              family="binomial", marginal_mu=TRUE)


####################################################################################
####################################################################################
## DAVOS ####
library(swCRTdesign)

S1<- swPwr(swDsn(c(3,4,3)), distn="binomial", n=38, mu0=0.33, mu1=0.22, tau=0.04, eta=0, rho=0, gamma = 0);S1 #.80


Nm<-matrix(c(7, 8, 3, 5, 8, 8, 9, 5, 12, 13,
             10, 7, 8, 4, 7, 7, 7,
             6, 11, 14, 9, 0, 8, 0, 4, 2, 7, 7, 5,
             11,1, 0, 0, 0, 0, 0, 0, 3, 0, 0),ncol=4)

Nm3<-cbind(Nm[,1:3],matrix(c(3,4,5),10,3,byrow=TRUE))


delta<-qnorm(.33)-qnorm(.22)  #=.3323

curve(pnorm)

construct_timeadjust(Cl=10,3)
construct_DesMat(Cl=5,N=)


S2<- swPwr( swDsn(c(3,4,3)), distn="gaussian", n=Nm, mu0=0, mu1=0.66, icc=.07,cac=.07,sigma=1);S2
