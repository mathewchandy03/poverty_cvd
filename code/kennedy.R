### INPUT: l is an n*p matrix, a and y are vectors of length n
### l = matrix of covariates
### a = vector of treatment values
### y = vector of observed outcomes
# set up evaluation points & matrices for predictions
a.min <- min(a); a.max <- max(a)
a.vals <- seq(a.min,a.max,length.out=100)
la.new <- rbind(cbind(l,a), cbind( l[rep(1:n,length(a.vals)),],
                                   a=rep(a.vals,rep(n,length(a.vals))) ))
l.new <- la.new[,-dim(la.new)[2]]
# fit super learner (other methods could be used here instead)
sl.lib <- c("SL.earth","SL.gam","SL.gbm","SL.glm","SL.glmnet")
pimod <- SuperLearner(Y=a, X=l, SL.library=sl.lib, newX=l.new)
pimod.vals <- pimod$SL.predict; sq.res <- (a-pimod.vals)ˆ2
pi2mod <- SuperLearner(Y=sq.res,X=l, SL.library=sl.lib, newX=l.new)
pi2mod.vals <- pi2mod$SL.predict
mumod <- SuperLearner(Y=y, X=cbind(l,a), SL.library=sl.lib,
                      newX=la.new,family=binomial); muhat.vals <- mumod$SL.predict
# construct estimated pi/varpi and mu/m values
approx.fn <- function(x,y,z){ predict(smooth.spline(x,y),x=x2)$y }
a.std <- (la.new$a-pimod.vals)/sqrt(pi2mod.vals)
pihat.vals <- approx.fn(density(a.std[1:n])$x, density(a.std[1:n])$y,
                        a.std); pihat <- pihat.vals[1:n]
pihat.mat <- matrix(pihat.vals[-(1:n)], nrow=n,ncol=length(a.vals))
varpihat <- approx.fn(a.vals, apply(pihat.mat,2,mean), a)
varpihat.mat <- matrix(rep(apply(pihat.mat,2,mean),n), byrow=T,nrow=n)
muhat <- muhat.vals[1:n]
muhat.mat <- matrix(muhat.vals[-(1:n)], nrow=n,ncol=length(a.vals))
mhat <- approx.fn(a.vals, apply(muhat.mat,2,mean), a)
mhat.mat <- matrix( rep(apply(muhat.mat,2,mean),n), byrow=T,nrow=n)
# form adjusted/pseudo outcome xi
pseudo.out <- (y-muhat)/(pihat/varpihat) + mhat
# leave-one-out cross-validation to select bandwidth
library(KernSmooth); kern <- function(x){ dnorm(x) }
w.fn <- function(bw){ w.avals <- NULL; for (a.val in a.vals){
  a.std <- (a-a.val)/bw; kern.std <- kern(a.std)/bw
  w.avals <- c(w.avals, mean(a.stdˆ2*kern.std)*(kern(0)/bw) /
                 (mean(kern.std)*mean(a.stdˆ2*kern.std)-mean(a.std*kern.std)ˆ2))
}; return(w.avals/n) }
hatvals <- function(bw){ approx(a.vals,w.fn(bw),xout=a)$y }
cts.eff <- function(out,bw){ approx(locpoly(a,out,bw),xout=a)$y }