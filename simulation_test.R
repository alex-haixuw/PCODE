rm(list=ls())

library(deSolve)
library(fda)
library(MASS)
library(R.matlab)
library(pracma)
library(tictoc)
source('par_cascade_functions.R')

#####
#Define ode model 
model.par   <- c(theta = c(0.1))
state       <- c(X     = 0.1)

ode.model <- function(t, state,parameters){
  with(as.list(c(state,parameters)),
       {
         dX <- theta*X*(1-X/10)
         return(list(dX))
       })
}

Dmodel <- function(state,parameters){
  with(as.list(c(state,parameters)),
       {
         dX <- theta*X*(1-X/10)
         return(list(dX))
       })
}

#Observation points
times <- seq(0,100,length.out=101)

desolve.mod <- ode(y=state,times=times,func=ode.model,parms = model.par)


#Plot solution
plot(desolve.mod,main='noisy data') #curve
points(times,desolve.mod[,2],col='red') #Points

theta.est.optim <- cbind(rep(NA,100),rep(NA,100))

for (jj in 1:100){
  print(jj)
#Simulate data
nobs  <- length(times)
scale <- 0.5
noise <- scale*rnorm(n = nobs, mean = 0 , sd = 1)

observ <- desolve.mod[,2] + noise
#plot simulated data against generating model
points(times, observ,pch='*',col='blue')


#--------------------------------------------------------
#Generating basis functions for interpolating observations


#knots located @ each observation time point
knots <- seq(0,100,length.out=21)
#order of basis functions
norder <- 5
#number of basis funtions
nbasis <- length(knots) + norder - 2
#creating basis
basis  <- create.bspline.basis(c(0,100),nbasis,norder,breaks = knots)


nls_temp = cascade_nls(0.3, ode.model = Dmodel, basis , observ, times,controls = list(tau = 0.1, tolx = 1e-8,tolg = 1e-8, maxeval = 500))

phi.mat <- eval.basis(times,basis)
temp <- phi.mat %*% nls_temp$nuisance.par 
points(times,temp,pch='^',col='green')




optim_temp = cascade_optim(0.1, ode.model = Dmodel, basis, observ, times,oned_range = c(0,1))




theta.est.optim[jj,1] <- nls_temp$structural.par
theta.est.optim[jj,2] <- optim_temp$structural.par

}

par(mfrow=c(1,2))
boxplot(theta.est.optim[,1],main='NLS')
abline(h=0.1,col='red')
boxplot(theta.est.optim[-39,2],main='OPTIM')
abline(h=0.1,col='red')
mtext("SD = 0.1", side = 3, line = -26, outer = TRUE)



#save(theta.est.optim, file = 'theta.est.optim_3.Rdata')
rm(list=ls())
theta.mat <- readMat('theta_est_nonlin_3.mat')
load('theta.est.optim_3.Rdata')
rm(list=ls())
theta.mat <- readMat('theta_est_lsqnonlin.mat')
load('theta.est.optim.Rdata')

theta_est_nonlin <- as.vector(theta.mat$theta.est)

type = c(rep('optim',1000),rep('lsqnonlin',1000))


temp <- c(theta.est.optim, -theta_est_nonlin)
data.mat = data.frame( theta.est = temp, est.type = type)
boxplot(data.mat$theta.est ~ data.mat$est.type)

