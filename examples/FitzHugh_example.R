rm(list=ls())

library(deSolve)
library(fda)
library(MASS)
library(pracma)
source('/Users/Inori/Desktop/PCODE/R/PCODE.R')

#####

#Define ode model
model.par   <- c(a=0.2,b=0.2,c=3)
state       <- c(V=-1,R=1)

ode.model <- function(t, state,parameters){
  with(as.list(c(state,parameters)),
       {
         dV <- c*(V - (V^3)/3 + R)
         dR <- -(1/c) * (V - a + b*R)
         return(list(c(dV,dR)))
       })
}

Dmodel <- function(state,parameters){
  with(as.list(c(state,parameters)),
       {
         dV <- c*(V - (V^3)/3 + R)
         dR <- -(1/c) * (V - a + b*R)
         return(list(c(dV,dR)))
       })
}
#Observation points
times <- seq(0,20,length.out=401)
#Generate ODE
desolve.mod <- ode(y=state,times=times,func=ode.model,parms = model.par)
#
par(mfrow=c(2,1))
par(mar = c(4,4,2,2) + 0.1)
plot(desolve.mod[,1],desolve.mod[,2],type='l',main='FitzHugh-Nagumo model',ylab='V(t)',xlab='',lwd=2)
plot(desolve.mod[,1],desolve.mod[,3],type='l',main='',ylab='R(t)',xlab='t',lwd=2)


#simulations
nobs  <- length(times)
scale <- 0.1
noise <- scale*rnorm(n = nobs, mean = 0 , sd = 1)

observ <- matrix(NA, nrow = length(times),ncol =3)
observ[,1] <- times
observ[,2] <- desolve.mod[,2] + noise
observ[,3] <- desolve.mod[,3] + noise
#------------------------------------------------
plot(desolve.mod[,1],desolve.mod[,2],type='l',main='V(t)',ylab='V(t)',xlab='t')
points(desolve.mod[,1],desolve.mod[,2],col='red')
points(desolve.mod[,1],observ[,2],col='blue')
#------------------------------------------------
plot(desolve.mod[,1],desolve.mod[,3],type='l',main='R(t)',ylab='R(t)',xlab='t')
points(desolve.mod[,1],desolve.mod[,3],col='red')
points(desolve.mod[,1],observ[,3],col='blue')
#------------------------------------------------


#Define basis
knots <- seq(0,20,length.out=101)
#order of basis functions
norder <- 4
#number of basis funtions
nbasis <- length(knots) + norder - 2
#creating basis
basis  <- create.bspline.basis(c(0,20),nbasis,norder,breaks = knots)

basis.list <- list(basis,basis)
data <- observ[,2:3]
colnames(data) <- c('V','R')


#1. Running parameter cascade
pcode.result <- pcode(data = data, time= observ[,1], ode.model = Dmodel, par.names = c('a','b','c'),
                     state.names = c('V','R'), par.initial = c(0.1,0.3,4),lambda = 1e2,basis.list = basis.list,
                     controls = list(smooth.lambda = 1e-3,verbal = 1,maxeval = 20))
pcode.result$structural.par



#2.
deltavar.res <- deltavar(data = observ[,2:3], time= observ[,1], ode.model = ode.model, par.names = c('a','b','c'),
                               state.names = c('V','R'), par.initial = c(0.1,0.3,4),lambda = 1e2,basis.list = basis.list,
                          controls = list(smooth.lambda = 10,verbal = 1,maxiter = 50, inner.maxeval = 30),stepsize = 0.001, y_stepsize =0.001)
deltavar.res

#3.
bootsvar.res <- bootsvar(data = observ[,2:3], time= observ[,1], ode.model = ode.model, par.names = c('a','b','c'),
                         state.names = c('V','R'), par.initial = rnorm(3),lambda = 1e2,basis.list = basis.list,
                         controls = list(smooth.lambda = 10,verbal = 1,maxiter = 50, inner.maxeval = 30),bootsrep = 20)
bootsvar.res
#4.
lambda.check <- tunelambda(data =  observ[,2:3], time = observ[,1], ode.model = ode.model, par.initial = c(a=0.2,b=0.2,c=3),
                    par.names = c('a','b','c'), state.names = c('V','R'),basis.list = basis.list, lambda_grid = 10^(-2:3),cv_portion = .05,
                    rep = 2, kfolds = 5, controls = list(smooth.lambda = 10,verbal = 1,maxiter = 50, inner.maxeval = 30))
lambda.check
