rm(list=ls())

library(deSolve)
library(fda)
library(MASS)
source('par_cascade_functions.R')

theta.true <- -5
x0.initial <- 10

#Derivative model
deriv_model <- function(t, y, parms){
  theta = parms
  list(theta*y)
}

#Observation points
times <- seq(0,1,length.out=21)

#Set up ODE system
ode.system <- ode(y = x0.initial, times = times, func = deriv_model,parms = theta.true, method='ode45')

#check values
head(ode.system, n = 5)
tru.obs <- ode.system[,2]


#Plot solution
plot(ode.system) #curve
points(times,ode.system[,2],col='red') #Points


#Simulate data
nobs  <- length(times)
scale <- 0.1
noise <- scale*rnorm(n = nobs, mean = 0 , sd = 1)

observ <- tru.obs + noise
#plot simulated data against generating model
points(times, observ,pch='*',col='blue')


#--------------------------------------------------------
#Generating basis functions for interpolating observations


#knots located @ each observation time point
knots <- times
#order of basis functions
norder <- 4
#number of basis funtions
nbasis <- length(knots) + norder - 2
#creating basis
basis  <- create.bspline.basis(c(0,1),nbasis,norder,breaks = knots)

cascade_est(-4, derive_model, basis, observ, times,oned_range = c(-6,-4))
