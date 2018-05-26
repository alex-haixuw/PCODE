rm(list=ls())

library(deSolve)
library(fda)
library(MASS)
library(R.matlab)
library(pracma)
library(tictoc)
#source('par_cascade_functions.R')

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
plot(desolve.mod[,1],desolve.mod[,2],type='l',main='V(t)',ylab='V(t)',xlab='t')
plot(desolve.mod[,1],desolve.mod[,3],type='l',main='R(t)',ylab='R(t)',xlab='t')

#simulations 

par.result <- matrix(NA, nrow = 100,ncol = 3)
comp.time  <- rep(NA, nrow = 100)

for (jj in 1:1){
  #Generate observations 
  
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
  knots <- seq(0,20,length.out=401)
  #order of basis functions
  norder <- 5
  #number of basis funtions
  nbasis <- length(knots) + norder - 2
  #creating basis
  basis  <- create.bspline.basis(c(0,20),nbasis,norder,breaks = knots)
  
  basis.list <- list(basis,basis)
  
  start_time <- Sys.time()
  temp <- nls_multidim(par.initial = c(-0.1,0.1,5),ode.model = Dmodel, basis.list = basis.list,
                       observations = observ[,2:3],observation.times = observ[,1], state.names=c('V','R'), model.names = c('a','b','c'))
  #PC_ODE()
  #different observation times 
   #observations : data = 
   #observation.times : times = 
   #data, time , ode.model, state.names, par.names, par.initial = NULL, basis.list = DEFAULT (order = 5, BSPLINE,knots = times)
  end_time <- Sys.time()
  
  par.result[jj,] <- temp$structural.par
  comp.time[jj]       <- end_time - start_time
  
}

boxplot(par.result)
boxplot(comp.time)

tic()
something <- outterobj_multi_nls(ode.parameter = par.initial, basis.initial = unlist(initial_coef), derivative.model = Dmodel, inner.input = inner.input)
toc()







