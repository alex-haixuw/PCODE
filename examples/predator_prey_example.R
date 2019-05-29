rm(list=ls())

library(deSolve)
library(fda)
library(MASS)
library(pracma)
source('/Users/Inori/Desktop/PCODE/R/PCODE.R')

model.par   <- c(
  birth.ch = 3.9081, birth.br = 1.9731, k.ch = 4.28, k.br = 15.6758,
  eps = 0.1119, alpha = 0.0143, m = 0.1515)
state       <- c(N = 0, C = 30, R = 0.2549, B = 1.1589)

ode.model <- function(t, state,parameters){
  with(as.list(c(state,parameters)),
       {
         dN <- 0.68 * (80 - N) - (birth.ch * N)/(k.ch + N) * C
         dC <- (birth.ch * N)/(k.ch + N) * C - (birth.br * C) / (k.br + C) * B / eps - 0.68 * C
         dR <- (birth.br * C) / (k.br + C) * R - (0.68 + m + alpha) * R
         dB <- (birth.br * C) / (k.br + C) * R - (0.68 + m ) * B
         return(list(c(dN,dC,dR,dB)))
       })
}


#Observation points
times <- seq(0,16,length.out=401)
#Generate ODE
desolve.mod <- ode(y=state,times=times,func=ode.model,parms = model.par)
#
par(mfrow=c(4,1))
par(mar = c(4,4,2,2) + 0.1)
plot(desolve.mod[,1],desolve.mod[,2],type='l',main='Predator-Prey model in ecology',ylab='N(t)',xlab='',lwd=2)
plot(desolve.mod[,1],desolve.mod[,3],type='l',main='',ylab='C(t)',xlab='t',lwd=2)
plot(desolve.mod[,1],desolve.mod[,4],type='l',main='',ylab='R(t)',xlab='t',lwd=2)
plot(desolve.mod[,1],desolve.mod[,5],type='l',main='',ylab='B(t)',xlab='t',lwd=2)


Dmodel <- function(state,parameters){
  with(as.list(c(state,parameters)),
       {
         dN <- .68 * (80 - N) - (birth.ch * N)/(k.ch + N) * C
         dC <- (birth.ch * N)/(k.ch + N) * C - (birth.br * C) / (k.br + C) * B / eps - .68 * C
         dR <- (birth.br * C) / (k.br + C) * R - (.68 + m + alpha) * R
         dB <- (birth.br * C) / (k.br + C) * R - (.68 + m ) * B
         return(list(c(dN,dC,dR,dB)))
       })
}

nobs  <- length(times)
scale <- 0.1

observ <- matrix(NA, nrow = length(times),ncol =length(state)+1)
observ[,1] <- times
observ[,2] <- desolve.mod[,2] + scale*rnorm(n = nobs, mean = 0 , sd = 1)
observ[,3] <- desolve.mod[,3] + scale*rnorm(n = nobs, mean = 0 , sd = 1)
observ[,4] <- desolve.mod[,4] + scale*rnorm(n = nobs, mean = 0 , sd = 1)
observ[,5] <- desolve.mod[,5] + scale*rnorm(n = nobs, mean = 0 , sd = 1)

#Define basis
knots <- seq(0,16,length.out=201)
#order of basis functions
norder <- 4
#number of basis funtions
nbasis <- length(knots) + norder - 2
#creating basis
basis  <- create.bspline.basis(c(0,16),nbasis,norder,breaks = knots)

basis.list <- list(basis,basis,basis,basis)
data <- observ[,c(3,5)]
colnames(data) <- c('C','B')

pcode.result <- pcode(data = data, time= observ[,1], ode.model = Dmodel, par.names = names(model.par),
                       state.names = names(state), par.initial = model.par,lambda = 1e-2,basis.list = basis.list,
                       controls = list(smooth.lambda = 10,verbal = 1,maxeval = 50))
pcode.result$structural.par
