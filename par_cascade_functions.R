#########################################################################
#
#
#
#
#
#
#
#
#########################################################################


# 1. 'innerobj'
#
# Description:
#
#   'innerobj' is the inner objective function to be optimzied for estimating nuisance
#    parameters, the basis coefficients, of the model when a new structural parameter
#    estimate is available.
#
#    This objective function contains two component:
#        1. 'SSE' corresponds to the sum of squared errors between the observations
#                and the estiamted trajectory by basis expansion.
#        2. 'ode_penalty' controls the size of the extend how the estimated trajectory
#                fails to satisfy the differential equation exactly given a smoothing
#                parameter lambda.
#
# Input:
#      basis_coef      :
#      ode.par         :
#      input           :
#      derive.model    :
#
# Output:
#
#     SSE + ode_penalty :
#     residual.vec      :
#
# Comment:
#
#
#
#

innerobj <- function(basis_coef, ode.par, input, derive.model,NLS=TRUE){
  yobs    <- input[[1]]
  Phi.mat <- input[[2]]
  lambda  <- input[[3]]
  Qmat    <- input[[4]]
  Q.D1mat <- input[[5]]
  quadts  <- input[[6]]
  quadwts <- input[[7]]

  xhat     <- Phi.mat %*% basis_coef
  #e_{ij}
  residual <- yobs - xhat
  #Sum of squared errors
  #SSE <- sum((yobs-xhat)^2)

  #PEN_{i}
  penalty_residual <- rep(NA, length(quadwts))
  #Use composite Simpson's rule to calculate the integral of residual
  Xt   <- Qmat    %*%  basis_coef
  dXdt <- Q.D1mat %*%  basis_coef

  #Evaluate f(X)
  rightside <- rep(NA, length(quadwts))

  tempfun <- function(x, temp = ode.par){
    return(derive.model(state = c(X=x), parameters= c(theta = temp)))
  }
  rightside <- unlist(lapply(Xt,tempfun))

  #for (jj in 1:length(quadwts)){
  #      rightside[jj] <- derive.model(state = c(X = Xt[jj]),parameters = c(theta = #ode.par))
  #}

  penalty_residual <- sqrt(lambda) * sqrt(quadwts) * (dXdt - rightside)

  #ode_penalty <- sum(lambda * quadwts * (penalty_residual)^2)

  residual.vec <- c(residual,penalty_residual)
  if (NLS){
    return(residual.vec)
  }else{
    return(sum(residual.vec^2))
  }
}


outterobj <- function(ode.parameter, basis.initial, derivative.model, inner.input,NLS){


  if (NLS){
    inner_coef  <- lsqnonlin(innerobj, basis.initial, options=list(maxeval = 50),ode.par = ode.parameter,derive.model = derivative.model, input = inner.input)$x
  }else{
    inner_coef <- optim(par=basis.initial, fn=innerobj,ode.par = ode.parameter, derive.model = derivative.model, input = inner.input,NLS = FALSE,
                        control = list(abstol=1e-10))$par
  }



  yobs    <- inner.input[[1]]
  Phi.mat <- inner.input[[2]]


  xhat     <- Phi.mat %*% inner_coef
  #e_{ij}
  residual <- yobs - xhat
  if (NLS){
    return(residual)
  }else{
    return(sum(residual^2))
  }

}


cascade_optim <- function(par.initial, ode.model, basis, observations, obs.times,ode_penalty = 1e6,nquadpts = NULL, controls = NULL, oned_range = NULL){

    #Input:
    #    par.initial : initial values for structural parameter
    #    ode.model   : ODE model
    #    basis       : basis object for approximating x(t)
    #    observations: observations including a vector indicating times of observations
    #    nquadpts    : number of quadrature points to use composite Simpson's rule for calculating integral


    #Output:



    #Comment:

    #number of parameters
    npar  <- length(par.initial)
    #time points of observations


    #If no 'quadpts' is given, then set it to 101
    if (is.null(nquadpts)){
        nquadpts = 101
    }



    #Evaluating basis functions at time points of observations
    #and stored as columns in Phi matrix
    Phi.mat <- eval.basis(obs.times, basis)
    #Evaluating 1st derivative of basis functions
    D1.mat <- eval.basis(obs.times, basis,1)
    #Evaluating 2nd derivative of basis functions
    D2.mat <- eval.basis(obs.times, basis,2)

    #Calculate L2 penalty
    quadts   <- seq(min(obs.times),max(obs.times),length.out=nquadpts)
    nquad    <- length(quadts)
    quadwts  <- rep(1,nquad)
    even.ind <- seq(2,(nquad-1),by=2)
    odd.ind  <- seq(3,(nquad-2),by=2)
    quadwts[even.ind] = 4
    quadwts[odd.ind] = 2
    h  <- quadts[2] - quadts[1]
    quadwts <- quadwts*(h/3)

    Qmat    <- eval.basis(quadts, basis)
    Q.D1mat <- eval.basis(quadts, basis, 1)
    Q.D2mat <- eval.basis(quadts, basis, 2)

    #Initial estimat of basis coefficients
    Rmat = t(Q.D2mat)%*%(Q.D2mat*(quadwts%*%t(rep(1,nbasis))))
    basismat2 = t(Phi.mat)%*%Phi.mat;
    smooth_lambda = 1e-6   # smoothing parameter
    Bmat    = basismat2 + smooth_lambda*Rmat;
    #Initial basis coefficients
    initial_coef = ginv(Bmat)%*%t(Phi.mat)%*%observations;

    #Passing to inner objective functions and obtain initial value for parameter cascading
    inner.input = list(observations, Phi.mat, ode_penalty, Qmat, Q.D1mat,quadts,quadwts)

    #Using optim for optimization
    new.ini.basiscoef <- optim(par=initial_coef, fn=innerobj,ode.par = par.initial,input = inner.input, derive.model = ode.model,NLS=FALSE,             control = list(abstol=1e-10))$par
      if (npar == 1){
         loweb  = oned_range[1]
         upperb = oned_range[2]
         theta.final = optim(par = par.initial, fn = outterobj, basis.initial = new.ini.basiscoef, derivative.model = ode.model , inner.input = inner.input,NLS=FALSE,method = ('Brent'),lower = loweb, upper = upperb)$par
      }else{
        theta.final = optim(par = par.initial, fn = outterobj, basis.initial = new.ini.basiscoef, derivative.model = ode.model , inner.input = inner.input,NLS =FALSE)$par
      }
    basiscoef <- optim(par=initial_coef, fn=innerobj,ode.par = theta.final,input = inner.input, derive.model = ode.model,NLS=FALSE,
                         control = list(abstol=1e-10))$par


    return(list(structural.par = theta.final, nuisance.par = basiscoef))


}

cascade_nls <- function(par.initial, ode.model, basis, observations, obs.times,ode_penalty = 1e6,nquadpts = NULL,controls = NULL){

    #Input:
    #    par.initial : initial values for structural parameter
    #    ode.model   : ODE model
    #    basis       : basis object for approximating x(t)
    #    observations: observations including a vector indicating times of observations
    #    nquadpts    : number of quadrature points to use composite Simpson's rule for calculating integral


    #Output:



    #Comment:

    #number of parameters
    npar  <- length(par.initial)
    #time points of observations


    #If no 'quadpts' is given, then set it to 101
    if (is.null(nquadpts)){
        nquadpts = 101
    }



    #Evaluating basis functions at time points of observations
    #and stored as columns in Phi matrix
    Phi.mat <- eval.basis(obs.times, basis)
    #Evaluating 1st derivative of basis functions
    D1.mat <- eval.basis(obs.times, basis,1)
    #Evaluating 2nd derivative of basis functions
    D2.mat <- eval.basis(obs.times, basis,2)

    #Calculate L2 penalty
    quadts   <- seq(min(obs.times),max(obs.times),length.out=nquadpts)
    nquad    <- length(quadts)
    quadwts  <- rep(1,nquad)
    even.ind <- seq(2,(nquad-1),by=2)
    odd.ind  <- seq(3,(nquad-2),by=2)
    quadwts[even.ind] = 4
    quadwts[odd.ind] = 2
    h  <- quadts[2] - quadts[1]
    quadwts <- quadwts*(h/3)

    Qmat    <- eval.basis(quadts, basis)
    Q.D1mat <- eval.basis(quadts, basis, 1)
    Q.D2mat <- eval.basis(quadts, basis, 2)

    #Initial estimat of basis coefficients
    Rmat = t(Q.D2mat)%*%(Q.D2mat*(quadwts%*%t(rep(1,nbasis))))
    basismat2 = t(Phi.mat)%*%Phi.mat;
    smooth_lambda = 1e-6   # smoothing parameter
    Bmat    = basismat2 + smooth_lambda*Rmat;
    #Initial basis coefficients
    initial_coef = ginv(Bmat)%*%t(Phi.mat)%*%observations

    #Passing to inner objective functions and obtain initial value for parameter cascading
    inner.input = list(observations, Phi.mat, 1e5, Qmat, Q.D1mat,quadts,quadwts)



      #Using nonlinear least square for optimization
      temp <- lsqnonlin(innerobj, initial_coef, ode.par = par.initial, derive.model = ode.model, input = inner.input,NLS = TRUE)
      new.ini.basiscoef <- matrix(temp$x,length(temp$x),1)
      #--------------------------------------------------------

      theta.final  <- lsqnonlin(outterobj, par.initial, options = controls,basis.initial = initial_coef, derivative.model = ode.model,inner.input = inner.input,NLS=TRUE)$x
      basiscoef <- lsqnonlin(innerobj, initial_coef, ode.par = theta.final, derive.model = ode.model, input = inner.input,NLS = TRUE)$x

    return(list(structural.par = theta.final, nuisance.par = basiscoef))
}

getlikelihood <- function(ode.model, basis,par_range){

}


cascade_est <- function(par.initial, derivative.model, basis, observations, obs.times,NLS = TRUE,nquadpts = NULL,controls = NULL,oned_range = NULL){
  #Input:
  #    par.initial : initial values for structural parameter
  #    ode.model   : ODE model
  #    basis       : basis object for approximating x(t)
  #    observations: observations including a vector indicating times of observations
  #    nquadpts    : number of quadrature points to use composite Simpson's rule for calculating integral


  #Output:



  #Comment:

  #number of parameters
  npar  <- length(par.initial)
  #time points of observations


  #If no 'quadpts' is given, then set it to 101
  if (is.null(nquadpts)){
      nquadpts = 101
  }



  #Evaluating basis functions at time points of observations
  #and stored as columns in Phi matrix
  Phi.mat <- eval.basis(times, basis)
  #Evaluating 1st derivative of basis functions
  D1.mat <- eval.basis(times, basis,1)
  #Evaluating 2nd derivative of basis functions
  D2.mat <- eval.basis(times, basis,2)

  #Calculate L2 penalty
  quadts   <- seq(0,1,length.out=nquadpts)
  nquad    <- length(quadts)
  quadwts  <- rep(1,nquad)
  even.ind <- seq(2,(nquad-1),by=2)
  odd.ind  <- seq(3,(nquad-2),by=2)
  quadwts[even.ind] = 4
  quadwts[odd.ind] = 2
  h  <- quadts[2] - quadts[1]
  quadwts <- quadwts*(h/3)

  Qmat    <- eval.basis(quadts, basis)
  Q.D1mat <- eval.basis(quadts, basis, 1)
  Q.D2mat <- eval.basis(quadts, basis, 2)

  #Initial estimat of basis coefficients
  Rmat = t(Q.D2mat)%*%(Q.D2mat*(quadwts%*%t(rep(1,nbasis))))
  basismat2 = t(Phi.mat)%*%Phi.mat;
  smooth_lambda = 1e-6   # smoothing parameter
  Bmat    = basismat2 + smooth_lambda*Rmat;
  #Initial basis coefficients
  initial_coef = ginv(Bmat)%*%t(Phi.mat)%*%observ;

  #Passing to inner objective functions and obtain initial value for parameter cascading
  inner.input = list(observ, Phi.mat, 1e5, Qmat, Q.D1mat,quadts,quadwts)

  if (NLS){

    #Using nonlinear least square for optimization
    temp <- lsqnonlin(innerobj, initial_coef, ode.par = par.initial, derive.model = derivative.model, input = inner.input,NLS = NLS)
    new.ini.basiscoef <- matrix(temp$x,length(temp$x),1)
    #--------------------------------------------------------

    theta.final  <- lsqnonlin(outterobj, par.initial, options = controls,basis.initial = initial_coef, derivative.model = derivative.model,inner.input = inner.input,NLS=NLS)$x
    basiscoef <- lsqnonlin(innerobj, initial_coef, ode.par = theta.final, derive.model = derivative.model, input = inner.input,NLS = NLS)$x

  }else{
    #Using optim for optimization
    new.ini.basiscoef <- optim(par=initial_coef, fn=innerobj,ode.par = par.initial,input = inner.input, derive.model = derivative.model,NLS=FALSE,
                                control = list(abstol=1e-10))$par
    if (npar == 1){
       loweb  = oned_range[1]
       upperb = oned_range[2]
       theta.final = optim(par = par.initial, fn = outterobj, basis.initial = new.ini.basiscoef, derivative.model = derivative.model , inner.input = inner.input,NLS=FALSE,method = ('Brent'),lower = loweb, upper = upperb)$par
    }else{
      theta.final = optim(par = par.initial, fn = outterobj, basis.initial = new.ini.basiscoef, derivative.model = derivative.model , inner.input = inner.input,NLS =FALSE)$par
    }
    basiscoef <- optim(par=initial_coef, fn=innerobj,ode.par = theta.final,input = inner.input, derive.model = derivative.model,NLS=FALSE,
                       control = list(abstol=1e-10))$par
  }

  return(list(structural.par = theta.final, nuisance.par = basiscoef))

}
