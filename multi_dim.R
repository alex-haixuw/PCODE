#############################################################################
#
#
#
#
#
#############################################################################


#Function list
#1. prepare_basis
#2. PC_ODE
#3. innerobj_multi
#4. outterobj_multi_nls
#5. innerobj
#6. outterobj
#7. cascade_nls
#
#
#



###############
#1. Function : 'prepare_basis'
#
#        Given a basis object, it evaluates the basis over both observation points and equally spaced quadtuare points.
#
#Input:
#           x: A basis object
#       times: A vector containing Observation points
#    nquadpts: A numeric value that defines the number of quadrature points
#
#Output:
#     Phi.mat: The matrix contains evaluations of all basis functions at every observation times, where evaluations of each basis
#                    function are  stored as the columns of the matrix.
#        Qmat: The matrix contains evaluations of all basis functions at every quadrature points in the support, where evaluations of
#                    each basis function are stored as the columns of the matrix.
#     Q.D1mat: The matrix contains evaluations of the first order derivative of all basis functions at every quadature points.
#     Q.D2mat: The matrix contains evaluations of the second order derivative of all basis functions at every quadature points.
#      quadts: A vector of quadrature points.
#     quadwts: A vector of quadrature weights for each quadrature points.
#
#Comment:
#
#
###############
prepare_basis <- function(x, times, nquadpts){

  #Evaluate basis functions over observation time points
  Phi.mat <- eval.basis(times, x)

  #Preparation to calculate L2 penalty
  #Evaluate basis function over quadrature points, and the number of
  #quadrature points, 'nquadpts', defined the density of points.
  quadts   <- seq(min(times),max(times),length.out=nquadpts)
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

  return(list(Phi.mat = Phi.mat, Qmat = Qmat, Q.D1mat = Q.D1mat, Q.D2mat = Q.D2mat, quadts= quadts, quadwts = quadwts))

}


###############
#2. Function: 'PC_ODE'
#
#
#Input:
#
#         data: A data frame or a matrix contain observations from each dimension of the ODE model.
#         time: A vector contain observation times or a matrix if time points are different between dimensions.
#    ode.model: A function that computes the time derivative of the ODE model given states variable evaluated at a given time.
#    par.names: A vector contains the names of structural parameters defined in the 'ode.model'.
#  state.names: A vector contains the names of state variables defined in the 'ode.model'.
#  par.initial: Initial value of structural parameters to be optimized.
#   basis.list: A list of basis objects for smoothing each dimension's observations. Can be the same or different across dimensions.
#       lambda: Penalty parameter.
#     controls: A list of control parameters. See ‘Details’.
#
#Output:
#
#   structural.par: The structural parameters of the ODE model.
#    nuissance.par: The nuissance parameters or the basis coefficients for interpolating observations.
#
#
#Comment:
#
###############
#data, time , ode.model, state.names, par.names, par.initial = NULL, basis.list = DEFAULT (order = 5,
PC_ODE <- function(data, time, ode.model,par.names,state.names,  par.initial, basis.list, lambda = NULL,controls = NULL){

      #Set up default controls for optimizations and quadrature evaluation
      con.default <- list(nquadpts = 101, smooth.lambda = 1e5, tau = 0.01, tolx = 1e-6,tolg = 1e-6, maxeval = 30)
      #Replace default with user's input
      con.default[(namc <- names(controls))] <- controls
      con.now  <- con.default


      #Evaluate basis functiosn for each state variable
      basis.eval.list <- lapply(basis.list, prepare_basis, times = time, nquadpts = con.now$nquadpts)

      #Number of dimensions of the ODE model
      ndim <- ncol(data)

      #Some initializations
      inner.input  <- list()
      initial_coef <- list()

      for (ii in 1:ndim){
        #For each dimension, obtain initial value for the nuissance parameters or the basis coefficients
        Rmat = t(basis.eval.list[[ii]]$Q.D2mat)%*%(basis.eval.list[[ii]]$Q.D2mat*(basis.eval.list[[ii]]$quadwts%*%t(rep(1,basis.list[[ii]]$nbasis))))
        basismat2 = t(basis.eval.list[[ii]]$Phi.mat)%*%basis.eval.list[[ii]]$Phi.mat;
        Bmat    = basismat2 + con.default$smooth_lambda*Rmat;
        #Initial basis coefficients
        initial_coef[[ii]] = ginv(Bmat)%*%t(basis.eval.list[[ii]]$Phi.mat)%*%data[,ii]
        inner.input[[ii]]  = list(data[,ii], basis.eval.list[[ii]]$Phi.mat, lambda,
                                  basis.eval.list[[ii]]$Qmat, basis.eval.list[[ii]]$Q.D1mat,
                                  basis.eval.list[[ii]]$quadts, basis.eval.list[[ii]]$quadwts,time,
                                  state.names,model.names)
      }


      #Searching for a better starting value for optimization
      #temp <- lsqnonlin(innerobj_multi, unlist(initial_coef), ode.par = par.initial, derive.model = ode.model, input = inner.input, NLS=TRUE,options = list(tau = 0.01))$x

      #Running optimization for outter objective function to obtain structural parameter
      par.final <- lsqnonlin(options = list(maxeval = con.default$maxeval,tau = con.default$tau),outterobj_multi_nls, par.initial,basis.initial = unlist(initial_coef), derivative.model = ode.model, inner.input = inner.input)$x
      #Condition on the obtained structural parameter, calculate the nuissance parameter or the basis coefficients to interpolate data
      basis.coef.final <- lsqnonlin(innerobj_multi, unlist(initial_coef), ode.par = par.final, derive.model = ode.model, input = inner.input, NLS=TRUE,options = list(tau =0.01))$x


      return(list(structural.par = par.final, nuissance.par = basis.coef.final))
}


#####
#'innerobj_multi':
#
#            Conditioning on the structural parameters of the ODE model, 'ode.par', this functions returns either the
#            inner objective function value (the sum of squared errors) or the residual between the observation and
#            and the basis expansion estimates. The form of returning value depends on the later called optimizers.
#
#
#Input:
#     basis_coef: A single vector containing the basis coefficients of all dimensions of the ODE model.
#        ode.par: A vector or a numeric value containing the structural parameter of the ODE model.
#   derive.model: The function defines the time derivative of of each state variable
#           NLS : A logic value determines whether the function should return a numeric value (NLS = FALSE)
#                                                                    or a vector of residuals (NLS = TRUE)
#Output
#
#         residual.vec : A vector containing the residuals from using basis expansion to estimate observations
#                        and the quarature evaluations to approximate the ODE penalty,
#   sum(residual.vec^2): A numeric value correspondes to the inner objective function value.
#
#Comment:
#
#
#
innerobj_multi  <- function(basis_coef, ode.par, input, derive.model,NLS=TRUE){


    #Get variables  from 'input'
    ndim     <- length(input)
    npoints  <- length(unlist(input[[1]][8]))
    nbasis   <- rep(NA, ndim+1)
    Xhat     <- matrix(NA, nrow = npoints, ncol = ndim)
    residual <- matrix(NA, nrow = npoints, ncol = ndim)

    state.names <- input[[1]][[9]]
    model.names <- input[[1]][[10]]
    Xt   <- matrix(NA, nrow = length(input[[1]][[7]]),ncol= ndim)
    dXdt <- matrix(NA, nrow = length(input[[1]][[7]]),ncol= ndim)

    #Turn a single vector 'basis_coef' into seperate locations
    coef.list <- list()
    nbasis[1] <- 0
    for (jj in 1:ndim){
       nbasis[jj+1]    <- ncol(input[[jj]][[2]])
       coef.list[[jj]] <- basis_coef[(nbasis[jj]+1):(nbasis[jj]+nbasis[jj+1])]
       #Basis expansion for j-th dimension
       Xhat[,jj]       <- input[[jj]][[2]] %*% coef.list[[jj]]
       #Residual from basis expansion
       residual[,jj]   <- input[[jj]][[1]] - Xhat[,jj]
       Xt[,jj]         <- input[[jj]][[4]] %*% coef.list[[jj]]
       dXdt[,jj]       <- input[[jj]][[5]] %*% coef.list[[jj]]
    }

   #Ode penalty


   names(ode.par) <- model.names

   temp_fun  <- function(x, temp = ode.par,names = state.names){
       names(x) <- state.names
       return(unlist(derive.model(state = x, parameters = temp)))
   }
   temp_eval <- t(apply(Xt,1,temp_fun))

   temp_list <- list(dXdt, temp_eval)
   temp_penalty_resid <- Reduce("-", temp_list)

   lambda <- input[[1]][[3]]
   penalty_residual <- matrix(NA, nrow = nrow(temp_eval),ncol=ncol(temp_eval))
   for (jj in 1:ndim){
     penalty_residual[,jj] <- sqrt(lambda) * sqrt(input[[jj]][[7]]) * temp_penalty_resid[,jj]
   }




   residual.vec <- c(as.vector(residual), as.vector(penalty_residual))
     if (NLS){
       return(residual.vec)
     }else{
       return(sum(residual.vec^2))
     }

}

#####
#
#####
outterobj_multi_nls <- function(ode.parameter, basis.initial, derivative.model, inner.input,NLS=TRUE){
  #Convergence of basis coefficients seems to happen before 'maxeval'.

  inner_coef <- lsqnonlin(innerobj_multi, basis.initial, ode.par = ode.parameter, derive.model = derivative.model, options = list(maxeval = 100,tolx=1e-6,tolg=1e-6), input = inner.input)$x
  #inner_coef <- optim(basis.initial,innerobj_multi,  ode.par = ode.parameter, derive.model = derivative.model, input = inner.input,NLS=FALSE)$par
  #Get variables  from 'input'
  ndim     <- length(inner.input)
  npoints  <- length(inner.input[[1]][[8]])
  basisnumber   <- rep(NA, ndim+1)
  Xhat     <- matrix(NA, nrow = npoints, ncol = ndim)
  residual <- matrix(NA, nrow = npoints, ncol = ndim)
  #Turn a single vector 'inner_coef' into seperate locations
  coef.list <- list()
  basisnumber[1] <- 0
  for (jj in 1:ndim){
    basisnumber[jj+1]    <- ncol(inner.input[[jj]][[2]])
     coef.list[[jj]] <- inner_coef[(basisnumber[jj]+1):(basisnumber[jj]+basisnumber[jj+1])]
     #Basis expansion for j-th dimension
     Xhat[,jj]       <- inner.input[[jj]][[2]] %*% coef.list[[jj]]
     #Residual from basis expansion
     residual[,jj]   <- inner.input[[jj]][[1]] - Xhat[,jj]
  }

  if(NLS){
    return(as.vector(residual))
  }else{
    return(sum(residual^2))
  }


}





#####
#  'innerobj'
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
#####
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
