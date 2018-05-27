#############################################################################
#
#
#
#
#
#############################################################################


#Function list
#1. prepare_basis
#2.
#
#
#
#
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
#         data:
#         time:
#    ode.model:
#  state.names:
#    par.names:
#  par.initial:
#   basis.list:
#       lambda:
#     controls:
#
#Output:
#
#   structural.par:
#    nuissance.par:
#
#
#Comment:
#
###############
#data, time , ode.model, state.names, par.names, par.initial = NULL, basis.list = DEFAULT (order = 5,



#PC_ODE <- function(par.initial, ode.model, basis.list,observations, observation.times, state.names, model.names,controls){
PC_ODE <- function(data, time, ode.model,state.names, par.names, par.initial, basis.list, lambda = NULL,controls = NULL){

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
