###########
#
#Copyright (C) <2019> <Haixu Wang>
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>
###########


#' @title Parameter Cascade Method for Ordinary Differential Equation Models
#' @description Obtain estimates of both structural and nuisance parameters of an ODE model by parameter cascade method.
#' @usage pcode(data, time, ode.model, par.names, state.names,
#'              likelihood.fun, par.initial, basis.list,lambda,controls)
#' @param        data A data frame or a matrix contain observations from each dimension of the ODE model.
#' @param         time A vector contain observation times or a matrix if time points are different between dimensions.
#' @param    ode.model An R function that computes the time derivative of the ODE model given observations of states variable and structural parameters.
#' @param    par.names The names of structural parameters defined in the 'ode.model'.
#' @param  state.names The names of state variables defined in the 'ode.model'.
#' @param likelihood.fun A likelihood function passed to PCODE in case of that the error terms do not have a Normal distribution.
#' @param  par.initial Initial value of structural parameters to be optimized.
#' @param   basis.list A list of basis objects for smoothing each dimension's observations. Can be the same or different across dimensions.
#' @param    lambda Penalty parameter for controling the fidelity of interpolation.
#' @param     controls A list of control parameters. See Details.
#'
#' @details The \code{controls} argument is a list providing addition inputs for the nonlinear least square optimizer or general optimizer \code{optim}:
#' \describe{
#'\item{\code{nquadpts}}{Determine the number of quadrature points for approximating an integral. Default is 101.}
#'\item{\code{smooth.lambda}}{Determine the smoothness penalty for obtaining initial value of nuisance parameters.}
#'\item{\code{tau}}{Initial value of Marquardt parameter. Small values indicate good initial values for structural parameters.}
#'\item{\code{tolx}}{Tolerance for parameters of objective functions. Default is set at 1e-6.}
#'\item{\code{tolg}}{Tolerance for the gradient of parameters of objective functions. Default is set at 1e-6.}
#'\item{\code{maxeval}}{The maximum number of evaluation of the outter optimizer. Default is set at 20.}
#'}
#'
#'
#' @return   \item{structural.par}{The structural parameters of the ODE model.}
#'
#' @return   \item{nuisance.par}{The nuisance parameters or the basis coefficients for interpolating observations.}
#' @import fda
#' @import MASS
#' @import stats
#' @importFrom deSolve ode
#' @importFrom pracma jacobian Norm eye
#' @examples \donttest{library(fda)
#'library(deSolve)
#'library(MASS)
#'library(pracma)
#'#Simple ode model example
#'#define model parameters
#'model.par   <- c(theta = c(0.1))
#'#define state initial value
#'state       <- c(X     = 0.1)
#'#Define model for function 'ode' to numerically solve the system
#'ode.model <- function(t, state,parameters){
#'  with(as.list(c(state,parameters)),
#'       {
#'         dX <- theta*X*(1-X/10)
#'         return(list(dX))
#'       })
#'}
#'#Observation time points
#'times <- seq(0,100,length.out=101)
#'#Solve the ode model
#'desolve.mod <- ode(y=state,times=times,func=ode.model,parms = model.par)
#'#Prepare for doing parameter cascading method
#'#Generate basis object for interpolation and as argument of pcode
#'#21 konts equally spaced within [0,100]
#'knots <- seq(0,100,length.out=21)
#'#order of basis functions
#'norder <- 4
#'#number of basis funtions
#'nbasis <- length(knots) + norder - 2
#'#creating Bspline basis
#'basis  <- create.bspline.basis(c(0,100),nbasis,norder,breaks = knots)

#'#Add random noise to ode solution for simulating data
#'nobs  <- length(times)
#'scale <- 0.1
#'noise <- scale*rnorm(n = nobs, mean = 0 , sd = 1)
#'observation <- desolve.mod[,2] + noise
#'#parameter estimation
#'pcode(data = observation, time = times, ode.model = ode.model,
#'                      par.initial = 0.1, par.names = 'theta',state.names = 'X',
#'                      basis.list = basis, lambda = 1e2)
#'}
#' @export
pcode <- function(data, time, ode.model, par.names, state.names, likelihood.fun = NULL, par.initial, basis.list, lambda, controls = list()) {
    # Set up default controls for optimizations and quadrature evaluation
    con.default <- list(nquadpts = 101, smooth.lambda = 100, tau = 0.01, tolx = 1e-06, tolg = 1e-06, maxeval = 20,verbal = 0)
    # Replace default with user's input
    con.default[(namc <- names(controls))] <- controls
    con.now <- con.default

    if (length(state.names) == 1) {
        if (!is.function(likelihood.fun)) {
            result <- pcode_1d(data = data, time = time, ode.model = ode.model, par.initial = par.initial, par.names = par.names,
                basis = basis.list, lambda = lambda, controls = con.now)
            return(list(structural.par = result$structural.par, nuisance.par = result$nuisance.par))
        } else {
            result <- pcode_lkh_1d(data = data, time = time, likelihood.fun = likelihood.fun, ode.model = ode.model,
                par.initial = par.initial, par.names = par.names, state.names = state.names, basis.list = basis.list,
                lambda = lambda, controls = con.now)
            return(list(structural.par = result$structural.par, nuisance.par = result$nuisance.par))
        }
    }else {

       #check whether dimension of lambda matches the
      if ((length(lambda) == 1)&&(length(state.names >1))){
          multi.lambda <- rep(lambda, length(state.names))
      }

      if (length(lambda) == length(state.names)){
        multi.lambda <- lambda
      }


      #Check if there are any missing variables from the ODE model
      #Calling fitting function for ODE model with missing variables
      #Multiple dimension case
      if(length(state.names) != ncol(data)){
              result <- pcode_missing(data = data, time = time, ode.model = ode.model,
                                       par.names = par.names, state.names = state.names, likelihood.fun = likelihood.fun,
                                       par.initial = par.initial , basis.list = basis.list, lambda = lambda, controls = con.now)
              return(list(structural.par = result$structural.par, nuisance.par = result$nuisance.par))
      }else{
        # Evaluate basis functiosn for each state variable
        basis.eval.list <- lapply(basis.list, prepare_basis, times = time, nquadpts = con.now$nquadpts)

        # Number of dimensions of the ODE model
        ndim <- ncol(data)

        # Some initializations
        inner.input <- list()
        initial_coef <- list()

        for (ii in 1:ndim) {
          # For each dimension, obtain initial value for the nuisance parameters or the basis coefficients
          Rmat = t(basis.eval.list[[ii]]$Q.D2mat) %*% (basis.eval.list[[ii]]$Q.D2mat * (basis.eval.list[[ii]]$quadwts %*%
                                                                                          t(rep(1, basis.list[[ii]]$nbasis))))
          basismat2 = t(basis.eval.list[[ii]]$Phi.mat) %*% basis.eval.list[[ii]]$Phi.mat
          Bmat = basismat2 + con.now$smooth.lambda * Rmat
          # Initial basis coefficients
          initial_coef[[ii]] = ginv(Bmat) %*% t(basis.eval.list[[ii]]$Phi.mat) %*% data[, ii]
          # // TODO: need to check colnames of data
          inner.input[[ii]] = list(data[, ii], basis.eval.list[[ii]]$Phi.mat, multi.lambda[ii], basis.eval.list[[ii]]$Qmat,
                                   basis.eval.list[[ii]]$Q.D1mat, basis.eval.list[[ii]]$quadts, basis.eval.list[[ii]]$quadwts, time, state.names,
                                   par.names)
        }


        # Searching for a better starting value for optimization temp <- nls_optimize(innerobj_multi, unlist(initial_coef),
        # ode.par = par.initial, derive.model = ode.model, input = inner.input, NLS=TRUE,options = list(tau = 0.01))$x

        # Running optimization for outter objective function to obtain structural parameter
        par.final <- nls_optimize(options = list(maxeval = con.now$maxeval, tau = con.now$tau), outterobj_multi_nls,
                                  par.initial, basis.initial = unlist(initial_coef), derivative.model = ode.model, inner.input = inner.input,
                                  verbal =con.now$verbal )$par
        # Condition on the obtained structural parameter, calculate the nuisance parameter or the basis coefficients to
        # interpolate data
        basis.coef.final <- nls_optimize.inner(innerobj_multi, unlist(initial_coef), ode.par = par.final, derive.model = ode.model,
                                               input = inner.input, NLS = TRUE, options = list(tau = 0.01))$par

        return(list(structural.par = par.final, nuisance.par = basis.coef.final))
      }

    }

}


#' @title Evaluate basis objects over observation times and quadrature points
#' @description Calculate all basis functions over observation time points and store them as columns in a single matrix for each dimension. Also include first and second order derivative. Repeat over quadrature points.
#' @usage prepare_basis(basis, times, nquadpts)
#' @param basis A basis object.
#' @param times The vector contain observation times for corresponding dimension.
#' @param nquadpts Number of quadrature points will be used later for approximating integrals.
#'
#' @return   \item{Phi.mat}{Evaluations of all basis functions stored as columns in the matrix.}
#' @return   \item{Qmat}{Evaluations of all basis functions over quadrature points stored as columns in the matrix.}
#' @return   \item{Q.D1mat}{Evaluations of first order derivative all basis functions over quadrature points stored as columns in the matrix.}
#' @return   \item{Q.D2mat}{Evaluations of second order derivative all basis functions over quadrature points stored as columns in the matrix.}
#' @return   \item{quadts}{Quadrature points.}
#' @return   \item{quadwts}{Quadrature weights.}
#'
prepare_basis <- function(basis, times, nquadpts) {

    # Evaluate basis functions over observation time points
    Phi.mat  <- eval.basis(times, basis)

    # Preparation to calculate L2 penalty Evaluate basis function over quadrature points, and the number of quadrature
    # points, 'nquadpts', defined the density of points.
    quadts   <- seq(min(times), max(times), length.out = nquadpts)
    nquad    <- length(quadts)
    quadwts  <- rep(1, nquad)
    even.ind <- seq(2, (nquad - 1), by = 2)
    odd.ind  <- seq(3, (nquad - 2), by = 2)
    quadwts[even.ind] = 4
    quadwts[odd.ind] = 2
    h        <- quadts[2] - quadts[1]
    quadwts  <- quadwts * (h/3)

    Qmat     <- eval.basis(quadts, basis)
    Q.D1mat  <- eval.basis(quadts, basis, 1)
    Q.D2mat  <- eval.basis(quadts, basis, 2)

    return(list(Phi.mat = Phi.mat, Qmat = Qmat, Q.D1mat = Q.D1mat, Q.D2mat = Q.D2mat, quadts = quadts, quadwts = quadwts))

}

#' @title Inner objective function (multiple dimension version)
#' @description An objective function combines the sum of squared error of basis expansion estimates and the penalty controls how those estimates fail to satisfies the ODE model
#' @usage innerobj_multi(basis_coef, ode.par, input, derive.model,NLS)
#' @param basis_coef Basis coefficients for interpolating observations given a basis object.
#' @param ode.par Structural parameters of the ODE model.
#' @param input Contains dependencies for the optimization, including observations, penalty parameter lambda, and etc..
#' @param derive.model The function defines the ODE model and is the same as the \code{ode.model} in \code{pcode}.
#' @param NLS Default is \code{TRUE} so the function returns vector of residuals, and otherwise returns sum of squared errors.
#'
#' @return   \item{residual.vec}{Vector of residuals and evaluation of penalty function on quadrature points for approximating the integral.}
innerobj_multi <- function(basis_coef, ode.par, input, derive.model, NLS = TRUE) {


    # Retrieve variables from 'input' get the dimesion of the ODE model
    ndim <- length(input)
    npoints <- length(unlist(input[[1]][8]))
    nbasis <- rep(NA, ndim + 1)
    Xhat <- matrix(NA, nrow = npoints, ncol = ndim)
    residual <- matrix(NA, nrow = npoints, ncol = ndim)

    state.names <- input[[1]][[9]]
    model.names <- input[[1]][[10]]
    Xt <- matrix(NA, nrow = length(input[[1]][[7]]), ncol = ndim)
    dXdt <- matrix(NA, nrow = length(input[[1]][[7]]), ncol = ndim)

    # Turn a single vector 'basis_coef' into seperate locations
    coef.list <- list()
    nbasis[1] <- 0
    for (jj in 1:ndim) {
        nbasis[jj + 1] <- ncol(input[[jj]][[2]])
        coef.list[[jj]] <- basis_coef[(nbasis[jj] + 1):(nbasis[jj] + nbasis[jj + 1])]
        # Basis expansion for j-th dimension
        Xhat[, jj] <- input[[jj]][[2]] %*% coef.list[[jj]]
        # Residual from basis expansion
        residual[, jj] <- input[[jj]][[1]] - Xhat[, jj]
        Xt[, jj] <- input[[jj]][[4]] %*% coef.list[[jj]]
        dXdt[, jj] <- input[[jj]][[5]] %*% coef.list[[jj]]
    }

    # Ode penalty


    names(ode.par) <- model.names

    temp_fun <- function(x, temp = ode.par, names = state.names) {
        names(x) <- state.names
        return(unlist(derive.model(state = x, parameters = temp)))
    }
    temp_eval <- t(apply(Xt, 1, temp_fun))

    temp_list <- list(dXdt, temp_eval)
    temp_penalty_resid <- Reduce("-", temp_list)

    lambda <- input[[1]][[3]]
    penalty_residual <- matrix(NA, nrow = nrow(temp_eval), ncol = ncol(temp_eval))
    for (jj in 1:ndim) {
        penalty_residual[, jj] <- sqrt(lambda) * sqrt(input[[jj]][[7]]) * temp_penalty_resid[, jj]
    }


    residual.mat <- rbind(residual,penalty_residual)
    residual.vec <- as.vector(residual.mat)
    if (NLS) {
      return(residual.vec)
    } else {
      return(sum(residual.vec^2))
    }

    # residual.vec <- c(as.vector(residual), as.vector(penalty_residual))
    # if (NLS) {
    #     return(residual.vec)
    # } else {
    #     return(sum(residual.vec^2))
    # }

}

#' @title Outter objective function (multiple dimension version)
#' @description An objective function of the structural parameter computes the measure of fit for the basis expansion.
#' @usage outterobj_multi_nls(ode.parameter, basis.initial, derivative.model, inner.input, NLS)
#' @param ode.parameter Structural parameters of the ODE model.
#' @param basis.initial Initial values of the basis coefficients for nonlinear least square optimization.
#' @param derivative.model The function defines the ODE model and is the same as the \code{ode.model} in \code{pcode}.
#' @param inner.input Input that will be passed to the inner objective function. Contains dependencies for the optimization, including observations, penalty parameter lambda, and etc..
#' @param NLS Default is \code{TRUE} so the function returns vector of residuals, and otherwise returns sum of squared errors.
#'
#' @return   \item{residual}{Vector of residuals and evaluation of penalty function on quadrature points for approximating the integral.}
#'
outterobj_multi_nls <- function(ode.parameter, basis.initial, derivative.model, inner.input, NLS = TRUE) {
    # Convergence of basis coefficients seems to happen before 'maxeval'.

    inner_coef <- nls_optimize.inner(innerobj_multi, basis.initial, ode.par = ode.parameter, derive.model = derivative.model,
        options = list(maxeval = 50, tolx = 1e-06, tolg = 1e-06), input = inner.input)$par
    ndim <- length(inner.input)
    npoints <- length(inner.input[[1]][[8]])
    basisnumber <- rep(NA, ndim + 1)
    Xhat <- matrix(NA, nrow = npoints, ncol = ndim)
    residual <- matrix(NA, nrow = npoints, ncol = ndim)
    # Turn a single vector 'inner_coef' into seperate locations
    coef.list <- list()
    basisnumber[1] <- 0
    for (jj in 1:ndim) {
        basisnumber[jj + 1] <- ncol(inner.input[[jj]][[2]])
        coef.list[[jj]] <- inner_coef[(basisnumber[jj] + 1):(basisnumber[jj] + basisnumber[jj + 1])]
        # Basis expansion for j-th dimension
        Xhat[, jj] <- inner.input[[jj]][[2]] %*% coef.list[[jj]]
        # Residual from basis expansion
        residual[, jj] <- inner.input[[jj]][[1]] - Xhat[, jj]
    }

    if (NLS) {
        return(as.vector(residual))
    } else {
        return(sum(residual^2))
    }


}





#' @title Inner objective function (Single dimension version)
#' @description An objective function combines the sum of squared error of basis expansion estimates and the penalty controls how those estimates fail to satisfies the ODE model
#' @usage innerobj(basis_coef, ode.par, input, derive.model,NLS)
#' @param basis_coef Basis coefficients for interpolating observations given a basis object.
#' @param ode.par Structural parameters of the ODE model.
#' @param input Contains dependencies for the optimization, including observations, penalty parameter lambda, and etc..
#' @param derive.model The function defines the ODE model and is the same as the ode.model in 'pcode'
#' @param NLS Default is \code{TRUE} so the function returns vector of residuals, and otherwise returns sum of squared errors.
#'
#' @return   \item{residual.vec}{Vector of residuals and evaluation of penalty function on quadrature points for approximating the integral.}
#'
innerobj <- function(basis_coef, ode.par, input, derive.model, NLS = TRUE) {
    yobs <- input[[1]]
    Phi.mat <- input[[2]]
    lambda <- input[[3]]
    Qmat <- input[[4]]
    Q.D1mat <- input[[5]]
    quadts <- input[[6]]
    quadwts <- input[[7]]

    xhat <- Phi.mat %*% basis_coef
    # e_{ij}
    residual <- yobs - xhat
    # Sum of squared errors SSE <- sum((yobs-xhat)^2)

    # PEN_{i}
    penalty_residual <- rep(NA, length(quadwts))
    # Use composite Simpson's rule to calculate the integral of residual
    Xt <- Qmat %*% basis_coef
    dXdt <- Q.D1mat %*% basis_coef

    # Evaluate f(X)
    rightside <- rep(NA, length(quadwts))

    tempfun <- function(x, temp = ode.par) {
        return(derive.model(state = c(X = x), parameters = temp))
    }
    rightside <- unlist(lapply(Xt, tempfun))


    penalty_residual <- sqrt(lambda) * sqrt(quadwts) * (dXdt - rightside)

    # ode_penalty <- sum(lambda * quadwts * (penalty_residual)^2)

    residual.vec <- c(residual, penalty_residual)
    if (NLS) {
        return(residual.vec)
    } else {
        return(sum(residual.vec^2))
    }
}


#' @title Outter objective function (Single dimension version)
#' @description An objective function of the structural parameter computes the measure of fit.
#' @usage outterobj(ode.parameter, basis.initial, derivative.model, inner.input, NLS)
#' @param ode.parameter Structural parameters of the ODE model.
#' @param basis.initial Initial values of the basis coefficients for nonlinear least square optimization.
#' @param derivative.model The function defines the ODE model and is the same as the ode.model in 'pcode'
#' @param inner.input Input that will be passed to the inner objective function. Contains dependencies for the optimization, including observations, penalty parameter lambda, and etc..
#' @param NLS Default is \code{TRUE} so the function returns vector of residuals, and otherwise returns sum of squared errors.
#'
#' @return   \item{residual}{Vector of residuals and evaluation of penalty function on quadrature points for approximating the integral.}
#'
outterobj <- function(ode.parameter, basis.initial, derivative.model, inner.input, NLS) {


    if (NLS) {
        inner_coef <- nls_optimize.inner(innerobj, basis.initial, options = list(maxeval = 50), ode.par = ode.parameter,
            derive.model = derivative.model, input = inner.input)$par
    } else {
        inner_coef <- optim(par = basis.initial, fn = innerobj, ode.par = ode.parameter, derive.model = derivative.model,
            input = inner.input, NLS = FALSE, control = list(abstol = 1e-10))$par
    }



    yobs <- inner.input[[1]]
    Phi.mat <- inner.input[[2]]


    xhat <- Phi.mat %*% inner_coef
    # e_{ij}
    residual <- yobs - xhat
    if (NLS) {
        return(residual)
    } else {
        return(sum(residual^2))
    }

}

#' @title Parameter Cascade Method for Ordinary Differential Equation Models (Single dimension version)
#' @description Obtain estiamtes of structural parameters of an ODE model by parameter cascade method.
#' @usage pcode_1d(data, time, ode.model, par.initial,par.names, basis,lambda,controls = list())
#' @param        data A data frame or a vector contains observations from the ODE model.
#' @param         time The vector contain observation times.
#' @param    ode.model Defined R function that computes the time derivative of the ODE model given observations of states variable.
#' @param  par.names The names of structural parameters defined in the 'ode.model'.
#' @param  par.initial Initial value of structural parameters to be optimized.
#' @param   basis A basis objects for smoothing observations.
#' @param       lambda Penalty parameter.
#' @param     controls A list of control parameters. See ‘Details’.
#'
#' @return   \item{structural.par}{The structural parameters of the ODE model.}
#' @return    \item{nuisance.par}{The nuisance parameters or the basis coefficients for interpolating observations.}
pcode_1d <- function(data, time, ode.model, par.initial, par.names, basis, lambda, controls = list()) {

    # number of parameters
    npar <- length(par.initial)

    nbasis <- basis$nbasis

    # Evaluating basis functions at time points of observations and stored as columns in Phi matrix
    Phi.mat <- eval.basis(time, basis)
    # Evaluating 1st derivative of basis functions
    D1.mat <- eval.basis(time, basis, 1)
    # Evaluating 2nd derivative of basis functions
    D2.mat <- eval.basis(time, basis, 2)

    # Calculate L2 penalty
    quadts <- seq(min(time), max(time), length.out = controls$nquadpts)
    nquad <- length(quadts)
    quadwts <- rep(1, nquad)
    even.ind <- seq(2, (nquad - 1), by = 2)
    odd.ind <- seq(3, (nquad - 2), by = 2)
    quadwts[even.ind] = 4
    quadwts[odd.ind] = 2
    h <- quadts[2] - quadts[1]
    quadwts <- quadwts * (h/3)

    Qmat <- eval.basis(quadts, basis)
    Q.D1mat <- eval.basis(quadts, basis, 1)
    Q.D2mat <- eval.basis(quadts, basis, 2)

    # Initial estimat of basis coefficients
    Rmat = t(Q.D2mat) %*% (Q.D2mat * (quadwts %*% t(rep(1, nbasis))))
    basismat2 = t(Phi.mat) %*% Phi.mat
    Bmat = basismat2 + controls$smooth.lambda * Rmat
    # Initial basis coefficients
    initial_coef = ginv(Bmat) %*% t(Phi.mat) %*% data

    # Passing to inner objective functions and obtain initial value for parameter cascading
    inner.input = list(data, Phi.mat, lambda, Qmat, Q.D1mat, quadts, quadwts)



    # Using nonlinear least square for optimization temp <- nls_optimize(innerobj, initial_coef, ode.par = par.initial,
    # derive.model = ode.model, input = inner.input,NLS = TRUE) new.ini.basiscoef <-
    # matrix(temp$par,length(temp$par),1)
    #--------------------------------------------------------
    names(par.initial) <- par.names
    theta.final <- nls_optimize(outterobj, par.initial, basis.initial = initial_coef, derivative.model = ode.model,
        inner.input = inner.input, NLS = TRUE, verbal = controls$verbal)$par

    basiscoef <- nls_optimize.inner(innerobj, initial_coef, ode.par = theta.final, derive.model = ode.model, input = inner.input,
        NLS = TRUE)$par

    return(list(structural.par = theta.final, nuisance.par = basiscoef))
}

#' @title       Optimizer for non-linear least square problems
#' @description Obtain the solution to minimize the sum of squared errors of the defined function \code{fun} by levenberg-marquardt method. Adapted from PRACMA package.
#' @usage       nls_optimize(fun, x0, ..., options,verbal)
#' @param       fun The function returns the vector of weighted residuals.
#' @param       x0  The initial value for optimization.
#' @param       verbal Default = 1 for printing iteration and other for suppressing
#' @param       ... Parameters to be passed for \code{fun}
#' @param       options Additional optimization controls.
#'
#' @return   \item{par}{The solution to the non-linear least square problem, the same size as \code{x0}}

nls_optimize <- function(fun, x0, ..., options = list(), verbal = 0) {
    stopifnot(is.numeric(x0))
    opts <- list(tau = 0.001, tolx = 1e-06, tolg = 1e-06, maxeval = 20)
    namedOpts <- match.arg(names(options), choices = names(opts), several.ok = TRUE)
    if (!is.null(names(options)))
        opts[namedOpts] <- options
    tau <- opts$tau
    tolx <- opts$tolx
    tolg <- opts$tolg
    maxeval <- opts$maxeval
    fct <- match.fun(fun)
    fun <- function(x) fct(x, ...)
    n <- length(x0)
    m <- length(fun(x0))
    x <- x0
    r <- fun(x)
    f <- 0.5 * sum(r^2)
    #Jacobian takes time to compute
    J <- jacobian(fun, x)
    g <- t(J) %*% r
    ng <- Norm(g, Inf)
    A <- t(J) %*% J
    mu <- tau * max(diag(A))
    nu <- 2
    nh <- 0
    errno <- 0
    k <- 1
    if (verbal == 1) {
        message("#######################################")
        message("Starting optimization")
        message("#######################################")
    }


    while (k < maxeval) {

        if (k == 1) {
            if (verbal == 1) {
                message(paste("Current iteration: ", k, collapse = ""))
                message(paste("Current loss value: ", round(2 * f, digits = 6), collapse = ""))
            }
        }

        k <- k + 1
        R <- chol(A + mu * eye(n))
        h <- c(-t(g) %*% chol2inv(R))
        nh <- Norm(h)
        nx <- tolx + Norm(x)
        if (nh <= tolx * nx) {
            errno <- 1
            break
        }
        xnew <- x + h
        h <- xnew - x
        dL <- sum(h * (mu * h - g))/2
        rn <- fun(xnew)
        fn <- 0.5 * sum(rn^2)
        if (verbal == 1) {
            message(paste("Current iteration: ", k, collapse = ""))
            message(paste("Current loss value: ", round(2 * fn, digits = 6), collapse = ""))
        }

        Jn <- jacobian(fun, xnew)
        if (length(rn) != length(r)) {
            df <- f - fn
        } else {
            df <- sum((r - rn) * (r + rn))/2
        }
        if (dL > 0 && df > 0) {
            x <- xnew
            f <- fn
            J <- Jn
            r <- rn
            A <- t(J) %*% J
            g <- t(J) %*% r
            ng <- Norm(g, Inf)
            mu <- mu * max(1/3, 1 - (2 * df/dL - 1)^3)
            nu <- 2
        } else {
            mu <- mu * nu
            nu <- 2 * nu
        }
        if (ng <= tolg) {
            errno <- 2
            break
        }
    }
    if (k >= maxeval)
        errno <- 3
    errmess <- switch(errno, "Stopped by small x-step.", "Stopped by small gradient.", "No. of function evaluations exceeded.")
    return(list(par = c(xnew), ssq = sum(fun(xnew)^2), ng = ng, nh = nh, mu = mu, neval = k, errno = errno, errmess = errmess))
}

#' @title       Optimizer for non-linear least square problems (for inner objective functions)
#' @description Obtain the solution to minimize the sum of squared errors of the defined function \code{fun} by levenberg-marquardt method. Adapted from PRACMA package.
#' @usage       nls_optimize.inner(fun, x0, ..., options)
#' @param       fun The function returns the vector of weighted residuals.
#' @param       x0  The initial value for optimization.
#' @param       ... Parameters to be passed for \code{fun}
#' @param       options Additional optimization controls.
#'
#' @return   \item{par}{The solution to the non-linear least square problem, the same size as \code{x0}}

nls_optimize.inner <- function(fun, x0,..., options = list()) {
    stopifnot(is.numeric(x0))
    opts <- list(tau = 0.001, tolx = 1e-06, tolg = 1e-06, maxeval = 30)
    namedOpts <- match.arg(names(options), choices = names(opts), several.ok = TRUE)
    if (!is.null(names(options)))
        opts[namedOpts] <- options
    tau <- opts$tau
    tolx <- opts$tolx
    tolg <- opts$tolg
    maxeval <- opts$maxeval
    fct <- match.fun(fun)
    fun <- function(x) fct(x, ...)
    n <- length(x0)
    m <- length(fun(x0))
    x <- x0
    r <- fun(x)
    f <- 0.5 * sum(r^2)
    J <- jacobian(fun, x)
    g <- t(J) %*% r
    ng <- Norm(g, Inf)
    A <- t(J) %*% J
    mu <- tau * max(diag(A))
    nu <- 2
    nh <- 0
    errno <- 0
    k <- 1


    while (k < maxeval) {


        k <- k + 1
        R <- chol(A + mu * eye(n))
        h <- c(-t(g) %*% chol2inv(R))
        nh <- Norm(h)
        nx <- tolx + Norm(x)
        if (nh <= tolx * nx) {
            errno <- 1
            break
        }
        xnew <- x + h
        h <- xnew - x
        dL <- sum(h * (mu * h - g))/2
        rn <- fun(xnew)
        fn <- 0.5 * sum(rn^2)

        Jn <- jacobian(fun, xnew)
        if (length(rn) != length(r)) {
            df <- f - fn
        } else {
            df <- sum((r - rn) * (r + rn))/2
        }
        if (dL > 0 && df > 0) {
            x <- xnew
            f <- fn
            J <- Jn
            r <- rn
            A <- t(J) %*% J
            g <- t(J) %*% r
            ng <- Norm(g, Inf)
            mu <- mu * max(1/3, 1 - (2 * df/dL - 1)^3)
            nu <- 2
        } else {
            mu <- mu * nu
            nu <- 2 * nu
        }
        if (ng <= tolg) {
            errno <- 2
            break
        }
    }
    if (k >= maxeval)
        errno <- 3
    errmess <- switch(errno, "Stopped by small x-step.", "Stopped by small gradient.", "No. of function evaluations exceeded.")
    return(list(par = c(xnew), ssq = sum(fun(xnew)^2), ng = ng, nh = nh, mu = mu, neval = k, errno = errno, errmess = errmess))
}

#' @title       Find optimial penalty parameter lambda by cross-validation.
#' @description Obtain the optimal sparsity parameter given a search grid based on cross validation score with replications.
#' @usage tunelambda(data, time, ode.model, par.names, state.names,
#'                   par.initial, basis.list,lambda_grid,cv_portion,kfolds, rep,controls)
#' @param data A data frame or matrix contrain observations from each dimension of the ODE model.
#' @param         time The vector contain observation times or a matrix if time points are different between dimensions.
#' @param    ode.model Defined R function that computes the time derivative of the ODE model given observations of states variable.
#' @param    par.names The names of structural parameters defined in the 'ode.model'.
#' @param  state.names The names of state variables defined in the 'ode.model'.
#' @param  par.initial Initial value of structural parameters to be optimized.
#' @param   basis.list A list of basis objects for smoothing each dimension's observations. Can be the same or different across dimensions.
#' @param lambda_grid A search grid for finding the optimial sparsity parameter lambda.
#' @param cv_portion A number indicating the proportion of data will be saved for doing cross validation. Default is set at 5 as minimum.
#' @param kfolds A number indicating the number of folds the data should be seprated into.
#' @param rep A integer controls the number of replication of doing cross-validation for each penalty parameter.
#' @param     controls A list of control parameters. See ‘Details’.

#'
#' @return \item{lambda_grid}{The original input vector of a search grid for the optimal lambda.}
#' @return \item{cv.score}{The matrix contains the cross validation score for each lambda of each replication}

#' @export
tunelambda <- function(data, time, ode.model, par.names, state.names, par.initial, basis.list, lambda_grid, cv_portion,
    kfolds, rep, controls = list()) { # nocov start

    # Determine the folding of the original data
    boundary.points <- seq(min(time), max(time), length.out = kfolds + 1)
    # Determine the number of points to keep in each fold
    numcvpoints <- ceiling(length(time) * cv_portion/kfolds)


    cv.score <- matrix(NA, nrow = length(lambda_grid), ncol = rep)

    for (jj in 1:length(lambda_grid)) {

        for (kk in 1:rep) {
            # Deter
            points.keep <- matrix(NA, nrow = kfolds, ncol = numcvpoints)
            for (jkjk in 1:kfolds) {
                points.keep[jkjk, ] <- sample(time[time > boundary.points[jkjk] & time < boundary.points[jkjk + 1]],
                  numcvpoints)
            }

            # Keep some observations for cross validation: Identifying time points
            time.index.keep <- which(time %in% as.vector(points.keep))
            # Seperate data into two parts
            if (length(state.names) == 1) {
                data.keep <- data[time.index.keep]
                data.fit <- data[-time.index.keep]
            } else {
                data.keep <- data[time.index.keep, ]
                data.fit <- data[-time.index.keep, ]
            }

            message(paste("Running on lambda = ", lambda_grid[jj], " for iteration ", kk, sep = ""))
            pcode.result <- pcode(data = data.fit, time = time[-time.index.keep], ode.model = ode.model, par.names = par.names,
                state.names = state.names, basis.list = basis.list, par.initial = par.initial, lambda = lambda_grid[jj],
                controls = controls)
            #
            par.res <- pcode.result$structural.par
            names(par.res) <- par.names
            nui.res <- pcode.result$nuisance.par
            index <- rep(NA, length(state.names) + 1)
            index[1] <- 0
            # Evaluate
            initial.val <- rep(NA, length(state.names))
            names(initial.val) <- state.names
            for (jk in 1:length(state.names)) {
                # evaluation of basis object
                if (length(state.names) == 1) {
                  phi.temp <- eval.basis(time[1], basis.list)
                  index[jk + 1] <- basis.list$nbasis
                  basis.coef <- nui.res[(index[jk] + 1):(index[jk] + index[jk + 1])]
                  initial.val[jk] <- phi.temp %*% basis.coef
                } else {
                  phi.temp <- eval.basis(time[1], basis.list[[jk]])
                  index[jk + 1] <- basis.list[[jk]]$nbasis
                  basis.coef <- nui.res[(index[jk] + 1):(index[jk] + index[jk + 1])]
                  initial.val[jk] <- phi.temp %*% basis.coef
                }

            }

            # simulate a dynamic system based on initial value and structural parameters
            sim.res <- ode(y = initial.val, times = sort(time.index.keep), func = ode.model, parms = par.res)
            residuals <- sim.res[, state.names] - data.keep
            if (length(state.names) == 1) {
                cv.score[jj, kk] <- mean(residuals^2)
            } else {
                cv.score[jj, kk] <- sum(apply(residuals, 2, function(x) {
                  mean(x^2)
                }))
            }

        }
    }



    return(list(cv.score = cv.score, lambda_grid = lambda_grid))
} # nocov end

#' @title Inner objective function (likelihood and multiple dimension version)
#' @description An objective function combines the likelihood or loglikelihood of errors from each dimension of state variables and the penalty controls how the state estimates fail to satisfy the ODE model.
#' @usage innerobj_lkh(basis_coef, ode.par, input, derive.model, likelihood.fun)
#' @param basis_coef Basis coefficients for interpolating observations given a basis boject.
#' @param ode.par Structural parameters of the ODD model.
#' @param input Contains dependencies for the optimization, including observations, ode penalty, and etc..
#' @param derive.model The function defines the ODE model and is the same as the ode.model in 'pcode'.
#' @param likelihood.fun The likelihood or loglikelihood function of the errors.
#'
#' @return obj.eval The evaluation of the inner objective function.

innerobj_lkh <- function(basis_coef, ode.par, input, derive.model, likelihood.fun) {
    # Retrieve variables from input which is passed through the outter objective function with 'inner.input' get the
    # dimesion of the ODE model
    ndim <- length(input)
    # get the smooth parameter
    lambda <- input[[1]][[3]]
    # get the number of time points of observations
    npoints <- length(unlist(input[[1]][8]))
    # get the names of the states
    state.names <- input[[1]][[9]]
    # get the names of the parameters
    model.names <- input[[1]][[10]]


    # Preallocate some spaces index for partition the long vector of basis coefficient
    nbasis <- rep(NA, ndim + 1)
    # Predictions of the states with their derivatives
    Xt <- matrix(NA, nrow = length(input[[1]][[7]]), ncol = ndim)
    dXdt <- matrix(NA, nrow = length(input[[1]][[7]]), ncol = ndim)

    # Turn a single vector 'basis_coef' into seperate vectors for each dimension
    coef.list <- list()
    nbasis[1] <- 0
    for (jj in 1:ndim) {
        nbasis[jj + 1] <- ncol(input[[jj]][[2]])
        coef.list[[jj]] <- basis_coef[(nbasis[jj] + 1):(nbasis[jj] + nbasis[jj + 1])]
        Xt[, jj] <- input[[jj]][[4]] %*% coef.list[[jj]]
        dXdt[, jj] <- input[[jj]][[5]] %*% coef.list[[jj]]
    }

    # Evaluate the likelihood function over estimation
    eval.lik <- sum(apply(Xt, 2, likelihood.fun))


    # Start of calculation for the ODE penalty

    names(ode.par) <- model.names

    temp_fun <- function(x, temp = ode.par, names = state.names) {
        names(x) <- state.names
        return(unlist(derive.model(state = x, parameters = temp)))
    }
    temp_eval <- t(apply(Xt, 1, temp_fun))

    temp_list <- list(dXdt, temp_eval)
    temp_penalty_resid <- Reduce("-", temp_list)

    lambda <- input[[1]][[3]]
    penalty_residual <- matrix(NA, nrow = nrow(temp_eval), ncol = ncol(temp_eval))
    for (jj in 1:ndim) {
        penalty_residual[, jj] <- sqrt(lambda) * sqrt(input[[jj]][[7]]) * temp_penalty_resid[, jj]
    }

    obj.eval <- -eval.lik + sum(as.vector(penalty_residual)^2)

    return(obj.eval)
}

#' @title Outter objective function (likelihood and multiple dimension version)
#' @description An objective function of the structural parameter computes the measure of fit.
#' @usage outterobj_lkh(ode.parameter, basis.initial, derivative.model, likelihood.fun, inner.input)
#' @param ode.parameter Structural parameters of the ODE model.
#' @param basis.initial Initial values of the basis coefficients for nonlinear least square optimization.
#' @param derivative.model The function defines the ODE model and is the same as the ode.model in 'pcode'
#' @param likelihood.fun The likelihood or loglikelihood function of the errors.
#' @param inner.input Input that will be passed to the inner objective function. Contains dependencies for the optimization, including observations, penalty parameter lambda, and etc..
#'
#' @return neglik The negative of the likelihood or the loglikelihood function that will be passed further to the 'optim' function.
outterobj_lkh <- function(ode.parameter, basis.initial, derivative.model, likelihood.fun, inner.input) {

    # Profiled estimation on the nuisance parameters
    inner_coef <- optim(basis.initial, innerobj_lkh, ode.par = ode.parameter, derive.model = derivative.model, likelihood.fun = likelihood.fun,
        control = list(maxit = 50))

    # Get estimation of states to calculate outter objective function
    ndim <- length(inner.input)
    npoints <- length(inner.input[[1]][[8]])
    basisnumber <- rep(NA, ndim + 1)
    Xhat <- matrix(NA, nrow = npoints, ncol = ndim)

    # Turn a single vector 'inner_coef' into seperate locations
    coef.list <- list()
    basisnumber[1] <- 0
    for (jj in 1:ndim) {
        basisnumber[jj + 1] <- ncol(inner.input[[jj]][[2]])
        coef.list[[jj]] <- inner_coef[(basisnumber[jj] + 1):(basisnumber[jj] + basisnumber[jj + 1])]
        # Basis expansion for j-th dimension
        Xhat[, jj] <- inner.input[[jj]][[2]] %*% coef.list[[jj]]
    }
    # Evaluate the likelihood function over estimation
    eval.lik <- sum(apply(Xhat, 2, likelihood.fun))

    return(eval.lik)
}


#' @title pcode_lkh (likelihood and multiple dimension version)
#' @description Obtain estimates of both structural and nuisance parameters of an ODE model by parameter cascade method.
#' @usage pcode_lkh(data, likelihood.fun, time, ode.model, par.names,
#'                  state.names, par.initial, basis.list, lambda, controls)
#' @param data A data frame or a matrix contain observations from each dimension of the ODE model.
#' @param likelihood.fun A function computes the likelihood or the loglikelihood of the errors.
#' @param time A vector contains observation ties or a matrix if time points are different between dimesion.
#' @param    ode.model An R function that computes the time derivative of the ODE model given observations of states variable and structural parameters.
#' @param    par.names The names of structural parameters defined in the 'ode.model'.
#' @param  state.names The names of state variables defined in the 'ode.model'.
#' @param  par.initial Initial value of structural parameters to be optimized.
#' @param   basis.list A list of basis objects for smoothing each dimension's observations. Can be the same or different across dimensions.
#' @param    lambda Penalty parameter.
#' @param     controls A list of control parameters. See ‘Details’.
#'
#' @details The \code{controls} argument is a list providing addition inputs for the nonlinear least square optimizer:
#' \itemize{
#'\item \code{nquadpts} Determine the number of quadrature points for approximating an integral. Default is 101.
#'\item \code{smooth.lambda} Determine the smoothness penalty for obtaining initial value of nuisance parameters.
#'\item \code{tau} Initial value of Marquardt parameter. Small values indicate good initial values for structural parameters.
#'\item \code{tolx} Tolerance for parameters of objective functions. Default is set at 1e-6.
#'\item \code{tolg} Tolerance for the gradient of parameters of objective functions. Default is set at 1e-6.
#'\item \code{maxeval} The maximum number of evaluation of the optimizer. Default is set at 20.
#'}
#'
#' @return   \item{structural.par}{The structural parameters of the ODE model.}
#'
#' @return    \item{nuisance.par}{The nuisance parameters or the basis coefficients for interpolating observations.}



pcode_lkh <- function(data, likelihood.fun, time, ode.model, par.names, state.names, par.initial, basis.list, lambda = NULL,
    controls = NULL) {
    # Set up default controls for optimizations and quadrature evaluation
    con.default <- list(nquadpts = 101, smooth.lambda = 100, tau = 0.01, tolx = 1e-06, tolg = 1e-06, maxeval = 20)
    # Replace default with user's input
    con.default[(namc <- names(controls))] <- controls
    con.now <- con.default

    # Evaluate basis functiosn for each state variable

    basis.eval.list <- lapply(basis.list, prepare_basis, times = time, nquadpts = con.now$nquadpts)

    # Number of dimensions of the ODE model
    ndim <- ncol(data)

    # Some initializations
    inner.input <- list()
    initial_coef <- list()

    for (ii in 1:ndim) {
        # For each dimension, obtain initial value for the nuisance parameters or the basis coefficients
        Rmat = t(basis.eval.list[[ii]]$Q.D2mat) %*% (basis.eval.list[[ii]]$Q.D2mat * (basis.eval.list[[ii]]$quadwts %*%
            t(rep(1, basis.list[[ii]]$nbasis))))
        basismat2 = t(basis.eval.list[[ii]]$Phi.mat) %*% basis.eval.list[[ii]]$Phi.mat
        Bmat = basismat2 + con.now$smooth.lambda * Rmat
        # Initial basis coefficients
        initial_coef[[ii]] = ginv(Bmat) %*% t(basis.eval.list[[ii]]$Phi.mat) %*% data[, ii]
        inner.input[[ii]] = list(data[, ii], basis.eval.list[[ii]]$Phi.mat, lambda, basis.eval.list[[ii]]$Qmat, basis.eval.list[[ii]]$Q.D1mat,
            basis.eval.list[[ii]]$quadts, basis.eval.list[[ii]]$quadwts, time, state.names, par.names)
    }


    par.final <- optim(par.initial, outterobj_lkh, basis.initial = unlist(initial_coef), derivative.model = ode.model,
        input = inner.input, likelihood.fun = likelihood.fun, control = list(maxit = con.now$maxeval))$par

    basis.coef.final <- optim(unlist(initial_coef), innerobj_lkh, ode.par = par.final, input = inner.input, derive.model = ode.model,
        likelihood.fun = likelihood.fun)

    return(list(structural.par = par.final, nuisance.par = basis.coef.final))

}



#' @title Parameter Cascade Method for Ordinary Differential Equation Models (likelihood and Single dimension version)
#' @description Obtain estimates of both structural and nuisance parameters of an ODE model by parameter cascade method.
#' @usage pcode_lkh_1d(data, likelihood.fun, time, ode.model, par.names,
#'                     state.names, par.initial, basis.list, lambda, controls)
#' @param data A data frame or a matrix contain observations from each dimension of the ODE model.
#' @param likelihood.fun A function computes the likelihood or the loglikelihood of the errors.
#' @param time A vector contains observation ties or a matrix if time points are different between dimesion.
#' @param    ode.model An R function that computes the time derivative of the ODE model given observations of states variable and structural parameters.
#' @param    par.names The names of structural parameters defined in the 'ode.model'.
#' @param  state.names The names of state variables defined in the 'ode.model'.
#' @param  par.initial Initial value of structural parameters to be optimized.
#' @param   basis.list A list of basis objects for smoothing each dimension's observations. Can be the same or different across dimensions.
#' @param    lambda Penalty parameter.
#' @param     controls A list of control parameters. See ‘Details’.
#'
#' @details The \code{controls} argument is a list providing addition inputs for the nonlinear least square optimizer:
#' \itemize{
#'\item \code{nquadpts} Determine the number of quadrature points for approximating an integral. Default is 101.
#'\item \code{smooth.lambda} Determine the smoothness penalty for obtaining initial value of nuisance parameters.
#'\item \code{tau} Initial value of Marquardt parameter. Small values indicate good initial values for structural parameters.
#'\item \code{tolx} Tolerance for parameters of objective functions. Default is set at 1e-6.
#'\item \code{tolg} Tolerance for the gradient of parameters of objective functions. Default is set at 1e-6.
#'\item \code{maxeval} The maximum number of evaluation of the optimizer. Default is set at 20.
#'}
#'
#' @return   \item{structural.par}{The structural parameters of the ODE model.}
#'
#' @return    \item{nuisance.par}{The nuisance parameters or the basis coefficients for interpolating observations.}
pcode_lkh_1d <- function(data,likelihood.fun,time, ode.model, par.names, state.names, par.initial, basis.list,
    lambda, controls = list()) {
    # Set up default controls for optimizations and quadrature evaluation
    con.default <- list(nquadpts = 101, smooth.lambda = 100, tau = 0.01, tolx = 1e-06, tolg = 1e-06, maxeval = 20,verbal=0)
    # Replace default with user's input
    con.default[(namc <- names(controls))] <- controls
    con.now <- con.default
    basis <- basis.list
    # number of parameters
    npar <- length(par.initial)
    nbasis <- basis.list$nbasis

    # Evaluating basis functions at time points of observations and stored as columns in Phi matrix
    Phi.mat <- eval.basis(time, basis)
    # Evaluating 1st derivative of basis functions
    D1.mat <- eval.basis(time, basis, 1)
    # Evaluating 2nd derivative of basis functions
    D2.mat <- eval.basis(time, basis, 2)

    # Calculate L2 penalty
    quadts <- seq(min(time), max(time), length.out = con.now$nquadpts)
    nquad <- length(quadts)
    quadwts <- rep(1, nquad)
    even.ind <- seq(2, (nquad - 1), by = 2)
    odd.ind <- seq(3, (nquad - 2), by = 2)
    quadwts[even.ind] = 4
    quadwts[odd.ind] = 2
    h <- quadts[2] - quadts[1]
    quadwts <- quadwts * (h/3)

    Qmat <- eval.basis(quadts, basis)
    Q.D1mat <- eval.basis(quadts, basis, 1)
    Q.D2mat <- eval.basis(quadts, basis, 2)

    # Initial estimat of basis coefficients
    Rmat = t(Q.D2mat) %*% (Q.D2mat * (quadwts %*% t(rep(1, nbasis))))
    basismat2 = t(Phi.mat) %*% Phi.mat
    Bmat = basismat2 + con.now$smooth.lambda * Rmat
    # Initial basis coefficients
    initial_coef = ginv(Bmat) %*% t(Phi.mat) %*% data

    # Passing to inner objective functions and obtain initial value for parameter cascading
    inner.input = list(data, Phi.mat, lambda, Qmat, Q.D1mat, quadts, quadwts, time, state.names, par.names)

    theta.final <- optim(par.initial, outterobj_lkh_1d, basis.initial = initial_coef, derivative.model = ode.model,
        control = list(trace = controls$verbal,maxit = con.now$maxeval), inner.input = inner.input, likelihood.fun = likelihood.fun)$par

    basiscoef <- optim(initial_coef, innerobj_lkh_1d, ode.par = theta.final, input = inner.input, derive.model = ode.model,
        likelihood.fun = likelihood.fun)$par


    return(list(structural.par = theta.final, nuisance.par = basiscoef))
}


#' @title Outter objective function (likelihood and single dimension version)
#' @description An objective function of the structural parameter computes the measure of fit.
#' @usage outterobj_lkh_1d(ode.parameter, basis.initial,
#'                         derivative.model, likelihood.fun, inner.input)
#' @param ode.parameter Structural parameters of the ODE model.
#' @param basis.initial Initial values of the basis coefficients for nonlinear least square optimization.
#' @param derivative.model The function defines the ODE model and is the same as the ode.model in 'pcode'
#' @param likelihood.fun The likelihood or loglikelihood function of the errors.
#' @param inner.input Input that will be passed to the inner objective function. Contains dependencies for the optimization, including observations, penalty parameter lambda, and etc..
#'
#' @return neglik The negative of the likelihood or the loglikelihood function that will be passed further to the 'optim' function.

outterobj_lkh_1d <- function(ode.parameter, basis.initial, derivative.model, likelihood.fun, inner.input) {

    # Profiled estimation on the nuisance parameters
    basis_coef <- optim(basis.initial, innerobj_lkh_1d, ode.par = ode.parameter, derive.model = derivative.model, likelihood.fun = likelihood.fun,
        input = inner.input)$par

    yobs <- inner.input[[1]]
    Phi.mat <- inner.input[[2]]
    xhat <- Phi.mat %*% basis_coef
    residual <- yobs - xhat
    if (myCount(likelihood.fun) == 1) {
        lik.eval <- likelihood.fun(residual)
    } else {
        if (myCount(likelihood.fun) == 2) {
            lik.eval <- likelihood.fun(yobs, xhat)
        }
    }

    return(neglik = -lik.eval)
}

#' @title Inner objective function (Likelihood and Single dimension version)
#' @description An objective function combines the likelihood or loglikelihood of errors from each dimension of state variables and the penalty controls how the state estimates fail to satisfy the ODE model.
#' @usage innerobj_lkh_1d(basis_coef, ode.par, input, derive.model, likelihood.fun)
#' @param basis_coef Basis coefficients for interpolating observations given a basis boject.
#' @param ode.par Structural parameters of the ODD model.
#' @param input Contains dependencies for the optimization, including observations, ode penalty, and etc..
#' @param derive.model The function defines the ODE model and is the same as the ode.model in 'pcode'.
#' @param likelihood.fun The likelihood or loglikelihood function of the errors.
#'
#' @return obj.eval The evaluation of the inner objective function.
innerobj_lkh_1d <- function(basis_coef, ode.par, input, derive.model, likelihood.fun) {

    yobs <- input[[1]]
    Phi.mat <- input[[2]]
    lambda <- input[[3]]
    Qmat <- input[[4]]
    Q.D1mat <- input[[5]]
    quadts <- input[[6]]
    quadwts <- input[[7]]
    state.names <- input[[9]]
    model.names <- input[[10]]
    # Part 1.
    xhat <- Phi.mat %*% basis_coef
    residual <- yobs - xhat
    if (myCount(likelihood.fun) == 1) {
        lik.eval <- likelihood.fun(residual)
    } else {
        if (myCount(likelihood.fun) == 2) {
            lik.eval <- likelihood.fun(yobs, xhat)
        }
    }

    # Part 2.
    penalty_residual <- rep(NA, length(quadwts))
    # Use composite Simpson's rule to calculate the integral of residual
    Xt <- Qmat %*% basis_coef
    dXdt <- Q.D1mat %*% basis_coef

    # Evaluate f(X)
    rightside <- rep(NA, length(quadwts))
    names(ode.par) <- model.names


    tempfun <- function(x, temp = ode.par, names = state.names) {
        names(x) <- state.names
        return(derive.model(state = x, parameters = temp))
    }
    rightside <- unlist(lapply(Xt, tempfun))

    penalty_residual <- sqrt(lambda) * sqrt(quadwts) * (dXdt - rightside)
    # ode_penalty <- sum(lambda * quadwts * (penalty_residual)^2)
    obj.eval <- sum(penalty_residual^2) - lik.eval
    return(obj.eval)
}

myCount <- function(...) {
  length(match.call())
}



#' @title Bootstrap variance estimator of structural parameters.
#' @description Obtaining an estimate of variance for structural parameters by bootstrap method.
#' @usage bootsvar(data, time, ode.model,par.names,state.names, likelihood.fun = NULL,
#'        par.initial, basis.list, lambda = NULL,bootsrep,controls = NULL)
#' @param        data A data frame or a matrix contain observations from each dimension of the ODE model.
#' @param         time A vector contain observation times or a matrix if time points are different between dimensions.
#' @param    ode.model An R function that computes the time derivative of the ODE model given observations of states variable and structural parameters.
#' @param    par.names The names of structural parameters defined in the 'ode.model'.
#' @param  state.names The names of state variables defined in the 'ode.model'.
#' @param  par.initial Initial value of structural parameters to be optimized.
#' @param likelihood.fun A likelihood function passed to PCODE in case of that the error termsdevtools::document()do not have a Normal distribution.
#' @param   basis.list A list of basis objects for smoothing each dimension's observations. Can be the same or different across dimensions.
#' @param    lambda Penalty parameter.
#' @param bootsrep Bootstrap sample to be used for estimating variance.
#' @param     controls A list of control parameters. Same as the controls in \code{pcode}.
#'
#' @return boots.var The bootstrap variance of each structural parameters.
#' @export
bootsvar <- function(data, time, ode.model, par.names, state.names, likelihood.fun = NULL, par.initial, basis.list,
    lambda = NULL, bootsrep, controls = NULL) {
    if (length(state.names) >= 2) {
        if (nrow(data) > ncol(data)) {
            colnames(data) <- state.names
        } else {
            rownames(data) <- state.names
        }
    }


    # Initial run
    result.ini <- pcode(data, time, ode.model, par.names, state.names, likelihood.fun, par.initial, basis.list, lambda = lambda,
        controls = controls)


    # 1D case for least square functions
    if (length(state.names) == 1 && length(par.names) == 1) {
        nuipar.ini <- result.ini$nuisance.par
        phi.ini <- eval.basis(time, basis.list)
        state.est <- phi.ini %*% nuipar.ini
        residual <- data - state.est
        var.est <- var(residual)

        names(state.est) <- state.names
        temp.initial <- result.ini$structural.par
        names(result.ini$structural.par) <- par.names
        # preallocate spae for structural parameters
        boots.res <- matrix(NA, nrow = bootsrep, ncol = length(par.names))
        # modify the arguments of the ode model
        tempmodel <- function(t, state, parameters) {
            return(ode.model(state = state, parameters = parameters))
        }
        for (iter in 1:bootsrep) {
            message(paste("Running on bootstrap iteration: ", iter, sep = ""))
            # data.boot <- state.est + rnorm(length(state.est),mean = 0 , sd = sqrt(var.est))
            data.boot <- ode(y = state.est[1], times = time, func = tempmodel, parms = result.ini$structural.par)[,
                -1] + rnorm(length(state.est), mean = 0, sd = sqrt(var.est))
            result <- pcode(data = data.boot, time = time, ode.model = ode.model, par.names = par.names, state.names = state.names,
                par.initial = temp.initial, basis.list = basis.list, lambda = lambda, controls = controls)
            boots.res[iter, ] <- result$structural.par
        }
        boots.var <- apply(boots.res, 2, var)
        names(boots.var) <- par.names

    }

    # MD case for least square

    if (length(state.names) >= 2 && length(par.names) >= 2) {
        # Get basis coefficients
        nuipar.ini <- result.ini$nuisance.par
        # allocate vector for storing state and variance estimate
        state.est <- matrix(NA, ncol = length(state.names), nrow = length(time))
        colnames(state.est) <- state.names
        var.est <- rep(NA, length(state.names))


        # Sort into each dimension
        basis.index <- length(state.names) + 1
        basis.index[1] <- 0
        coef.list <- list()

        for (jj in 1:length(state.names)) {
            basis.index[jj + 1] <- basis.list[[jj]]$nbasis
            coef.list[[jj]] <- nuipar.ini[(basis.index[jj] + 1):(basis.index[jj] + basis.index[jj + 1])]
            phi.ini <- eval.basis(time, basis.list[[jj]])
            state.est[, jj] <- phi.ini %*% coef.list[[jj]]
            residual <- (data[, state.names[jj]] - state.est[, state.names[jj]])
            var.est[jj] <- var(residual)
        }

        temp.initial <- result.ini$structural.par
        names(result.ini$structural.par) <- par.names
        # preallocate spae for structural parameters
        boots.res <- matrix(NA, nrow = bootsrep, ncol = length(par.names))
        # modify the arguments of the ode model
        tempmodel <- function(t, state, parameters) {
            return(ode.model(state = state, parameters = parameters))
        }
        base.est <- ode(y = state.est[1, ], times = time, func = tempmodel, parms = result.ini$structural.par)[, -1]

        for (iter in 1:bootsrep) {
            data.boot <- matrix(NA,nrow = length(time),ncol = length(state.names))
            print(paste("Running on bootstrap iteration: ", iter, sep = ""))
            # data.boot <- state.est + rnorm(length(state.est),mean = 0 , sd = sqrt(var.est))
            for (jj in 1:length(state.names)){
              data.boot[,jj] <- base.est[,jj] + rnorm(length(time), mean = 0, sd = sqrt(var.est[jj]))
            }
            result <- pcode(data = data.boot, time = time, ode.model = ode.model, par.names = par.names, state.names = state.names,
                par.initial = temp.initial, basis.list = basis.list, lambda = lambda, controls = controls)
            boots.res[iter, ] <- result$structural.par
        }
        boots.var <- apply(boots.res, 2, var)
        names(boots.var) <- par.names

    }

    return(boots.var)
}



#' @title Numeric estimation of variance of structural parameters by Delta method.
#' @description Obtaining variance of structural parameters by Delta method.
#' @usage  deltavar(data, time, ode.model,par.names,state.names,
#'         likelihood.fun, par.initial, basis.list, lambda,stepsize,y_stepsize,controls)
#' @param        data A data frame or a matrix contain observations from each dimension of the ODE model.
#' @param         time A vector contain observation times or a matrix if time points are different between dimensions.
#' @param    ode.model An R function that computes the time derivative of the ODE model given observations of states variable and structural parameters.
#' @param    par.names The names of structural parameters defined in the 'ode.model'.
#' @param  state.names The names of state variables defined in the 'ode.model'.
#' @param  par.initial Initial value of structural parameters to be optimized.
#' @param likelihood.fun A likelihood function passed to PCODE in case of that the error termsdevtools::document()do not have a Normal distribution.
#' @param   basis.list A list of basis objects for smoothing each dimension's observations. Can be the same or different across dimensions.
#' @param    lambda Penalty parameter.
#' @param stepsize Stepsize used in estimating partial derivatives with respect to structural parameters for the Delta method.
#' @param y_stepsize Stepsize used in estimating partial derivatives with respect to observations for the Delta method.
#' @param     controls A list of control parameters. Same as the controls in \code{pcode}.
#'
#' @return par.var The variance of structural parameters obtained by Delta method.
#' @export
deltavar <- function(data, time, ode.model, par.names, state.names, likelihood.fun = NULL, par.initial, basis.list,
    lambda = NULL, stepsize, y_stepsize, controls = NULL) {
    if (length(state.names) >= 2) {
        if (nrow(data) > ncol(data)) {
            colnames(data) <- state.names
        } else {
            rownames(data) <- state.names
        }
    }


    # Set up default controls for optimizations and quadrature evaluation
    con.default <- list(nquadpts = 101, smooth.lambda = 100, tau = 0.01, tolx = 1e-06, tolg = 1e-06, maxeval = 20,verbal=0)
    # Replace default with user's input
    con.default[(namc <- names(controls))] <- controls
    con.now <- con.default



    # Initial run
    result <- pcode(data, time, ode.model, par.names, state.names, likelihood.fun, par.initial, basis.list, lambda = lambda,
        controls = con.now)

    struc.res <- result$structural.par
    names(struc.res) <- par.names
    nuis.res <- result$nuisance.par
    if (length(state.names) == 1) {
        basis.eval.list <- prepare_basis(basis.list, times = time, nquadpts = con.now$nquadpts)
        inner.input <- list(data, basis.eval.list$Phi.mat, lambda, basis.eval.list$Qmat, basis.eval.list$Q.D1mat, basis.eval.list$quadts,
            basis.eval.list$quadwts, time, state.names, par.names)
        upper.eval <- sum(outterobj(struc.res + stepsize, basis.initial = nuis.res, derivative.model = ode.model, inner.input = inner.input,
            NLS = TRUE)^2)
        center.eval <- sum((inner.input[[2]] %*% nuis.res - data)^2)
        lower.eval <- sum(outterobj(struc.res - stepsize, basis.initial = nuis.res, derivative.model = ode.model, inner.input = inner.input,
            NLS = TRUE)^2)

        d_H2_theta2 <- (upper.eval - 2 * center.eval + lower.eval)/(stepsize^2)

        #-------------------------------
        H_deriv_wrt_y <- rep(NA, length(time))
        inner_coef_upper <- nls_optimize.inner(innerobj, nuis.res, ode.par = struc.res + stepsize, derive.model = ode.model,
            input = inner.input, options = list(maxeval = 50))$par
        inner_coef_lower <- nls_optimize.inner(innerobj, nuis.res, ode.par = struc.res - stepsize, derive.model = ode.model,
            input = inner.input, options = list(maxeval = 50))$par
        obs_at_upper <- inner.input[[2]] %*% inner_coef_upper
        obs_at_lower <- inner.input[[2]] %*% inner_coef_lower
        for (index in 1:length(time)) {
            a = (data[index] + y_stepsize - obs_at_upper[index])^2
            b = (data[index] - y_stepsize - obs_at_upper[index])^2
            c = (data[index] + y_stepsize - obs_at_lower[index])^2
            d = (data[index] - y_stepsize - obs_at_lower[index])^2
            H_deriv_wrt_y[index] <- (a - b - c + d)/(4 * stepsize * y_stepsize)
        }

        if (length(par.names) == 1) {
            theta_deriv_wrt_y <- -(1/d_H2_theta2) * H_deriv_wrt_y
            par.var <- sum(theta_deriv_wrt_y^2) * var(data)
        } else {

        }


        # Multiple dimension and parameters
    } else {
        if (length(stepsize) == 1) {
            stepsize <- rep(stepsize, length(par.names))
        }
        names(stepsize) <- par.names
        # evaluation of basis objects for each dimension over time points and quadrature points for further use
        basis.eval.list <- lapply(basis.list, prepare_basis, times = time, nquadpts = con.now$nquadpts)
        # Sort nuisance parameters into each dimension and prepare inner input for inner objective function
        basis.index <- length(state.names) + 1
        basis.index[1] <- 0
        coef.list <- list()
        inner.input <- list()

        state.est <- matrix(NA, nrow = length(time), ncol = length(state.names))
        state.residual <- matrix(NA, nrow = length(time), ncol = length(state.names))
        colnames(state.est) <- state.names

        for (hh in 1:length(state.names)) {
            inner.input[[hh]] <- list(data[, state.names[hh]], basis.eval.list[[hh]]$Phi.mat, lambda, basis.eval.list[[hh]]$Qmat,
                basis.eval.list[[hh]]$Q.D1mat, basis.eval.list[[hh]]$quadts, basis.eval.list[[hh]]$quadwts, time, state.names,
                par.names)
            basis.index[hh + 1] <- basis.list[[hh]]$nbasis
            coef.list[[hh]] <- nuis.res[(basis.index[hh] + 1):(basis.index[hh] + basis.index[hh + 1])]
            state.est[, state.names[hh]] <- inner.input[[hh]][[2]] %*% coef.list[[hh]]
            state.residual[, hh] <- state.est[, state.names[hh]] - data[, state.names[hh]]
        }



        # preallocate space for evaluations
        d_H2_theta2 <- matrix(NA, nrow = length(par.names), ncol = length(par.names))
        upper.eval <- rep(NA, length(par.names))
        center.eval <- rep(NA, length(par.names))
        lower.eval <- rep(NA, length(par.names))


        for (par.index in 1:length(par.names)) {
            # Adding or subtrating stepsize
            par.moveup <- struc.res
            par.movedown <- struc.res
            par.moveup[par.names[par.index]] <- struc.res[par.names[par.index]] + stepsize[par.names[par.index]]
            par.movedown[par.names[par.index]] <- struc.res[par.names[par.index]] - stepsize[par.names[par.index]]

            upper.eval[par.index] <- outterobj_multi_nls(ode.parameter = par.moveup, basis.initial = nuis.res, derivative.model = ode.model,
                inner.input = inner.input, NLS = FALSE)
            center.eval[par.index] <- sum(as.vector(state.residual)^2)
            lower.eval[par.index] <- outterobj_multi_nls(ode.parameter = par.movedown, basis.initial = nuis.res, derivative.model = ode.model,
                inner.input = inner.input, NLS = FALSE)

            d_H2_theta2[par.index, par.index] <- (upper.eval[par.index] - 2 * center.eval[par.index] + lower.eval[par.index])/(stepsize[par.index]^2)
        }


        for (par.index in 1:(length(par.names) - 1)) {
            for (ii in par.index:(length(par.names) - 1)) {
                par.a <- par.b <- par.c <- par.d <- struc.res
                par.a[c(par.index, (ii + 1))] <- struc.res[c(par.index, (ii + 1))] + stepsize[c(par.index, (ii + 1))]
                par.b[par.index] <- struc.res[par.index] + stepsize[par.index]
                par.b[ii + 1] <- struc.res[ii + 1] - stepsize[ii + 1]
                par.c[par.index] <- struc.res[par.index] - stepsize[par.index]
                par.c[ii + 1] <- struc.res[ii + 1] + stepsize[ii + 1]
                par.d[c(par.index, (ii + 1))] <- struc.res[c(par.index, (ii + 1))] - stepsize[c(par.index, (ii + 1))]
                eval.a = outterobj_multi_nls(ode.parameter = par.a, basis.initial = nuis.res, derivative.model = ode.model,
                  inner.input = inner.input, NLS = FALSE)
                eval.b = outterobj_multi_nls(ode.parameter = par.b, basis.initial = nuis.res, derivative.model = ode.model,
                  inner.input = inner.input, NLS = FALSE)
                eval.c = outterobj_multi_nls(ode.parameter = par.c, basis.initial = nuis.res, derivative.model = ode.model,
                  inner.input = inner.input, NLS = FALSE)
                eval.d = outterobj_multi_nls(ode.parameter = par.d, basis.initial = nuis.res, derivative.model = ode.model,
                  inner.input = inner.input, NLS = FALSE)
                d_H2_theta2[par.index, ii + 1] <- (eval.a - eval.b - eval.c + eval.d)/(4 * stepsize[par.index] * stepsize[ii +
                  1])
            }
        }
        d_H2_theta2[lower.tri(d_H2_theta2)] <- t(d_H2_theta2)[lower.tri(d_H2_theta2)]
        #---------

        # second term (MD)
        H_deriv_wrt_y <- matrix(NA, nrow = length(time), ncol = length(state.names))

        inner.coef.upper <- list()
        inner.coef.lower <- list()

        for (par.index in 1:length(par.names)) {
            par.moveup <- struc.res
            par.movedown <- struc.res
            par.moveup[par.names[par.index]] <- struc.res[par.names[par.index]] + stepsize[par.names[par.index]]
            par.movedown[par.names[par.index]] <- struc.res[par.names[par.index]] - stepsize[par.names[par.index]]

            inner.coef.upper[[par.index]] <- nls_optimize.inner(innerobj_multi, nuis.res, ode.par = par.moveup, derive.model = ode.model,
                options = list(maxeval = 50, tolx = 1e-06, tolg = 1e-06), input = inner.input)$par
            inner.coef.lower[[par.index]] <- nls_optimize.inner(innerobj_multi, nuis.res, ode.par = par.movedown, derive.model = ode.model,
                options = list(maxeval = 50, tolx = 1e-06, tolg = 1e-06), input = inner.input)$par
        }


        coef.upper <- list()
        coef.lower <- list()

        basis.index <- length(state.names) + 1
        basis.index[1] <- 0

        for (ii in 1:length(par.names)) {
            coef.upper[[ii]] <- list()
            coef.lower[[ii]] <- list()
            for (hh in 1:length(state.names)) {
                basis.index[hh + 1] <- basis.list[[hh]]$nbasis
                coef.upper[[ii]][[state.names[hh]]] <- inner.coef.upper[[ii]][(basis.index[hh] + 1):(basis.index[hh] +
                  basis.index[hh + 1])]
                coef.lower[[ii]][[state.names[hh]]] <- inner.coef.lower[[ii]][(basis.index[hh] + 1):(basis.index[hh] +
                  basis.index[hh + 1])]
            }
        }


        obs_at_upper <- list()
        for (ii in 1:length(par.names)) {
            obs_at_upper[[ii]] <- matrix(NA, nrow = length(time), ncol = length(state.names))
            colnames(obs_at_upper[[ii]]) <- state.names
            for (jj in 1:length(state.names)) {
                obs_at_upper[[ii]][, state.names[jj]] <- inner.input[[jj]][[2]] %*% coef.upper[[ii]][[state.names[jj]]]
            }
        }
        obs_at_lower <- list()
        for (ii in 1:length(par.names)) {
            obs_at_lower[[ii]] <- matrix(NA, nrow = length(time), ncol = length(state.names))
            colnames(obs_at_lower[[ii]]) <- state.names
            for (jj in 1:length(state.names)) {
                obs_at_lower[[ii]][, state.names[jj]] <- inner.input[[jj]][[2]] %*% coef.lower[[ii]][[state.names[jj]]]
            }
        }

        H_deriv_wrt_y <- list()

        for (state.index in 1:length(state.names)) {
            H_deriv_wrt_y[[state.index]] <- matrix(NA, nrow = length(time), ncol = length(par.names))
            for (par.index in 1:length(par.names)) {

                for (time.index in 1:length(time)) {
                  a <- (data[time.index, state.names[state.index]] + y_stepsize - obs_at_upper[[par.index]][time.index,
                    state.names[state.index]])^2
                  b <- (data[time.index, state.names[state.index]] - y_stepsize - obs_at_upper[[par.index]][time.index,
                    state.names[state.index]])^2
                  c <- (data[time.index, state.names[state.index]] + y_stepsize - obs_at_lower[[par.index]][time.index,
                    state.names[state.index]])^2
                  d <- (data[time.index, state.names[state.index]] - y_stepsize - obs_at_lower[[par.index]][time.index,
                    state.names[state.index]])^2
                  H_deriv_wrt_y[[state.index]][time.index, par.index] <- (a - b - c + d)/(4 * stepsize[par.index] *
                    y_stepsize)

                }

            }
        }
        #

        second_H <- do.call(rbind, H_deriv_wrt_y)  #
        dtheta_dy <- t(-solve(d_H2_theta2) %*% t(second_H))

        var.vec <- rep(NA, length(state.names))
        for (jj in 1:length(state.names)) {
            var.vec[jj] <- var(state.residual[, jj])
        }
        var.rep <- rep(var.vec, each = length(time))
        Sigma.mat <- diag(var.rep)

        par.var <- t(dtheta_dy) %*% Sigma.mat %*% dtheta_dy

        rownames(par.var) <- colnames(par.var) <- par.names

    }

    return(par.var)
}

#' @title Parameter Cascade Method for Ordinary Differential Equation Models with missing state variable
#' @description Obtain estiamtes of both structural and nuisance parameters of an ODE model by parameter cascade method when the dynamics are partially observed.
#' @usage pcode_missing(data, time, ode.model, par.names, state.names,
#'                      likelihood.fun,par.initial, basis.list,lambda,controls)
#' @param        data A data frame or a matrix contain observations from each dimension of the ODE model.
#' @param         time A vector contain observation times or a matrix if time points are different between dimensions.
#' @param    ode.model An R function that computes the time derivative of the ODE model given observations of states variable and structural parameters.
#' @param    par.names The names of structural parameters defined in the 'ode.model'.
#' @param  state.names The names of state variables defined in the 'ode.model'.
#' @param  par.initial Initial value of structural parameters to be optimized.
#' @param likelihood.fun A likelihood function passed to PCODE in case of that the error termsdevtools::document()do not have a Normal distribution.
#' @param   basis.list A list of basis objects for smoothing each dimension's observations. Can be the same or different across dimensions.
#' @param    lambda Penalty parameter.
#' @param     controls A list of control parameters. See Details.
#'
#' @details The \code{controls} argument is a list providing addition inputs for the nonlinear least square optimizer or general optimizer \code{optim}:
#' \itemize{
#'\item \code{nquadpts} Determine the number of quadrature points for approximating an integral. Default is 101.
#'\item \code{smooth.lambda} Determine the smoothness penalty for obtaining initial value of nuisance parameters.
#'\item \code{tau} Initial value of Marquardt parameter. Small values indicate good initial values for structural parameters.
#'\item \code{tolx} Tolerance for parameters of objective functions. Default is set at 1e-6.
#'\item \code{tolg} Tolerance for the gradient of parameters of objective functions. Default is set at 1e-6.
#'\item \code{maxeval} The maximum number of evaluation of the optimizer. Default is set at 20.
#'}
#'
#'
#' @return   \item{structural.par}{The structural parameters of the ODE model.}
#'
#' @return    \item{nuisance.par}{The nuisance parameters or the basis coefficients for interpolating observations.}
pcode_missing <- function(data,time, ode.model,par.names, state.names, likelihood.fun,
                           par.initial, basis.list, lambda, controls){ # nocov start

  # Evaluate basis functiosn for each state variable
  basis.eval.list <- lapply(basis.list, prepare_basis, times = time, nquadpts = controls$nquadpts)

  # Number of dimensions of the ODE model
  ndim <- ncol(data)

  # Some initializations
  inner.input <- list()
  initial_coef <- list()

  #indicate which state variables are observed and stored in data matrix as columns
  observ.index <- which(state.names %in% colnames(data))

  for (ii in 1:length(observ.index)) {
    index <- observ.index[ii]
    # For each dimension, obtain initial value for the nuisance parameters or the basis coefficients
    Rmat = t(basis.eval.list[[index]]$Q.D2mat) %*% (basis.eval.list[[index]]$Q.D2mat * (basis.eval.list[[index]]$quadwts %*%
                                                                                          t(rep(1, basis.list[[index]]$nbasis))))
    basismat2 = t(basis.eval.list[[index]]$Phi.mat) %*% basis.eval.list[[index]]$Phi.mat
    Bmat = basismat2 + controls$smooth.lambda * Rmat
    # Initial basis coefficients
    initial_coef[[index]] = ginv(Bmat) %*% t(basis.eval.list[[index]]$Phi.mat) %*% data[, state.names[index]]
    # // TODO: need to check colnames of data
    inner.input[[index]] = list(data[, state.names[index]], basis.eval.list[[index]]$Phi.mat, lambda, basis.eval.list[[index]]$Qmat,
                                basis.eval.list[[index]]$Q.D1mat, basis.eval.list[[index]]$quadts, basis.eval.list[[index]]$quadwts, time, state.names,
                                par.names)
  }

  noobserv.index <- (1:length(state.names))[-observ.index]
  for (ii in 1:length(noobserv.index)){
    index <- noobserv.index[ii]
    # For each dimension, obtain initial value for the nuisance parameters or the basis coefficients
    Rmat = t(basis.eval.list[[index]]$Q.D2mat) %*% (basis.eval.list[[index]]$Q.D2mat * (basis.eval.list[[index]]$quadwts %*%
                                                                                          t(rep(1, basis.list[[index]]$nbasis))))
    basismat2 = t(basis.eval.list[[index]]$Phi.mat) %*% basis.eval.list[[index]]$Phi.mat
    Bmat = basismat2 + controls$smooth.lambda * Rmat
    # Initial basis coefficients
    initial_coef[[index]] = rnorm(ncol(Rmat))
    # // TODO: need to check colnames of data
    inner.input[[index]] = list(rep(NA,length(time)), basis.eval.list[[index]]$Phi.mat, lambda, basis.eval.list[[index]]$Qmat,
                                basis.eval.list[[index]]$Q.D1mat, basis.eval.list[[index]]$quadts, basis.eval.list[[index]]$quadwts, time, state.names,
                                par.names)
  }

  par.final <- nls_optimize(options = list(maxeval = controls$maxeval, tau = controls$tau), outterobj_multi_missing,
                            par.initial, basis.initial = unlist(initial_coef), derivative.model = ode.model, inner.input = inner.input,
                            verbal = 1)$par
  # Condition on the obtained structural parameter, calculate the nuisance parameter or the basis coefficients to
  # interpolate data
  basis.coef.final <- nls_optimize.inner(innerobj_multi_missing, unlist(initial_coef), ode.par = par.final, derive.model = ode.model,
                                         input = inner.input, NLS = TRUE, options = list(tau = 0.01))$par


  return(list(structural.par = par.final, nuisance.par = basis.coef.final))

} # nocov end


#' @title Outter objective function (multiple dimension version with unobserved state variables)
#' @description An objective function of the structural parameter computes the measure of fit for the basis expansion.
#' @usage outterobj_multi_missing(ode.parameter, basis.initial, derivative.model, inner.input, NLS)
#' @param ode.parameter Structural parameters of the ODE model.
#' @param basis.initial Initial values of the basis coefficients for nonlinear least square optimization.
#' @param derivative.model The function defines the ODE model and is the same as the ode.model in 'pcode'
#' @param inner.input Input that will be passed to the inner objective function. Contains dependencies for the optimization, including observations, penalty parameter lambda, and etc..
#' @param NLS Default is \code{TRUE} so the function returns vector of residuals, and otherwise returns sum of squared errors.
#'
#' @return   \item{residual}{Vector of residuals and evaluation of penalty function on quadrature points for approximating the integral.}
#'
outterobj_multi_missing <- function(ode.parameter, basis.initial, derivative.model, inner.input, NLS = TRUE) {# nocov start
  # Convergence of basis coefficients seems to happen before 'maxeval'.

  inner_coef <- nls_optimize.inner(innerobj_multi_missing, basis.initial, ode.par = ode.parameter, derive.model = derivative.model,
                                   options = list(maxeval = 100, tolx = 1e-06, tolg = 1e-06), input = inner.input)$par
  ndim <- length(inner.input)
  npoints <- length(inner.input[[1]][[8]])
  basisnumber <- rep(NA, ndim + 1)
  Xhat <- matrix(NA, nrow = npoints, ncol = ndim)
  residual <- matrix(NA, nrow = npoints, ncol = ndim)
  # Turn a single vector 'inner_coef' into seperate locations
  coef.list <- list()
  basisnumber[1] <- 0
  for (jj in 1:ndim) {
    basisnumber[jj + 1] <- ncol(inner.input[[jj]][[2]])
    coef.list[[jj]] <- inner_coef[(basisnumber[jj] + 1):(basisnumber[jj] + basisnumber[jj + 1])]
    # Basis expansion for j-th dimension
    Xhat[, jj] <- inner.input[[jj]][[2]] %*% coef.list[[jj]]
    # Residual from basis expansion
    residual[, jj] <- inner.input[[jj]][[1]] - Xhat[, jj]
  }

  if (NLS) {
    return(as.vector(residual)[!is.na(as.vector(residual))])
  } else {
    return(sum(residual^2))
  }

}# nocov end

#' @title Inner objective function (multiple dimension version with unobserved state variables)
#' @description An objective function combines the sum of squared error of basis expansion estimates and the penalty controls how those estimates fail to satisfies the ODE model
#' @usage innerobj_multi_missing(basis_coef, ode.par, input, derive.model,NLS)
#' @param basis_coef Basis coefficients for interpolating observations given a basis object.
#' @param ode.par Structural parameters of the ODE model.
#' @param input Contains dependencies for the optimization, including observations, penalty parameter lambda, and etc..
#' @param derive.model The function defines the ODE model and is the same as the ode.model in 'pcode'
#' @param NLS Default is \code{TRUE} so the function returns vector of residuals, and otherwise returns sum of squared errors.
#'
#' @return   \item{residual.vec}{Vector of residuals and evaluation of penalty function on quadrature points for approximating the integral.}
innerobj_multi_missing <- function(basis_coef, ode.par, input, derive.model, NLS = TRUE) { # nocov start

  # Retrieve variables from 'input' get the dimesion of the ODE model
  ndim <- length(input)
  npoints <- length(unlist(input[[1]][8]))
  nbasis <- rep(NA, ndim + 1)
  Xhat <- matrix(NA, nrow = npoints, ncol = ndim)
  residual <- matrix(NA, nrow = npoints, ncol = ndim)

  state.names <- input[[1]][[9]]
  model.names <- input[[1]][[10]]
  Xt <- matrix(NA, nrow = length(input[[1]][[7]]), ncol = ndim)
  dXdt <- matrix(NA, nrow = length(input[[1]][[7]]), ncol = ndim)

  # Turn a single vector 'basis_coef' into seperate locations
  coef.list <- list()
  nbasis[1] <- 0
  for (jj in 1:ndim) {
    nbasis[jj + 1] <- ncol(input[[jj]][[2]])
    coef.list[[jj]] <- basis_coef[(nbasis[jj] + 1):(nbasis[jj] + nbasis[jj + 1])]
    # Basis expansion for j-th dimension
    Xhat[, jj] <- input[[jj]][[2]] %*% coef.list[[jj]]
    # Residual from basis expansion
    residual[, jj] <- input[[jj]][[1]] - Xhat[, jj]
    Xt[, jj] <- input[[jj]][[4]] %*% coef.list[[jj]]
    dXdt[, jj] <- input[[jj]][[5]] %*% coef.list[[jj]]
  }

  # Ode penalty

  names(ode.par) <- model.names

  temp_fun <- function(x, temp = ode.par, names = state.names) {
    names(x) <- state.names
    return(unlist(derive.model(state = x, parameters = temp)))
  }
  temp_eval <- t(apply(Xt, 1, temp_fun))

  temp_list <- list(dXdt, temp_eval)
  temp_penalty_resid <- Reduce("-", temp_list)

  lambda <- input[[1]][[3]]
  penalty_residual <- matrix(NA, nrow = nrow(temp_eval), ncol = ncol(temp_eval))
  for (jj in 1:ndim) {
    penalty_residual[, jj] <- sqrt(lambda) * sqrt(input[[jj]][[7]]) * temp_penalty_resid[, jj]
  }
  #Ignore na residuals
  residual.vec <- c(as.vector(residual)[!is.na(as.vector(residual))], as.vector(penalty_residual))
  if (NLS) {
    return(residual.vec)
  } else {
    return(sum(residual.vec^2))
  }

} # nocov end
