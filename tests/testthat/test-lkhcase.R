context('pcode_lkh')
test_that("pcode works in parameter estimation (likelihood case)",{
  set.seed(123)
  #Define model parameters
  model.par   <- c(theta = c(0.1))
  #Define initial value of state variable
  state       <- c(X     = 0.1)
  #Define ODE model

  ode.model <- function(t, state,parameters){
    with(as.list(c(state,parameters)),
         {
           dX <- theta*X*(1-X/10)
           return(list(dX))
         })
  }
  #Observation points
  times <- seq(0,100,length.out=101)
  #Obain numeric solution of the ode model
  desolve.mod <- ode(y=state,times=times,func=ode.model,parms = model.par)
  #Generate noisay observations
  nobs  <- length(times)
  noise <- rnorm(n = nobs, mean = 0 , sd = 0.1)
  observ <- exp(log(desolve.mod[,2]) + noise)

  #parameter estimation
  #Define basis for doing interpolation
  #knots location
  knots <- seq(0,100,length.out=21)
  #order of basis functions
  norder <- 4
  #number of basis funtions
  nbasis <- length(knots) + norder - 2
  #creating basis
  basis  <- create.bspline.basis(c(0,100),nbasis,norder,breaks = knots)


  #Define the likelihood function for evaluating the fitness of the model
  likfun <- function(y,x){
    residuals <- log(y)-log(x)
    res <- lapply(residuals,function(t){ dnorm(t, mean= 0, sd= 0.1,log=TRUE)})
    return(sum(unlist(res)))
  }


  lkh.result <- pcode(data = observ,time = times, likelihood.fun = likfun, par.initial = 0.3,
                      ode.model = ode.model, par.names = 'theta',state.names ='X',basis.list = basis, lambda = 1e2,controls=list(verbal = 0))

  expect_equal( (abs(unname(lkh.result$structural.par) - unname(model.par))/abs(unname(model.par))) <= 1e-1, TRUE)
})
