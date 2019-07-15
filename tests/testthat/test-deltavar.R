context('deltavar')
test_that('Delta variance estimation works',{
  set.seed(123)
  ode.model <- function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      dX <- theta*X*(1-X/10)
      return(list(dX))})}
  #define model parameters
  model.par   <- c(theta = c(0.1))
  #define state initial value
  state       <- c(X     = 0.1)
  times  <- seq(0,100,length.out=101)
  mod    <- ode(y=state,times=times,func=ode.model,parms = model.par)
  nobs   <- length(times)
  scale  <- 0.5
  noise  <- scale * rnorm(n = nobs, mean = 0 , sd = 1)
  observ <- mod[,2] + noise
  #Generate basis object for interpolation and as argument of pcode
  #21 konts equally spaced within [0,100]
  knots <- seq(0,100,length.out=21)
  #order of basis functions
  norder <- 4
  #number of basis funtions
  nbasis <- length(knots) + norder - 2
  #creating Bspline basis
  basis  <- create.bspline.basis(c(0,100),nbasis,norder,breaks = knots)

  delta_var = deltavar(data = observ, time = times, ode.model = ode.model,
           par.initial = 0.3, par.names = 'theta',state.names = 'X',
           basis.list = basis, lambda = 1e2,
           stepsize = 1e-5,y_stepsize = 1e-5)

expect_equal( (abs(delta_var - 1e-5)/1e-5) < 1,TRUE)
})
