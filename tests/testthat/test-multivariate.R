context('pcode_multivariate')
test_that('pcode pameter estimation works (multivariate)',{
  set.seed(123)
  #Define model parameters
  model.par   <- c(a=0.2,b=0.2,c=3)
  #Define initial value of state variables
  state       <- c(V=-1,R=1)
  #Define ode model for obtaining numeric solution
  ode.model <- function(t, state,parameters){
    with(as.list(c(state,parameters)),
         {
           dV <- c*(V - (V^3)/3 + R)
           dR <- -(1/c) * (V - a + b*R)
           return(list(c(dV,dR)))
         })}


  #Define observation points
  times <- seq(0,20,length.out=401)
  #Obtain ode solution
  desolve.mod <- ode(y=state,times=times,func=ode.model,parms = model.par)


  #Generate measurement noises
  nobs  <- length(times)
  scale <- 0.1
  noise_v <- scale*rnorm(n = nobs, mean = 0 , sd = 1)
  noise_r <- scale*rnorm(n = nobs, mean = 0 , sd = 1)
  #Store observations
  observ <- matrix(NA, nrow = length(times),ncol =3)
  observ[,1] <- times
  observ[,2] <- desolve.mod[,2] + noise_v
  observ[,3] <- desolve.mod[,3] + noise_r

  #Define basis for interpolating observation of both state variables
  #can use same basis for both dimensions
  knots <- seq(0,20,length.out=101)
  #order of basis functions
  norder <- 4
  #number of basis funtions
  nbasis <- length(knots) + norder - 2
  #creating basis
  basis  <- create.bspline.basis(c(0,20),nbasis,norder,breaks = knots)
  #Make a list of basis object for interpolation
  basis.list <- list(basis, basis)

  data <- observ[,2:3]
  colnames(data) <- c('V','R')

  #parameter estimation
  pcode.result <- pcode(data = observ[,2:3], time= observ[,1], ode.model = ode.model,
                        par.names = c('a','b','c'),state.names = c('V','R'), par.initial = c(0.1,0.3,4),
                        lambda = c(1e2,1e2),basis.list = basis.list,
                        controls = list(smooth.lambda = 10,verbal = 1,maxeval = 150))
  pcode.result$structural.par

  expect_equal( max(abs(pcode.result$structural.par - model.par) / model.par) < 1,TRUE)
})
