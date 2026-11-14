# Main DGP from dgp.R, Section 4.2; not the override in vq_sim_qhete.R.
# normal/ar/smooth use base R; garch requires fGarch.
# Retain unused N normal draws to preserve original random-number stream.

dgp_hete<-function(beta, N, type, delta, kstar){
  # under H1
  y=c(rep(NA,N))
  v = matrix(rnorm(N), N, 1)
  sigma = c(rep(NA,N))

  kN <- floor(kstar * N) 

  e = c(rep(NA,N))
  
  omegas = c(3, 0.5)

  # length of sub-intervals
  size1 <- round(N * 0.45)
  size2 <- length( (round(N*0.45)+1):N ) 
  switch(type,
         normal = {
            e[1:round(N*0.45)] = rnorm( size1 ,0,omegas[1])
            e[ (round(N*0.45)+1) : N] = rnorm( size2, 0,omegas[2])
         },
         ar = {
            e[1:round(N*0.45)] = arima.sim(n = size1, list(order=c(1,0,0), ar=.3), rand.gen = rnorm, sd = omegas[1] )
            e[(round(N*0.45)+1):N] = arima.sim(n = size2, list(order=c(1,0,0), ar=.3), rand.gen = rnorm, sd = omegas[2] )
         },
         garch = {
            e[1:round(N*0.45)] = as.vector(fGarch::garchSim(fGarch::garchSpec(model = list(omega = omegas[1], alpha = 0.01, beta = 0.8)), n = size1 ))
            e[(round(N*0.45)+1):N] = as.vector(fGarch::garchSim(fGarch::garchSpec(model = list(omega = omegas[2], alpha = 0.01, beta = 0.8)), n = size2 ))
         },
         stop("Enter something to switch me!")
         )

  x = matrix(1,N,2)
  x[,2] = c( rnorm(size1,0,3), rnorm(size2,0,0.5) )

  y[1:kN] = x[1:kN,]%*%beta+e[1:kN]
  y[(kN+1):N] = x[(kN+1):N,]%*%(beta+delta )+e[(kN+1):N]

  return(list(x = x,y = y,e = e))
}


dgp_hetero5<-function(beta, N, delta, thetastar){
  # smooth change - log growth smooth function  
  # under H1
  y=c(rep(NA,N))
  v = matrix(rnorm(N), N, 1)
  sigma = c(rep(NA,N))
  kN = floor(thetastar*N)

  omegas = 1

  iid = rnorm( N, mean = 0, sd = omegas )
  smooth = 1+ 2/(1+exp(-10*( (1:N)/N-0.5 ) ) )
  e = smooth*iid

  size1 <- round(N * 0.45)
  size2 <- length( (round(N*0.45)+1):N ) 
  
  x=matrix(1,N,2)
  x[,2] = c( rnorm(size1,0,3), rnorm(size2,0,0.5) )

  y[1:kN] = x[1:kN,]%*%beta+e[1:kN]
  y[(kN+1):N] = x[(kN+1):N,]%*%(beta+delta)+e[(kN+1):N]

  return(list(x = x,y = y,e = e))
}



# Convenience interface: scalar delta shifts both intercept and slope.
dgp_main <- function(N=250, delta=0, kstar=0.5,
                     type=c("normal", "ar", "garch", "smooth"), beta=c(1,1)) {
  type <- match.arg(type)
  stopifnot(length(N)==1, is.finite(N), N >= 30, N==as.integer(N),
            length(beta)==2, all(is.finite(beta)),
            length(delta)==1, is.finite(delta),
            length(kstar)==1, kstar>0, kstar<1,
            floor(kstar*N)>=1, floor(kstar*N)<N)
  if (type == "garch" && !requireNamespace("fGarch", quietly=TRUE))
    stop("Install fGarch to use GARCH errors.")
  ans <- if (type == "smooth") dgp_hetero5(beta,N,delta,kstar) else
    dgp_hete(beta,N,type,delta,kstar)
  ans$change_after <- floor(N*kstar)
  ans
}
