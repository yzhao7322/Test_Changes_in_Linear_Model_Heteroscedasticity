# Core extracted from vq_sim_qhete.R. Main setting: d=2, kappa=1/2.
# MnK = V_HET; RnK = original Q implementation (see README).
# Only cointReg is required. Sourcing this file runs no simulations.

ols_e <- function(dv,indv){
    dv = as.matrix(dv); indv = as.matrix(indv)
    beta_hat = solve(t(indv)%*%indv) %*% (t(indv)%*%dv)
    return(beta_hat)
}

## residuals
resid <- function(dv,indv){
    N = length(dv)
    beta_hat = ols_e(dv,indv)
    residual = c()
    for(i in 1:N){
       residual[i] = dv[i]-as.matrix(indv)[i,]%*%beta_hat
    }
    return(residual)
}

## weighted residuals
weight_resid <- function(xs,uind){
    xs = as.matrix(xs)
    N = dim(xs)[1]
    if (floor(uind*(N+1))<1){k = 1}else{k = floor(uind*(N+1))}
    zh = N^(-1/2)*(colSums(as.matrix(xs[1:k,]))-floor(uind*(N+1))/N*colSums(xs))
    return(zh)
}



mylongcov_G <-function(data, demean, kw){
  N = nrow(data)
  d = ncol(data)
  times <- c(seq(1/100,1,1/100))
  ts = 6
  if(demean == 1){
    dmd = as.matrix(apply(data,2,function(x){x-mean(x)}))
  }else{dmd = as.matrix(data)}

  longv = array(0, c(d, d, 100))
  bandwidth = "and"; kervar = "ba"
  nlag <- floor(cointReg::getBandwidth(dmd, bandwidth = bandwidth, kernel = kervar))
  #nlag = floor(log(N)) #floor(N^(1/5)) 
  
  # if (nlag >= ts-1){nlag = ts-2}

  if(nlag == 0){
    # sample covariance
    for (m in 4:100){
      tm = floor(times[m]*N)
      if(tm==1){tm=2} # in small sample, dmd[1:tm,] could be a scalar
      dm = as.matrix(dmd[1:tm,])
      # dm = as.matrix(apply(dm,2,function(x){x-mean(x)}))
      V = (t(dm)%*%dm)/ N
      longv[,,m] = V
    }

  }else{
    # long-run covariance
    # kernel weights
    w = 1 - c(0, seq(1/nlag, 1, 1/nlag))

    for (m in 4:100){
      tm = floor(times[m]*N)
      if(tm <= (nlag+1)){tm = nlag+2} # in small sample, make sure dmd[1:tm,] & dm[(nlag+1):tm,] is not a scalar

      dm = as.matrix(dmd[1:tm,]) 
      V = (t(dm)%*%dm)/ N

      for (i in 1:nlag){  
        gammai = (t(dm[(i+1):tm,])%*%dm[1:(tm-i),]) / (N-i)
        Gprime = gammai + t(gammai)

        V =  V + w[i+1] * Gprime
      }

      longv[,,m] = V
    }
    
  }
  return(longv)
}

############# covariance

Gbar_matrix <-function(y, x, wdat){
  N=length(y)
  times <- c(seq(1/100,1,1/100)) 
  # longrun covariance estimator
  de = ncol(as.matrix(wdat))

  Gbar = array(0, c(de, de, 100))
  GNt = mylongcov_G(wdat, 1, 'ba')

  for (t in 4:100){
    Nt = floor(N*times[t])
    uNt1 = matrix(rowSums(apply(x[1:Nt,],1,function(xi){return(xi%*%t(xi))} )),2,2)
    uNt2 = matrix(rowSums(apply(x[1:N,],1,function(xi){return(xi%*%t(xi))} )),2,2)
    uNt = (uNt1 - floor((N+1)*times[t])/N*uNt2)%*%solve(t(x)%*%x)
        
    firstline = GNt[,,t] - 2*times[t]*GNt[,,t] + times[t]^2*GNt[,,100]
    secondline =  - uNt%*%GNt[,,t] - GNt[,,t]%*%t(uNt) + (times[t]*uNt)%*%GNt[,,t]
    thirdline =  + times[t]*GNt[,,t]%*%t(uNt) + uNt%*%GNt[,,100]%*%t(uNt)
    Gbar[,,t] = firstline + secondline + thirdline
  }

  return(Gbar)
}


###################
######  statistics
###################
MnK<-function(y,x,kappa=1/2){
  N = length(y)
  d = ncol(x)

  res = resid(y,x)
  xs = cbind(x,res)
  wdat = as.matrix(apply(xs,1, function(xi){return(xi[1:(length(xi)-1)]*xi[length(xi)])} ) )
  if (ncol(wdat) > nrow(wdat)){wdat = t(wdat)}

  # longrun covariance estimator
  Gbar = Gbar_matrix(y, x, wdat)

  # statistics
  decums = c()
  for(i in 8:(N-1)){
    t = round(i/N*100)
    if(t<4){t=4}
    if(t==100){t=99}

    zh = weight_resid(wdat,i/(N+1))
    decums[i] = (zh%*%solve(Gbar[,,t])%*%as.matrix(zh))^(1/2) 
  }

  stat2 = max(decums[8:(N-d-1 - ceiling(N/100) )], na.rm=TRUE) # -N/100, i.e., t*N
    
  cp = which(decums==stat2)
  return(list(stat2, cp))
}


RnK<-function(y,x,kappa=1/2){
  N = length(y)
  d = ncol(x)

  res = resid(y,x)
  xs = cbind(x,res)
  wdat = as.matrix(apply(xs,1, function(xi){return(xi[1:(length(xi)-1)]*xi[length(xi)])} ) )
  if (ncol(wdat) > nrow(wdat)){wdat = t(wdat)}

  # longrun covariance estimator
  Gbar = Gbar_matrix(y, x, wdat)
  # statistics
  decums = matrix(NA,(N-1),d)  
  for(i in 8:(N-1)){
    t = round(i/N*100)
    if(t<4){t=4}
    if(t==100){t=99}

    zh = weight_resid(wdat,i/N)
    Gsq = (solve(Gbar[,,t]))^(1/2); Gsq[is.nan(Gsq)]=0
    decums[i,] = Gsq%*%zh
  }

  decums[is.na(decums)] <- 0
  mnorm = abs(decums[,which(max(colSums(abs(decums)))==colSums(abs(decums)))])
  wmnorm = mnorm
 
  stat2 = max(wmnorm[8:(N-d-1- ceiling(N/100) )], na.rm=TRUE) #
  cp = which(wmnorm==stat2)

  return(list(stat2,cp))
}


vos_f1 = function(x,d,N){
  1- x^d * exp(-x^2/2) / (2^(d/2)* gamma(d/2)  )* (log(N^2) - d/x^2 * log(N^2) + 4/x^2) #1- x^d * exp(-x^2/2) / (2^(d/2)* (pi/(2*d)))* (log(N^2) - d/x^2 * log(N^2) + 4/x^2)
  }
# curve(vos_f1(x,d=2,N=250), from = 2, to= 10)

vos_cv1 = function(d,N){
   cv=as.vector(c(NA,NA,NA))
   discv = seq(2,7,by=0.001)
   percent= vos_f1(discv,d,N)
   ind90 = which.min(abs(percent - 0.90))
   ind95 = which.min(abs(percent - 0.95))
   ind99 = which.min(abs(percent - 0.99))
   cv[1] = discv[ind90]
   cv[2] = discv[ind95]
   cv[3] = discv[ind99]
   return(cv)
}


# Q using symmetric matrix inverse square root and max over all coordinates.
# Keeps the original covariance estimator, grid, and endpoint trimming.
Q_matrix <- function(y, x) {
  N <- length(y); d <- ncol(x)
  wdat <- x * as.numeric(resid(y,x))
  G <- Gbar_matrix(y,x,wdat)
  index <- 8:(N-d-1-ceiling(N/100))
  path <- vapply(index, function(i) {
    m <- min(99,max(4,round(i/N*100)))
    eg <- eigen((G[,,m]+t(G[,,m]))/2, symmetric=TRUE)
    if (min(eg$values) <= max(eg$values)*1e-12)
      stop("Covariance is not positive definite; no silent regularization applied.")
    root <- eg$vectors %*% diag(1/sqrt(eg$values),d) %*% t(eg$vectors)
    z <- weight_resid(wdat,i/(N+1))
    max(abs(root %*% z))
  }, numeric(1))
  list(max(path),index[which.max(path)])
}

# Default reproduces the supplied code. Use q_method="matrix" for the
# paper's matrix-root/infinity-norm definition with the same numerical grid.
vq_test <- function(y, x, q_method=c("original", "matrix")) {
  q_method <- match.arg(q_method)
  x <- as.matrix(x); y <- as.numeric(y)
  stopifnot(nrow(x)==length(y), ncol(x)==2, length(y)>=30,
            all(is.finite(x)), all(is.finite(y)))
  if (qr(x)$rank < ncol(x)) stop("x must have full column rank.")
  if (!requireNamespace("cointReg", quietly=TRUE))
    stop("Install cointReg for Andrews automatic bandwidth selection.")
  N <- length(y); d <- ncol(x)
  # Check the original partial-covariance indexing before calling it.
  h <- floor(cointReg::getBandwidth(x*as.numeric(resid(y,x)),
                                    bandwidth="and",kernel="ba"))
  if (!is.finite(h) || h<0 || h>N-2) stop("Unsupported automatic bandwidth.")
  v <- MnK(y,x)
  q <- if (q_method=="original") suppressWarnings(RnK(y,x)) else Q_matrix(y,x)
  cv <- c(vos_cv1(d,N)[2], (qn(d)+b_d(log(N)))/a(log(N)))
  data.frame(test=c("V_HET", "Q_HET"),
             statistic=c(v[[1]],q[[1]]), critical_5pct=cv,
             reject=c(v[[1]],q[[1]])>cv,
             peak_index=c(v[[2]][1],q[[2]][1]),
             q_method=q_method, row.names=NULL)
}

# Darling--Erdos critical value for Q at alpha=0.05.
a <- function(x) sqrt(2*log(x))
b_d <- function(x,d=1) 2*log(x)+(d/2)*log(log(x))-log(gamma(d/2))
qn <- function(d) log(2*d/(-log(0.95)))
