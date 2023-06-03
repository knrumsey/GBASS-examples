library(GBASS)
library(BASS)
library(hetGP)
library(BART)
library(laGP)
library(skewt)
library(lhs)
#library(duqling)
library(tictoc)
library(stargazer)

piston <- function(x, scale01=FALSE){
  if(scale01){
    RR <- cbind(c(30, 0.005, 0.002, 1000, 90000, 290, 340),
                c(60, 0.020, 0.010, 5000, 110000, 296, 360))
    x[1:7] <- x[1:7]*(RR[,2] - RR[,1]) + RR[,1]
  }

  M  <- x[1]
  S  <- x[2]
  V0 <- x[3]
  k  <- x[4]
  P0 <- x[5]
  Ta <- x[6]
  T0 <- x[7]

  Aterm1 <- P0 * S
  Aterm2 <- 19.62 * M
  Aterm3 <- -k*V0 / S
  A <- Aterm1 + Aterm2 + Aterm3

  Vfact1 <- S / (2*k)
  Vfact2 <- sqrt(A^2 + 4*k*(P0*V0/T0)*Ta)
  V <- Vfact1 * (Vfact2 - A)

  fact1 <- M
  fact2 <- k + (S^2)*(P0*V0/T0)*(Ta/(V^2))

  C <- 2 * pi * sqrt(fact1/fact2)
  return(C)
}


skew <- function(x){
  zz <- (x - mean(x))/sd(x)
  mean(zz^3)
}
kurt <- function(x){
  zz <- (x - mean(x))/sd(x)
  mean(zz^4)
}
rmse <- function(x, y){
  sqrt(mean((x-y)^2))
}
nwbass_error_mean <- function(mod1, nn=100000){
  ww <- mean(mod1$w)
  cc <- mod1$scale
  beta <- mean(mod1$beta)
  xi <- rnorm(nn, 0, sqrt(cc))
  gamma <- mean(mod1$gamma)
  vv <- rep(NA, nn)
  for(i in 1:nn){
    vv[i] <- rgig2(-1/2, gamma^2, 1)
  }
  epsilon <- sqrt(ww)*(beta*vv + sqrt(vv)*xi)
}

nwbass_error <- function(mod1, nn=100000){
  epsilon <- rep(NA, nn)
  ww <- mod1$w
  cc <- mod1$scale
  beta <- mod1$beta
  gamma <- mod1$gamma
  xi <- rnorm(nn, 0, sqrt(cc))
  mm <- length(ww)
  for(i in 1:nn){
    ii <- sample(mm, 1)
    vv <- rgig2(-1/2, gamma[ii], 1)
    epsilon[i] <- sqrt(ww[ii])*(beta[ii]*vv + sqrt(vv)*xi[ii])
  }
  return(epsilon)
}

nwbass_error_i <- function(mod1, ii=1, nn=10){
  epsilon <- rep(NA, nn)
  ww <- mod1$w
  cc <- mod1$scale
  beta <- mod1$beta
  gamma <- mod1$gamma
  xi <- rnorm(nn, 0, sqrt(cc))
  mm <- length(ww)
  for(i in 1:nn){
    vv <- rgig2(-1/2, gamma[ii], 1)
    epsilon[i] <- sqrt(ww[ii])*(beta[ii]*vv + sqrt(vv)*xi[i])
  }
  return(epsilon)
}

theta <- c(0.02, 11, 84)
stochastic_piston <- function(xx, scale01=TRUE){
  x <- c(xx[1:4], 5*(rbeta(1, 10, 2)-0.88), rbeta(1, 1.4, 1.5), xx[5])
  if(scale01){
    RR <- cbind(c(30, 0.005, 0.002, 1000, 90000, 290, 340),
                c(60, 0.020, 0.010, 5000, 110000, 296, 360))
    x[1:7] <- x[1:7]*(RR[,2] - RR[,1]) + RR[,1]
  }

  M  <- x[1]
  S  <- x[2]
  V0 <- x[3]
  k  <- x[4]
  P0 <- x[5]
  Ta <- x[6]
  T0 <- x[7]

  Aterm1 <- P0 * S
  Aterm2 <- 19.62 * M
  Aterm3 <- -k*V0 / S
  A <- Aterm1 + Aterm2 + Aterm3

  Vfact1 <- S / (2*k)
  Vfact2 <- sqrt(A^2 + 4*k*(P0*V0/T0)*Ta)
  V <- Vfact1 * (Vfact2 - A)

  fact1 <- M
  fact2 <- k + (S^2)*(P0*V0/T0)*(Ta/(V^2))

  C <- 2 * pi * sqrt(fact1/fact2)
  return(C)
}

myLHS <- function(n, k){
  if(n > 2000){
    d <- randomLHS(n, k)
  }else{
    d <- maximinLHS(n, k)
  }
  return(d)
}

# Generates draws from asymmetric laplace distribution
# with skewness kappa (mean 0 and variance 1)
ralap_mod <- function(n, kappa=1){
  lambda <- sqrt((1+kappa^4)/kappa^2)
  mu <- (1-kappa^2)/(lambda*kappa)
  #u <- runif(n, -kappa, 1/kappa)
  u <- rbeta(n, 1/0.8, 1/0.8)*(1/kappa + kappa) - kappa
  s <- sign(u)
  z <- s*kappa^s
  x <- -1/(lambda*z)*log(1 - u*z) - mu
}

ralap <- function(n, kappa=1){
  lambda <- sqrt((1+kappa^4)/kappa^2)
  mu <- (1-kappa^2)/(lambda*kappa)
  u <- runif(n, -kappa, 1/kappa)
  s <- sign(u)
  z <- s*kappa^s
  x <- -1/(lambda*z)*log(1 - u*z) - mu
}


hpd_quantile <- function(x, alpha){
  xx <- sort(x)
  n <- length(xx)
  k <- ceiling(n*(1-alpha))
  best <- Inf
  for(i in 1:(n-k+1)){
    curr <- x[i+k-1] - x[i]
    if(curr < best){
      best <- curr
      best_i <- i
    }
  }

  probs <- c(best_i, best_i+k-1)/n
  return(quantile(x, probs))
}
