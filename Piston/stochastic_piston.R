source("Piston/piston_setup.R")
library(quantkriging)
library(qrsvm) #install.packges("qrsvm/", type="source", repos=NULL)

set.seed(11149127)

# SIMULATE DATA
n <- 500
X <- if(n > 2000){randomLHS(n, 5)}else{maximinLHS(n, 5)}
X <- rbind(X, X)
y <- apply(X, 1, stochastic_piston, scale01=TRUE)
Q <- 0.75

# Get true quantiles using Monte Carlo
N_mc <- 10000
qtrue <- rep(NA, n)
for(i in 1:n){
  XX <- matrix(rep(X[i,], each=N_mc), nrow=N_mc)
  yy <- apply(XX, 1, stochastic_piston)
  qtrue[i] <- quantile(yy, Q)
}


# Fit quantile regression models
mod1 <- qbass(X, y, q=Q, maxInt=3,
              w_prior=list(type="GIG",
                           p=-.001, a=0, b=0.001,
                           prop_sigma=0.15),
              a_lambda=.03, b_lambda=.03)

mod2 <- quantKrig(X, y, Q, lower=c(1e-4, 1e-4), upper=c(10, 10),
                  rs=FALSE, nm=FALSE)

# For SVM do 5-fold CV to choose cost parameter
cost_vec <- 10^seq(-2, 1, by=0.25)
nfolds <- 5
best <- Inf
foo <- rep(NA, length(cost_vec))
for(i in seq_along(cost_vec)){
  qhat <- NULL
  for(k in 1:nfolds){
    holdout <- (1 + (k-1)*n/nfolds):(k*n/nfolds)
    mod3_curr <- qrsvm(X[-holdout,], y[-holdout],
                       cost=cost_vec[i], tau=Q)
    qhat <- c(qhat,
              qrsvm:::predict.qrsvm(mod3_curr, X[holdout,]))
  }
  curr <- rmse(qtrue, qhat)
  if(curr < best){
    istar <- i
    best <- curr
  }
  foo[i] <- curr
  print(i)
}
mod3 <- qrsvm(X, y, cost=cost_vec[istar], tau=Q)


# GENERATE TESTING DATA
n2 <- 500
X2 <- if(n2 > 2000){randomLHS(n2, 5)}else{maximinLHS(n2, 5)}
y2 <- apply(X, 1, stochastic_piston, scale01=TRUE)
Q <- 0.75

# Get true quantiles using Monte Carlo
N_mc <- 10000
qtrue2 <- rep(NA, n2)
for(i in 1:n2){
  XX <- matrix(rep(X[i,], each=N_mc), nrow=N_mc)
  yy <- apply(XX, 1, stochastic_piston)
  qtrue2[i] <- quantile(yy, Q)
}

qhat1 <- apply(predict(mod1, X2), 2, mean)
qhat2 <- as.numeric(predict(mod2, X2))
qhat3 <- qrsvm:::predict.qrsvm(mod3, X2)

rmse(qtrue2, qhat1)
rmse(qtrue2, qhat2)
rmse(qtrue2, qhat3)
