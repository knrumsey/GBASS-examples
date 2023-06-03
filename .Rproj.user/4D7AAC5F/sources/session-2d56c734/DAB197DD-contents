source("Fish/fish_setup.R")
Q <- 0.9 #target quantile

data <- read.csv("https://raw.githubusercontent.com/jhuang672/fish/master/data/GridData_2.csv", header=TRUE)

X <- matrix(data$population_size, ncol=1)
xx <- matrix(seq(0, 1, length.out=20))
X01 <- BASS:::scale.range(X)
y <- data$recaptures
plot(X01+rnorm(X01, 0, 0.002), y, pch=16, col='gray', cex=0.5)

# Get "true" quantiles
qtrue <- rep(NA, length(unique(X01)))
for(i in 1:length(unique(X01))){
  indx <- which(X01 == unique(X01)[i])
  qtrue[i] <- quantile(y[indx], Q)
}
points(xx, qtrue, pch=21 , bg="dodgerblue")

Nreps_vec <- 2^(0:8)  # Number of replications
Nsims <- 20 # Number of simulations

# COMPARISON
res1 <- res2 <- res3 <- matrix(NA, nrow=Nsims, ncol=length(Nreps_vec))
for(k in 1:length(Nreps_vec)){
  cat("Starting Nreps = ", Nreps_vec[k], "\n")
  Nreps <- Nreps_vec[k]
  for(i in 1:Nsims){
    indx <- sample(500, Nreps)
    curr <- unlist(lapply(indx, function(ii) (1 + 20*(ii-1)):(20*ii)))
    Xi <- X01[curr]
    yi <- y[curr]

    # GBASS
    mod1 <- qbass(Xi, yi, q=Q, maxInt=1,
                         w_prior=list(type="GIG",
                                      p=-.001, a=0, b=0.001,
                                      prop_sigma=0.15),
                         a_lambda=.03, b_lambda=.03)
    yhat1 <- apply(predict(mod1, xx), 2, mean)

    # Quantile Kriging
    mod2 <- quantKrig(Xi, yi, Q, lower=c(1e-4, 1e-4), upper=c(10, 10),
                      rs=FALSE, nm=FALSE)
    yhat2 <- predict(mod2, xx)

    # QRSVM
    mod3 <- qrsvm(matrix(Xi), yi, tau=Q, cost=1000)
    yhat3 <- qrsvm:::predict.qrsvm(mod3, xx)

    res1[i,k] <- rmse(yhat1, qtrue)
    res2[i,k] <- rmse(yhat2, qtrue)
    res3[i,k] <- rmse(yhat3, qtrue)
  }
  save(res1, res2, res3, file="Fish/data/fish_sims.Rda")
}


