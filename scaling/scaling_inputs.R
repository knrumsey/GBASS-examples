library(GBASS)
library(BASS)
library(lhs)
#library(duqling)
library(tictoc)
library(stargazer)

# This is copied over from the duqling R package
# github.com/knrumsey/duqling
micwicz <- function(x, scale01=TRUE, m=10, active_dim=length(x)){
  p <- min(active_dim, length(x))
  if(scale01){
    RR <- cbind(rep(0, p), rep(2*pi, p))
    x <- 2*(base::pi)*x
  }

  y <- 0
  for(i in 1:p){
    y <- y - sin(x[i])*sin(i*x[i]^2/(base::pi))^(2*m)
  }
  return(y)
}

dms_simple <- function(x, scale01=TRUE) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)

# Different p_vec we want to try
p_vec <- 2^seq(2, 10, by=2)
n     <- 5000

r1 <- r2 <- rep(NA, length(p_vec))
t1 <- t2 <- rep(NA, length(p_vec))
for(i in 1:length(p_vec)){
  p <- p_vec[i]

  X <- randomLHS(n, p)
  y1 <- apply(X, 1, micwicz, m=2)
  y2 <- apply(X, 1, dms_simple)

  # Add stochasticity (SNR = 6)
  y1 <- y1 + rnorm(n, 0, sqrt(var(y1)/6))
  y2 <- y2 + rnorm(n, 0, sqrt(var(y2)/6))

  tic()
  mod1 <- tbass(X, y1, df=5)
  tf1 <- toc()
  t1[i] <- tf1$toc - tf1$tic

  tic()
  mod2 <- tbass(X, y2, df=5)
  tf2 <- toc()
  t2[i] <- tf2$toc - tf2$tic

  # Check predictions
  Xt <- randomLHS(100, p)
  y1t <- apply(Xt, 1, micwicz, m=2)
  y2t <- apply(Xt, 1, dms_simple)

  y1hat <- colMeans(predict(mod1, Xt))
  y2hat <- colMeans(predict(mod2, Xt))

  r1[i] <- sqrt(mean((y1hat-y1t)^2))
  r2[i] <- sqrt(mean((y2hat-y2t)^2))
}

sd_vec <- rep(NA, length(p_vec))
for(i in seq_along(p_vec)){
  p <- p_vec[i]
  Xtmp <- matrix(runif(10000*p), ncol=p)
  ytmp <- apply(Xtmp, 1, micwicz, m=2)
  sd_vec[i] <- sd(ytmp)
}

fve1 <- 1 - (r1/sd_vec)^2
fve2 <- 1 - (r2/sd(y2))^2

png("scaling/figs/FVE.png", height=5, width=5, units="in", res=300)
plot(fve1, type='o', lwd=2, pch=15, lty=1, cex=1.2, col='orange',
     ylim=c(0, 1), xaxt='n', xlab='p', ylab='Fraction of Variance Explained')
lines(fve2, lwd=2, col='dodgerblue', lty=2)
points(fve2, cex=1.2, pch=16, col='dodgerblue')
axis(1, at=1:5, p_vec)
legend("bottomleft", c("Micwicz Function", "DMS Function"),
       col=c("orange", "dodgerblue"), pch=15:16, lwd=2, lty=1:2,
       cex=1, bty='n')
dev.off()

png("scaling/figs/timing.png", height=5, width=5, units="in", res=300)
plot(t1/60, type='o', lwd=2, pch=15, lty=1, cex=1.2, col='orange',
     xaxt='n', xlab='p', ylab='Time to Fit (m)')
lines(t2/60, lwd=2, col='dodgerblue', lty=2)
points(t2/60, cex=1.2, pch=16, col='dodgerblue')
axis(1, at=1:5, p_vec)
legend("topright", c("Micwicz Function", "DMS Function"),
       col=c("orange", "dodgerblue"), pch=15:16, lwd=2, lty=1:2,
       cex=1, bty='n')
axis(1, at=1:5, p_vec)
dev.off()
