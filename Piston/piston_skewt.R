source("Piston/piston_setup.R")

VARY <- 0.0198
SNR  <- 3
N    <- 1000
M    <- 1
Nt   <- 1000
nt   <- 30 # predictive draws per posterior sample.
alpha<- 1 - c(0.5, 0.8, 0.9, 0.95, 0.99)

set.seed(123109811)
df_vec <- c(8,   12,  20,     40,    100,   4)
lg_vec <- c(0.6,  0, -0.792, log(2), 0.25, -0.2)
sk_vec <- rep(NA, length(df_vec))
for(k in 1:length(df_vec)){
  df <- df_vec[k]
  lg <- lg_vec[k]
  tmp <- rskt(10000, df, exp(lg))
  mu <- mean(tmp)
  sd <- sd(tmp)
  sk <- skew(tmp)
  ku <- kurt(tmp)

  tmp2 <- (rskt(10000000, df, exp(lg)) - mu)/sd*sqrt(VARY)/sqrt(SNR)
  sk_vec[k] <- skew(tmp2)
}


# Models:
# mod1: GBASS
# mod2: BASS
# mod3: BART
# mod4: laGP
# mod5: TBASS
# mod6: QBASS
TAB_list <- list()
mod1_list <- mod2_list <- mod3_list <- mod4_list <- mod5_list <- mod6_list <- list()
for(k in 1:(length(df_vec)+1)){
  if(k <= length(df_vec)){
    # GET INFO FOR ERROR DISTRIBUTION
    df <- df_vec[k]
    lg <- lg_vec[k]
    tmp <- rskt(10000000, df, exp(lg))
    mu <- mean(tmp)
    sd <- sd(tmp)
    sk <- skew(tmp)
    ku <- kurt(tmp)

    eps <- (rskt(5*n, df, exp(lg)) - mu)/sd*sqrt(VARY)/sqrt(SNR) # Sample 5*n for better ecdf later on
  }else{
    df <- 10
    eps <- ralap(5*n, 2)*sqrt(VARY)/sqrt(SNR)
  }


  # GENERATE TRAINING DATA
  n <- N/M
  X <- myLHS(n, 7)
  if(M > 1){
    X0 <- X
    for(m in 1:(M-1)){
      X <- rbind(X, X0)
    }
  }
  y <- apply(X, 1, piston, scale01=TRUE)

  Feps <- ecdf(eps)
  y <- y + eps[1:n]

  # GENERATE TESTING DATA
  X2 <- myLHS(Nt, 7)
  y2 <- apply(X2, 1, piston, scale01=TRUE)
  #y3 <- y2 + (rskt(Nt, df, exp(lg)) - mu)/sd*sqrt(VARY)/sqrt(SNR)
  #y3 <- y2 + ralap(Nt, 2)*sqrt(VARY)/sqrt(SNR)
  y3 <- y2 + eps[(n+1):(n+Nt)]

  # 1. GBASS
  mod1 <- nwbass2(X, y, m_gamma=90, s_gamma=25,
                          m_beta=0, s_beta=10,
                          scale=1, lag_beta=50)
  pred1 <- predict(mod1, X2)
  yhat1 <- apply(pred1, 2, mean)
  rmse1a <- rmse(y2, yhat1 + mean(sqrt(mod1$w)*mod1$beta/mod1$gamma))
  rmse1b <- rmse(y3, yhat1 + mean(sqrt(mod1$w)*mod1$beta/mod1$gamma))
  intr1 <- matrix(NA, nrow=Nt, ncol=5)
  covr1 <- matrix(NA, nrow=Nt, ncol=5)
  KS1   <- rep(NA, Nt)
  for(i in 1:Nt){
    tmp <- NULL
    for(j in seq(1, length(mod1$w), by=1)){
      tmp <- c(tmp, pred1[j,i] + nwbass_error_i(mod1, j, nt))
    }
    yy <- y2[i] + eps
    F0 <- ecdf(yy)
    F1 <- ecdf(tmp)
    xx <- seq(min(yy), max(yy), length=1001)
    KS1[i] <- max(abs(F0(xx)-F1(xx)))
    bounds <- list()
    for(aa in seq_along(alpha)){
      #bounds[[aa]] <- hpd_quantile(tmp, alpha[aa])
      bounds[[aa]] <- quantile(tmp, c(alpha[aa]/2, 1- alpha[aa]/2))
      covr1[i,aa] <- as.numeric((bounds[[aa]][1] <= y3[i]) & (bounds[[aa]][2] >= y3[i]))
      intr1[i,aa] <- diff(bounds[[aa]])
    }
    if( i %% 10 == 0) print(i)
  }
  cvrg1 <- apply(covr1, 2, mean)
  diff1 <- apply(intr1, 2, median)

  # 2. BASS
  mod2 <- bass(X, y)

  pred2 <- predict(mod2, X2)
  yhat2 <- apply(pred2, 2, mean)
  rmse2a <- rmse(y2, yhat2)
  rmse2b <- rmse(y3, yhat2)
  intr2 <- matrix(NA, nrow=Nt, ncol=5)
  covr2 <- matrix(NA, nrow=Nt, ncol=5)
  KS2   <- rep(NA, Nt)
  for(i in 1:Nt){
    tmp <- NULL
    for(j in seq(1, length(mod2$s2), by=1)){
      tmp <- c(tmp, pred2[j,i] + rnorm(nt, 0, sqrt(mod2$s2[j])))
    }
    yy <- y2[i] + eps
    F0 <- ecdf(yy)
    F2 <- ecdf(tmp)
    xx <- seq(min(yy), max(yy), length=1001)
    KS2[i] <- max(abs(F0(xx)-F2(xx)))
    bounds <- list()
    for(aa in seq_along(alpha)){
      bounds[[aa]] <- quantile(tmp, c(alpha[aa]/2, 1-alpha[aa]/2))
      covr2[i,aa] <- as.numeric((bounds[[aa]][1] <= y3[i]) & (bounds[[aa]][2] >= y3[i]))
      intr2[i,aa] <- diff(bounds[[aa]])
    }
    if( i %% 10 == 0) print(i)
  }
  cvrg2 <- apply(covr2, 2, mean)
  diff2 <- apply(intr2, 2, median)

  par(mfrow=c(1,2))
  plot(1-alpha, cvrg1, pch=16, cex=2, xlim=c(0.4, 1), ylim=c(0.4, 1))
  points(1-alpha, cvrg2, pch=15, cex=1, col='dodgerblue')
  abline(0,1, lty=3, col='red')

  plot(1-alpha, diff1, pch=16, cex=2, xlim=c(0.4, 1), ylim=c(0.4, 1))
  points(1-alpha, diff2, pch=15, cex=2, col='dodgerblue')
  abline(0,1, lty=3, col='red')


  # 3. BART
  mod3 <- wbart(X, y)

  pred3 <- predict(mod3, X2)
  yhat3 <- apply(pred3, 2, mean)
  rmse3a <- rmse(y2, yhat3)
  rmse3b <- rmse(y3, yhat3)
  intr3 <- matrix(NA, nrow=Nt, ncol=5)
  covr3 <- matrix(NA, nrow=Nt, ncol=5)
  KS3   <- rep(NA, Nt)
  for(i in 1:Nt){
    tmp <- NULL
    for(j in seq(1, nrow(pred3), by=1)){
      tmp <- c(tmp, pred3[j,i] + rnorm(nt, 0,  mod3$sigma[100+j]))
    }
    yy <- y2[i] + eps
    F0 <- ecdf(yy)
    F3 <- ecdf(tmp)
    xx <- seq(min(yy), max(yy), length=1001)
    KS3[i] <- max(abs(F0(xx)-F3(xx)))
    bounds <- list()
    for(aa in seq_along(alpha)){
      bounds[[aa]] <- quantile(tmp, c(alpha[aa]/2, 1-alpha[aa]/2))
      covr3[i,aa] <- as.numeric((bounds[[aa]][1] <= y3[i]) & (bounds[[aa]][2] >= y3[i]))
      intr3[i,aa] <- diff(bounds[[aa]])
    }
    if( i %% 10 == 0) print(i)
  }
  cvrg3 <- apply(covr3, 2, mean)
  diff3 <- apply(intr3, 2, median)

  # 4. laGP
  yhat4 <- covr4 <- rep(NA, Nt)
  intr4 <- matrix(NA, nrow=Nt, ncol=5)
  covr4 <- matrix(NA, nrow=Nt, ncol=5)
  KS4   <- rep(NA, Nt)
  for(i in 1:Nt){
    fit4 <- laGP(t(X2[i,]), start=10, end=round(sqrt(N)),
                 X, y, g=garg(g=list(mle=TRUE), y=y3),
                 method="alc")
    yhat4[i] <- fit4$mean
    yy <- y2[i] + eps
    F0 <- ecdf(yy)
    F4 <- function(xx) pnorm(xx, fit4$mean, sqrt(fit4$s2))
    xx <- seq(min(yy), max(yy), length=1001)
    KS4[i] <- max(abs(F0(xx)-F4(xx)))
    bounds <- list()
    for(aa in seq_along(alpha)){
      bounds[[aa]] <- fit4$mean + c(-1, 1)*qt(1 - alpha[aa]/2, fit4$df)*sqrt(fit4$s2)
      covr4[i,aa] <- as.numeric((bounds[[aa]][1] <= y3[i]) & (bounds[[aa]][2] >= y3[i]))
      intr4[i,aa] <- diff(bounds[[aa]])
    }
    if( i %% 10 == 0) print(i)
  }
  cvrg4 <- apply(covr4, 2, mean)
  diff4 <- apply(intr4, 2, median)
  rmse4a <- rmse(y2, yhat4)
  rmse4b <- rmse(y3, yhat4)

  # 5. TBASS
  df0 <- df
  mod5 <- tbass(X, y, df=df0)
  pred5 <- predict(mod5, X2)
  yhat5 <- apply(pred5, 2, mean)
  rmse5a <- rmse(y2, yhat5)
  rmse5b <- rmse(y3, yhat5)
  intr5 <- matrix(NA, nrow=Nt, ncol=5)
  covr5 <- matrix(NA, nrow=Nt, ncol=5)
  KS5   <- rep(NA, Nt)
  for(i in 1:Nt){
    tmp <- NULL
    for(j in seq(1, length(mod5$w), by=1)){
      tmp <- c(tmp, pred5[j,i] + sqrt(mod5$w[j])*rt(nt, df=df0))
    }
    yy <- y2[i] + eps
    F0 <- ecdf(yy)
    F5 <- ecdf(tmp)
    xx <- seq(min(yy), max(yy), length=1001)
    KS5[i] <- max(abs(F0(xx)-F5(xx)))
    bounds <- list()
    for(aa in seq_along(alpha)){
      bounds[[aa]] <- quantile(tmp, c(alpha[aa]/2, 1-alpha[aa]/2))
      covr5[i,aa] <- as.numeric((bounds[[aa]][1] <= y3[i]) & (bounds[[aa]][2] >= y3[i]))
      intr5[i,aa] <- diff(bounds[[aa]])
    }
    if( i %% 10 == 0) print(i)
  }
  cvrg5 <- apply(covr5, 2, mean)
  diff5 <- apply(intr5, 2, median)

  # 6. QBASS(0.5)
  mod6 <- qbass(X, y, q=0.5)
  pred6 <- predict(mod6, X2)
  yhat6 <- apply(pred6, 2, mean)
  rmse6a <- rmse(y2, yhat6)
  rmse6b <- rmse(y3, yhat6)
  intr6 <- matrix(NA, nrow=Nt, ncol=5)
  covr6 <- matrix(NA, nrow=Nt, ncol=5)
  KS6   <- rep(NA, Nt)
  for(i in 1:Nt){
    tmp <- NULL
    for(j in seq(1, length(mod6$w), by=1)){
      tmp <- c(tmp, pred6[j,i] + sqrt(mod6$w[j])*ralap(nt, 1))
    }
    yy <- y2[i] + eps
    F0 <- ecdf(yy)
    F6 <- ecdf(tmp)
    xx <- seq(min(yy), max(yy), length=1001)
    KS6[i] <- max(abs(F0(xx)-F6(xx)))
    bounds <- list()
    for(aa in seq_along(alpha)){
      bounds[[aa]] <- quantile(tmp, c(alpha[aa]/2, 1-alpha[aa]/2))
      covr6[i,aa] <- as.numeric((bounds[[aa]][1] <= y3[i]) & (bounds[[aa]][2] >= y3[i]))
      intr6[i,aa] <- diff(bounds[[aa]])
    }
    if( i %% 10 == 0) print(i)
  }
  cvrg6 <- apply(covr6, 2, mean)
  diff6 <- apply(intr6, 2, median)

  # MAKE TABLE
  TAB <- matrix(NA, nrow=6, ncol=13)
  rownames(TAB) <- c("GBASS", "TBASS", "QBASS(0.5)", "BASS", "BART", "laGP")
  colnames(TAB) <- c("RMSE*", "RMSE", "KS-dist",
                     "Cov 0.5", "Cov 0.8", "Cov 0.9", "Cov 0.95", "Cov 0.99",
                     "Width 0.5", "Width 0.8", "Width 0.9", "Width 0.95", "Width 0.99")

  TAB[1,] <- c(rmse1a, rmse1b, median(KS1),  cvrg1, diff1)
  TAB[2,] <- c(rmse5a, rmse5b, median(KS5),  cvrg5, diff5)
  TAB[3,] <- c(rmse6a, rmse6b, median(KS6),  cvrg6, diff6)
  TAB[4,] <- c(rmse2a, rmse2b, median(KS2),  cvrg2, diff2)
  TAB[5,] <- c(rmse3a, rmse3b, median(KS3),  cvrg3, diff3)
  TAB[6,] <- c(rmse4a, rmse4b, median(KS4),  cvrg4, diff4)

  TAB_list[[k]] <- TAB

  mod1_list[[k]] <- mod1
  mod2_list[[k]] <- mod2
  mod3_list[[k]] <- mod3
  #mod4_list[[k]] <- mod4 # No model for laGP
  mod5_list[[k]] <- mod5
  mod6_list[[k]] <- mod6
}
save(TAB_list, mod1_list, mod2_list,
     mod3_list, mod4_list, mod5_list, mod6_list, file="Piston/data/piston_skewt_results_final3.rda")


bob <- RColorBrewer::brewer.pal(5, "Dark2")
nw_triangle(mod1_list[[1]], pch=3, col=bob[1])
for(k in 2:4){
  nw_triangle(mod1_list[[k]], add=TRUE, col=bob[k])
}


png("Piston/figs/coverage0.png", height=5, width=5, units="in", res=300)
par(mfrow=c(1,1))
bob <- RColorBrewer::brewer.pal(6, "Set2")
plot(NULL, xlim=c(0.5, 1), ylim=c(-0.105, 0.075),
     xlab="Desired Coverage", ylab="Difference in Nominal Coverage")
abline(h=0, lwd=2, lty=3)
for(i in c(1:6)[-3]){
  lines(1-alpha, TAB_list[[7]][i,4:8]-(1-alpha), lwd=2, col=bob[i])
  points(1-alpha, TAB_list[[7]][i,4:8]-(1-alpha), pch=i, col=bob[i], cex=1.5, lwd=2)
}
legend(0.65, -0.025, c("Normal-Wald (BMARS)", "t(10) (BMARS)", "Laplace (BMARS)", "Gaussian (BASS)", "BART", "laGP")[-3],
       lwd=2, col=bob[-3], pch=(1:6)[-3], cex=0.75)
dev.off()

# Make tables
for(i in 1:6){
  stargazer(TAB_list[[i]][-3,c(1,3:8)], digits=3, summary=FALSE)
}
tab1 <- TAB_list[[7]][,c(1,3:8)]
stargazer(tab1, digits=3, summary=FALSE)


tab2 <- matrix(NA, nrow=6, ncol=12)
rownames(tab2) <- rownames(TAB_list[[1]])
colnames(tab2) <- rep(c("RMSPE", "KS"), 6)
for(i in 1:6){
  tmp <- TAB_list[[i]]
  tab2[,2*(i-1) + 1:2] <- tmp[,c(1,3)]
}
tab2[,seq(1,11,by=2)] <- round(tab2[,seq(1,11,by=2)] , 3)
tab2[,seq(2,12,by=2)] <- round(tab2[,seq(2,12,by=2)] , 3)
stargazer(t(tab2), digits=NA, summary=F)


