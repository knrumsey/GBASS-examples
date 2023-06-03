source("Piston/piston_setup.R")

n_vec <- c(100, 500, 1000, 2000, 5000, 10000)
VARY <- 0.0198
SNR  <- 3

t1 <- t2 <- t3 <- t4 <- t5 <- rep(NA, length(n_vec))
rmse1 <- rmse2 <- rmse3 <- rmse4 <- rmse5 <- rep(NA, length(n_vec))
for(i in seq_along(n_vec)){
  # Generate training data
  n <- n_vec[i]
  X <- myLHS(n, 7)
  y <- apply(X, 1, piston) + rt(n, 12)*10/12*sqrt(VARY)/sqrt(SNR)

  # Generate testing data
  X2 <- myLHS(1000, 7)
  y2 <- apply(X2, 1, piston)

  # Fit and time models
  # GBASS
  tic()
  mod1 <- nwbass2(X, y, m_gamma=90, s_gamma=25,
                  m_beta=0, s_beta=10,
                  scale=1, lag_beta=50,
                  a_lambda=10, b_lambda=10)
  tt <- toc()
  t1[i] <- tt$toc - tt$tic
  yhat1 <- apply(predict(mod1, X2), 2, mean)
  rmse1[i] <- rmse(y2, yhat1 + mean(sqrt(mod1$w)*mod1$beta/mod1$gamma))

  # BASS
  tic()
  mod2 <- bass(X, y)
  tt <- toc()
  t2[i] <- tt$toc - tt$tic
  yhat2 <- apply(predict(mod2, X2), 2, mean)
  rmse2[i] <- rmse(y2, yhat2)

  # BART
  tic()
  mod3 <- wbart(X, y)
  tt <- toc()
  t3[i] <- tt$toc - tt$tic
  yhat3 <- apply(predict(mod3, X2), 2, mean)
  rmse3[i] <- rmse(y2, yhat3)

  # GP (Full)
  if(i <= 3){
    tic()
    mod4 <- mleHomGP(X, y)
    tt <- toc()
    t4[i] <- tt$toc - tt$tic
    yhat4 <- predict(mod4, X2)
    rmse4[i] <- rmse(y2, yhat4$mean)
  }

  # TBASS
  tic()
  mod5 <- tbass(X, y, 12, a_lambda=1, b_lambda=1)
  tt <- toc()
  t5[i] <- tt$toc - tt$tic
  yhat5 <- apply(predict(mod5, X2), 2, mean)
  rmse5[i] <- rmse(y2, yhat5)
}

save(t1, t2, t3, t4, t5,
     rmse1, rmse2, rmse3, rmse4, rmse5,
     file="Piston/data/scaling.Rda")

##Make plot
lm4 <- lm(t4[1:4]~+I(n_vec[1:4]^3))
plot(n_vec[1:3], t4)
t4[4:6] <- lm4$coefficients[1] + lm4$coefficients[2]*n_vec[4:6]^3

png("Piston/figs/scaling.png", height=5, width=8, units="in", res=300)
bob <- RColorBrewer::brewer.pal(6, "Dark2")
ord <- c(1, 4, 5, 6, 2)
ln_vec <- log10(n_vec)
plot(ln_vec, log10(t1), type='o',
     ylim=log10(range(c(t4, t1))),
     lwd=2, col=bob[ord[1]], pch=ord[1], cex=1.5,
     xaxt='n', yaxt='n',
     xlab="n", ylab="Time (s)")
for(i in 2:5){
  ti <- get(paste0("t", i))
  if(i == 4){
    lines(ln_vec[1:4], log10(ti[1:4]), col=bob[ord[i]], lwd=2)
    lines(ln_vec[4:6], log10(ti[4:6]), col=bob[ord[i]], lwd=2, lty=3)
  }else{
    lines(ln_vec, log10(ti), col=bob[ord[i]], lwd=2)
  }
  points(ln_vec, log10(ti), pch=ord[i], col=bob[ord[i]], cex=1.5, lwd=2)
}
axis(1, ln_vec, n_vec)
axis(2, 0:4, 10^(0:4))
legend("topleft", c("Normal-Wald (BMARS)", "t(12) (BMARS)" , "Gaussian (BASS)", "BART", "Full GP"),
       lwd=2, col=bob[c(1,2,4,5,6)], pch=c(1,2,4,5,6), cex=1)
dev.off()




