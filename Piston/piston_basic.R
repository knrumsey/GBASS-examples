source("Piston/piston_setup.R")

df <- 20
lg <- -0.792
tmp <- rskt(1000000, df, exp(lg))
#hist(tmp, breaks=50)
mu <- mean(tmp)
sd <- sd(tmp)
sk <- skew(tmp)
ku <- kurt(tmp)

# SIMULATE DATA
SNR <- 3
n <- 2000
X <- if(n > 2000){randomLHS(n, 7)}else{maximinLHS(n, 7)}
y <- apply(X, 1, piston, scale01=TRUE)
sdy <- sd(y)
#y <- y + (rskt(n, df, exp(lg)) - mu)/sd*sdy/sqrt(SNR)
y <- y + ralap(n, 2)*sdy/sqrt(SNR)
alpha <- 0.05

mod1 <- nwbass2(X, y, m_gamma=90, s_gamma=25, scale=10, m_beta=0, s_beta=10, lag_beta=50)
mod2 <- bass(X, y)

# Make plot
par(mfrow=c(1,1))
png("Piston/figs/error_distribution.png", height=5, width=5, units="in", res=300)
epsilon_esti <- nwbass_error(mod1)
epsilon_true <- (rskt(1000000, df, exp(lg)) - mu)/sd*sdy/sqrt(SNR)
epsilon_true <- ralap(100000, 2)*sdy/sqrt(SNR)
plot(density(epsilon_true), lwd=2,
     xlim=c(-0.4, 0.3), ylim=c(0, 9.5), main="", xlab="Error")
lines(density(epsilon_esti - mean(epsilon_esti)), col=bob[1], lwd=2, lty=2)
curve(dnorm(x, 0, mean(sqrt(mod2$s2))), add=TRUE, col=bob[4], lwd=2, lty=3)
legend("topleft", c("True Distribution", "Normal Wald BMARS", "Gaussian BMARS"),
       lty=1:3, lwd=2, col=c("black", bob[1], bob[4]), cex=0.8, bty='y')
dev.off()


