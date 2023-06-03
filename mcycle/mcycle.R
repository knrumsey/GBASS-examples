library(MASS)
library(BASS)
library(parallel)
library(RColorBrewer)
mcycle$x <- BASS:::scale.range(mcycle$times)
xx <- matrix(seq(0,1,by=.01), ncol=1)

mods <- mclapply(c(.1,.25,.5,.75,.9),
                 function(qq) qbass(matrix(mcycle$x),
                                    mcycle$accel,q=qq,
                                    maxInt=1,
                                    w_prior=list(type="GIG", p=-0.0001, a=0, b=0.0001, prop_sigma=0.2),
                                    a_lambda=.03, b_lambda=.03,
                                    nmcmc=10000, nburn=8001, thin=4),
                 mc.cores = 5, mc.preschedule = F)


png("mcycle/figs/mcycle.png", height=5, width=8, units="in", res=300)
bob.ross <- brewer.pal(6, "Blues")
plot(mcycle$x, mcycle$accel, pch=16, xlab="Time", ylab='Acceleration', cex.lab=1.25, ylim=c(-150, 90))
grid()
for(i in 1:5){
  yhat <- colMeans(predict(mods[[i]], xx))
  lines(xx, yhat, lwd=2, col=bob.ross[i+1])
}
points(mcycle$x, mcycle$accel, pch=16)
legend("bottomright", c("Percentile", c("10%", "25%", "50%", "75%", "90%")),
       lwd=c(0, rep(1.5, 5)), col=bob.ross)
dev.off()
