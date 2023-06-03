source("Fish/fish_setup.R")
load("Fish/data/fish_sims.Rda")

# QUick figure (not in paper)
plot(apply(res1, 2, mean), type='o', ylim=c(0, 8))
for(i in 2:3){
  resi <- get(paste0("res", i))
  lines(apply(resi, 2, mean), col=i)
  points(apply(resi, 2, mean), col=i)
}


# Make figure
bob <- RColorBrewer::brewer.pal(8, "YlOrBr")[c(2,5,8)]
png("Fish/figs/boxplots.png", height=5, width=8, units="in", res=300)
boxplot(res1, width=rep(0.25, 9), at=seq(1, 33, by=4), range=0, ylim=c(0, max(c(res1, res2, res3))), xaxt='n',
        ylab="RMSE", xlab="Number of Replications", xlim=c(0, 35),
        col=bob[1])
boxplot(res2, add=TRUE, width=rep(0.25, 9), at=seq(2, 34, by=4), col=bob[2], range=0, xaxt='n')
boxplot(res3, add=TRUE, width=rep(0.25, 9), at=seq(3, 35, by=4), col=bob[3], range=0, xaxt='n')
axis(1, seq(2, 34, by=4), c("1", "2", "4", "8", "16", "32", "64", "128", "256"))
legend("topright",
       c("GBMARS", "quantKrig", "qrsvm"),
       fill=bob, cex=1.2)
dev.off()


# Make table
cnt <- matrix(0, nrow=3, ncol=9)
for(i in 1:20){
  for(j in 1:9){
    k <- which.min(c(res1[i,j], res2[i,j], res3[i,j]))
    cnt[k, j] <- cnt[k,j] + 1
  }
}
tab <- 5*cnt
colnames(tab) <- c("1", "2", "4", "8", "16", "32", "64", "128", "256")
rownames(tab) <- c("GBMARS", "QuantKrig", "qrsvm")
stargazer(tab, summary=FALSE)
