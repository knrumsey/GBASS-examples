source("Epi/sir_setup.R")
bob <- brewer.pal(10, "Set3")
bob[2] <- "gold1"
bob2 <- brewer.pal(4, "Set1")

# Fit models in parallel
q_vec <- c(0.05, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, .95)
mods <- mclapply(q_vec,
                 function(qq) qbass(X,
                                    y,q=qq,
                                    maxInt=3,
                                    w_prior=list(type="GIG", p=-0.001, a=0, b=0.001, prop_sigma=0.1),
                                    a_lambda=.01, b_lambda=.01,
                                    nmcmc=10000, nburn=8001, thin=2,
                                    Iw0=c(10, 30, 5),
                                    Zw0=c(50, 100, 20, 20)),
                 mc.cores = 5, mc.preschedule = F)

#save(mods, file="Epi/data/mods_list.Rda")

# Get sobol decompositions
sob <- list()
for(i in 1:length(q_vec)){
  mod <- mods[[i]]
  mod <- gm2bm(mod)
  mod$xx.des <- X
  sob[[i]] <- sobol(mod)
  plot(sob[[i]], main=i)
}

for(i in 1:length(q_vec)){
  print(names(sob[[i]]$S))
}


# Make legend keys
leg1 <- c(expression(paste("x"[1])),
          expression(paste("x"[2])),
          expression(paste("x"[3])),
          expression(paste("x"[4])))

leg2 <- c(expression(paste("x"[1], "x"[2])),
          expression(paste("x"[1], "x"[3])),
          expression(paste("x"[1], "x"[4])),
          expression(paste("x"[2], "x"[3])),
          expression(paste("x"[2], "x"[4])),
          expression(paste("x"[3], "x"[4])))

leg3 <- c(expression(paste("x"[1], "x"[2], "x"[3])),
          expression(paste("x"[1], "x"[2], "x"[4])),
          expression(paste("x"[1], "x"[3], "x"[4])),
          expression(paste("x"[2], "x"[3], "x"[4])))


png("Epi/figs/sobol.png", height=6, width=8, units="in", res=300)
par(mfrow=c(2,2))
par(mar=c(5, 4, 4, 2) + 0.1 - c(0.6, 0, 0.6, 0))
# Make Total Sensitivity Plot
Q <- array(NA, dim=c(length(q_vec), 4, 3))
for(i in 1:length(q_vec)){
  A <- sob[[i]]$T
  Q[i,,] <- t(apply(A, 2, function(zz) quantile(zz, c(0.1, 0.5, 0.9))))
}
plot(NULL, xlim=c(0, 1), ylim=c(0, 1), main="Total Sensitivity",
     xlab="Quantile", ylab="Proportion of Variance")
for(j in 1:4){
  polygon(c(q_vec, rev(q_vec)), c(Q[,j,1], rev(Q[,j,3])), border=bob[j], col=adjustcolor(bob[j], alpha.f=0.2))
  points(q_vec, Q[,j,2], pch=j, col=bob[j], lwd=2)
}
legend("right", leg1, col=bob, pch=1:4, horiz = TRUE, pt.lwd=2, cex=0.8)



# Make Sensitivity Plot (Main Effects)
Q <- array(NA, dim=c(length(q_vec), 4, 3))
for(i in 1:length(q_vec)){
  A <- sob[[i]]$S[,1:4]
  Q[i,,] <- t(apply(A, 2, function(zz) quantile(zz, c(0.1, 0.5, 0.9))))
}
plot(NULL, xlim=c(0, 1), ylim=c(0, 1), main="Sensitivity (Main Effects)",
     xlab="Quantile", ylab="Proportion of Variance")
for(j in 1:4){
  polygon(c(q_vec, rev(q_vec)), c(Q[,j,1], rev(Q[,j,3])), border=bob[j], col=adjustcolor(bob[j], alpha.f=0.2))
  points(q_vec, Q[,j,2], pch=j, col=bob[j], lwd=2)
}
legend("right", leg1, col=bob, pch=1:4, horiz = TRUE, pt.lwd=2, cex=0.8)


# Make Sensitivity Plot (2 Way Interactions)
Q <- array(NA, dim=c(length(q_vec), 6, 3))
for(i in 1:length(q_vec)){
  A <- sob[[i]]$S[,5:10]
  Q[i,,] <- t(apply(A, 2, function(zz) quantile(zz, c(0.1, 0.5, 0.9))))
}
plot(NULL, xlim=c(0, 1), ylim=c(0, 1), main="Sensitivity (2 Way Interactions)",
     xlab="Quantile", ylab="Proportion of Variance")
for(j in 1:6){
  polygon(c(q_vec, rev(q_vec)), c(Q[,j,1], rev(Q[,j,3])), border=bob[4+j], col=adjustcolor(bob[4+j], alpha.f=0.2))
  points(q_vec, Q[,j,2], pch=j, col=bob[4+j], lwd=2)
}
legend("topright", leg2, col=bob[5:10], pch=1:6, horiz = TRUE, pt.lwd=2, cex=0.8)


# Make Sensitivity Plot (3 Way Interactions)
Q <- array(NA, dim=c(length(q_vec), 4, 3))
for(i in 1:length(q_vec)){
  A <- sob[[i]]$S
  cnt <- 1
  for(j in 1:2){
    for(k in (j+1):3){
      for(l in (k+1):4){
        ind <- which(names(A) == paste0(j, "x", k, "x", l))
        if(length(ind) == 0){
          Q[i,cnt,] <- rep(0, 3)
        }else{
          Q[i,cnt,] <- quantile(A[,ind], c(0.1, 0.5, 0.9))
        }
        cnt <- cnt + 1
      }
    }
  }
}
plot(NULL, xlim=c(0, 1), ylim=c(0, 1), main="Sensitivity (3 Way Interactions)",
     xlab="Quantile", ylab="Proportion of Variance")
for(j in 1:4){
  polygon(c(q_vec, rev(q_vec)), c(Q[,j,1], rev(Q[,j,3])), border=bob2[j], col=adjustcolor(bob2[j], alpha.f=0.2))
  points(q_vec, Q[,j,2], pch=j, col=bob2[j], lwd=2)
}
legend("topright", leg3, col=bob2, pch=1:4, horiz = TRUE, pt.lwd=2, cex=0.9)

dev.off()



