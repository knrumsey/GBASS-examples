source("Epi/stochastic_sir.R")
#X <- maximinLHS(2000, 4)
tic()
y <- apply(X, 1, stochastic_sir)
toc()
save(X, y, file="Epi/data/stochastic_sir_data2.Rda")

hist(y, breaks=30)
