#' Inputs
#' x1 inf.prob  (0.25, 0.75)
#' x2 R0  (0.95, 2.0)
#' x3 rec.rate  (0.1, 0.2)
#' x4 inter.eff (0.2, 0.8)
#' x5 Nothing

stochastic_sir <- function(x, scale01=TRUE, n=1){
  if(scale01){
    x[1] <- x[1]*(0.7 - 0.4)  + 0.4
    x[2] <- x[2]*(4.0 - 1) + 1
    x[3] <- x[3]*(0.2 - 0.08) + 0.08
    x[4] <- x[4]*(1.0 - 0.5) + 0.5
  }

  N_tot <- 5000 # population
  N_inf <- 1
  N_t   <- 21
  N_i   <- 14

  param <- param.icm(inf.prob = x[1], act.rate = x[2], rec.rate=x[3],
                     inter.start=N_i, inter.eff=x[4]) # rescale inputs...
  init <- init.icm(s.num = N_tot, i.num = N_inf, r.num=0) # 500 suceptible, 1 infected
  control <- control.icm(type = "SIR", nsims = n, nsteps = N_t) # single simulation
  mod <- icm(param, init, control)
  res <- mod$epi$i.num[N_t,] + mod$epi$r.num[N_t,]
  return(res)
}



