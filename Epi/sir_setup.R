library(GBASS)
library(BASS)
library(lhs)
library(tictoc)
library(stargazer)
library(EpiModel)
library(parallel)
library(RColorBrewer)
source("Epi/stochastic_sir.R")

rmse <- function(x, y){
  sqrt(mean((x-y)^2))
}
