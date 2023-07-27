library(GBASS)
library(BASS)
library(lhs)
library(tictoc)
library(stargazer)
library(EpiModel)
library(parallel)
library(RColorBrewer)
source("Epi/stochastic_sir.R")

load("Epi/data/stohcastic_sir_data2.Rda")

rmse <- function(x, y){
  sqrt(mean((x-y)^2))
}
