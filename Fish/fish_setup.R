library(GBASS)
library(BASS)
library(lhs)
library(tictoc)
library(stargazer)
library(quantkriging)
# qrsvm package is no longer on CRAN so use an archived local version
# May need to install some dependencies first
# install.packages("qrsvm/", type="source", repos=NULL)
library(qrsvm)

rmse <- function(x, y){
  sqrt(mean((x-y)^2))
}
