#' Predict an Object of class "qrsvm"
#'
#' @param model An object of class "qrsvm"
#' @param newdata The predictors of the predictable data in an n X m Matrix
#' @return A numeric vector of predicted values
#' @import kernlab
predict.qrsvm <- function(model, newdata) {
  requireNamespace("kernlab", quietly = TRUE)



  if (class(model)=="qrsvm"){
  xold <- model$xtrain

  if (ncol(newdata) != ncol(xold)) {
    cat("Newdata has different number of columns than xtrain please check consistency!",
        fill = TRUE)
  }
  alpha <- model$alpha
  kern <- model$kernel
  b <- model$b0
  prediction <- kernelMult(kern, newdata, xold, alpha) + b
  }

  return(prediction)
}
