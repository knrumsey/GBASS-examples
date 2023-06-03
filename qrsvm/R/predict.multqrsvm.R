#' Predict an Object of class "multqrsvm"
#'
#' @param model An object of class "multqrsvm"
#' @param newdata The predictors of the predictable data in an n X m Matrix
#' @return A list of predicted values
#' @import kernlab
predict.multqrsvm <- function(model, newdata) {
  requireNamespace("kernlab", quietly = TRUE)



  if (class(model)=="multqrsvm"){
    prediction<-list()

    if (ncol(newdata) != ncol(model[[1]]$xtrain)) {
      cat("Newdata has different number of columns than xtrain please check consistency!",
          fill = TRUE)
    }
    for (i in 1:length(model)){
      xold <- model[[i]]$xtrain
      alpha <- model[[i]]$alpha
      kern <- model[[i]]$kernel
      b <- model[[i]]$b0

      prediction[[i]] <- kernelMult(kern, newdata, xold, alpha) + b
    }
  }

  return(prediction)
}
