#' Fits a quantile regression SVM based on the Pinball Loss
#'
#' @param x An n X m matrix containing the predictors (n= number of observatiosn, m = number of predictors) .
#' @param y The Response onto which the qrsvm shall be fitted
#' @param kernel a string giving the type of kernels from package kernlab to use f.e. "rbfdot" for Radial Basis Function Kernel. All Kernels except "stringdot" supported.
#' @param cost The Cost parameter see f.e. package "e1071" and "kernlab"
#' @param tau The Quantile that shall be estimated. 0<=tau<=1
#' @param sigma A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param degree A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param scale A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param offset A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param order A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @details There is no preimplemented scaling of the input variables which should be considered beforehand. Also optimization is based on "quadprog:solve.QP" function which can be considerably slow compared to other SVM implementations.
#' @return An object of class "qrsvm"
#' @references "Nonparametric Quantile Regression" by I.Takeuchi, Q.V. Le, T. Sears, A.J. Smola (2004)
#' @import kernlab
#' @import Matrix
#' @import quadprog
#' @export
#' @examples
#' n<-200
#'
#'x<-as.matrix(seq(-1.5,1.5,length.out = n))
#'y<-rnorm(n)*(0.3+abs(sin(x)))
#'
#'plot(x,y)
#'
#'mod005<-qrsvm(x,y, tau=0.05)
#'mod095<-qrsvm(x,y, tau=0.95)
#'lines(x, mod005$fitted, col="red")
#'lines(x, mod095$fitted, col="red")
qrsvm <- function(x, y, kernel = "rbfdot", cost = 1,    tau = 0.95,
                  sigma = 5, degree = 2, scale = 1, offset = 1, order = 1) {

  requireNamespace("kernlab")
  requireNamespace("quadprog")
  requireNamespace("Matrix")

  if(class(kernel)=="character"){
  if (kernel == "rbfdot") {
    kern <- rbfdot(sigma = sigma)
    kernmat <- kernelMat(kern, x)
  }
  if (kernel == "polydot") {
    kern <- polydot(degree = degree, scale = scale, offset = offset)
    kernmat <- kernelMat(kern, x)
  }
  if (kernel == "vanilladot") {
    kern <- vanilladot()
    kernmat <- kernelMat(kern, x)
  }
  if (kernel == "tanhdot") {
    kern <- tanhdot(scale = scale, offset = offset)
    kernmat <- kernelMat(kern, x)
  }
  if (kernel == "laplacedot") {
    kern <- laplacedot(sigma = sigma)
    kernmat <- kernelMat(kern, x)
  }
  if (kernel == "besseldot") {
    kern <- besseldot(sigma = sigma, order = order, degree = degree)
    kernmat <- kernelMat(kern, x)
  }
  if (kernel == "anovadot") {
    kern <- anovadot(sigma = sigma, degree = degree)
    kernmat <- kernelMat(kern, x)
  }
  if (nrow(kernmat) < nrow(x)) {
    print("kernelMat not valid! Check if valid kernel type stated!")
  }
  pdcalc <- nearPD(kernmat)
  pdmat <- pdcalc$mat
}



  if(class(kernel)=="matrix"){
    pdmat<- kernel
    kernmat<- kernel
    if(nrow(kernel)!=ncol(kernel)){
      cat("Given kernelMat has different col- rownumbers!! Check!", fill=TRUE)
    }
  }

    Amat <- rbind(rep(1, nrow(x)), diag(x = 1, nrow = nrow(x)),
        diag(x = -1, nrow = nrow(x)))

    Dmat <- pdmat
    dvec <- y
    b0 <- c(0, rep((cost * (tau - 1)), nrow(x)), rep(-cost *
        tau, nrow(x)))
    erg <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat),
        bvec = b0, meq = 1, factorized = FALSE)
    alpha <- erg$solution
      f <- alpha %*% kernmat
    offshift <- which.min((round(alpha, 3) -(cost * tau))^2 + (round(alpha, 3) - (cost * (tau - 1)))^2)
    b <- y[offshift] - f[offshift]
    fnew <- alpha %*% kernmat + b



    if(class(kernel)=="character"){
      model <- list(alpha = alpha, xtrain = x, kernel = kern,
                    sigma = sigma, cost = cost, b0 = b, fitted = as.numeric(fnew),
                    tau = tau, scale = scale, offset = offset, order = order,
                    kernstring = kernel, y = y)
    }



    if(class(kernel)=="matrix"){
      model <- list(alpha = alpha, xtrain = x, kernel = 0,
                    sigma = sigma, cost = cost, b0 = b, fitted = as.numeric(fnew),
                    tau = tau, scale = scale, offset = offset, order = order,
                    kernstring = kernel, y = y)
  }


    class(model) <- "qrsvm"
    return(model)
}







