#' Fits multiple Quantile Regression SVM
#'
#' @param x An n X m matrix containing the predictors (n= number of observatiosn, m = number of predictors) .
#' @param y The Response onto which the qrsvm shall be fitted
#' @param kernel a string giving the type of kernels from package kernlab to use f.e. "rbfdot" for Radial Basis Function Kernel. All Kernels except "stringdot" supported.
#' @param cost The Cost parameter see f.e. package "e1071" and "kernlab"
#' @param tau The Quantile that shall be estimated. A Vector of values (0<=tau<=1)
#' @param sigma A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param degree A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param scale A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param offset A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param order A possible tuning parameter for specific Kernelfunctions, see package kernlab.
#' @param doPar Should a parallel backend be used. Logical.
#' @param clustnum The number of parallel tasks to use given doPar==TRUE. Default = 2.
#' @details There is no preimplemented scaling of the input variables which should be considered beforehand. Also optimization is based on "quadprog:solve.QP" function which can be considerably slow compared to other SVM implementations.
#' @return An object of class "qrsvm"
#' @references "Nonparametric Quantile Regression" by I.Takeuchi, Q.V. Le, T. Sears, A.J. Smola (2004)
#' @import kernlab
#' @import Matrix
#' @import quadprog
#' @import doParallel
#' @import foreach
#' @export
#' @examples
#'n<-200
#'
#'x<-as.matrix(seq(-2,2,length.out = n))
#'y<-rnorm(n)*(0.3+abs(sin(x)))
#'
#'plot(x,y)
#'
#'models<-list()
#'quant<-c(0.01,0.25,0.5,0.75,0.99)
#'models<-multqrsvm(x,y,tau = quant, doPar=FALSE, sigma = 1)
#'for(i in 1:length(models)){
#'  lines(x, models[[i]]$fitted, col="red")
#'}
multqrsvm<- function(x, y, kernel = "rbfdot", cost = 1,    tau = c(0.05,0.25,0.5,0.75,0.95),
                     sigma = 5, degree = 2, scale = 1, offset = 1, order = 1, doPar=FALSE, clustnum=2){

  requireNamespace("kernlab", quietly = TRUE)
  requireNamespace("Matrix", quietly = TRUE)


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
  pdmat <- as.matrix(pdcalc$mat)





  models<-list()
  if(doPar==FALSE){
    for(i in 1:length(tau)){
      models[[i]]<-qrsvm(x=x, y=y, kernel = pdmat, cost=cost,
                         tau=tau[i], sigma = sigma, degree = degree,
                         scale = scale, offset = offset, order = order)

    }
  }
  else{
    requireNamespace("doParallel", quietly = TRUE)

     registerDoParallel(clustnum)

    models<-foreach(i=1:length(tau), .multicombine = TRUE, .combine = list, .packages = c("qrsvm","kernlab")) %dopar% {
      qrsvm(x=x, y=y, kernel = pdmat, cost=cost,
            tau=tau[i], sigma = sigma, degree = degree,
            scale = scale, offset = offset, order = order)
    }
  }

  for(i in 1:length(models)){
    models[[i]]$kernel<-kern
  }

  class(models)<-"multqrsvm"
  return(models)
}
