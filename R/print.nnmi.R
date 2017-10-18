#' Print function for object of 'nnmi' class.
#' 
#' @description print basic information for an object of class 'nnmi'
#'
#' @param x a 'nnmi' object
#' @param ... further arguments passed to function
#' 
#' 
#' @examples
#' # load required packages
#' library(NNMIS)
#' library(survival)
#'
#' # load data set - stanford2 in package 'survival'
#' data("stanford2")
#' head(stanford2)
#' attach(stanford2)
#'
#' # performance multiple imputation on missing covariate t5
#' imp.dat <- NNMIS(t5, xa=age, xb=age, time=time, event=status, Seed = 1234, mc.core=1)
#' 
#' # print
#' imp.dat
#'
#' @export

print.nnmi <- function(x, ...)
{
  cat("\n------------\n")
  cat("  NNMIS\n")
  cat("------------\n")

  cat("Parameters:\n")
  cat("\tNumber of imputation:",x$MI,"\n")
  cat("\tNumber of nearst neighbours considered for imputing missing covariate:", x$NN, "\n")
  cat("\tWeight in the working model1:", x$w1, "\n")
  cat("\tWeight in the working model2:", x$w2, "\n")
  cat("\tDistribution family used in working model2:", x$mfamily, "\n")

  if(x$imputeCT){
    cat("\tNumber of nearst neighbours considered for imputing censoring time:", x$NN.t, "\n")
  }

  cat("\tNumber of cpu cores were used:", x$mc.cores, "\n")

  cat("\n")

}
