#' A 'nnmi' class
#' 
#' @description A list-based S3 class for storing imputation results from nearest neighbor based multiple imputation
#'
#' @param dat a list of imputed data set. 
#' @param ... further arguments passed to or from other methods.
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
#' # S3 method for class 'nnmi'
#' print(imp.dat)
#' 
#' 
#' @export

nnmi <- function(dat, ...) UseMethod("nnmi")

