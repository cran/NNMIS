#' Estimate Cox regression model pooling over the imputed datasets
#'
#' @description
#' This function estimates Cox regression model, taking into account
#' the additional uncertainty that arises due to a finite number of
#' imputations of the missing data.
#'
#' @param obj A 'nnmi' object, that contains a finite number of imputations
#'            of the missing data.
#' @param time A vector contains the observed time.
#' @param status A vector contains the event indicator.
#' @param Z A vector or matrix that contains other covariates.
#' @param forceNumeric Logical, if it is True, the class of imputed variable will force to be numeric.
#'                      The default is FALSE.
#' @param setRef Optional, a reference group can be set for binary or categorical variable.
#'
#' @return A data frame contains pooled estimation of Cox regression model.
#'
#' @examples
#'
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
#' imp.dat <- NNMIS(t5, xa=age, xb=age, time=time, event=status, Seed = 2016)
#'
#' # this program can impute censoring time based on the imputed missing covariate
#' # imp.dat <- NNMIS(t5, xa=age, xb=age, time=time, event=status, imputeCT=T, Seed = 2016)
#' # check imputation results
#' # head(imp.dat$dat.NNMI)    #> missing covariates
#' # head(imp.dat$dat.T.NNMI)  #> censoring time
#' # head(imp.dat$dat.Id.NNMI) #> censoring indicator
#'
#' # check imputation results
#' head(imp.dat$dat.NNMI)
#'
#' # combine inference from imputed data sets by using Rubin's rules
#' # estimates in Cox regression
#' coxph.pool(imp.dat, time, status, age)
#'
#' @export
#'


coxph.pool <- function(obj, time, status, Z, forceNumeric=FALSE, setRef=NULL)
{
  X.NNMI <- obj$dat.NNMI

  MI <- obj$MI

  coefs <- vars <- NULL

  for(i in 1:MI){
    temp <- X.NNMI[,i]
    if((obj$mfamily != 'gaussian') && (forceNumeric==FALSE)){
      temp <- as.factor(temp)
      if(!is.null(setRef)){
        temp_ref <- try(stats::relevel(temp, ref=setRef), silent = T)
        if(!inherits(temp_ref, 'try-error')){
          X <- data.frame(cbind(X=temp_ref, Z))
        }else{
          stop('Could not find this reference category.')
        }
      }else{
        X <- data.frame(cbind(X=temp, Z))
      }
    } else {
      X <- data.frame(cbind(X=as.numeric(temp), Z))
    }
    fit.nnmi <- survival::coxph(survival::Surv(time, status) ~ ., data = X)
    coefs <- cbind(coefs,fit.nnmi$coef)
    vars <- cbind(vars, (summary(fit.nnmi)$coef[,3])^2)
  }

  coef.imp <- rowMeans(coefs)

  w <- rowMeans(vars)

  B <- apply(coefs,1,Bfun)

  Timp <- w + (1+1/MI)*B

  vm <- (MI-1)*(1+w/B/(1+1/MI))^2

  t <- (coef.imp-mean(coef.imp))/(Timp)^(0.5)

  rslt <- data.frame(cbind(coef=coef.imp,
                           coef.exp=exp(coef.imp),
                           coef.se=sqrt(Timp),
                           t=t,
                           df=vm,
                           Pvalue=2*stats::pt(abs(t),df=vm, lower.tail = F)))

  return(rslt)

}


Bfun <- function(x){
  xh <- mean(x)
  xl <- length(x)
  B <- (x-xh) %*% (x-xh)
  return(B/(xl-1))
}
