#' Perform Kaplan-Meier estmation over the multiply imputed survival data sets
#'
#' @description
#' This function estimates Kaplan-Meier estimates based on Rubin's rules (multiple imputation algorithms) (Rubin, 2004).
#'
#' @param obj A 'nnmi' object, that contains imputed data for the missing covariate and the censored observations.
#' @param time A vector contains the observed time.
#' @param status A vector contains the event indicator.
#'
#' @return A data frame contains pooled Kaplan-Meier estimates.
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
#' # performance multiple imputation on missing covariate t5 and
#' # censored observations based on the imputed missing covariates
#' imp.dat <- NNMIS(t5, xa=age, xb=age, time=time, event=status, imputeCT=TRUE, Seed = 1234)
#'
#' # check imputation results
#' head(imp.dat$dat.T.NNMI)
#'
#' # combine inference from imputed data sets using Rubin's rules
#' # Kaplan-Meier estimates
#' kmfit <- km.pool(imp.dat, time, status)
#' plotKM(kmfit)
#'
#' @references
#' Rubin DB. Multiple imputation for nonresponse in surveys. New York: John Wiley and Sons; 2004.
#'
#' @export
#'


km.pool <- function(obj, time, status) {
  imp.time <- obj$dat.T.NNMI
  imp.status <- obj$dat.Id.NNMI
  MI <- obj$MI

  fit.org <- survival::survfit(Surv(time, status) ~ 1)
  time.point <- fit.org$time

  coefs <- vars <- matrix(NA, nrow=length(time.point), ncol=MI)

  for(i in 1:MI) {
    kmfit <- survival::survfit(survival::Surv(as.vector(imp.time[,i]), as.vector(imp.status[,i])) ~ 1)
    for(j in 1:length(time.point)) {
      if(time.point[j] %in% kmfit$time) {
        coefs[j,i] <- kmfit$surv[kmfit$time==time.point[j]]
        vars[j,i] <- (kmfit$std.err[kmfit$time==time.point[j]])^2
      } else {
        if(j==1) {
          coefs[j,i] <- vars[j,i] <- NA
        } else {
          coefs[j,i] <- coefs[j-1,i]
          vars[j,i] <- vars[j-1,i]
        }
      }
    }
    if(is.na(coefs[1,i])){
      coefs[,i][is.na(coefs[,i])] <- coefs[,i][!is.na(coefs[,i])][1]
      vars[,i][is.na(vars[,i])] <- vars[,i][!is.na(vars[,i])][1]
    }
  }


  coef.imp <- rowMeans(coefs)

  w <- rowMeans(vars)

  B <- apply(coefs,1,Bfun)

  Timp <- w + (1+1/MI)*B

  res <- data.frame(cbind(time=time.point, survival=coef.imp, std.err=sqrt(Timp)))

  res$lower.95CI <- apply(res[,2:3],1,function(x){lci <- x[1]-stats::pnorm(0.975)*x[2]; return(max(0,lci))})
  res$upper.95CI <- apply(res[,2:3],1,function(x){uci <- x[1]+stats::pnorm(0.975)*x[2]; return(min(1,uci))})

  names(res) <- c("time","survival", "std.err", "lower 95% CI", "upper 95% CI")

  return(res)
}


Bfun <- function(x){
  xh <- mean(x)
  xl <- length(x)
  B <- (x-xh) %*% (x-xh)
  return(B/(xl-1))
}
