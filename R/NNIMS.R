###################################################################################################
#                                                                                                 #
#   NNMIS (Nearest Neighbor Based Multiple Imputation for Survival Data with Missing Covariates)  #
#                                                                                                 #
###################################################################################################



#' @title Nearest Neighbor Based Multiple Imputation for Survival Data with Mssing Covariates (NNMIS)
#'
#' @description
#' This function performs the nearest neighbor based multiple imputation approach proposed by
#' Hsu et al. (2006), Long et al. (2012), Hsu et al. (2014) and Hsu and Yu (2017, 2018) to impute for missing covariates
#' and censored observations (optional). To perform imputation for missing covariates, the approach
#' requires one to fit two working models: one for predicting the missing covariate values and the other
#' for predicting the missing probabilities based on the observed data. The distribution of the working
#' model for predicting the missing covariate values will be automatically decided by the data type of
#' the missing covariate. A logistic regression model will be fitted to predict the missing probabilities.
#' The estimation results of the two working models are then used to select a nearest neighborhood for
#' each missing covariate observation. Once the nearest neighborhood is chosen, multiple impuation is then
#' performed on the neighborhood non-parametrically. The detailed procedures can be found in Long et al. (2012), Hsu et al. (2014),
#' and Hsu and Yu (2017, 2018). Similarily, to perform imputation for censored observations, one has to fit two
#' working models first: one for predicting the survival time and the other for predicting the censoring time.
#' These two working models are derived using Cox regression. The estimation results of the two working models are then
#' used to select a nearest neighborhood for each censored observation. Once the nearest
#' neighborhood is chosen, multiple impuation is then performed on the neighborhood non-parametrically.
#' The detailed procedures can be found in Hsu et al. (2006).
#'
#' Note that the current version can only perform imputation for a situation with only one missing covariate.
#' Before you use this package, please check the input covariates matrix to see if there is more than one
#' missing covariate.
#'
#' @param y Can be any vector of covariate, which contains missing values to be imputed. Missing values are coded as NA.
#' @param xa Can be any vector or matrix, which will be used as the covariates along with the estimated cumulative baseline hazard
#' and the observed censoring indicator for the working model of predicting the missing covariate values.
#'            Note that no missing values are allowed for this.
#' @param xb Can be any vector or matrix, which will be used as the covariates along with the estimated cumulative baseline hazard
#' and the observed censoring indicator for the working model of predicting the missing probabilities.
#'            Note that no missing values are allowed for this.
#' @param time This is the observed time.
#' @param event This is the censoring indicator, i.e. 0:censored; 1: event.
#' @param MI Number of imputation. The default is MI=10.
#' @param NN Size of the nearest neighborhood considered for imputing missing covariate. Default is NN=5.
#' @param w1 Weight will be used in the working model of predicting the missing covariate values. The default is w1=0.8.
#' @param w2 Weight will be used in the working model of predicting the missing probabilities. The default is w1=0.2.
#' @param imputeCT Logical. If TRUE, survival times for censored observations will be imputed and exported as part of output. (optional)
#' @param NN.t Size of the nearest neighborhood considered for imputing survival times for each censored observation. Default is NN.t=10.
#' @param mc.cores Number of cpu cores to be used. This option depends on package "parallel". The default is mc.core=1.
#' @param Seed An integer that is used as argument by the set.seed() for
#'  offsetting the random number generator. Default is to leave the random number generator alone.
#' @param verbose If True, print messages.
#'
#' @return An object of class "nnmi" is a list containing parameters used in multiple imputation and all outputs.
#' \item{N}{ Number of observations.}
#' \item{MI}{ Number of imputation.}
#' \item{NN}{ Size of the nearest neighborhood considered for imputing missing covariate.}
#' \item{w1}{ Weight in the working model for predicting the missing covariate values/survival times.}
#' \item{w2}{ Weight in the working model for predicting the missing probabilities/censoring times.}
#' \item{mfamily}{ Distribution family used in the working model for predicting the missing covariate values.}
#' \item{imputeCT}{ Logical, whether to impute survival times for censored observations or not.}
#' \item{dat.NNMI}{ data frame containing imputed missing covariate values.}
#' \item{dat.T.NNMI}{ data frame containing imputed survival times.}
#' \item{dat.Id.NNMI}{ data frame containing censoring indicator.}
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
#' imp.dat <- NNMIS(t5, xa=age, xb=age, time=time, event=status, Seed = 2016, mc.core=1)
#'
#' # check imputation results
#' head(imp.dat$dat.NNMI)
#'
#' # this program can impute survival times for censored observations based on 
#' # the imputed missing covariate values
#' # imp.dat <- NNMIS(t5, xa=age, xb=age, time=time, event=status, imputeCT=TRUE, Seed = 2016)
#' # check imputation results
#' # head(imp.dat$dat.NNMI)    # imputed missing covariate values
#' # head(imp.dat$dat.T.NNMI)  # imputed survival times
#' # head(imp.dat$dat.Id.NNMI) # censoring indicator
#'
#' @references
#'
#' Hsu CH, Taylor JM, Murray S, Commenges D. Survival analysis using auxiliary variables via nonparametric multiple
#' imputation. Statistics in Medicine 2006; 25: 3503-17.
#'
#' Hsu CH, Long Q, Li Y, Jacobs E. A Nonparametric Multiple Imputation Approach for Data with Missing Covariate Values with Application to Colorectal Adenoma Data. Journal of Biopharmaceutical Statistics 2014; 24: 634-648.
#' 
#' Hsu CH, Yu M. Cox regression analysis with missing covariates via nonparametric multiple imputation. arXiv 2017; 1710.04721.
#' 
#' Hsu CH, Yu M. Cox regression analysis with missing covariates via nonparametric multiple imputation. Statistical Methods in Medical Research 2018; doi: 10.1177/0962280218772592.
#' 
#' Long Q, Hsu CH, Li Y. Doubly robust nonparametric multiple imputation for ignorable missing data. Statistica Sinica 2012; 22: 149-172.
#' 
#' @export
#'




NNMIS <- function (y, xa=NULL, xb=NULL, time, event, MI = 10, NN = 5, w1 = 0.8, w2 = 0.2, Seed=NA, imputeCT = FALSE, NN.t=10, mc.cores=1, verbose=TRUE)
{
  #require(survival)

  #------------------------------------------------------------------------------------------------------------------------
  # check the input data

  y <- as.vector(y)

  # number of observations
  N <- length(time)

  xa <- as.matrix(xa)
  xb <- as.matrix(xb)

  # check missingness

  # 1. must have NA's in target vector
  if(!anyNA(y)) {stop("No missing values found.")}

  # 2. no missing values allowed in xa, xb, and xc
  if(is.null(xa)) {stop("Need to specify a vector/matrix used for the working model of predicting the missing covariate values.")}
  if(is.null(xb)) {stop("Need to specify a vector/matrix used for the working model of predicting the missing probabilities.)")}
  if(anyNA(xa)) {stop("No missing values allowed in xa.")}
  if(anyNA(xb)) {stop("No missing values allowed in xb.")}

  # check consistency
  if(N != length(y)) {stop("Length of y is not consistent with the length of survival time.")}
  if(N != length(event)) {stop("Length of y is not consistent with the length of censoring indicator.")}
  if(N != nrow(xa)) {stop("Row number of xa is not consistent with the length of y.")}
  if(N != nrow(xb)) {stop("Row number of xb is not consistent with the length of y.")}

  p <- ncol(xa)
  q <- ncol(xb)

  # if want to impute missing covariate, all unique colums will be used in working model
  if(imputeCT) {
    xc.temp <- as.matrix(cbind(xa,xb))
    xc <- as.matrix(xc.temp[,!duplicated(t(xc.temp))])
  }

  ## check data type
  # 1. four data types are allowed for target vector, numeric, character, factor, and logical
  # 2. need to assign mfamily to each type of data
  #         gaussian   - numeric (#unique value > 5)
  #         binomial   - numeric (#unique value = 2)
  #                    - logical
  #                    - factor/character (#unique level = 2)
  #         multinomial- factor/character (#unique level > 2)
  #                    - numeric (2 < #unique value <= 5)

  datatype <- mode(y)
  count <- pmatch(datatype,c('numeric','character','factor','logical'),nomatch = 0L) #legal data type

  if(count==0L) {stop('Can not recognize the type of input data.')}

  # distinguish binomial, multinomial and gaussian
  if(datatype=='character') {y <- as.factor(y)}

  datalevels <- levels(as.factor(y))

  if(length(datalevels)==1) stop('Input data is not valid.')

  if(datatype != 'numeric'){
    y <- as.factor(y)
    if(length(datalevels)==2) {
      mfamily = 'binomial'
      } else {
      mfamily = 'multinomial'
      }
  } else {
    if(length(datalevels) >5) {
      mfamily = 'gaussian'
      } else {
      y <- as.factor(y)
      if(length(datalevels) == 2) {
        mfamily = 'binomial'
        } else {
          mfamily = 'multinomial'
        }
    }
  }
  #-------------------------------------------------------------------------------------------------------------------------


  if(!is.na(Seed)) {set.seed(Seed)}

  # need #(MI) random number as seed for each imputation
  saveSeed <- floor(stats::runif(MI, min=1, max=2001))

  # how many threads are requested
  if(verbose) {message("Request for ", mc.cores, " threads.")}

  if(mc.cores > 1) { # request for multiple threads
    #require(parallel)              # depends on parallel package
    no_cores <- parallel::detectCores() - 1   # probably not going to occupy all threads
    if(mc.cores > no_cores) {
      mc.cores <- no_cores
      if(verbose) message("Only ", mc.cores, " threads can use.")
    } else {
      if(verbose) message(mc.cores, " threads are using.")
    }

    if(.Platform$OS.type == "unix") {
      Dat <- parallel::mclapply(saveSeed, ff, y, xa, xb, time, event, 1, NN, w1, w2, imputeCT, NN.t, mfamily, datalevels, xc)
    }else{
      type <- if (exists("mcfork", mode="function")) "FORK" else "PSOCK"
      cl <- parallel::makeCluster(mc.cores, type=type)
      parallel::setDefaultCluster(cl)
      parallel::clusterEvalQ(NULL, library(NNMIS))
      parallel::clusterEvalQ(NULL, library(survival))
      Dat <- parallel::parLapply(NULL, saveSeed, ff, y, xa, xb, time, event, 1, NN, w1, w2, imputeCT, NN.t, mfamily, datalevels, xc)
      parallel::stopCluster(cl)
    }

    dat.NNMI <- dat.T.NNMI <- dat.Id.NNMI <- as.data.frame(matrix(NA, nrow=N, ncol=MI))
    colnames(dat.NNMI) <- colnames(dat.T.NNMI) <- colnames(dat.Id.NNMI) <- paste("M",1:MI,sep="")
    if(imputeCT==FALSE){
      for(i in 1:MI){
        temp <- Dat[[i]][[1]][,1]
        if(mfamily != 'gaussian') temp <- factor(Dat[[i]][[1]][,1],labels=datalevels)
        dat.NNMI[,i] <- temp
      }
      dat.T.NNMI <- NA
      dat.Id.NNMI <- NA
    }else{
      for(i in 1:MI){
        temp <- Dat[[i]][[1]][,1]
        if(mfamily != 'gaussian') temp <- factor(Dat[[i]][[1]][,1],labels=datalevels)
        dat.NNMI[,i] <- temp
        dat.T.NNMI[,i] <- Dat[[i]][[2]][,1]
        dat.Id.NNMI[,i] <- Dat[[i]][[3]][,1]

        dat.T.NNMI[is.na(dat.T.NNMI[,i]),i] <- time[is.na(dat.T.NNMI[,i])]
        dat.Id.NNMI[is.na(dat.Id.NNMI[,i]),i] <- event[is.na(dat.Id.NNMI[,i])]
      }
    }

  } else {
    if (verbose) message("Only 1 thread is used.")
    Dat <- sapply(saveSeed, ff, y, xa, xb, time, event, 1, NN, w1, w2, imputeCT, NN.t, mfamily, datalevels, xc)

    dat.NNMI <- dat.T.NNMI <- dat.Id.NNMI <- as.data.frame(matrix(NA, nrow=N, ncol=MI))
    colnames(dat.NNMI) <- colnames(dat.T.NNMI) <- colnames(dat.Id.NNMI) <- paste("M",1:MI,sep="")
    if(imputeCT==FALSE){
      for(i in 1:MI){
        temp <- Dat[[3*i-2]][,1]
        if(mfamily != 'gaussian') temp <- factor(Dat[[3*i-2]][,1],labels=datalevels)
        dat.NNMI[,i] <- temp
      }
      dat.T.NNMI <- NA
      dat.Id.NNMI <- NA
    }else{
      for(i in 1:MI){
        temp <- Dat[[3*i-2]][,1]
        if(mfamily != 'gaussian') temp <- factor(Dat[[3*i-2]][,1],labels=datalevels)
        dat.NNMI[,i] <- temp
        dat.T.NNMI[,i] <- Dat[[3*i-1]][,1]
        dat.Id.NNMI[,i] <- Dat[[3*i]][,1]

        dat.T.NNMI[is.na(dat.T.NNMI[,i]),i] <- time[is.na(dat.T.NNMI[,i])]
        dat.Id.NNMI[is.na(dat.Id.NNMI[,i]),i] <- event[is.na(dat.Id.NNMI[,i])]
      }
    }
  }

  rslt <- list(N=N,
               MI=MI,
               mfamily=mfamily,
               NN=NN,
               NN.t=NN.t,
               w1=w1,
               w2=w2,
               imputeCT=imputeCT,
               mc.cores=mc.cores,
               dat.NNMI=dat.NNMI,
               dat.T.NNMI=dat.T.NNMI,
               dat.Id.NNMI=dat.Id.NNMI
               )

  class(rslt) <- 'nnmi'

  return(rslt)

}



ff <- function(X, ...){
  return(NNMI(Seed = X,...))
}
