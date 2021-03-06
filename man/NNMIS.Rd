% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NNIMS.R
\name{NNMIS}
\alias{NNMIS}
\title{Nearest Neighbor Based Multiple Imputation for Survival Data with Mssing Covariates (NNMIS)}
\usage{
NNMIS(y, xa = NULL, xb = NULL, time, event, MI = 10, NN = 5,
  w1 = 0.8, w2 = 0.2, Seed = NA, imputeCT = FALSE, NN.t = 10,
  mc.cores = 1, verbose = TRUE)
}
\arguments{
\item{y}{Can be any vector of covariate, which contains missing values to be imputed. Missing values are coded as NA.}

\item{xa}{Can be any vector or matrix, which will be used as the covariates along with the estimated cumulative baseline hazard
and the observed censoring indicator for the working model of predicting the missing covariate values.
           Note that no missing values are allowed for this.}

\item{xb}{Can be any vector or matrix, which will be used as the covariates along with the estimated cumulative baseline hazard
and the observed censoring indicator for the working model of predicting the missing probabilities.
           Note that no missing values are allowed for this.}

\item{time}{This is the observed time.}

\item{event}{This is the censoring indicator, i.e. 0:censored; 1: event.}

\item{MI}{Number of imputation. The default is MI=10.}

\item{NN}{Size of the nearest neighborhood considered for imputing missing covariate. Default is NN=5.}

\item{w1}{Weight will be used in the working model of predicting the missing covariate values. The default is w1=0.8.}

\item{w2}{Weight will be used in the working model of predicting the missing probabilities. The default is w1=0.2.}

\item{Seed}{An integer that is used as argument by the set.seed() for
offsetting the random number generator. Default is to leave the random number generator alone.}

\item{imputeCT}{Logical. If TRUE, survival times for censored observations will be imputed and exported as part of output. (optional)}

\item{NN.t}{Size of the nearest neighborhood considered for imputing survival times for each censored observation. Default is NN.t=10.}

\item{mc.cores}{Number of cpu cores to be used. This option depends on package "parallel". The default is mc.core=1.}

\item{verbose}{If True, print messages.}
}
\value{
An object of class "nnmi" is a list containing parameters used in multiple imputation and all outputs.
\item{N}{ Number of observations.}
\item{MI}{ Number of imputation.}
\item{NN}{ Size of the nearest neighborhood considered for imputing missing covariate.}
\item{w1}{ Weight in the working model for predicting the missing covariate values/survival times.}
\item{w2}{ Weight in the working model for predicting the missing probabilities/censoring times.}
\item{mfamily}{ Distribution family used in the working model for predicting the missing covariate values.}
\item{imputeCT}{ Logical, whether to impute survival times for censored observations or not.}
\item{dat.NNMI}{ data frame containing imputed missing covariate values.}
\item{dat.T.NNMI}{ data frame containing imputed survival times.}
\item{dat.Id.NNMI}{ data frame containing censoring indicator.}
}
\description{
This function performs the nearest neighbor based multiple imputation approach proposed by
Hsu et al. (2006), Long et al. (2012), Hsu et al. (2014) and Hsu and Yu (2017, 2018) to impute for missing covariates
and censored observations (optional). To perform imputation for missing covariates, the approach
requires one to fit two working models: one for predicting the missing covariate values and the other
for predicting the missing probabilities based on the observed data. The distribution of the working
model for predicting the missing covariate values will be automatically decided by the data type of
the missing covariate. A logistic regression model will be fitted to predict the missing probabilities.
The estimation results of the two working models are then used to select a nearest neighborhood for
each missing covariate observation. Once the nearest neighborhood is chosen, multiple impuation is then
performed on the neighborhood non-parametrically. The detailed procedures can be found in Long et al. (2012), Hsu et al. (2014),
and Hsu and Yu (2017, 2018). Similarily, to perform imputation for censored observations, one has to fit two
working models first: one for predicting the survival time and the other for predicting the censoring time.
These two working models are derived using Cox regression. The estimation results of the two working models are then
used to select a nearest neighborhood for each censored observation. Once the nearest
neighborhood is chosen, multiple impuation is then performed on the neighborhood non-parametrically.
The detailed procedures can be found in Hsu et al. (2006).

Note that the current version can only perform imputation for a situation with only one missing covariate.
Before you use this package, please check the input covariates matrix to see if there is more than one
missing covariate.
}
\examples{
# load required packages
library(NNMIS)
library(survival)

# load data set - stanford2 in package 'survival'
data("stanford2")
head(stanford2)
attach(stanford2)

# performance multiple imputation on missing covariate t5
imp.dat <- NNMIS(t5, xa=age, xb=age, time=time, event=status, Seed = 2016, mc.core=1)

# check imputation results
head(imp.dat$dat.NNMI)

# this program can impute survival times for censored observations based on 
# the imputed missing covariate values
# imp.dat <- NNMIS(t5, xa=age, xb=age, time=time, event=status, imputeCT=TRUE, Seed = 2016)
# check imputation results
# head(imp.dat$dat.NNMI)    # imputed missing covariate values
# head(imp.dat$dat.T.NNMI)  # imputed survival times
# head(imp.dat$dat.Id.NNMI) # censoring indicator

}
\references{
Hsu CH, Taylor JM, Murray S, Commenges D. Survival analysis using auxiliary variables via nonparametric multiple
imputation. Statistics in Medicine 2006; 25: 3503-17.

Hsu CH, Long Q, Li Y, Jacobs E. A Nonparametric Multiple Imputation Approach for Data with Missing Covariate Values with Application to Colorectal Adenoma Data. Journal of Biopharmaceutical Statistics 2014; 24: 634-648.

Hsu CH, Yu M. Cox regression analysis with missing covariates via nonparametric multiple imputation. arXiv 2017; 1710.04721.

Hsu CH, Yu M. Cox regression analysis with missing covariates via nonparametric multiple imputation. Statistical Methods in Medical Research 2018; doi: 10.1177/0962280218772592.

Long Q, Hsu CH, Li Y. Doubly robust nonparametric multiple imputation for ignorable missing data. Statistica Sinica 2012; 22: 149-172.
}
