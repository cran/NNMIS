#' Plot function for pooled Kaplan-Meier estimates
#'
#' A plot of survival curves is produced.
#'
#' @param x a data.frame contains pooled estimates of survival function generated from function 'km.pool'.
#'
#' @seealso km.pool
#'
#' @export
#'


plotKM <- function(x){
  #require(survival)
  class(x) <- 'survfit'
  graphics::plot(x,ylab="Survival",xlab="Time")
}
