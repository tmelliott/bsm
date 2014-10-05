##' Returns the BUGS posterior mean for a chosen parameter(s)
##'
##' Just makes life easy ...
##' @title Get the posterior mean of a parameter
##' @param fit a bsmfit object
##' @param par the parameter(s) to obtain means of
##' @return numeric vector
##' @author Tom Elliott
##' @export
postMean <- function(fit, par) {
    fit$fit$BUGSoutput$mean[[par]]
}
