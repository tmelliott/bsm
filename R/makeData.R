##' Creates the required data list for JAGS for a given model.
##'
##' Not much to know ...
##' @title Make JAGS Data from bsmdata Object
##' @param x an object of class \code{bsmdata}, from the \code{bsmData} function
##' @param family \code{"logistic"} (default) or \code{"poisson"}, the family for the likelihood
##' @param combine logical, if \code{TRUE}, then the "combined hauls" approach is used, otherwise a
##' hierarchical approach is used.
##' @param random vector of parameters to have hierarchical or random effects
##' @param L50 the formula for L50, default is L50 = ~haul if random = TRUE
##' @param SR the formula for SR, default is SR = ~haul
##' @param phi the formula for phi, default is phi = ~haul
##' @param delta the formula for delta, default is delta = ~1
##' @param ... additional arguments
##' @return A list of data for JAGS
##' @author Tom Elliott
##' @export
makeData <- function(x, family = "binomial", combine = is.null(random), random = NULL,
                     L50 = NULL, SR = NULL, phi = NULL, delta = NULL, ...) {
    
    if (class(x) != "bsmdata")
        stop("x must be a bsmdata object from the bsmData function")
    
    dat <- list(N = attr(x, "nlength"), M = attr(x, "nhaul"), 
                x = sort(unique(x$length)),
                q1 = tapply(x$q1, x$haul, mean), q2 = tapply(x$q2, x$haul, mean))

    if (!is.null(L50)) {
       # dat$PL50 <- ncol(L50)
        dat$L50des <- L50
    }
    if (!is.null(SR)) {
       # dat$PSR <- ncol(SR)
        dat$SRdes <- SR
    }
    if (!is.null(phi)) {
       # dat$Pphi <- ncol(phi)
        dat$phides <- phi
    }
    if (!is.null(delta)) {
       # dat$Pdelta <- ncol(delta)
        dat$deltades <- delta
    }
    
    ## The count information:
    switch(family,
           "binomial" = {
               dat$y <- toMatrix(x, "y1", dat$x)
               dat$n <- toMatrix(x, "n", dat$x)
           },
           "poisson" = {
               dat$y1 <- toMatrix(x, "y1", dat$x)
               dat$y2 <- toMatrix(x, "y2", dat$x)
           })

    dat
}


toMatrix <- function(x, par, all.len) {
    ## I hope to rewrite this using something from the `reshape` package, but not right now ...
    t(do.call(rbind, as.list(tapply(1:length(x[, par]), x$haul, function(i) {
        df <- as.data.frame(x[i, c(par, "length")])
        y <- merge(df, data.frame(length = all.len), by = "length", all = TRUE)[, par]
        y[is.na(y)] <- 0
        y
    }))))
}
