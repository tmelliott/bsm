##' Converts the glm coefficients from intercept and slope to the more interpretable L50 and SR.
##'
##' 
##' @title Reparameterize Coefficients
##' @param beta vector of length 2, \code{c("intercept", "slope")}
##' @param L50 numeric, length of 50\% retention
##' @param SR numeric, selection range
##' @param delta numeric, additional Richards' parameter, default = 1
##' @return vector of length 2, \code{c("L50", "SR")} or \code{beta} if L50 and SR specified
##' @author Tom Elliott
##' @export
repar <- function(beta, L50, SR, delta = 1) {
    if (missing(beta)) {
        if (missing(L50) | missing(SR)) 
            stop("You must specify L50 and SR")
        
        beta1 <- (delta * log(3) - log(4^delta - 3^delta) + log(4^delta - 1)) / SR
        beta0 <- - (L50 / SR) * (delta * log(3) - log(4^delta - 3^delta) +
                                 log(4^delta -1)) - log(2^delta - 1)

        out <- c(beta0, beta1)
    } else {
        alpha <- beta[1]
        beta <- beta[2]
        L50 <- - (alpha + log(2^delta - 1)) / beta
        SR <- (delta * log(3) - log(4^delta - 3^delta) +
           log(4^delta - 1)) / beta
        
        out <- structure(c(L50, SR), .Names = c("L50", "SR"))
    }

    out
}
