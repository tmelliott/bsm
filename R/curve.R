##' Draws a selection curve for a given set of parameter values.
##'
##' If \code{obj} if specified, then the appropriate axes are drawn using that. Alternatively,
##' generic x-limits can be specified by \code{xlim}. If these are both empty, but \code{add =
##' TRUE}, then the values are taken from \code{par()}.
##' @title Draw Selection Curves
##' @param x vector of L50, SR, and delta (optional)
##' @param obj a \code{bsmdata} object that includes the data, which will be drawn if \code{add = FALSE}
##' @param xlim numeric vector of length 2, specifying the range of x values to draw the curve over
##' @param add logical, if \code{TRUE}, the curve is added to an existing plot, otherwise a new one
##'     is drawn.
##' @param ... additional arguments
##' @return NULL
##' @author Tom Elliott
##' @export
bsmCurve <- function(x, obj, xlim, add = TRUE, ...) {
    if (missing(xlim)) {
        if (missing(obj)) {
            xlim <- par()$usr[1:2]
        } else {
            xlim <- range(obj$x, finite = TRUE)
        }
    }
    
    if (length(x) < 2)
        stop("x must be a vector of values for L50, SR, and delta (= 1 by default)")

    if (!is.null(names(x))) {
        names(x) <- gsub("mu_", "", names(x))
        if (any(!names(x) %in% c("L50", "SR", "delta", "phi")))
            stop("Names of x must be `L50`, `SR`, `delta` and `phi`")
    }
    
    if (length(x) == 2) {
        delta <- 1
        phi <- 1
    } else if (length(x) == 3) {
        if (is.null(names(x))) {
            warning("Assuming richards curve (not paired).")
            
            phi <- 1
            delta <- ifelse(is.null(names(x)), x[3], x["delta"])
        } else {
            if ("phi" %in% names(x)) {
                phi <- ilogit(x["phi"])
                delta <- 1
            } else {
                delta <- x["delta"]
                phi <- 1
            }
        }
    } else if (length(x) == 4) {
        if (is.null(names(x))) {
            delta <- x[3]
            phi <- ilogit(x[4])
        } else {
            delta <- x["delta"]
            phi <- ilogit(x["phi"])
        }
    }
    
    if (is.null(names(x))) {
        L50 <- x[1]
        SR <- x[2]
    } else {
        L50 <- x["L50"]
        SR <- x["SR"]
    }

    xx <- seq(xlim[1], xlim[2], length = 1001)
    eta <- richards(xx, L50, SR, delta)
    yy <- (1 / (1 + exp(-eta))) ^ (1 / delta) * phi

    if (add)
        lines(xx, yy, ...)
    else
        plot(xx, yy, type = "l", ...)

    return(invisible(NULL))
}


richards <- function(x, L50, SR, delta = 1)
    (delta * log(3) - log(4^delta - 3^delta) +
     log(4^delta - 1)) / SR * (x - L50) - log(2^delta - 1)


ilogit <- function(x) {
    exp(x) / (1 + exp(x))
}
