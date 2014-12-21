##' This function handles manipulation of a provided data set (in the form of a data
##' frame, from \code{read.csv} or similar) to a form that can be used by other
##' \code{\link{bsmData}} functions.
##'
##' The simplest use is where a data frame is provided, and the required column names are
##' specified as function arguments.
##'
##' It is also possible to supply vectors and no data frame.
##' @title Import Data for BSM
##' @param y1 the count for the codend or experimental trawl
##' @param y2 the count fo the cover or control trawl
##' @param n the total fish caught in the trawl
##' @param length the length of fish
##' @param length.unit the unit of length measurements
##' @param haul the haul id
##' @param q1 sampling fraction for codend/experimental trawl
##' @param q2 sampling fraction for cover/experimental trawl
##' @param paired logical, is the data from a paired trawl?
##' @param data the name of the data frame. If \code{NULL}, then the other supplied argument
##' must be vectors. Otherwise if \code{df} is supplied, then the other arguments are the
##' variable names.
##' @param ... additional variables can be passed to the object, such as weights, mesh sizes, and
##' any other covariates to be used in later modelling
##' @return An object of class \code{bsmdf}. This is a list with necessary information,
##' and is what should be supplied as the \code{data} argument to other
##' \code{\link{bsmData}} functions.
##' @author Tom Elliott
##' @export
bsmData <- function(y1, y2, n, length, length.unit = NULL, haul, q1, q2, paired = FALSE,
                    data = NULL, ...) {
    ## This allows users to specify data = `something` but avoids other complicated scenarios
    if (!is.null(data)) {
        call <- match.call()
        call$data <- NULL
        return(with(data, eval(call)))
    }
    
    ## Function proper starts here ----------------------------------------
    if (missing(y1) + missing(y2) + missing(n) != 1)
        stop("Please specify two of: `y1`, `y2`, `n`")
    
    if (missing(n)) {
        n <- y1 + y2
    } else if (missing(y2)) {
        y2 <- n - y1
    } else if (missing(y1)) {
        y1 <- n - y2
    }
    
    if (missing(q1))
        q1 <- rep(1, length(y1))
    if (missing(q2))
        q2 <- rep(1, length(y1))
    
    if (missing(haul))
        haul <- rep(1, length(y1))
    
    df <- list(y1 = y1, y2 = y2, n = n, q1 = q1, q2 = q2,
               length = length, haul = haul)
    dots <- list(...)
    df <- c(df, dots)
    
    class(df) <- "bsmdata"
    attr(df, "paired") <- paired
    attr(df, "length.unit") <- length.unit
    attr(df, "nhaul") <- length(unique(haul))
    attr(df, "nlength") <- length(unique(length))
    
    df
}



##' @param x a bsmData object
##' @describeIn bsmData Print information about the data object
##' @export
print.bsmdata <-
    function(x, ...) {
       # print(do.call(data.frame, x))
        cat("\nData set for ", ifelse(attr(x, "paired"), "paired", "covered-codend"),
            " experimental trawl.\n\n", sep = "")
        cat("Total of ", length(x$y1), " observations from ",
            attr(x, "nhaul"), " hauls.\n", sep = "")
        cat("Lengths range from ", min(x$length), " to ", max(x$length),
            attr(x, "length.unit"), ".\n", sep = "")
        cat(" \n")
        invisible(NULL)
    }


##' Return the first few rows of the data set.
##'
##' Acts just like head.data.frame
##' 
##' @title Head of BSM Data
##' @param x a bsmdata object
##' @param ... extra arguments
##' @return data.frame
##' @author Tom Elliott
##' @export head.bsmdata
head.bsmdata <-
    function(x, ...) {
        print(head(as.data.frame(x)))
        invisible(NULL)
    }


##' Return the last few rows of the data set.
##'
##' Acts just like tail.data.frame
##' 
##' @title Tail of BSM Data
##' @param x a bsmdata object
##' @param ... extra arguments
##' @return data.frame
##' @author Tom Elliott
##' @export head.bsmdata
tail.bsmdata <-
    function(x, ...) {
        print(tail(as.data.frame(x)))
        invisible(NULL)
    }

##' @export 
as.data.frame.bsmdata <- function(x, ...) {
    do.call(data.frame, x)
}

##' @export
`[.bsmdata` <- function(x, i, j) {
    df <- `[.data.frame`(as.data.frame(x), i, j)
    att <- attributes(x)
    att$names <- attr(df, "names")
    attributes(df) <- att
    df
}
