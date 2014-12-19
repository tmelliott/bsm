##' .. content for description{} (no empty lines) ..
##'
##' .. content for details{} ..
##' @title Generate a mathematical expression for a parameter
##' @param par parameter being formulated
##' @param expression an expression
##' @param data data from fitted object
##' @param est parameter estimates for formula
##' @param ... additional arguments
##' @return a bsmFormula object for a given parameter
##' @author Tom Elliott
generateFormula <- function(par, expression, data, est, ...) {
    vars <- attr(terms(expression), "term.labels")
    #vars <- vars[vars != "haul"]
    varL <- lapply(vars, function(x) eval(parse(text = x), data))
    levels <- lapply(varL, levels)

    ## remove "as.factor(  )" from variable names etc
    vars <- sapply(vars, function(v) {
        s1 <- strsplit(v, "[(]")[[1]]
        s <- s1[length(s1)]
        strsplit(s, "[)]")[[1]][1]
    })
    names(levels) <- vars

    parname <- paste0("mu_", par)
    is.int <- parname %in% rownames(est)
    int <- ifelse(is.int, est[parname, ], FALSE)
    
    # drop "baseline" IF no intercept
    levs <- lapply(levels, function(l)
        if (length(l) == 0) {
            NULL
        } else if (is.int) {
            l[-1]
        } else {
            l
        })
    
    levels <- lapply(levels, function(l) if (is.null(l)) 0 else l)

    
    coef <- switch(par, "L50" = "beta", "SR" = "gamma",
                   "phi" = "omega", "delta" = "zeta")

    ri <- grep(coef, rownames(est))
    
    
    out <- list(parameter   = par,
                intercept   = int,
                coef.names  = rownames(est)[ri],
                coef.vals   = est[ri],
                variables   = levels,
                used.levels = levs,
                formula     = expression,
                factors     = levels[!sapply(levels, is.null)])
    
    class(out) <- "bsmFormula"
    out    
}


##' @param x object of class bsmFormula
##' @param use.values if \code{TRUE}, the parameter estimates are used
##' @describeIn generateFormula
##' @export
print.bsmFormula <- function(x, use.values = FALSE, ...) {
    par    <- x$parameter
    int    <- signif(x$intercept, 4)
    coefs  <- x$coef.names
    values <- signif(x$coef.vals, 4)
    levels <- x$variables
    levs   <- x$used.levels
    vars   <- names(levs)

    form.pars <- unlist(lapply(1:length(levs), function(i) {
        if (is.null(levs[[i]])) {
            vars[i]
        } else {
            paste0("(", vars[i], " = ", levs[[i]], ")")
        }
    }))
    
    values.str <- paste(ifelse(values>=0, "+", "-"), c(abs(values)))

    print.coefs  <- paste0(coefs, " * ", form.pars, collapse = " + ")
    print.values <- paste0(values.str, " * ", form.pars, collapse = " ")

    if (use.values) {
        txt <- paste0(par, " = ", int, " ", print.values)
    } else {
        txt <- paste0(par, " = ", if (int != 0) paste0("mu_", par, " + "), print.coefs)
    }

    cat("  ", txt, "\n")
}




##' .. content for description{} (no empty lines) ..
##'
##' .. content for details{} ..
##' @title Predict the values for fixed factor levels
##' @param x bsmFormula object
##' @param ... extra arguments
##' @return Estimates
##' @author Tom Elliott
##' @export
predictPar <- function(x, ...) {
    df <- do.call(expand.grid, x$variables)
    
    mat <- model.matrix(x$formula, df)
    if ("haul" %in% colnames(mat))
        mat <- mat[, colnames(mat) != "haul"]
    
    pred.mat <- cbind(df, mat %*%
                          c(if ("(Intercept)" %in% colnames(mat)) x$intercept, x$coef.vals))
    colnames(pred.mat) <- c(colnames(df), x$parameter)

    attr(pred.mat, "factor.names") <- colnames(pred.mat)
    pred.mat
}
