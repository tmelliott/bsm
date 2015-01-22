##' Computes the posterior probability of the null hypothesis being true.
##'
##' The parameter has a posterior distribution, and this function computes
##' the posterior probability that the parameter is consistently larger or
##' smaller than the hypothesised value, which ever is smaller,
##'
##' min(Pr(theta > H0), Pr(theta < H0))
##'
##' In the case of a one-sided test, the probability is simply
##' the posterior probability that the parameter is greater than,
##' or smaller than, the hypothesised value.
##' 
##' @title Generate a Bayesian p-value for a Hypothesis Test
##' @param fit a bsmfit object
##' @param parameter the name of the parameter to test
##' @param H0 the null hypothesis value
##' @param alternative the alternative hypotheis, one of "!=", "<" or ">"
##' @return numeric
##' @author Tom Elliott
##' @export
testH0 <- function(fit, parameter, H0, alternative = "!=") {
    post <- fit$fit$BUGSoutput$sims.matrix[, parameter]
    p <- mean(post < H0)

    switch(alternative,
           "!=" = min(p, 1 - p),
           "<" = 1 - p,
           ">" = p)
}
