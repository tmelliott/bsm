##' Returns lists of parameter names for summaries, initial values, etc.
##'
##' Not much to know.
##' @title Parameter Names for the Model
##' @param family \code{"logistic"} (default) or \code{"poisson"}, the family for the likelihood
##' @param curve \code{"logistic"} or \code{"richards"}, the type of selection curve to be fitted
##' (will likely allow other options, for example "Bspline", future)
##' @param check.od logical, if \code{TRUE}, then overdispersion estimates will be produced
##' @param od logical, if \code{TRUE}, the model will be fit allowing for overdispersion
##' @param pred.check logical, if \code{TRUE}, then random observations will be drawn from the
##' posterior predictive distribution.
##' @param combine logical, if \code{TRUE}, then the "combined hauls" approach is used, otherwise a
##' hierarchical approach is used.
##' @param random vector of parameters to have hierarchical or random effects
##' @param L50 the formula for L50, default is L50 = ~haul
##' @param SR the formula for SR, default is SR = ~haul
##' @param phi the formula for phi, default is phi = ~haul
##' @param delta the formula for delta, default is delta = ~1
##' @param length.dist "iid" or "multinomial", the type of length distribution for the lambda
##' parameters. See details for more information.
##' @param paired whether or not the data are paired
##' @param ... additional arguments
##' @return lists of parameters
##' @author Tom Elliott
##' @export
getPars <- function(family = "binomial", curve = "logistic", check.od = FALSE, od = FALSE,
                    pred.check = FALSE,
                    combine = is.null(random), random = NULL, L50 = NULL, SR = NULL, phi = NULL,
                    delta = NULL, length.dist = "iid", paired = FALSE, ...) {

    allpars <- c()
    
    ## top level parameters:
    allpars <- c(allpars, "mu_L50", "mu_SR")
    if (curve == "richards")
        allpars <- c(allpars, "mu_delta", "log_mu_delta")
    if (paired)
        allpars <- c(allpars, "mu_phi")

    ## Random effects variances
    if (!is.null(random))
        allpars <- c(allpars,
                     paste0("sig2_", random))

    if (family == "poisson") {
        allpars <-
            switch(length.dist,
                   "iid" = c(allpars, "log_lambda", "lambda"),
                   "multinomial" = c(allpars, "Lambda", "pi", "alpha"))  # not correct probably.
    }

    if (check.od)
        allpars <- c(allpars, "od_obs", "od_exp", "p_od")
    if (od)
        allpars <- c(allpars, "sig2_od")
    if (pred.check & family == "binomial")
        allpars <- c(allpars, "yrep")

    ## Extra covariates
    if (!is.null(L50))
        allpars <- c(allpars, "beta")
    if (!is.null(SR))
        allpars <- c(allpars, "gamma")
    if (!is.null(phi))
        allpars <- c(allpars, "omega")
    if (!is.null(delta))
        allpars <- c(allpars, "zeta")

    
    

    all.summary.pars <-
        c("mu_L50", "mu_SR", "mu_delta", "mu_phi",
          "sig2_L50", "sig2_SR", "sig2_delta", "sig2_phi",
          "beta", "omega", "gamma", "zeta",
          "pi",
          "od_obs", "od_exp", "p_od", "sig2_od", "yrep")
    
    summary.pars <- all.summary.pars[all.summary.pars %in% allpars]

    list(summary.parameters = summary.pars,
         all.parameters = allpars)
}
