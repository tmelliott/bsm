##' Creates and fits a selection curve to the supplied \code{bsmdata} object. The correct data list
##' is created, the JAGS model constructed, and the model is run with default settings. It is also
##' possible to suppress the running of the JAGS and simply output the model and data.
##'
##' There are going to be a lot of details ...
##' 
##' @title Fit a Bayesian Selectivity Model
##' 
##' @param x an object of class \code{bsmdata}, from the \code{bsmData} function
##' 
##' @param family \code{"logistic"} (default) or \code{"poisson"}, the family for the likelihood
##' 
##' @param curve \code{"logistic"} or \code{"richards"}, the type of selection curve to be fitted
##' (will likely allow other options, for example "Bspline", future)
##' 
##' @param check.od logical, if \code{TRUE}, then overdispersion estimates will be produced
##' 
##' @param od logical, if \code{TRUE}, the model will be fit allowing for overdispersion
##' 
##' @param combine logical, if \code{TRUE}, then the 'combined hauls' approach is used, otherwise a
##' hierarchical approach is used.
##' 
##' @param random vector of parameters to have hierarchical or random effects
##' 
##' @param L50 the formula for L50, default is L50 = ~1
##' 
##' @param SR the formula for SR, default is SR = ~1
##' 
##' @param phi the formula for phi, default is phi = ~1
##' 
##' @param delta the formula for delta, default is delta = ~1
##' 
##' @param priors the prior distributions for specified parameters. See details.
##' 
##' @param inits initial values for parameters
##' 
##' @param length.dist 'iid' or "multinomial", the type of length distribution for the lambda
##' parameters. See details for more information.
##' 
##' @param parameters the parameters to save in the JAGS output. NOTE: only use this if you really
##' really only want these parameters - may cause errors in summary output and plots depending on
##' parameters selected.
##' 
##' @param extra.pars additional parameters to be tracked as well as \code{parameters}. Lets you
##' keep the default parameters.
##' 
##' @param file the file name to save the JAGS model in. If \code{NULL}, then a temporary file is
##' created
##' 
##' @param fit logical, if \code{TRUE} then the JAGS model is fit, otherwise only the model and data
##' are returned.
##' 
##' @param include.lambda logical, if \code{TRUE}, the lambda parameter will be saved as well
##' 
##' @param n.samples the total number of samples to obtain (excludes burnin and thin)
##' 
##' @param n.thin the thinning interval
##' 
##' @param n.burn the burn in period
##' 
##' @param n.chains the number of chains
##' 
##' @param max.attempts the maximum number of times bsm will attempt to fit the JAGS model (often
##' due to bad initial values, the model cannot be initiated and fails)
##' 
##' @param quiet suppress all output messages from JAGS?
##' 
##' @param progress.bar one of "none", "text" or "gui" for the JAGS progress bar
##' 
##' @param ... additional arguments
##' 
##' @return an object of class \code{bsmfit}
##' 
##' @author Tom Elliott
##' 
##' @export
bsm <- function(x, family = "binomial", curve = "logistic", check.od = TRUE, od = FALSE,
                combine = FALSE, random = NULL, L50 = ~1, SR = ~1,
                phi = if (attr(x, "paired")) ~1 else NULL,
                delta = if (curve == "richards") ~1 else NULL,
                priors = NULL, inits = NULL, length.dist = "iid",
                parameters = par$summary.parameters, extra.pars = NULL, file = NULL, fit = TRUE,
                include.lambda = FALSE,
                n.samples = 1000, n.thin = 1, n.burn = n.samples * n.thin, n.chains = 3,
                max.attempts = 10, quiet = FALSE, progress.bar = "text", ...) {

    if (class(x) != "bsmdata")
        stop("x must be a bsmdata object from the bsmData function")

    default.priors <- list(L50    = "dnorm(0, 1E-6)T(0,)",
                           SR     = "dlnorm(0, 1E-5)",
                           delta  = "dlnorm(0, 1E-3)",
                           phi    = "dnorm(0, 1E-3)")
#                           sig_od = "dgamma(1E-6, 1E-6)")
    if (!is.null(priors)) {
        priors <- c(priors, default.priors[!names(default.priors) %in% names(priors)])
    } else priors <- default.priors


    if (combine) {
        ## If we are to combine, we need to account for sampling fractions:
        ## We will recreate 'x' to remove haul etc etc etc

        df <- as.data.frame(x)
        y1 <- df$y1 / df$q1
        y2 <- df$y2 / df$q2

        lens <- unique(df$length)

        x <- bsmData(y1 = round(tapply(y1, df$length, sum)),
                     y2 = round(tapply(y2, df$length, sum)),
                     length = unique(df$length),
                     length.unit = attr(x, "length.unit"),
                     paired = attr(x, "paired"))
        combine <- FALSE
        ## Can't check overdispersion ...?
        ## check.od <- FALSE  # (for now ...)
    }

    
    
    
    ## Check for "haul" in the design

    ## Create design matrices for parameters:
    L50des <- SRdes <- phides <- deltades <- NULL
    if (!is.null(L50)) {
        if ("haul" %in% attr(terms(L50), "term.labels")) {
            if (!"L50" %in% random)
                random <- c(random, "L50")

            L50 <- update.formula(L50, ~.-haul)
        }
        L50des <- makeDesign(L50, as.data.frame(x))
    }
    if (!is.null(SR)) {
        if ("haul" %in% attr(terms(SR), "term.labels")) {
            if (!"SR" %in% random)
                random <- c(random, "SR")

            SR <- update.formula(SR, ~.-haul)
        }
        SRdes <- makeDesign(SR, as.data.frame(x))
    }
    if (!is.null(phi) & attr(x, "paired")) {
        if ("haul" %in% attr(terms(phi), "term.labels") & !"phi" %in% random)
            random <- c(random, "phi")      
        phides <- makeDesign(phi, as.data.frame(x))
    }
    if (!is.null(delta)) {
        if ("haul" %in% attr(terms(delta), "term.labels") & !"delta" %in% random)
            random <- c(random, "delta")      
        deltades <- makeDesign(delta, as.data.frame(x))
    }
    
    mod <- makeModel(family, curve, check.od, od, combine, random,
                     L50des, SRdes, phides, deltades, priors, length.dist,
                     paired = attr(x, "paired"), file)
    dat <- makeData(x, family, combine, random, L50des, SRdes, phides, deltades)
    par <- getPars(family, curve, check.od, od, combine, random,
                   L50des, SRdes, phides, deltades, length.dist,
                   paired = attr(x, "paired"))


    ## Some more convenient ways of including additional parameters:
    if (include.lambda & family == "poisson") {
        parameters <- c(parameters, "lambda")
    }
    if (!is.null(extra.pars)) {
        parameters <- c(parameters, extra.pars)
    }

    default.inits <- function()
        list(mu_L50 = rnorm(1, median(dat$x), diff(range(dat$x))/2),
             mu_SR = rlnorm(1, log(diff(range(dat$x))/3)))

    if (is.null(inits)) {
        final.inits <- default.inits
    } else if (is.function(inits)) {
        ## Combine output of two functions ......??
        final.inits <- function() {
            ini <- inits()
            def.ini <- default.inits()
            c(ini, def.ini[!names(def.ini) %in% names(ini)])
        }
    } else {
        if (!all(sapply(inits, class) == "list"))
            stop("Inits needs to be a list of lists (one per chain), or a function returning a list. See ?jags")

        if (length(inits) < n.chains)
            stop(sprintf("Inits needs to be a list of %s lists. See ?jags", n.chains))
        
        final.inits <- lapply(1:n.chains, function(i) {
            ini <- inits[[1]]
            def.ini <- default.inits()
            c(ini, def.ini[!names(def.ini) %in% names(ini)])
        })
    }
    
    ## Fit the model!!
    out <- list()
    if (fit) {
        for (i in 1:max.attempts) {
            if (quiet)
                capture.output(suppressMessages(
                    j <- try(R2jags::jags(dat, final.inits, parameters, mod,
                                          n.chains = n.chains, n.iter = n.burn + n.samples * n.thin,
                                          n.burnin = n.burn, n.thin = n.thin,
                                          progress.bar = progress.bar),
                             TRUE)),
                               file = tempfile()
                               )
            else
                j <- try(R2jags::jags(dat, final.inits, parameters, mod,
                                      n.chains = n.chains, n.iter = n.burn + n.samples * n.thin,
                                      n.burnin = n.burn, n.thin = n.thin,
                                      progress.bar = progress.bar),
                         TRUE)
            
            if (!inherits(j, "try-error"))
                break
        }
        if (inherits(j, "try-error")) {
            cat("\n", paste0(rep("=", options()$width), collapse=""), "\n")
            warning(paste0(
                "The following error occured while attempting to fit the model:\n\n",
                as.character(j),
                "  You may need to specify initial values. See details in ?bsm for help."))
        } else {
            out$fit <- j
        }
    }

    bugs <- out$fit$BUGSoutput
    ordL <- bugs$long.short
    names(ordL) <- pp <- bugs$root.short
    indL <- lapply(bugs$indexes.short, function(i) if (is.null(i)) 0 else unlist(i))
    pp <- rep(pp, sapply(ordL, length))
    ord <- unlist(sapply(parameters, function(x) ordL[[x]]))
    ind <- unlist(indL, use.names = FALSE)[ord]
    indo <- ifelse(ind == 0, "", paste0("[", ind, "]"))
    all.parameters <- paste0(pp[ord], indo)
    
    out <- c(out, list(data = dat, file = mod,
                       summary.parameters = parameters,
                       all.parameters = all.parameters,
                       object = x, combine = combine,
                       L50 = L50des, SR = SRdes, phi = phides, delta = deltades,
                       L50f = L50, SRf = SR, phif = phi, deltaf = delta,
                       family = family, curve = curve, check.od = check.od, od = od))
    
    class(out) <- "bsmfit"
    out
}


##' @param model logical, if \code{TRUE}, then the JAGS model is printed, otherwise summary
##' information is provided
##' @param coda logical, if \code{TRUE}, then the summary from the \code{coda} package is used
##' 
##' @describeIn bsm Print the summary of the model, or the model code
##' @export
print.bsmfit <- function(x, model = FALSE, coda = FALSE, ...) {
    if (model) {
        cat(readLines(x$file), sep = "\n")
    } else if ("fit" %in% names(x)) {
        if (coda) {
            ## Maintain the ordering of the parameters
            pp <- rownames(x$fit$BUGSoutput$summary)
            print(summary(coda::as.mcmc(x$fit)[, pp]), ...)
        } else {
            print(x$fit, ...)
        }
    } else {
        cat("\n")
        cat(" JAGS model written to \"", x$file, "\"\n", sep = "")
        cat(" Data can be accessed by `$data`\n\n")
    }

    return(invisible(NULL))
}


##' @param object a bsmfit object
##' @param p.values logical, include significance tests for parameters != 0?
##' @param formula.values logical, display values in formulae?
##' @param predict.values values to be included in the prediction
##' 
##' @describeIn bsm Generate summary output for a bsmfit object
##' @export
summary.bsmfit <- function(object, p.values = FALSE, formula.values = FALSE,
                           predict.values = NULL, ...) {
    x <- object
    fit <- x$fit

    ## In several parts - first, create summary "list" with a class - then give this a print method


    ## 1: Summary information of the model.
    out <- list()
    out$type <- ifelse(attr(x$object, "paired"), "paired", "covered-codend")
    out$family <- x$family
    out$curve <- x$curve

    bugs <- fit$BUGSoutput
    ordL <- bugs$long.short
    names(ordL) <- pp <- bugs$root.short
    indL <- lapply(bugs$indexes.short, function(i) if (is.null(i)) 0 else unlist(i))
    pp <- rep(pp, sapply(ordL, length))
    
    ord <- c(ordL$mu_L50, ordL$mu_SR)

    if ("mu_phi" %in% pp)
        ord <- c(ord, ordL$mu_phi)
    if ("mu_delta" %in% pp)
        ord <- c(ord, ordL$mu_delta)
    
    if ("sig2_L50" %in% pp)
        ord <- c(ord, ordL$sig2_L50)
    if ("sig2_SR" %in% pp)
        ord <- c(ord, ordL$sig2_SR)
    if ("sig2_phi" %in% pp)
        ord <- c(ord, ordL$sig2_phi)
    if ("sig2_delta" %in% pp)
        ord <- c(ord, ordL$sig2_delta)
    
    if ("beta" %in% pp)
        ord <- c(ord, ordL$beta)
    if ("gamma" %in% pp)
        ord <- c(ord, ordL$gamma)
    if ("omega" %in% pp)
        ord <- c(ord, ordL$omega)
    if ("zeta" %in% pp)
        ord <- c(ord, ordL$zeta)

    mat <- bugs$summary[ord, c("mean", "sd", "2.5%", "50%", "97.5%", "Rhat")]

    ind <- unlist(indL, use.names = FALSE)[ord]
    indo <- ifelse(ind == 0, "", paste0("[", ind, "]"))
    ppo <- pp[ord]

    rownames(mat) <- paste(ppo, indo, sep = "")
    out$summary.matrix <- mat

    if (!is.null(x$L50))
        out$L50f <- generateFormula("L50", x$L50f, x$object,
                                    mat[, "mean", drop = FALSE], getMean(x, "L50"))
    if (!is.null(x$SR))
        out$SRf <- generateFormula("SR", x$SRf, x$object,
                                   mat[, "mean", drop = FALSE], getMean(x, "SR"))
    if (!is.null(x$phi))
        out$phif <- generateFormula("phi", x$phif, x$object,
                                    mat[, "mean", drop = FALSE], getMean(x, "phi"))
    if (!is.null(x$delta))
        out$deltaf <- generateFormula("delta", x$deltaf, x$object,
                                      mat[, "mean", drop = FALSE], getMean(x, "delta"))


    ## Generate predictions for each 'curve'

    ## Have to update variables to take away center
    if (!is.null(predict.values)) {
        predict.values <- list(original = predict.values,
                               centers =
                                   lapply(pred.names <- names(predict.values), function(var) {
                                       obj.df <- as.data.frame(x$object)
                                       if (var %in% colnames(obj.df)) {
                                           center <- mean(obj.df[, var])
                                       } else {
                                           warning(paste0(var, " is not a known variable."))
                                           center <- 0
                                       }
                                       center
                                   }))
        names(predict.values$centers) <- pred.names
    }
    
    if (!is.null(x$L50))
        L50mat <- predictPar(out$L50f, predict.values = predict.values)
    if (!is.null(x$SR))
        SRmat <- predictPar(out$SRf, predict.values = predict.values)
    if (!is.null(x$phi))
        phimat <- predictPar(out$phif, predict.values = predict.values)
    if (!is.null(x$delta))
        deltamat <- predictPar(out$deltaf, predict.values = predict.values)

    ## Prediction matrix:
    if (!is.null(x$L50)) {
        pred.mat <- L50mat
    } else {
        pred.mat <- matrix(mat["mu_L50", "mean"], dimnames = list(1, "L50"))
    }
    if (!is.null(x$SR)) {
        pred.mat <-
            if (is.null(pred.mat)) SRmat
            else merge(pred.mat, SRmat)
    } else {
        pred.mat <-
            if (is.null(pred.mat)) matrix(mat["mu_SR", "mean"], dimnames = list(1, "SR"))
            else cbind(pred.mat, SR = mat["mu_SR", "mean"])
    }
    if (!is.null(x$phi)) {
        pred.mat <-
            if (is.null(pred.mat)) phimat
            else merge(pred.mat, phimat)
    } else if (attr(x$object, "paired")) {
        pred.mat <-
            if (is.null(pred.mat)) matrix(mat["mu_phi", "mean"], dimnames = list(1, "phi"))
            else cbind(pred.mat, phi = mat["mu_phi", "mean"])
    }
    if (!is.null(x$delta)) {
        pred.mat <-
            if (is.null(pred.mat)) deltamat
            else merge(pred.mat, deltamat)
    } else if (x$curve == "richards") {
        pred.mat <-
            if (is.null(pred.mat)) matrix(mat["mu_delta", "mean"], dimnames = list(1, "delta"))
            else cbind(pred.mat, delta = mat["mu_delta", "mean"])
    }
    
    which <- colnames(pred.mat) %in% c("L50", "SR", "phi", "delta")
    pred.response <- pred.mat[, which, drop = FALSE]
    pred.explanatory <- pred.mat[, !which, drop = FALSE]
    
    predict <- cbind(pred.explanatory, pred.response)
    
    predict.df <- as.data.frame(predict)
    predict.df <- do.call(data.frame, lapply(predict.df, function(x) {
        tmp <- try(as.numeric(x))
        if (inherits(tmp, "try-error"))
            x
        else
            tmp
    }))
    
    dropCols <- apply(predict.df, 2, function(x) if (is.numeric(x)) all(x == 0) else FALSE)               
    out$predict <- as.data.frame(predict[, !dropCols, drop = FALSE])

    

    if (x$check.od) {
        odmat <- bugs$summary[c("od_exp", "od_obs", "p_od"), "mean"]
        
        out$overdispersion <- list(
            check = list(
                expected.chisquare = odmat[1],
                observed.chisquare = odmat[2],
                p.overdispersion   = odmat[3]
                )
            )
    }
    if (x$od) {
        if (!x$check.od)
            out$overdispersion <- list()
        out$overdispersion$variance = bugs$summary[c("sig2_od"), "mean"]
    }

    
    ## "significance" tests
    if (p.values)
        out$p.values <- apply(bugs$sims.matrix[, rownames(mat)], 2,
                              function(x) {
                                  m <- mean(x < 0)
                                  min(m, 1 - m)
                              })
        
    ## MCMC info
    out$mcmc <- with(x$fit$BUGSoutput,
                     list(
                         chains     = n.chains,
                         samples    = n.keep,
                         thin       = n.thin,
                         burn       = n.burnin,
                         iterations = n.iter,
                         total      = n.sims,
                         DIC        = DIC,
                         pD         = pD
                         ))

    
    
    out$gelman <- coda::gelman.diag(as.mcmc(x))
    out$geweke <- coda::geweke.diag(as.mcmc(x))
    out$heidel <- coda::heidel.diag(as.mcmc(x))
    out$formula.values <- formula.values
    
    class(out) <- "summary.bsmfit"
    out
}


##' @export
print.summary.bsmfit <- function(x, ...) {
    cat("\n")
    cat(sprintf("Summary of model fitted to %s experimental trawl using JAGS.\n\n", x$type))

    cat("  Likelihood: ", x$family, "\n", sep = "")
    cat("  Curve: ", x$curve, "\n\n", sep = "")
    

    cat("Summary of posterior distribution for the model:\n\n")
    
    cat(apply(cbind("",
                    format(c("", rownames(x$summary.matrix)), justify = "right"),
                    do.call(cbind, lapply(colnames(x$summary.matrix), function(v) {
                        format(c(v, format(x$summary.matrix[, v], digits = 4)), justify = "right")
                    })),
                    if (!is.null(x$p.values)) {
                        p <- x$p.values
                        format(c("Pr(x=0)",
                                 ifelse(p < 2e-16, "<2e-16", signif(p, 2))),
                               justify = "right")
                    } else { "" },
                    "\n"),
              1, paste, collapse = "   "), sep = "")
    

    if (any(NL <- sapply(x[c("L50f", "SRf", "phif", "deltaf")], function(x) !is.null(x))))
        cat(sprintf("\nFormula%s for selection curve parameter%s:\n",
                    ifelse(sum(NL) > 1, "e", ""), ifelse(sum(NL) > 1, "s", "")))
    
    if (!is.null(x$L50f))
        print(x$L50f, use.values = x$formula.values)
    if (!is.null(x$SRf))
        print(x$SRf, use.values = x$formula.values)
    if (!is.null(x$phif))
        print(x$phif, use.values = x$formula.values)
    if (!is.null(x$deltaf))
        print(x$deltaf, use.values = x$formula.values)
    
    cat("\n")

    if (!is.null(x$overdispersion)) {
        cat("Overdispersion information:\n")

        od <- x$overdispersion

        if (exists("variance", od)) {
            cat("The estimated (mean) overdispersion parameter is:",
                format(od$variance, digits = 3), "\n")
        }
        
        if (exists("check", od)) {
            cat("\n")
            pval <- unlist(od$check)[3]
            cat(apply(cbind("  ",
                            format(c(paste(format(c("Expected", "Observed"), justify = "right"),
                                           "Pearson Chi-square:"),
                                     "Pr[expected < observed]:"), justify = "right"),
                            format(c(format(unlist(od$check)[1:2], digits = 4),
                                     ifelse(pval < 1e-12, "< 1e-12", signif(pval, 2))),
                                   justify = "right"),
                            "\n"),
                      1, paste, collapse = " "), sep = "")
        }
        cat("\n")
    }

    cat("MCMC Information:\n")
    mcmc <- x$mcmc
    cat(apply(cbind("",
                    format(c("Chains", mcmc$chains), justify = "right"),
                    format(c("Iterations/chain", mcmc$iterations), justify = "right"),
                    format(c("Burn-in", mcmc$burn), justify = "right"),
                    format(c("Thin", mcmc$thin), justify = "right"),
                    format(c("Total", mcmc$total), justify = "right"),
                    "\n"),
              1, paste, collapse = "  "), sep = "")

    cat(sprintf("\n  Model DIC = %s, pD = %s\n",
                format(mcmc$DIC, digits = 5), format(mcmc$pD, digits = 5)))
    cat(sprintf("  Gelman and Rubin Convergence Statistic (multivariate) = %s\n",
                format(x$gelman$mpsrf, digits = 4)))
    ## check if Heidelberger says it's stationary:
    heidel <- x$heidel
    stest <- all(sapply(heidel, function(x) {
        X <- x[, "stest"]
        all(X[!is.na(X)] == 1)        
    }))
    htest <- all(sapply(heidel, function(x) {
        X <- x[, "htest"]
        all(X[!is.na(X)] == 1)        
    }))

    if (stest & htest) {
        cat("  Heidelberger & Welch's Stationarity and Half-width Tests passed.\n")
    } else {
        if (!stest & !htest)
            txt <- "Heidelberger & Welch's Stationarity and Half-width Tests failed."
        else if (!stest)
            txt <- "Heidelberger & Welch's Stationarity Test failed."
        else if (!htest)
            txt <- "Heidelberger & Welch's Half-width Test failed."

        cat("\nWarning message:\n  ", txt, "\n  You should try increasing the chain length.", sep = "")
    }


    cat("\n")
}


##' .. content for description{} (no empty lines) ..
##'
##' .. content for details{} ..
##' @title Predict parameter values for selection curve
##' @param object bsmfit or summary.bsmfit object
##' @param sort character vector of variables to sort by
##' @param predict.values values to be included in the prediction
##' @param ... extra arguments
##' @return matrix
##' @author Tom Elliott
##' @export
predict.bsmfit <- function(object, sort = NULL, predict.values = NULL, ...) {
    predict(summary(object, predict.values = predict.values), sort = sort, ...)
}

##' @method predict summary.bsmfit
##' @rdname predict.bsmfit
##' @export
predict.summary.bsmfit <- function(object, sort = NULL, ...) {
    df <- object$predict
    if (!is.null(sort))
        df <- df[do.call(order, df[, sort, drop = FALSE]), ]
    
    class(df) <- c("bsm.predict", "data.frame")
    df
}

##' @export
print.bsm.predict <- function(x, ...) {
    m <- do.call(cbind, lapply(colnames(x), function(i) {
        if (i %in% c("L50", "SR", "phi", "delta"))
            format(c(i, format(x[, i], digits = 4)), justify = "right")
        else
            format(c(i, as.character(x[, i])), justify = "right")
    }))
    apply(cbind("", m, "\n"), 1, cat, sep = "   ")
}










##' @export
as.mcmc.bsmfit <- function(x) {
    drop <- c("p_od")
    coda::as.mcmc(x$fit)[, c(x$all.parameters[!x$all.parameters %in% drop], "deviance")]
}



##' Extract the DIC value from an object
##'
##' Generic Method
##' @title Extract DIC Value
##' @param x an object, such as bsmfit 
##' @param ... extra arguments
##' @return the DIC summary information
##' @author Tom Elliott
##' @export
dic <- function(x, ...)
    UseMethod("dic")


##' @describeIn dic Extract the DIC information from the JAGS object
##' @export
dic.bsmfit <- function(x, ...) {
    dd <- list(...)
    if (length(dd) > 0)
        xL <- c(list(x), dd)
    else
        xL <- list(x)

    res <- lapply(xL, function(y)
                  with(summary(y)$mcmc, structure(c(DIC, pD), .Names = c("DIC", "pD"))))

    names(res) <- as.character(match.call())[-1]

    class(res) <- "dic.summary"
    res
}


##' @export
print.dic.summary <- function(x, ...) {
    m <- do.call(rbind, x)
    apply(cbind("",
                format(c("Model", rownames(m)), justify = "right"),
                format(c("DIC", format(m[, "DIC"], digits = 4)), justify = "right"),
                format(c("pD", format(m[, "pD"], digits = 2)), justify = "right"),
                "\n"), 1, cat, sep = "   ")
}

##' @export
`[.dic.summary` <- function(x, i, j) {
    mat <- do.call(rbind, x)
    mat[i, j]
}

##' @export
`$.dic.summary` <- function(x, name) {
    x[, name]
}
