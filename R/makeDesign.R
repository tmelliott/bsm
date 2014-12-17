makeDesign <- function(par, dat) {
    # par <- update.formula(par, ~. - 1)  # remove the intercept
    mf <- model.matrix(par, dat)[, -1, drop = FALSE]
    
    des <- as.matrix(mf[tapply(1:nrow(mf), dat$haul,
                               function(x) x[1]), colnames(mf) != "haul"])

    ## Center covariates?
    terms <- attr(terms(par), "term.labels")
    terms <- terms[terms != "haul"]
    if (length(terms) > 0) {
        isnum <- sapply(terms, function(p) {
            is.numeric(eval(parse(text = p), dat))
        }, USE.NAMES = FALSE)
                
        for (i in which(isnum))
            des[, i] <- des[, i, drop = FALSE] - mean(des[, i])
    }
    
    if (all(dim(des)) > 0)
        des
    else
        NULL
}
