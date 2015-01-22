##' A plot method to easily obtain a graphical representation of the data, with several options to
##' enable coding of some variables. For example, hauls can be colour coded, and weighting variables
##' used to control point size.
##'
##' A lot of the models used are hierarchical, with random effects on hauls. Therefore, to see the
##' variation in hauls graphically, you can pass \code{col = "haul"} and the proportions will be
##' coloured differently for each haul.
##'
##' The same idea is applied to plotting character (\code{pch}). To use a different plotting
##' character for each covariate (e.g., mesh size), simply pass \code{pch = mesh} into the
##' function. However, note that using \code{pch} and \code{col} will not be easy to see trends
##' unless they are very clear. It is recommended that you use both symbol and colour only if there
##' are individual points you wish to investigate, rather than for finding trends, which is best
##' done using multiple plots and colour.
##' @title Plot Method for bsmdata Object
##' @param x a \code{bsmdata} object
##' @param scale logical, if \code{TRUE}, the data is scaled by sampling fractions
##' @param col the color of points. This can also be used to colour points by a factor variable such
##' as \code{haul}. See details.
##' @param pch the plotting symbol. This can also be the name of a (factor) variable in the data set
##' @param legend logical, if \code{TRUE} then a legend is drawn if points are coloured by some
##' variable (eg. haul)
##' @param weight the weights that will decide bubble sizes. Can be a variable or a vector
##' @param leg.posx position of legend. This can be text such as "topleft", "bottomright" etc, or
##' specific \code{x}-coordinates
##' @param leg.posy if \code{leg.posx} is an x-value, then this is the y-value
##' @param leg.cex size of legend
##' @param leg.bty the box type for the legend. Defaults to \code{"n"} for no box
##' @param col.palette the color palette to use by default; can be \code{rainbow} or a vector of colors
##' @param ... additional arguments that will be passed to the \code{plot.default} function
##' @return NULL
##' @author Tom Elliott
##' @export
plot.bsmdata <- function(x, scale = TRUE, col, pch, legend = FALSE, weight,
                         leg.posx = "topleft", leg.posy = NULL, leg.cex = 0.7, leg.bty = "n",
                         col.palette = "rainbow",
                         ...) {
    ## grab some variables from x (bsmdata object)
    paired <- attr(x, "paired")
    length.unit <- attr(x, "length.unit")
    nhaul <- attr(x, "nhaul")

    dots <- list(...)

    ## look into the call and save the arguments+values
    mc <- as.list(match.call()[-1])

    ## examine the `col` argument and decide if it is missing, a colour, or a variable to colour by
    if (missing(col))
        col <- "#000000"
    else
        col <- eval(mc$col, x)
    
    if (length(col) == 1) {
        bycol <- col %in% names(x)
        
        if (bycol)
            colid <- as.factor(x[[col]])      
    } else {
        bycol <- length(col) == length(x$y1)
        colid <- as.factor(col)
    }

    ## examine the `pch` argument and do the same as for `col`
    if (missing(pch))
        pch <- 1
    else
        pch <- eval(mc$pch, x)

    if (length(pch) == 1) {
        bypch <- pch %in% names(x)
        
        if (bypch)
            pchid <- as.factor(x[[pch]])      
    } else {
        bypch <- length(pch) == length(x$y1)
        pchid <- as.factor(pch)
    }

    ## figure out the name for the colour and pch variables
    colname <- pchname <- ""
    if (bycol)
        colname <- names(x)[sapply(x, function(x) {
            if (is.factor(x)) {
                if (length(levels(x)) == length(levels(colid)))
                    if (all(levels(x) == levels(colid)))
                        all(x == as.factor(colid))
                    else
                        FALSE
                else
                    FALSE
            } else {
                all(x == colid)
            }
        })]
    if (bypch)
        pchname <- names(x)[sapply(x, function(x) all(x == pchid))] 

    ## implement appropriate scaling
    Y1 <- with(x, if (scale) y1 / q1 else y1)
    Y2 <- with(x, if (scale) y2 / q2 else y2)
    N <- Y1 + Y2

    ## generate the colours and plotting symbols 
    if (bycol) {
        ncol <- length(levels(colid))
        if (length(col.palette == 1) & col.palette[1] == "rainbow")
            cols <- rainbow(ncol, start = 0/6, end = 5/6,
                            s = 0.9, v = 0.9, alpha = 0.6)
        else
            cols <- rep(col.palette, ncol)
    } else {
        cols <- col
        colid <- 1
    }

    if (bypch) {
        pchs <-
            if (bycol & length(levels(pchid)) <= 4) c(19, 17, 15, 18)
            else 1:length(levels(pchid))
    } else {
        pchs <- if (bycol & pch == 1) 19 else pch
        pchid <- 1
    }
    npch <- length(unique(pchs))

    ## figure out point sizes
    if (missing(weight)) {
        size <- 1
    } else {
        w <- eval(mc$weight, x)
        ## if specifying N as the weighting variable, sum over lengths
        if (all(w == N))
            w <- (tapply(N, x$haul, sum))[x$haul]
        
        size <- 4 * (w - min(w)) / diff(range(w)) + 0.5
        if (length(pchs == 1))
            if (pchs == 1)
                pchs <- 19
    }
    if ("cex" %in% names(dots))
        size <- size * dots$cex

    ## draw the plot
    with(x, plot(length, Y1 / N,
                 xlab = paste0("Length",
                     ifelse(is.null(length.unit), "", paste0(" (", length.unit, ")"))),
                 ylab = paste0(
                     ifelse(scale, "Scaled p", "P"),
                     "roportion caught in ",
                     ifelse(paired, "experimental trawl", "codend")
                     ),
                 ylim = 0:1, col = cols[as.numeric(colid)], pch = pchs[as.numeric(pchid)],
                 cex = size))

    ## draw the legend if it's asked for
    if (legend && (bycol | bypch)) {
        LEG <- COL <- character()
        PCH <- numeric()
        
        if (bycol) {
            LEG <- c(LEG, paste(colname, levels(colid)))
            COL <- c(COL, cols)
            PCH <- c(PCH, rep(19, ncol))
        }
        if (bypch) {
            LEG <- c(LEG, paste(pchname, levels(pchid)))
            COL <- c(COL, rep("#000000", npch))
            PCH <- c(PCH, pchs)
        }
        
        legend(leg.posx, leg.posy, LEG, pch = PCH, col = COL,
               bty = leg.bty, cex = leg.cex)
    }

    return(invisible(NULL))
}




##' Various types of plots of the fitted BSM object
##'
##' Lots ??
##' @title Plot a BSM Fit
##' @param x a bsmfit object
##' @param which which plot to draw, c("posterior", "curve")
##' @param parameters if \code{which = "posterior"}, \code{"pairs"}, then these parameters are used
##' @param estimate "mean" or "median", the estimate to be used when drawing the curve
##' @param cred.int logical, if \code{TRUE}, then a credible interval is calcualted from the
##' posterior samples of L50 and SR (and delta)
##' @param cred.alpha the significance level for the credible interval, default is 0.95
##' @param cred.col the colour of the credible interval
##' @param new logical, if \code{TRUE}, then a new plot of the data is drawn, otherwise the curve is
##' plot over any existing one
##' @param col colour of lines
##' @param lty line type
##' @param lwd line width
##' @param legend.order order of variables for legend purposes (color, line type, line width)
##' @param predict.values extra values to predict
##' @param n.posterior.rows the number of rows in the posterior plot
##' @param leg.posx x-position of legend, can be "topleft" etc
##' @param leg.posy NULL, or y-position of legend
##' @param leg.cex size of legend
##' @param leg.bty box type of legend
##' @param chain.cols colours of lines for traceplots and density plots, NOTE: not used for pairs
##' and curve
##' @param interactive logical, if \code{TRUE}, use clicking on the plot to navigate through plots
##' @param ... additional parameters to \code{plot}
##' @return NULL
##' @author Tom Elliott
##' @export
plot.bsmfit <- function(x, which = "posterior", parameters = NULL,
                        estimate = "mean", cred.int = FALSE, cred.alpha = 0.95, cred.col = "#99999950",
                        new = TRUE, col = NULL, lty = 1, lwd = 2, legend.order = NULL,
                        predict.values = NULL, n.posterior.rows = 3,
                        leg.posx = "topleft", leg.posy = NULL, leg.cex = 0.7, leg.bty = "n",
                        chain.cols = rainbow(x$fit$BUGSoutput$n.chains + 1, s = 0.8, v = 0.8),
                        interactive = FALSE,
                        ...) {
    if (length(which) > 1)
        warning("`which` can only be a single value, only the first being used")

    if (is.null(parameters)) {
        all.par <- c("mu_L50", "sig2_L50", "mu_SR", "sig2_SR",
                     "mu_phi", "sig2_phi", "mu_delta", "sig2_delta")
        
        parameters <- all.par[all.par %in% (ap <- x$fit$BUGSoutput$root.short)]
        mn <- colnames(x$fit$BUGSoutput$sims.matrix)
        for (coef in c("beta", "gamma", "omega", "zeta"))
            if (coef %in% ap)
                parameters <- c(parameters, mn[grep(coef, mn)])
        if (x$check.od)
            parameters <- c(parameters, "od_exp", "od_obs")
        if (x$od)
            parameters <- c(parameters, "sig2_od")

        parameters <- c(parameters, "deviance")
##        parameters <- x$all.parameters
    }

    dots <- list(...)

    switch(which,
           "posterior" = {
               m <- coda::as.mcmc(x$fit)[, parameters]
               nr <- np <- length(parameters)
               Nit <- x$fit$BUGSoutput$n.keep
               if (np > n.posterior.rows & !interactive)                   
                   devAskNewPage(TRUE)
               
               nr <- n.posterior.rows
               op <- par(mfrow = c(nr, 2))
               theplot <- function(i) {
                   coda::traceplot(mi <- m[, i], col = chain.cols, 
                                   main = paste0("Trace of ", parameters[i])) 
                   dens <- lapply(mi, density)
                   plot.new()
                   plot.window(ylim = c(0, max(sapply(dens, function(x) max(x$y)))),
                               xlim = range(sapply(dens, function(x) range(x$x))))
                   box(); axis(1); axis(2)
                   for (j in 1:length(dens))
                       lines(dens[[j]], col = chain.cols[j])
                   title(main = paste0("Density of ", parameters[i]),
                         xlab = paste0("N = ", Nit),
                         ylab = "Density")
               }
               
               if (interactive) {
                   cat("\n")
                   cat("  Click top-right for next and top-left for previous.\n",
                       "  Right-click to exit, or click the bottom of the graphics window.\n\n", sep = "")
                   I <- 1
                   while (I > 0) {
                       ow <- par(mfrow = c(nr, 2))
                       
                       ii <- 3 * (I - 1) + (1:3)
                           
                       for (i in ii[ii <= np]) {
                           dev.hold()
                           theplot(i)
                           dev.flush()
                       }
                       
                       ou <- par(usr = c(0, 1, 0, 1))
                       
                       xy <- locator(1)
                       if (is.null(xy))
                           I <- 0
                       else if (xy$y < 0.2)
                           I <- 0
                       else if (xy$x > 0.7)
                           I <- ifelse(I == ceiling(np / 3), I, I + 1)
                       else if (xy$x < 0.3)
                           I <- ifelse(I == 1, I, I - 1)
                       par(ow)
                   }
               } else {
                   for (i in 1:np)
                       theplot(i)
               }
               par(op)
               devAskNewPage(FALSE)
           },
           "pairs" = {
               s <- x$fit$BUGSoutput$sims.matrix
               s <- s[, parameters]# %in% colnames(s)]
               pairs(s, panel = function(x, y, ...) points(x, y, col = "#0055bb30", pch = 19))
           },
           "curve" = {
               if (new)
                   plot(x$object)
               
               if (estimate == "median") estimate <- "50%"
               s <- x$fit$BUGSoutput$summary[, estimate]
               s <- s[names(s) %in% c("mu_L50", "mu_SR", "mu_delta", "mu_phi")]
               
               predmat <- predict(x, predict.values = predict.values,
                                  sort = legend.order[legend.order != "null"])

               if (nrow(predmat) == 1) {
                   ## If they ask for a confidence interval:
                   if (cred.int) {
                       if (cred.alpha < 0 | cred.alpha > 1) {
                           warning("cred.alpha must be between 0 and 1. Using default = 0.95")
                           cred.alpha <- 0.95
                       }
                       
                       post <- x$fit$BUGSoutput$sims.matrix[, names(s)]
                       if (!"mu_delta" %in% (cc <- colnames(post))) {
                           post <- cbind(post, 1)
                           colnames(post) <- c(cc, "mu_delta")
                       }
                       if (!"mu_phi" %in% (cc <- colnames(post))) {
                           post <- cbind(post, 100)
                           colnames(post) <- c(cc, "mu_phi")
                       }
                       
                       xl <- par()$usr[1:2]
                       xx <- seq(xl[1], xl[2], length = length(x$data$x) * 2)
                       curves <- apply(post, 1, function(par) {
                           eta <- richards(xx, par["mu_L50"], par["mu_SR"], par["mu_delta"])
                           yy <- (1 / (1 + exp(-eta))) ^ (1 / par["mu_delta"])
                           yy * ilogit(par["mu_phi"])
                       })
                       
                       ci <- 0.5 + cred.alpha / 2 * c(-1, 1)
                       qx <- apply(curves, 1, quantile, probs = ci)
                       polygon(c(xx, rev(xx)), c(qx[1, ], rev(qx[2, ])),
                               col = cred.col, lty = 3)
                   }

                   names(s) <- gsub("mu_", "", names(s))
                   if (is.null(col))
                       col <- "black"
                   bsmCurve(s, col = col, lwd = lwd, lty = lty, ...)
               } else {           
                   p <- c("L50", "SR", "delta", "phi")
                   resp <- p[p %in% colnames(predmat)]

                   expl <- colnames(predmat)[!colnames(predmat) %in% resp]
                   if (!is.null(legend.order))
                       expl <- legend.order
                   
                   nr <- nrow(predmat)
                   
                   LEG.LAB = character()
                   LEG.COL = character()
                   LEG.LTY = numeric()
                   LEG.LWD = numeric()
                   
                   COLS <- rep(ifelse(is.null(col), "#000000", col[1]), nr)
                   if (length(expl) > 0) {
                       if (expl[1] != "null") {
                           v1 <- as.integer(as.factor(predmat[[expl[1]]]))
                           
                           if (is.null(col) | length(col) < length(unique(v1))) {
                               COL <- rainbow(length(unique(v1)), start = 0/6, end = 5/6,
                                              s = 0.8, v = 0.8)

                               LEG.LAB <- c(LEG.LAB, paste(expl[1], "=",
                                                           levels(as.factor(predmat[[expl[1]]]))))
                               LEG.COL <- c(LEG.COL, COL)
                               LEG.LTY <- c(LEG.LTY, rep(1, length(COL)))
                               LEG.LWD <- c(LEG.LWD, rep(2, length(COL)))
                           } else {
                               COL <- col
                           }
                           
                           COLS <- COL[v1]
                       }
                   }

                   
                   if (length(expl) > 1) {
                       v2 <- as.integer(as.factor(predmat[[expl[2]]]))
                       LTY <- v2
                       
                       LEG.LAB <- c(LEG.LAB, paste(expl[2], "=",
                                                   levels(as.factor(predmat[[expl[2]]]))))
                       LEG.COL <- c(LEG.COL, rep("#000000", length(unique(LTY))))
                       LEG.LTY <- c(LEG.LTY, unique(LTY))
                       LEG.LWD <- c(LEG.LWD, rep(1, length(unique(LTY))))
                   } else {
                       LTY <- rep(lty[1], nr)
                   }
                   
                   if (length(expl) > 2) {
                       v3 <- as.integer(as.factor(predmat[[expl[3]]]))
                       LWD <- v3
                       
                       LEG.LAB <- c(LEG.LAB, paste(expl[3], "=",
                                                   levels(as.factor(predmat[[expl[3]]]))))
                       LEG.COL <- c(LEG.COL, rep("#000000", length(unique(LWD))))
                       LEG.LTY <- c(LEG.LTY, rep(1, length(unique(LWD))))
                       LEG.LWD <- c(LEG.LWD, unique(LWD))
                   } else {
                       LWD <- rep(lwd[1], nr)
                   }
                   
                   for (i in 1:nr) {
                       c <- as.numeric(predmat[i, resp])
                       names(c) <- resp
                       bsmCurve(c, col = COLS[i], lty = LTY[i], lwd = LWD[i], ...)
                   }
                   if (length(LEG.LAB) > 0)
                       legend(leg.posx, leg.posy, LEG.LAB, col = LEG.COL, lty = LEG.LTY, lwd = LEG.LWD,
                              bty = leg.bty, cex = leg.cex)
               }
           })
    
    return(invisible(NULL))
}

