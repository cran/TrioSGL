
TrioSGL <- function(X, index = NULL, maxit = 10000, thresh = 0.0001, min.frac = 0.01, nlam = 25, lambdas = NULL, alpha = 0.95, gamma = 0.8, step = 1, reset = 20, standardize = FALSE, verbose = FALSE) {

    ## sanity check ##
    if (nrow(X) %% 4 != 0) {
        stop("Number of rows is not a multiple of 4.")
    }

    if (is.null(index)) {
        message("Note: No groups specified. Putting each element into its own group.")
        index <- 1:ncol(X)
    }

    ## Aus dem original SGL übernommen ##
    X.transform <- NULL
    if (standardize == TRUE) {
        means <- apply(X, 2, mean)
        X <- t(t(X) - means)
        var <- apply(X, 2, function(x)(sqrt(sum(x^2))))
        X <- t(t(X) / var)
        X.transform <- list(X.scale = var, X.means = means)
    }

    Sol <- trioFit(X, index, inner.iter = maxit, outer.iter = maxit, thresh = thresh, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, alpha = alpha, gamma = gamma, step = step, reset = reset, verbose = verbose)

    Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, X.transform = X.transform)

    class(Sol) = "TrioSGL"
    return(Sol)
}


trioFit <- function(X, index, inner.iter = 10000, outer.iter = 10000, thresh = 0.0001, outer.thresh = 0.0001, min.frac = 0.01, nlam = 25, lambdas = NULL, alpha = 0.95, gamma = 0.8, step = 1, reset = 20, verbose = FALSE) {

    ## Setting up group lasso stuff ##

    ord <- order(index)
    index <- index[ord]
    X <- X[, ord]
    unOrd <- match(1:length(ord), ord)

    ## Coming up with other C++ info ##

    groups <- unique(index)
    num.groups <- length(groups)
    range.group.ind <- rep(0, (num.groups + 1))
    for (i in 1:num.groups) {
        range.group.ind[i] <- min(which(index == groups[i])) - 1
    }
    range.group.ind[num.groups + 1] <- ncol(X)

    group.length <- diff(range.group.ind)

    ## DONE SETTING UP C STUFF ##

    ## Setting up lambdas, if not given ##

    if (is.null(lambdas)) {
        # Calculating componentwise upper bounds for the function

        minusSumX_i <- -colSums(X[seq(1, nrow(X), 4), ])

        sumOfMax <- rep(0, ncol(X))
        sumOfMin <- rep(0, ncol(X))
        for (i in seq(1, nrow(X), 4)) {
            sumOfMax <- sumOfMax + apply(X[i:(i + 3), ], 2, max)
            sumOfMin <- sumOfMin + apply(X[i:(i + 3), ], 2, min)
        }
        bounds1 <- (sumOfMax + minusSumX_i) * 4 / nrow(X)
        bounds2 <- (sumOfMin + minusSumX_i) * 4 / nrow(X)

        # groupwise estimation not necessary, because for group k: sqrt(groupsize(k)) * norm_2(group-k-part of beta) >= norm_1(group-k-part of beta)
        # so max.lam is already a bound for groupwise lasso

        max.lam <- max(max(abs(bounds1)), max(abs(bounds2)))

        min.lam <- min.frac * max.lam
        lambdas <- exp(seq(log(max.lam), log(min.lam), (log(min.lam) - log(max.lam)) / (nlam - 1)))
        #lambdas <- seq(max.lam ^ 0.25, min.lam ^ 0.25, (min.lam ^ 0.25 - max.lam ^ 0.25) / (nlam - 1)) ^ 4
        #lambdas <- log(seq(exp(max.lam), exp(min.lam), (exp(min.lam) - exp(max.lam)) / (nlam - 1)))
        #lambdas <- seq(max.lam, min.lam, (min.lam - max.lam) / (nlam - 1))
    }

    ## time to start the fitting ##

    nlam = length(lambdas)
    beta.old <- rep(0, ncol(X))
    beta.is.zero <- rep(1, num.groups)
    beta <- matrix(0, nrow = ncol(X), ncol = nlam)

    eta <- rep(0, nrow(X))
    for (i in 1:nlam) {
        if (verbose == TRUE) {
            write(paste("***Lambda", i, "=", lambdas[i], "***"),"")
        }

        junk <- .C("triofit", X = as.double(as.vector(X)), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(lambdas[i] * alpha), lambda2 = as.double((1 - alpha) * lambdas[i]), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), step = as.double(step), reset = as.integer(reset))

        beta.new <- junk$beta
        beta[, i] <- beta.new
        beta.is.zero <- junk$betaIsZero
        eta <- junk$eta
        beta.old <- beta.new
        gc(FALSE)
    }

    return(list(beta = beta[unOrd, ], X = X, lambdas = lambdas))
}

