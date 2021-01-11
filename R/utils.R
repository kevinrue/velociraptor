#' @importFrom DelayedArray is_sparse
.make_np_friendly <- function(x) {
    if (is_sparse(x)) {
        as(x, "dgCMatrix")
    } else {
        as.matrix(x)
    }
}

#' @importFrom DelayedArray t
.extractor_python_dict <- function(thing, names, single=FALSE, transpose=FALSE) {
    if (single) {
        values <- lapply(names, function(x) {
            if (transpose) t(thing[x])
            else thing[x]
        })
    } else {
        values <- lapply(names, function(x) {
            if (transpose) t(thing[[x]])
            else thing[[x]]
        })
    }
    names(values) <- names
    values
}

# ggplot2-like colour scale in HCL space
.gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

# check if x is a valid R colour
#' @author Michael Stadler
.isValidColour <- function (x) {
    vapply(X = x, FUN = function(y) tryCatch(is.matrix(grDevices::col2rgb(y)), 
                                             error = function(e) FALSE),
           FUN.VALUE = logical(1), USE.NAMES = FALSE)
}

# map numeric values to colours
#' @author Michael Stadler
.valueToColour <- function (x, rng = range(x, na.rm = TRUE),
                           col = c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
                                   "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
                                   "#F46D43", "#D53E4F", "#9E0142"),
                           NA.col = "lightgray", alpha = NULL) {
    stopifnot(exprs = {
        is.numeric(x)
        is.numeric(rng) && length(rng) == 2L && rng[1] < rng[2]
        all(.isValidColour(col))
        .isValidColour(NA.col) && length(NA.col) == 1L
        is.null(alpha) || (is.numeric(alpha) && length(alpha) == 1L &&
                               alpha >= 0 && alpha <= 255)
    })
    colMatrix <- grDevices::col2rgb(col = col, alpha = TRUE)
    if (!is.null(alpha)) {
        if (all(colMatrix[4, ] == 255)) 
            colMatrix[4, ] <- as.integer(alpha)
        else warning("ignoring 'alpha', as 'col' already specifies alpha channels")
    }
    col <- grDevices::rgb(colMatrix[1, ], colMatrix[2, ], colMatrix[3,],
                          colMatrix[4, ], maxColorValue = 255)
    i <- !is.na(x)
    x[i & x < rng[1]] <- rng[1]
    x[i & x > rng[2]] <- rng[2]
    colfunc <- grDevices::colorRamp(col, alpha = TRUE)
    xcolMatrix <- colfunc((x[i] - rng[1])/(rng[2] - rng[1]))
    xcol <- rep(NA.col, length(x))
    xcol[i] <- grDevices::rgb(xcolMatrix[, 1], xcolMatrix[, 2], 
                              xcolMatrix[, 3], xcolMatrix[, 4],
                              maxColorValue = 255)
    return(xcol)
}

# add a color bar to the bottom right corner of a plot
#' @author Michael Stadler
.colbar <- function(rng, cols, width = 0.3, n = 32) {
    pusr <- graphics::par("usr")
    x0 <- pusr[1] + (1 - width) * diff(pusr[1:2])
    x1 <- pusr[2] - graphics::par("cxy")[1]
    xs <- seq(x0, x1, length.out = n + 1)
    y0 <- pusr[3]
    y1 <- y0 + 0.5 * graphics::par("cxy")[2]
    lcols <- grDevices::colorRampPalette(cols)(n)
    graphics::rect(xs[-(n+1)], y0, xs[-1], y1, col = lcols, border = NA)
    xl <- c(xs[1], mean(xs), xs[n+1])
    graphics::text(x = xl, y = y0, adj = c(0.5, 1.1), xpd = NA,
                   signif(c(rng[1], mean(rng), rng[2]), digits = 2))
}

# helper functions to calculate (un)spliced counts from starting point and model parameters
# see https://github.com/theislab/scvelo/blob/95d90de3d0935ce58a01218c9f179c9494ff593e/scvelo/tools/dynamical_model_utils.py#L109
.unspliced <- function(tau, u0, alpha, beta) {
    expu <- exp(-beta * tau)
    return(u0 * expu + alpha / beta * (1 - expu))
}
# see https://github.com/theislab/scvelo/blob/95d90de3d0935ce58a01218c9f179c9494ff593e/scvelo/tools/dynamical_model_utils.py#L114
.spliced <- function(tau, s0, u0, alpha, beta, gamma) {
    cn <- (alpha - u0 * beta) * 1 / (gamma - beta)
    expu <- exp(-beta * tau)
    exps <- exp(-gamma * tau)
    return(s0 * exps + alpha / gamma * (1 - exps) + cn * (exps - expu))
}
