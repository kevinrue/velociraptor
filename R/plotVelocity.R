#' Phase and velocity graphs for a set of genes
#'
#' For a each gene in a set of genes, show the phase graph (spliced versus
#' unspliced counts and fitted model) and reduced dimension graphs with
#' cell colored by velocity and (spliced) expression.
#'
#' @param x A \linkS4class{SingleCellExperiment} object with RNA velocity results
#'   as returned by \code{\link{scvelo}}, and low-dimensional coordinates, e.g.,
#'   after t-SNE, in its \code{\link{reducedDims}}.
#' @param genes A character vector with one or several genes for which to plot
#'   phase and velocity graphs. \code{genes} have to be in \code{rownames(x)}.
#' @param use.dimred String or integer scalar specifying the reduced dimensions
#'   to retrieve from \code{x}.
#' @param assay.splicedM An integer scalar or string specifying the assay of
#'   \code{x} containing the moments of spliced abundances.
#' @param assay.unsplicedM An integer scalar or string specifying the assay of
#'   \code{x} containing the moments unspliced abundances.
#' @param which.plots A character vector specifying which plots to create for
#'   each gene. Possible values are \code{"phase", "velocity", "expression"} and
#'   correspond to the phase graph or reduced dimension graphs with cells
#'   colored by velocity or (spliced) expression.
#' @param genes.per.row An integer scalar with the numbers of genes to visualize
#'   per row of plots. For example, if \code{which.plots =  c("phase","expression")}
#'   and \code{genes.per.row = 2}, the resulting figure will have four plot
#'   panels per row.
#' @param color_by A character scalar specifying a column in \code{colData(x)}
#'   to color cells in the phase graph. Alternatively, \code{color_by} can be
#'   set to vector of valid R colors, either of length one (recycled for all
#'   cells) or of length \code{ncol(x)}, which will then be used to color cells
#'   in the phase graph.
#' @param color.alpha An integer scalar giving the transparency of colored
#'   cells. Possible values are between 0 (fully transparent) and 1.0 (opaque).
#' @param colors.velocity,colors.expression Character vectors specifying the
#'   color ranges used for mapping velocities and expression values. The
#'   defaults are \code{RColorBrewer::brewer.pal(11, "RdYlGn")} for the
#'   velocities and \code{viridisLite::viridis(11)} for the expression values.
#' @param max.abs.velo A numeric scalar greater than zero giving the maximum
#'   absolute velocity to limit the color scale for the \code{"velocity"} graph.
#'
#' @details Please note that \code{plotVelocity} will modify parameters of
#'   the current graphics device using \code{\link{layout}} and \code{\link{par}},
#'   in order to create the layout for the generated graph panels.
#' 
#' @return A patchwork object with the plots selected by \code{which.plot} for
#'   the genes in \code{genes}, arranged in a grid according to \code{genes.per.row}.
#'
#' @author Michael Stadler
#' 
#' @examples
#' library(scuttle)
#' set.seed(42)
#' sce1 <- mockSCE(ncells = 100, ngenes = 500)
#' sce2 <- mockSCE(ncells = 100, ngenes = 500)
#'
#' datlist <- list(X=counts(sce1), spliced=counts(sce1), unspliced=counts(sce2))
#'
#' out1 <- scvelo(datlist, mode = "steady_state")
#' out2 <- scvelo(datlist, mode = "dynamical")
#' 
#' plotVelocity(out1, c("Gene_0031","Gene_0268"))
#' plotVelocity(out2, c("Gene_0031","Gene_0268"))
#'
#' @seealso 
#' \code{\link{scvelo}}, to generate \code{x},
#' \code{\link[RColorBrewer]{brewer.pal}} and \code{\link[viridisLite]{viridis}}
#' for creation of color palettes, packages \pkg{ggplot2} and \pkg{patchwork}
#' used to generate and arrange the plots.
#' 
#' @export
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom SingleCellExperiment reducedDim reducedDims reducedDimNames
plotVelocity <- function(x, genes, use.dimred = 1,
                         assay.splicedM = "Ms", assay.unsplicedM = "Mu",
                         which.plots = c("phase", "velocity", "expression"),
                         genes.per.row = 1,
                         color_by = "#222222",
                         color.alpha = 0.4,
                         colors.velocity = c("#A50026", "#D73027", "#F46D43",
                                              "#FDAE61", "#FEE08B", "#FFFFBF",
                                              "#D9EF8B", "#A6D96A", "#66BD63",
                                              "#1A9850", "#006837"),
                         colors.expression = c("#440154", "#482576", "#414487",
                                                "#35608D", "#2A788E", "#21908C",
                                                "#22A884", "#43BF71", "#7AD151",
                                                "#BBDF27", "#FDE725"),
                         max.abs.velo = 0.001) {
    # check arguments
    genes <- unique(genes)
    if (!all(genes %in% rownames(x))) {
        stop("not all 'genes' are not found in 'x'")
    }
    stopifnot(exprs = {
        all(which.plots %in% c("phase", "velocity", "expression"))
        is.numeric(genes.per.row)
        length(genes.per.row) == 1L
        is.numeric(max.abs.velo)
        length(max.abs.velo) == 1L
        max.abs.velo >= 0.0
    })
    if (is.numeric(use.dimred)) {
        stopifnot(exprs = {
            length(use.dimred) == 1L
            use.dimred <= length(reducedDims(x))
        })
        use.dimred <- reducedDimNames(x)[use.dimred]
    } else if (is.character(use.dimred)) {
        stopifnot(exprs = {
            length(use.dimred) == 1L
            use.dimred %in% reducedDimNames(x)
        })
    } else {
        stop("'use.dimred' is not a valid value for use in reducedDim(x, use.dimred)")
    }

    # extract data from x
    S <- assay(x, assay.splicedM)
    U <- assay(x, assay.unsplicedM)
    V <- assay(x, "velocity")
    rd <- rowData(x)
    xy <- reducedDim(x, use.dimred)
    df2 <- data.frame(x = xy[, 1],
                      y = xy[, 2])

    
    # iterate over genes and create graphs
    pL <- list()
    for (gene in genes) {
        df1 <- data.frame(s = as.vector(S[gene, ]),
                          u = as.vector(U[gene, ]))

        # spliced/unspliced phase portrait with model estimates
        if ("phase" %in% which.plots) {
            if (is.character(color_by) && length(color_by) == 1L && color_by %in% colnames(colData(x))) {
                df1$col <- factor(colData(x)[, color_by])
                colByFeat <- TRUE
            } else if (length(color_by) == 1L || length(color_by) == ncol(x)) {
                df1$col <- color_by
                colByFeat <- FALSE
            } else {
                stop("invalid 'colour_by' (not in colData(x), nor of the right length)")
            }
            p1 <- ggplot2::ggplot(df1, ggplot2::aes(x = !!ggplot2::sym("s"),
                                                    y = !!ggplot2::sym("u"))) +
                ggplot2::xlim(0, max(df1$s, na.rm = TRUE)) +
                ggplot2::ylim(0, max(df1$u, na.rm = TRUE)) +
                ggplot2::labs(title = gene, x = "spliced", y = "unspliced") +
                ggplot2::theme_minimal() +
                ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(),
                               axis.ticks = ggplot2::element_line(color = "grey20"),
                               panel.border = ggplot2::element_rect(fill = NA, color = "grey20"))
            if (colByFeat) {
                p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = !!ggplot2::sym("col")),
                                               alpha = color.alpha, show.legend = FALSE)
            } else {
                p1 <- p1 + ggplot2::geom_point(color = df1$col, alpha = color.alpha)
            }

            # gene specific model parameters
            #   scvelo(..., mode = "dynamical") creates "fit_" prefixes,
            #   all other modes only have "velocity_gamma"
            if ("fit_gamma" %in% colnames(rd)) {
                alpha <- rd[gene, "fit_alpha"]
                beta <- rd[gene, "fit_beta"] * rd[gene, "fit_scaling"]
                gamma <- rd[gene, "fit_gamma"]
                scaling <- rd[gene, "fit_scaling"]
                # t_ is the pseudotime separating the increase/decrease phases
                t_ <- rd[gene, "fit_t_"]
            } else {
                alpha <- beta <- scaling <- t_ <- NULL
                if ("velocity_gamma" %in% colnames(rd)) {
                    gamma <- rd[gene, "velocity_gamma"]
                } else {
                    gamma <- NULL
                }
            }
            
            # steady-state ratio
            if (!is.null(gamma) && is.finite(gamma)) {
                if (!is.null(beta) && is.finite(beta)) {
                    # in the dynamical model
                    p1 <- p1 + ggplot2::geom_abline(intercept = 0,
                                                    slope = gamma / beta * scaling,
                                                    linetype = "dashed")
                } else {
                    # in a static model
                    p1 <- p1 + ggplot2::geom_abline(intercept = 0,
                                                    slope = gamma,
                                                    linetype = "dashed")
                }
            }

            # dynamic model curve
            if (!is.null(beta) && is.finite(beta)) {
                # calculate spliced and unspliced abundances from initial
                #   conditions and fitted model parameters
                # see https://github.com/theislab/scvelo/blob/95d90de3d0935ce58a01218c9f179c9494ff593e/scvelo/plotting/simulation.py#L26
                if (all(c("fit_u0","fit_s0") %in% colnames(rd))) {
                    u0_offset <- rd[gene, "fit_u0"]
                    s0_offset <- rd[gene, "fit_s0"]
                } else {
                    u0_offset <- s0_offset <- 0
                }
                # see https://github.com/theislab/scvelo/blob/95d90de3d0935ce58a01218c9f179c9494ff593e/scvelo/tools/dynamical_model_utils.py#L602
                ft <- as.vector(assay(x, "fit_t")[gene, ])
                o <- as.numeric(ft < t_) # 1 or 0 for cells increasing or decreasing this gene's expression
                tau <- ft * o + (ft - t_) * (1 - o)
                u0 <- s0 <- alpha0 <- 0
                u0_ <- .unspliced(t_, u0, alpha, beta)
                s0_ <- .spliced(t_, s0, u0, alpha, beta, gamma)
                u0 <- u0 * o + u0_ * (1 - o)
                s0 <- s0 * o + s0_ * (1 - o)
                alpha <- alpha * o + alpha0 * (1 - o)
                # see https://github.com/theislab/scvelo/blob/95d90de3d0935ce58a01218c9f179c9494ff593e/scvelo/plotting/simulation.py#L54
                ut <- .unspliced(tau, u0, alpha, beta) * scaling + u0_offset
                st <- .spliced(tau, s0, u0, alpha, beta, gamma) + s0_offset
                p1 <- p1 + ggplot2::geom_point(data = data.frame(x = st, y = ut),
                                               mapping = ggplot2::aes(x = !!ggplot2::sym("x"),
                                                                      y = !!ggplot2::sym("y")),
                                               color = "purple", shape = 46)
            }
            pL <- c(pL, list(p1))
        }
        
        # velocity plot
        if ("velocity" %in% which.plots) {
            df2$z <- as.vector(V[gene, ])
            if (all(!is.finite(df2$z))) {
                warning("velocity estimates for ", gene, " are all non-finite")
            }
            mx <- max(c(max.abs.velo, abs(df2$z)), na.rm = TRUE)
            p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = !!ggplot2::sym("x"),
                                                    y = !!ggplot2::sym("y"),
                                                    color = !!ggplot2::sym("z"))) +
                ggplot2::geom_point(alpha = color.alpha) +
                ggplot2::scale_color_gradientn(colours = colors.velocity,
                                               limits = c(-mx, mx)) +
                ggplot2::labs(title = "velocity", x = paste0(use.dimred, " 1"),
                              y = paste0(use.dimred, " 2"), color = ggplot2::element_blank()) +
                ggplot2::theme_minimal() +
                ggplot2::theme(axis.text = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank())
            pL <- c(pL, list(p2))
        }
        
        # expression plot
        if ("expression" %in% which.plots) {
            df2$z <- as.vector(S[gene, ])
            p3 <- ggplot2::ggplot(df2, ggplot2::aes(x = !!ggplot2::sym("x"),
                                                    y = !!ggplot2::sym("y"),
                                                    color = !!ggplot2::sym("z"))) +
                ggplot2::geom_point(alpha = color.alpha) +
                ggplot2::scale_color_gradientn(colours = colors.expression) +
                ggplot2::labs(title = "expression", x = paste0(use.dimred, " 1"),
                              y = paste0(use.dimred, " 2"), color = ggplot2::element_blank()) +
                ggplot2::theme_minimal() +
                ggplot2::theme(axis.text = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank())
            pL <- c(pL, list(p3))
        }
    }
    
    # create plot panel layout
    nplts <- length(unique(which.plots))
    nrows <- ceiling(length(genes) / genes.per.row)
    p <- do.call(patchwork::wrap_plots, c(pL, list(nrow = nrows)))

    return(p)
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
