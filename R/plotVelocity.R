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
#'   coloured by velocity or (spliced) expression.
#' @param genes.per.row An integer scalar with the numbers of genes to visualize
#'   per row of plots. For example, if \code{which.plots =  c("phase","expression")}
#'   and \code{genes.per.row = 2}, the resulting figure will have four plot
#'   panels per row.
#' @param colour_by A character scalar specifying a column in \code{colData(x)}
#'   to colour cells in the phase graph. Alternatively, \code{colour_by} can be
#'   set to vector of valid R colours, either of length one (recycled for all
#'   cells) or of length \code{ncol(x)}, which will then be used to colour cells
#'   in the phase graph.
#' @param colour.alpha An integer scalar giving the transparency of the
#'   cell colours in all graphs. Possible values between 0 (fully transparent)
#'   and 255 (opaque).
#' @param colours.velocity,colours.expression Character vectors specifying the
#'   colour ranges used for mapping velocities and expression values. The
#'   defaults are \code{RColorBrewer::brewer.pal(11, "RdYlGn")} for the
#'   velocities and \code{viridisLite::viridis(11)} for the expression values.
#' @param max.abs.velo A numeric scalar greater than zero giving the maximum
#'   absolute velocity to limit the color scale for the \code{"velocity"} graph.
#'
#' @details Please note that \code{plotVelocity} will modify parameters of
#'   the current graphics device using \code{par()}, in order to create the
#'   layout for the generated graph panels.
#' 
#' @return Invisible a two element vector with the numbers of rows and columns
#'   of the generated graph panel layout.
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
#' for creation of color palettes.
#' 
#' @export
#' @importFrom grDevices col2rgb
#' @importFrom graphics layout
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom SingleCellExperiment reducedDim reducedDims reducedDimNames
plotVelocity <- function(x, genes, use.dimred = 1,
                         assay.splicedM = "Ms", assay.unsplicedM = "Mu",
                         which.plots = c("phase", "velocity", "expression"),
                         genes.per.row = 1,
                         colour_by = "#222222",
                         colour.alpha = 100,
                         colours.velocity = c("#A50026", "#D73027", "#F46D43",
                                              "#FDAE61", "#FEE08B", "#FFFFBF",
                                              "#D9EF8B", "#A6D96A", "#66BD63",
                                              "#1A9850", "#006837"),
                         colours.expression = c("#440154", "#482576", "#414487",
                                                "#35608D", "#2A788E", "#21908C",
                                                "#22A884", "#43BF71", "#7AD151",
                                                "#BBDF27", "#FDE725"),
                         max.abs.velo = 0.001) {
    # check arguments
    genes <- unique(genes)
    stopifnot(exprs = {
        is.character(genes)
        all(genes %in% rownames(x))
        assay.splicedM %in% assayNames(x)
        assay.unsplicedM %in% assayNames(x)
        "velocity" %in% assayNames(x)
        is.character(which.plots)
        all(which.plots %in% c("phase", "velocity", "expression"))
        is.numeric(genes.per.row)
        length(genes.per.row) == 1L
        (.isValidColour(colour_by) && (length(colour_by) == 1L || length(colour_by) == ncol(x))) ||
            (is.character(colour_by) && length(colour_by) == 1L && colour_by %in% colnames(colData(x)))
        is.numeric(colour.alpha)
        length(colour.alpha) == 1L
        colour.alpha >= 0 && colour.alpha <= 255
        .isValidColour(colours.velocity)
        .isValidColour(colours.expression)
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
    }
    else if (is.character(use.dimred)) {
        stopifnot(exprs = {
            length(use.dimred) == 1L
            use.dimred %in% reducedDimNames(x)
        })
    }
    else {
        stop("'use.dimred' is not a valid value for use in reducedDim(x, use.dimred)")
    }

    # create plot panel layout
    which.plots <- unique(which.plots)
    nplts <- length(which.plots)
    nrows <- ceiling(length(genes) / genes.per.row)
    ncols <- genes.per.row * nplts
    
    layout(mat = matrix(seq.int(length(genes) * nplts),
                        nrow = nrows, ncol = ncols, byrow = TRUE))
    
    # extract data from x
    S <- assay(x, assay.splicedM)
    U <- assay(x, assay.unsplicedM)
    V <- assay(x, "velocity")
    rd <- rowData(x)
    xy <- reducedDim(x, use.dimred)
    
    # iterate over genes and create graphs
    for (gene in genes) {
        s <- as.vector(S[gene, ])
        u <- as.vector(U[gene, ])

        # spliced/unspliced phase portrait with model estimates
        if ("phase" %in% which.plots) {
            par(mar = c(3,3,2,0), mgp = c(1.75, 0.75, 0))
            if (!.isValidColour(colour_by)) {
                coli <- factor(x[[colour_by]])
                colour_by <- .gg_color_hue(nlevels(coli))[as.numeric(coli)]
            }
            colMatrix <- grDevices::col2rgb(col = colour_by, alpha = TRUE)
            if (any(colMatrix[4, ] != 255)) {
                warning("ignoring 'colour.alpha' in phase plot, ",
                        "as 'colour_by' already specifies alpha channels")
            } else {
                cols <- paste0(colour_by, as.hexmode(colour.alpha))
            }
            plot(s, u, xlab = "spliced", ylab = "unspliced", pch = "*",
                 col = cols, main = gene,
                 xlim = range(c(0, s), na.rm = TRUE),
                 ylim = range(c(0, u), na.rm = TRUE))

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
                    abline(a = 0, b = (gamma / beta), lty = 3)
                } else {
                    # in a static model
                    abline(a = 0, b = gamma, lty = 3)
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
                points(st, ut, pch = ".", col = "purple")
            }
        }
        
        # velocity plot
        if ("velocity" %in% which.plots) {
            z <- as.vector(V[gene, ])
            if (all(!is.finite(z))) {
                warning("velocity estimates for ", gene, " are all non-finite")
            }
            mx <- max(c(max.abs.velo, abs(z)), na.rm = TRUE)
            cols <- .valueToColour(z, rng = c(-mx, mx), col = colours.velocity,
                                   alpha = colour.alpha)
            par(mar = c(3,3,2,0), mgp = c(1.75, 0.75, 0))
            plot(xy[,1], xy[,2],
                 xlab = paste0(use.dimred, " 1"), ylab = paste0(use.dimred, " 2"),
                 main = "velocity", pch = "*", col = cols, axes = FALSE)
            .colbar(c(-mx, mx), colours.velocity)
        }
        
        # expression plot
        if ("expression" %in% which.plots) {
            z <- s
            cols <- .valueToColour(z, col = colours.expression,
                                   alpha = colour.alpha)
            par(mar = c(3,3,2,0), mgp = c(1.75, 0.75, 0))
            plot(xy[,1], xy[,2],
                 xlab = paste0(use.dimred, " 1"), ylab = paste0(use.dimred, " 2"),
                 main = "expression", pch = "*", col = cols, axes = FALSE)
            .colbar(range(z, na.rm = TRUE), colours.expression)
        }
    }
    
    return(invisible(c(nrows, ncols)))
}
