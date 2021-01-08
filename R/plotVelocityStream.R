#' Velocity stream plot in low-dimensional space
#'
#' Plot velocities embedded into low-dimensional space as a stream plot. Stream
#' lines are lines that follow the gradient in the velocity field and illustrate
#' paths that cells could follow based on observed RNA velocities.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object containing
#'   low-dimensional coordinates, e.g., after t-SNE, in its
#'   \code{\link{reducedDims}}.
#' @param embedded A low-dimensional projection of the velocity vectors into the
#'   embedding of \code{sce}. This should be of the same dimensions as \code{sce}
#'   and is typically produced by \code{\link{embedVelocity}}.
#' @param use.dimred String or integer scalar specifying the reduced dimensions
#'   to retrieve from \code{sce}.
#' @param colour_by A character scalar specifying a column in \code{colData(sce)}
#'   to colour cells in the phase graph. Alternatively, \code{colour_by} can be
#'   set to a valid R colour to be used to colour cells.
#' @param colour.alpha An integer scalar giving the transparency of coloured
#'   cells. Possible values are between 0 (fully transparent) and 1.0 (opaque).
#' @param grid.resolution Integer scalar specifying the resolution of the grid,
#'   in terms of the number of grid intervals along each axis.
#' @param scale Logical scalar indicating whether the averaged vectors should be
#'   scaled by the grid resolution.
#' @param stream.L Integer scalar giving the typical length of a streamline
#'   low-dimensional space units.
#' @param stream.min.L A numeric scalar with the minimum length of segments to be shown.
#' @param stream.res Numeric scalar specifying the resolution of estimated
#'   streamlines (higher numbers increase smoothness of lines but also the time
#'   for computation).
#' @param stream.width A numeric scalar controlling the width of streamlines.
#' @param colour.streamlines Logical scalar. If \code{TRUE} streamlines will
#'   be colored by local velocity. Arrows cannot be shown in that case.
#' @param colour.streamlines.map A character vector specifying the
#'   colour range used for mapping local velocities to streamline colors. The
#'   default is \code{viridisLite::viridis(11)}.
#' @param arrow.angle,arrow.length Numeric scalars giving the \code{angle} and
#'   \code{length} of arrowheads.
#'   
#' @details \code{grid.resolution} and \code{scale} are passed to
#'   \code{\link{gridVectors}}, which is used to summarized the velocity vectors
#'   into an initial grid. A full regular grid is computed from that and used
#'   in \code{\link[metR]{geom_streamline}} to calculate streamlines. The
#'   following arguments are passed to the arguments given in parenthesis of
#'   \code{\link[metR]{geom_streamline}}:
#'   \code{stream.L} (\code{L}), \code{stream.res} (\code{res}),
#'   \code{stream.min.L} (\code{min.L}), \code{arrow.angle} (\code{arrow.angle})
#'   and \code{arrow.length} (\code{arrow.length}).
#'   Streamlines are computed by simple integration with a forward Euler method,
#'   and \code{stream.L} and \code{stream.res} are used to compute the number of
#'   steps and the time interval between steps for the integration.
#'   \code{stream.width} is multiplied with \code{..step..} estimated by
#'   \code{\link[metR]{geom_streamline}} to control the width of streamlines.
#'   
#' @return A \code{ggplot2} object with the streamline plot.
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
#' out <- scvelo(datlist, mode = "dynamical")
#' 
#' em <- embedVelocity(reducedDim(out, 1), out)[,1:2]
#' 
#' plotVelocityStream(out, em)
#' plotVelocityStream(out, em, colour.streamlines = TRUE)
#' 
#' @seealso \code{\link{gridVectors}} used to summarize velocity vectors into
#'   a grid (velocity field), the \pkg{ggplot2} package used for plotting,
#'   \code{\link[metR]{geom_streamline}} in package \pkg{metR} used to
#'   calculate and add streamlines from the RNA velocity field to the plot,
#'   \code{\link[viridisLite]{viridis}} for creation of color palettes.
#' 
#' @importFrom S4Vectors DataFrame
#' @importFrom ggplot2 ggplot aes geom_point labs scale_colour_gradientn scale_alpha_continuous coord_munch theme_minimal theme element_blank
#' @importFrom metR geom_streamline
plotVelocityStream <- function(sce, embedded, use.dimred = 1,
                               colour_by = "#444444", colour.alpha = 0.2,
                               grid.resolution = 60, scale = TRUE,
                               stream.L = 10, stream.min.L = 0, stream.res = 4,
                               stream.width = 8,
                               colour.streamlines = FALSE,
                               colour.streamlines.map = c("#440154", "#482576", "#414487",
                                                          "#35608D", "#2A788E", "#21908C",
                                                          "#22A884", "#43BF71", "#7AD151",
                                                          "#BBDF27", "#FDE725"),
                               arrow.angle = 8, arrow.length = 0.8) {
    stopifnot(exprs = {
        is(sce, "SingleCellExperiment")
        is.matrix(embedded)
        ncol(embedded) == 2L
        ncol(sce) == nrow(embedded)
        (.isValidColour(colour_by) && (length(colour_by) == 1L || length(colour_by) == ncol(sce))) ||
            (is.character(colour_by) && length(colour_by) == 1L && colour_by %in% colnames(colData(sce)))
        is.numeric(colour.alpha)
        length(colour.alpha) == 1L
        colour.alpha >= 0 && colour.alpha <= 1.0
        is.numeric(stream.L)
        length(stream.L) == 1L
        is.numeric(stream.min.L)
        length(stream.min.L) == 1L
        is.numeric(stream.res)
        length(stream.res) == 1L
        is.numeric(stream.width)
        length(stream.width) == 1L
        is.logical(colour.streamlines)
        length(colour.streamlines) == 1L
        .isValidColour(colour.streamlines.map)
    })
    if (is.numeric(use.dimred)) {
        stopifnot(exprs = {
            length(use.dimred) == 1L
            use.dimred <= length(reducedDims(sce))
        })
        use.dimred <- reducedDimNames(sce)[use.dimred]
    }
    else if (is.character(use.dimred)) {
        stopifnot(exprs = {
            length(use.dimred) == 1L
            use.dimred %in% reducedDimNames(sce)
        })
    }
    else {
        stop("'use.dimred' is not a valid value for use in reducedDim(sce, use.dimred)")
    }
    
    # get coordinates in reduced dimensional space
    xy <- reducedDim(sce, use.dimred)[, 1:2]
    
    # summarize velocities in a grid
    gr <- gridVectors(x = xy, embedded = embedded,
                      resolution = grid.resolution, scale = scale,
                      as.data.frame = FALSE,
                      return.intermediates = TRUE)
    
    # now make it a regular grid needed for metR::geom_streamline
    xbreaks <- seq(gr$limits[1,1], gr$limits[2,1], by = gr$delta[1])
    ybreaks <- seq(gr$limits[1,2], gr$limits[2,2], by = gr$delta[2])
    plotdat2 <- expand.grid(x = xbreaks + gr$delta[1] / 2,
                            y = ybreaks + gr$delta[2] / 2,
                            dx = 0, dy = 0)
    allcategories <- DataFrame(expand.grid(V1 = seq(0, grid.resolution),
                                           V2 = seq(0, grid.resolution)))
    ivec <- match(gr$categories[sort(unique(gr$grp)), ], allcategories)
    plotdat2[ivec, c("dx", "dy")] <- gr$vec
    
    
    # plot it using ggplot2 and metR::geom_streamline
    plotdat1 <- data.frame(xy)
    colnames(plotdat1) <- c("x", "y")
    if (!.isValidColour(colour_by)) {
        plotdat1 <- cbind(plotdat1, col = colData(sce)[, colour_by])
    }
    p <- ggplot(plotdat1, aes(x = x, y = y)) +
        labs(x = paste(use.dimred, "1"), y = paste(use.dimred, "2"))
    if (.isValidColour(colour_by)) {
        colMatrix <- grDevices::col2rgb(col = colour_by, alpha = TRUE)
        if (any(colMatrix[4, ] != 255)) {
            warning("ignoring 'colour.alpha' as 'colour_by' already specifies alpha channels")
            colour.alpha <- colMatrix[4, ] / 255
        }
        p <- p + geom_point(colour = colour_by, alpha = colour.alpha)
    } else {
        p <- p + geom_point(aes(colour = col), alpha = colour.alpha) + labs(colour = colour_by)
    }
    if (colour.streamlines) {
        # remark: when coloring streamlines, we currently cannot have any arrows
        p <- p + geom_streamline(mapping = aes(x = x, y = y, dx = dx, dy = dy,
                                               size = stream.width * ..step.., alpha = ..step..,
                                               colour = sqrt(..dx..^2 + ..dy..^2)),
                                 arrow = NULL, lineend = "round",
                                 data = plotdat2, size = 0.6, jitter = 2, L = stream.L,
                                 min.L = stream.min.L, res = stream.res,
                                 inherit.aes = FALSE) +
            scale_colour_gradientn(colours = colour.streamlines.map, guide = "none") +
            scale_alpha_continuous(guide = "none") +
            theme_minimal() + theme(axis.text = element_blank(),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank())
    } else {
        p <- p + geom_streamline(mapping = aes(x = x, y = y, dx = dx, dy = dy,
                                               size = stream.width * ..step..),
                                 data = plotdat2, size = 0.3, jitter = 2, L = stream.L,
                                 min.L = stream.min.L, res = stream.res,
                                 arrow.angle = arrow.angle, arrow.length = arrow.length,
                                 inherit.aes = FALSE) +
            theme_minimal() + theme(axis.text = element_blank(),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank())
    }
    
    return(p)
}

