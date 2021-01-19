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
#' @param color_by A character scalar specifying a column in \code{colData(sce)}
#'   to color cells in the phase graph. Alternatively, \code{color_by} can be
#'   set to a valid R color to be used to color cells.
#' @param color.alpha An integer scalar giving the transparency of colored
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
#' @param color.streamlines Logical scalar. If \code{TRUE} streamlines will
#'   be colored by local velocity. Arrows cannot be shown in that case.
#' @param color.streamlines.map A character vector specifying the
#'   color range used for mapping local velocities to streamline colors. The
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
#' plotVelocityStream(out, em, color.streamlines = TRUE)
#' 
#' @seealso \code{\link{gridVectors}} used to summarize velocity vectors into
#'   a grid (velocity field), the \pkg{ggplot2} package used for plotting,
#'   \code{\link[metR]{geom_streamline}} in package \pkg{metR} used to
#'   calculate and add streamlines from the RNA velocity field to the plot,
#'   \code{\link[viridisLite]{viridis}} for creation of color palettes.
#' 
#' @export
#' @importFrom S4Vectors DataFrame
plotVelocityStream <- function(sce, embedded, use.dimred = 1,
                               color_by = "#444444", color.alpha = 0.2,
                               grid.resolution = 60, scale = TRUE,
                               stream.L = 10, stream.min.L = 0, stream.res = 4,
                               stream.width = 8,
                               color.streamlines = FALSE,
                               color.streamlines.map = c("#440154", "#482576", "#414487",
                                                          "#35608D", "#2A788E", "#21908C",
                                                          "#22A884", "#43BF71", "#7AD151",
                                                          "#BBDF27", "#FDE725"),
                               arrow.angle = 8, arrow.length = 0.8) {
    stopifnot(exprs = {
        is(sce, "SingleCellExperiment")
        is.matrix(embedded)
        ncol(embedded) == 2L
        ncol(sce) == nrow(embedded)
        (.isValidColor(color_by) && (length(color_by) == 1L || length(color_by) == ncol(sce))) ||
            (is.character(color_by) && length(color_by) == 1L && color_by %in% colnames(colData(sce)))
        is.numeric(color.alpha)
        length(color.alpha) == 1L
        color.alpha >= 0 && color.alpha <= 1.0
        is.numeric(stream.L)
        length(stream.L) == 1L
        is.numeric(stream.min.L)
        length(stream.min.L) == 1L
        is.numeric(stream.res)
        length(stream.res) == 1L
        is.numeric(stream.width)
        length(stream.width) == 1L
        is.logical(color.streamlines)
        length(color.streamlines) == 1L
        .isValidColor(color.streamlines.map)
        is.numeric(arrow.angle)
        length(arrow.angle) == 1L
        is.numeric(arrow.length)
        length(arrow.length) == 1L
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
    if (!require(ggplot2)) {
        stop("'plotVelocityStream' requires the package 'ggplot2'.")
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
    if (!.isValidColor(color_by)) {
        plotdat1 <- cbind(plotdat1, col = colData(sce)[, color_by])
    }
    p <- ggplot2::ggplot(plotdat1, ggplot2::aes(x = x, y = y)) +
        ggplot2::labs(x = paste(use.dimred, "1"), y = paste(use.dimred, "2"))
    if (.isValidColor(color_by)) {
        colMatrix <- grDevices::col2rgb(col = color_by, alpha = TRUE)
        if (any(colMatrix[4, ] != 255)) {
            warning("ignoring 'color.alpha' as 'color_by' already specifies alpha channels")
            color.alpha <- colMatrix[4, ] / 255
        }
        p <- p + ggplot2::geom_point(color = color_by, alpha = color.alpha)
    } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(color = col), alpha = color.alpha) +
            ggplot2::labs(color = color_by)
    }
    if (color.streamlines) {
        # remark: when coloring streamlines, we currently cannot have any arrows
        p <- p +
            metR::geom_streamline(mapping = ggplot2::aes(x = x, y = y, dx = dx, dy = dy,
                                                         size = stream.width * ..step..,
                                                         alpha = ..step..,
                                                         color = sqrt(..dx..^2 + ..dy..^2)),
                                  arrow = NULL, lineend = "round",
                                  data = plotdat2, size = 0.6, jitter = 2,
                                  L = stream.L, min.L = stream.min.L,
                                  res = stream.res, inherit.aes = FALSE) +
            ggplot2::scale_color_gradientn(colors = color.streamlines.map,
                                           guide = "none") +
            ggplot2::scale_alpha_continuous(guide = "none") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank())
    } else {
        p <- p +
            metR::geom_streamline(mapping = ggplot2::aes(x = x, y = y, dx = dx, dy = dy,
                                                         size = stream.width * ..step..),
                                  data = plotdat2, size = 0.3, jitter = 2,
                                  L = stream.L, min.L = stream.min.L,
                                  res = stream.res, arrow.angle = arrow.angle,
                                  arrow.length = arrow.length, inherit.aes = FALSE) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank())
    }
    
    return(p)
}

