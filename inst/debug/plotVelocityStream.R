library(velociraptor)

###

# man page

# library(scuttle)
# sce1 <- mockSCE()
# sce2 <- mockSCE()
# 
# spliced <- counts(sce1)
# unspliced <- counts(sce2)
# 
# out <- scvelo(list(X=spliced, spliced=spliced, unspliced=unspliced), mode = "dynamical")

# vignette

library(scRNAseq)
sce <- HermannSpermatogenesisData()
sce

sce <- sce[, 1:500]

library(scuttle)
sce <- logNormCounts(sce, assay.type=1)

library(scran)
dec <- modelGeneVar(sce)
top.hvgs <- getTopHVGs(dec, n=2000)

library(velociraptor)
velo.out <- scvelo(
  sce, subset.row=top.hvgs, assay.X="spliced",
  scvelo.params=list(neighbors=list(n_neighbors=30L))
)
out <- velo.out

###

em <- embedVelocity(reducedDim(out, 1), out)[,1:2]

###

sce <- out
embedded <- em


use.dimred = 1
color_by = "#444444"
color.alpha = 0.2
grid.resolution = 60
scale = TRUE
stream.L = 10
stream.min.L = 0
stream.res = 4
stream.width = 8
color.streamlines = FALSE
color.streamlines.map = c("#440154", "#482576", "#414487",
    "#35608D", "#2A788E", "#21908C",
    "#22A884", "#43BF71", "#7AD151",
    "#BBDF27", "#FDE725")
arrow.angle = 8
arrow.length = 0.8

###
if (!identical(ncol(sce), nrow(embedded))) {
  stop("'sce' and 'embedded' do not have consistent dimensions.")
}
if (is.numeric(use.dimred)) {
  stopifnot(exprs = {
    identical(length(use.dimred), 1L)
    use.dimred <= length(reducedDims(sce))
  })
  use.dimred <- reducedDimNames(sce)[use.dimred]
} else if (is.character(use.dimred)) {
  stopifnot(exprs = {
    length(use.dimred) == 1L
    use.dimred %in% reducedDimNames(sce)
  })
} else {
  stop("'use.dimred' is not a valid value for use in reducedDim(sce, use.dimred)")
}
if (!requireNamespace("ggplot2")) {
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
if (is.character(color_by) && length(color_by) == 1L && color_by %in% colnames(colData(sce))) {
  plotdat1 <- cbind(plotdat1, col = colData(sce)[, color_by])
  colByFeat <- TRUE
} else {
  colByFeat <- FALSE
}
p <- ggplot2::ggplot(plotdat1, ggplot2::aes(x = !!ggplot2::sym("x"), y = !!ggplot2::sym("y"))) +
  ggplot2::labs(x = paste(use.dimred, "1"), y = paste(use.dimred, "2"))
if (!colByFeat) {
  colMatrix <- grDevices::col2rgb(col = color_by, alpha = TRUE)
  if (any(colMatrix[4, ] != 255)) {
    warning("ignoring 'color.alpha' as 'color_by' already specifies alpha channels")
    color.alpha <- colMatrix[4, ] / 255
  }
  p <- p + ggplot2::geom_point(color = color_by, alpha = color.alpha)
} else {
  p <- p + ggplot2::geom_point(ggplot2::aes(color = !!ggplot2::sym("col")), alpha = color.alpha) +
    ggplot2::labs(color = color_by)
}

plotdat2.backup <- plotdat2
# plotdat2 <- subset(plotdat2.backup, dx != 0 & dy != 0)
plotdat2 <- plotdat2.backup

if (color.streamlines) {
  # remark: when coloring streamlines, we currently cannot have any arrows
  # remark: ..dx.., ..dy.. and ..step.. are calculated by metR::geom_streamline
  p <- p +
    metR::geom_streamline(
      mapping = ggplot2::aes(
        x = !!ggplot2::sym("x"),
        y = !!ggplot2::sym("y"),
        dx = !!ggplot2::sym("dx"),
        dy = !!ggplot2::sym("dy")
      ),
      arrow = NULL, lineend = "round",
      data = plotdat2, linewidth = 0.6, jitter = 2,
      L = stream.L, min.L = stream.min.L,
      res = stream.res, inherit.aes = FALSE
    ) +
    ggplot2::scale_color_gradientn(colors = color.streamlines.map,
      guide = "none") +
    ggplot2::scale_alpha_continuous(guide = "none") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())
} else {
  p <- p +
    metR::geom_streamline(
      mapping = ggplot2::aes(
        x = !!ggplot2::sym("x"),
        y = !!ggplot2::sym("y"),
        dx = !!ggplot2::sym("dx"),
        dy = !!ggplot2::sym("dy")
      ),
      data = plotdat2, size = 0.3, jitter = 2,
      L = stream.L, min.L = stream.min.L,
      res = stream.res, arrow.angle = arrow.angle,
      arrow.length = arrow.length, inherit.aes = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())
}

p

