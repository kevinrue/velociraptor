# This tests the plotting functionality.
# library(testthat); library(velociraptor); source("test-plotting.R")

library(scuttle)
#' set.seed(42)
sce1 <- mockSCE()
sce2 <- mockSCE()

spliced <- counts(sce1)
unspliced <- counts(sce2)
datlist <- list(X=counts(sce1), spliced=counts(sce1), unspliced=counts(sce2))

set.seed(100001)
out1 <- scvelo(datlist, mode = "steady_state")
out2 <- scvelo(datlist, mode = "dynamical")

genes1 <- rownames(out2)[!is.finite(rowData(out2)$fit_alpha)][1:2]
genes2 <- rownames(out2)[is.finite(rowData(out2)$fit_alpha)][1:2]

tsne.results <- matrix(rnorm(2*ncol(out1)), ncol=2)

em1 <- embedVelocity(reducedDim(out1, 1), out1)[,1:2]
em2 <- embedVelocity(reducedDim(out2, 1), out2)[,1:2]

out3 <- out1
colData(out3) <- cbind(colData(out3), type = sample(letters, ncol(out3), replace = TRUE))
rowData(out3) <- rowData(out3)[, -grep("gamma", colnames(rowData(out3)))]

out4 <- out2
rowData(out4) <- rowData(out4)[, -grep("fit_[us]0", colnames(rowData(out4)))]

test_that("plotVelocity runs", {
    
    expect_error(plotVelocity("error", genes2))
    expect_error(plotVelocity(out1, "error"))
    expect_error(plotVelocity(out1, genes2, use.dimred = "error"))
    expect_error(plotVelocity(out1, genes2, use.dimred = FALSE))
    expect_error(plotVelocity(out1, genes2, 1, "error"))
    expect_error(plotVelocity(out1, genes2, 1, which.plots = "error"))
    expect_error(plotVelocity(out1, genes2, 1, genes.per.row = "error"))
    expect_error(plotVelocity(out1, genes2, 1, color_by = "error"))
    expect_error(plotVelocity(out1, genes2, 1, color.alpha = "error"))
    expect_error(plotVelocity(out1, genes2, 1, colors.velocity = "error"))
    expect_error(plotVelocity(out1, genes2, 1, colors.expression = "error"))
    expect_error(plotVelocity(out1, genes2, 1, max.abs.velo = "error"))
    expect_error(plotVelocity(out1, genes2, 1, max.abs.velo = -1))
    

    tf <- tempfile(fileext = ".png")
    png(tf)
    expect_warning(res0 <- plotVelocity(out1, genes1, which.plots = "phase", color_by = "#44444422"))
    expect_identical(res0, c(2, 1))
    res1 <- plotVelocity(out1, genes1)
    expect_identical(res1, c(2, 3))
    expect_warning(res2 <- plotVelocity(out2, genes1))
    expect_identical(res1, res2)
    res3 <- plotVelocity(out4, genes2)
    expect_identical(res1, res3)
    res4 <- plotVelocity(out3, genes2, which.plots = "phase", color_by = "type")
    expect_identical(res0, res4)
    res5 <- plotVelocity(out2, genes2)
    expect_identical(res1, res5)
    dev.off()
    expect_true(file.exists(tf))
    unlink(tf)
})

test_that("plotVelocityStream runs", {

    skip_if_not_installed("ggplot2")
    skip_if_not_installed("metR")
    
    expect_error(plotVelocityStream("error", em2))
    expect_error(plotVelocityStream(out2, "error"))
    expect_error(plotVelocityStream(out2, em2, use.dimred = "error"))
    expect_error(plotVelocityStream(out2, em2, use.dimred = FALSE))
    expect_error(plotVelocityStream(out2, em2, color_by = "error"))
    expect_error(plotVelocityStream(out2, em2, color.alpha = 2))
    expect_error(plotVelocityStream(out2, em2, grid.resolution = "error"))
    expect_error(plotVelocityStream(out2, em2, scale = "error"))
    expect_error(plotVelocityStream(out2, em2, stream.L = "error"))
    expect_error(plotVelocityStream(out2, em2, stream.min.L = "error"))
    expect_error(plotVelocityStream(out2, em2, stream.res = "error"))
    expect_error(plotVelocityStream(out2, em2, stream.width = "error"))
    expect_error(plotVelocityStream(out2, em2, color.streamlines = "error"))
    expect_error(plotVelocityStream(out2, em2, color.streamlines.map = "error"))
    expect_error(plotVelocityStream(out2, em2, arrow.angle = "error"))
    expect_error(plotVelocityStream(out2, em2, arrow.length = "error"))
    
    tf <- tempfile(fileext = ".png")
    png(tf)
    expect_warning(print(plotVelocityStream(out2, em2, color_by = "#44444422")))
    print(plotVelocityStream(out2, em2))
    print(plotVelocityStream(out3, em2, color_by = "type"))
    print(plotVelocityStream(out2, em2, color.streamlines = TRUE))
    dev.off()
    expect_true(file.exists(tf))
    unlink(tf)
})
