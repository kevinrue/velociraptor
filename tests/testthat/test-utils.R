# This tests the utils functionality.
# library(testthat); library(velociraptor); source("test-utils.R")

test_that(".gg_color_hue works correctly", {
    res <- lapply(1:10, .gg_color_hue)
    expect_identical(lengths(res), 1:10)
    expect_true(all(unlist(lapply(res, is.character))))
    expect_true(all(unlist(lapply(res, .isValidColour))))
})

test_that(".isValidColour works correctly", {
    expect_false(.isValidColour("error"))
    expect_true(all(.isValidColour(1:10)))
    expect_true(all(.isValidColour(c("red","green","blue"))))
    expect_true(all(.isValidColour(c("#2244FF", "#44556677"))))
})

test_that(".valueToColour() works properly", {
    expect_error(.valueToColour(x = "error"))
    expect_error(.valueToColour(x = 1:2, rng = "error"))
    expect_error(.valueToColour(x = 1:2, rng = 1:3))
    expect_error(.valueToColour(x = 1:2, col = "error"))
    expect_error(.valueToColour(x = 1:2, NA.col = "error"))
    expect_error(.valueToColour(x = 1:2, alpha = 256))
    expect_warning(.valueToColour(x = 1:2, col = "#10101020", alpha = 128))
    
    res1 <- .valueToColour(x = 1:10)
    res2 <- .valueToColour(x = 1:10, rng = 4:5)
    res3 <- .valueToColour(x = 1:2, col = c("white","black"), alpha = 128)
    expect_is(res1, "character")
    expect_length(res1, 10L)
    expect_length(res2, 10L)
    expect_length(unique(res2), 2L)
    expect_identical(res3, c("#FFFFFF80", "#00000080"))
})
