# This tests the utils functionality.
# library(testthat); library(velociraptor); source("test-utils.R")

test_that(".isValidColor works correctly", {
    expect_false(.isValidColor("error"))
    expect_true(all(.isValidColor(1:10)))
    expect_true(all(.isValidColor(c("red","green","blue"))))
    expect_true(all(.isValidColor(c("#2244FF", "#44556677"))))
})
