

test_that("norm_type_or_null works", {

  x <- norm_type_or_null(NULL)
  expect_null(x)

  x <- norm_type_or_null(1L, integer())
  expect_identical(x, 1L)

  x <- norm_type(1.2, integer())
  expect_identical(x, 1L)

  x <- norm_type(c(1, 2), numeric(), length = NULL)
  expect_identical(x, c(1, 2))

})

# test_that("norm_type throws appropriate errors", {
#
#   x <- 1
#   z <- substitute(norm_type(x, numeric(), length = 2))
#   expect_error(
#     eval(z),
#     "`x` has length 1, expected length is 2."
#   )
#
#   x <- c(1, NA, 2)
#   z <- substitute(norm_type(x, numeric(), length = 3))
#   expect_error(
#     eval(z),
#     "`x` contains NAs"
#   )
# })
