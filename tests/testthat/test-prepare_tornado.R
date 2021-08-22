
test_that("scales can be resolved", {

  # continuous scale should be accepted
  scale <- ggplot2::scale_fill_distiller()

  out <- resolve_scale(c(0, 1), scale, 1, 0, title = "dummy")
  expect_s3_class(out, "ScaleContinuous")
  expect_s3_class(out$guide, "colorbar")
  expect_equal(out$guide$title, "dummy")

  # legend guides should be accepted
  scale <- ggplot2::scale_fill_distiller(guide = ggplot2::guide_legend())

  out <- resolve_scale(c(0, 1), scale, 1, 0, title = "dummy")
  expect_s3_class(out$guide, "legend")

  # nonsense guides should be rejected
  scale <- ggplot2::scale_fill_distiller(guide = "nonsense")
  out <- substitute(resolve_scale(c(0, 1), scale, 1, 0, title = "dummy"))

  expect_error({out <- eval(out)}, "object 'guide_nonsense'")

  # discrete scale should be rejected
  scale <- ggplot2::scale_fill_brewer()

  expect_error(resolve_scale(c(0, 1), scale, upper = 1, lower = 1),
               "continuous colour scale.")

})

test_that("limits can be chosen", {
  data <- 1:100
  lim <- choose_limits(data = data)
  expect_equal(lim, c(0, 99), tolerance = 0.02)

  lim <- choose_limits(data = data, upper = "p99")
  expect_equal(lim, c(0, 99), tolerance = 0.02)

  lim <- choose_limits(lower = "q0.01", upper = 100, data = data)
  expect_equal(lim, c(2, 100), tolerance = 0.02)

  lim <- choose_limits(lower = "p1", upper = 100, data = data)
  expect_equal(lim, c(2, 100), tolerance = 0.02)

  lim <- choose_limits(data = data, oldlim = c(20, 40))
  expect_equal(lim, c(20, 40))

  lim <- substitute(choose_limits(lower = "z1", upper = "z99"))
  warn <- capture_warnings({lim <- eval(lim)})
  expect_match(warn[1], "Could not interpret 'upper'.")
  expect_match(warn[2], "Could not interpret 'lower'.")
})


test_that("tornado_df can be rbind'ed", {
  tor <- dummy_tornado()

  df1 <- prep_tornado(
    tor, scale = ggplot2::scale_fill_distiller(palette = "Blues")
  )
  df2 <- prep_tornado(
    tor, scale = ggplot2::scale_fill_distiller(palette = "Reds")
  )
  df2$sample_name <- paste0(df2$sample_name, "2")

  df3 <- rbind(df1, df2)
  expect_s3_class(df3, "tornado_df")

  tlist <- df3$tornado
  expect_length(scale_list(tlist), 2)

  expect_equal(
    format(tlist),
    rep("[2 x 50] colour matrix", length(tlist))
  )

  expect_output(
    print(tlist), "<list_of_tornados\\[8\\]>"
  )

  p <- autoplot(df3)
  expect_s3_class(p, "ggplot")


  grDevices::png({tmp <- tempfile(fileext = ".png")})
  gt <- ggplot2::ggplotGrob(p)
  grDevices::dev.off()

  expect_s3_class(gt, "gtable")
})


test_that("tornado_list works with imperfect input", {

  tlist <- tornado_list(matrix(LETTERS[1:4], 2, 2))

  expect_null(scale_list(tlist)[[1]])
  expect_null(levels(tlist))

  expect_equal(vctrs::vec_ptype_abbr(tlist), "lot")

  expect_null(vctrs::obj_print_data(tlist[0]))

})
