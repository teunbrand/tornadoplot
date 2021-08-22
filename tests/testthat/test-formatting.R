


test_that("sample data can be formatted", {

  d <- dummy_bigwig(pattern = "dmmy")

  o <- format_sample_data(d, "d", n = 1)

  expect_true(startsWith(o$sample_name, "dmmy"))

  dd <- rtracklayer::BigWigFileList(
    setNames(rep(path(d), 2), c("X", "Y"))
  )

  o <- format_sample_data(dd, "d", n = 2)

  expect_equal(o$sample_name, c("X", "Y"))

  d <- path(d)

  o <- format_sample_data(d, "d", n = 1)

  expect_match(o$sample_name, "test bigwig")

  unlink(d)

  d <- dummy_tabix("dmmy")

  o <- format_sample_data(d, "d", barcode_groups = list("A" = "A", "B" = "B"),
                          n = 2)
  expect_equal(o$sample_name, c("A", "B"))

  nms <- 10

  nms <- format_sample_name(names = nms, n = 2)
  expect_equal(nms, c("10_1", "10_2"))

})
