context("TFactSR: an R package for flexible TFactS analysis")


test_that("TFactS calculation", {
  data(DEGs)
  data(catalog)

  set.seed(123456789)

  tftg <- extractTFTG(DEGs, catalog)
  TFs <- tftg$TFs
  all.targets <- tftg$all.targets

  ## calculation of TFacts by DEMO datasets
  res <- calculateTFactS(DEGs, catalog, TFs, all.targets)

  ## test for calculated data.frame
  expect_identical(class(res), "data.frame")
  ## data frame size == 29
  expect_equal(dim(res)[1], 29)  ## Total number of tests

  ## TFactS values
  expect_equal(res[1, ]$m, 78)
  expect_equal(res[1, ]$n, 18)
  expect_equal(res[1, ]$N, 6838)
  expect_equal(res[1, ]$k, 7)
  expect_equal(res[1, ]$RC, 1)
  expect_equal(res[2, ]$RC, 1)
  expect_equal(res[3, ]$RC, 0)


  data(DEGs39)
  set.seed(1234)

  tftg <- extractTFTG(DEGs39, catalog)
  TFs <- tftg$TFs
  all.targets <- tftg$all.targets

  res39 <- calculateTFactS(DEGs39, catalog, TFs, all.targets)
  ## TFactS values
  expect_equal(res39[1, ]$m, 78)
  expect_equal(res39[1, ]$n, 39)
  expect_equal(res39[1, ]$N, 6838)
  expect_equal(res39[1, ]$k, 10)
  expect_equal(dim(res39)[1], 73)  ## Total number of tests
})

test_that("fast RC calculation", {
  data(DEGs)
  data(catalog)
  data(example.df)

  tftg <- extractTFTG(DEGs, catalog)
  TFs <- tftg$TFs
  all.targets <- tftg$all.targets

  ## calculation of TFacts by DEMO datasets
  set.seed(1234)
  res <- calculateRC(example.df, DEGs, catalog, TFs, all.targets)
  set.seed(1234)
  resRC <- FASTcalculateRC(example.df, DEGs, catalog, TFs, all.targets)

  ## TFactS values
  expect_equal(res[1, ]$RC, 0)
  expect_equal(res[2, ]$RC, 2)
  expect_equal(res[7, ]$RC, 6)
  ##
  expect_equal(resRC[1, ]$RC, 0)
  expect_equal(resRC[2, ]$RC, 2)
  expect_equal(resRC[7, ]$RC, 6)
})
