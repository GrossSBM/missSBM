context("test repair_empty_classes()")

repair_empty_classes <- missSBM:::repair_empty_classes

test_that("repair_empty_classes() fills a single absent label without emptying any other class", {
  set.seed(1)
  labels <- c(1, 1, 1, 1, 2, 2, 3, 3, 5) # label 4 is absent, K = 5
  counts_before <- tabulate(labels, nbins = 5)
  out <- repair_empty_classes(labels, 5)

  expect_length(out, length(labels))
  expect_true(all(tabulate(out, nbins = 5) > 0))
  ## every class that had >= 1 member before still has >= 1 member after
  expect_true(all(tabulate(out, nbins = 5)[counts_before > 0] > 0))
})

test_that("repair_empty_classes() fills several absent labels at once", {
  set.seed(2)
  labels <- c(1, 1, 1, 1, 1, 1, 2, 2) # labels 3, 4 absent, K = 4
  counts_before <- tabulate(labels, nbins = 4)
  out <- repair_empty_classes(labels, 4)

  expect_length(out, length(labels))
  expect_true(all(tabulate(out, nbins = 4) > 0))
  expect_true(all(tabulate(out, nbins = 4)[counts_before > 0] > 0))
})

test_that("repair_empty_classes() is a no-op when every class is already occupied", {
  labels <- c(1, 2, 3, 1, 2, 3)
  out <- repair_empty_classes(labels, 3)
  expect_equal(tabulate(out, nbins = 3), tabulate(labels, nbins = 3))
})

test_that("repair_empty_classes() never empties an already-valid class, even repeatedly", {
  set.seed(3)
  for (i in 1:20) {
    labels <- c(rep(1, 10), 2, 4) # labels 3, 5 absent, K = 5
    out <- repair_empty_classes(labels, 5)
    expect_true(all(tabulate(out, nbins = 5) > 0))
  }
})
