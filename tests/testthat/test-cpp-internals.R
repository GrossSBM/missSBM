context("test internal C++ helpers (nlopt wrapper, packing/unpacking)")

test_that("nlopt wrapper internals are consistent", {
  expect_true(missSBM:::cpp_test_nlopt())
})

test_that("packing/unpacking internals are consistent", {
  expect_true(missSBM:::cpp_test_packing())
})
