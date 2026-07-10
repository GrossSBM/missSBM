context("test internal C++ helpers (nlopt wrapper, packing/unpacking)")

test_that("nlopt wrapper internals are consistent", {
  expect_true(missSBM:::cpp_test_nlopt())
})

test_that("an exception thrown inside the nlopt objective is caught, not crashing the session", {
  expect_true(missSBM:::cpp_test_nlopt_exception_safety())
})

test_that("packing/unpacking internals are consistent", {
  expect_true(missSBM:::cpp_test_packing())
})
