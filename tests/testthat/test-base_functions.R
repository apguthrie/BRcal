
#############################################
#  LLO() Tests                              #
#############################################

test_that("LLO() only accepts single numeric inputs > 0 for delta", {
  # Set up
  set.seed(47)
  n <- 100
  g <- 2

  # delta = 0 - error
  d2 <- 0
  expect_error(LLO(x, d2, g))

  # delta < 0 - error
  d3 <- -5
  expect_error(LLO(x, d3, g))

  # very large delta - no error
  d4 <- 10000
  expect_no_error(LLO(x, d4, g))

  # very large delta & vector of 0s and 1s - no error
  x8 <- rbinom(n, 1, prob=x)
  expect_no_error(LLO(x8, d4, g))

  # character input for delta - error
  d5 <- "hello"
  expect_error(LLO(x, d5, g))
  d7 <- c("a", "b")
  expect_error(LLO(x, d7, g))

  # vector length > 1 for delta - error
  d6 <- c(1, 1)
  expect_error(LLO(x, d6, g))


})

test_that("LLO() only accepts single numeric inputs for gamma", {
  # Set up
  set.seed(47)
  n <- 100
  d8 <- 1

  # gamma = Inf - warning
  g2 <- Inf
  expect_warning(LLO(x, d8, g2))

  # Non-real number input for gamma
  g3 <- -1+2i
  expect_error(LLO(x, d8, g3))

  # very large gamma - no error
  g4 <- 10000
  expect_no_error(LLO(x, d8, g4))
  g5 <- -10000
  expect_no_error(LLO(x, d8, g5))

  # character input for gamma - error
  g6 <- "hello"
  expect_error(LLO(x, d8, g6))
  g8 <- c("a", "b")
  expect_error(LLO(x, d8, g8))

  # vector length > 1 for gamma
  g7 <- c(1, 1)
  expect_error(LLO(x, d8, g7))

})

test_that("LLO() only accepts p in correct format", {
  # Set up
  set.seed(47)
  n <- 100
  d <- 2
  g <- 2

  # Numeric vector input - no error
  x <- runif(n)
  expect_no_error(p <- LLO(x, d, g))

  # Numeric vector input with non [0,1] values - warning
  x2 <- rnorm(n)
  expect_warning(p2 <- LLO(x2, d, g))

  # Single numeric value - no error
  x3 <- 0.5
  expect_no_error(p3 <- LLO(x3, d, g))

  # Single character value - error
  x4 <- "0.5"
  expect_error(LLO(x4, d, g))

  # Single numeric value outside [0,1] - warning
  x5 <- 3
  expect_warning(p5 <- LLO(x5, d, g))
  x6 <- -3
  expect_warning(p6 <- LLO(x6, d, g))

  # Character vector input - error
  x7 <- c("a", "b")
  expect_error(LLO(x7, d, g))

  # Vector of 0s and 1s - no error
  x8 <- rbinom(n, 1, prob=x)
  expect_no_error(LLO(x8, d, g))


})

test_that("LLO() returns vector of correct size", {
  # Set up
  set.seed(47)
  n <- 100
  d <- 2
  g <- 2

  # Numeric vector input - vector of size n=100
  x <- runif(n)
  expect_vector(LLO(x, d, g), ptype=numeric(), size=n)

  # Numeric vector input with non [0,1] values - warning
  x2 <- rnorm(n)
  expect_vector(LLO(x2, d, g), ptype=numeric(), size=n)

  # Single numeric value - single element return
  x3 <- 0.5
  expect_vector(LLO(x3, d, g), ptype=numeric(), size=1)

  # Single numeric value outside [0,1] - warning, single element return
  x5 <- 3
  expect_vector(LLO(x5, d, g), ptype = numeric(), size=1)
  x6 <- -3
  expect_vector(LLO(x6, d, g), ptype = numeric(), size=1)

  # very large delta - vector of size n=100
  d4 <- 10000
  expect_vector(LLO(x, d4, g), ptype=numeric(), size=n)

  # very large delta & vector of 0s and 1s - no error
  x8 <- rbinom(n, 1, prob=x)
  expect_vector(LLO(x8, d4, g), ptype=numeric(), size=n)

  # gamma = Inf - warning
  g2 <- Inf
  expect_vector(LLO(x, d, g2), ptype=numeric(), size=n)

  # very large gamma - no error
  g4 <- 10000
  expect_vector(LLO(x, d, g4), ptype=numeric(), size=n)
  g5 <- -10000
  expect_vector(LLO(x, d, g5), ptype=numeric(), size=n)

})



#############################################
#  Tests                               #
#############################################

# test_that("LRT gives correct results",{
#
#   lrt_538 <- LLO_LRT(hockey$x, hockey$y)
#
#   # check that LLO_LRT gives correct p-value for fivethirtyeight
#   expect_equal(lrt_538$pval, 0.118396594)
#
#   # check test stat is consistent for both
#   expect_equal(lrt_538$test_stat, 4.2674306)
#
#   # Check params
#   expect_equal(lrt_538$est_params[1], 0.9453853)
#   expect_equal(lrt_538$est_params[2], 1.4014034)
#
# })
#
# test_that("LRT gives correct test stat",{
#
#   #lrt_538 <- LLO_LRT(hockey$x, hockey$y)
#
#   # check that LLO_LRT gives correct p-value for fivethirtyeight
#   #expect_equal(lrt_538$pval, 0.118396594)
#
#   # check test stat is consistent for both
#   #expect_equal(lrt_538$test_stat, 4.2674306)
#
#   # Check params
#   #expect_equal(lrt_538$est_params[1], 0.9453853)
#   #expect_equal(lrt_538$est_params[2], 1.4014034)
#
# })



