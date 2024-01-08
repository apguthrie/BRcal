
#############################################
#  LLO() Tests                              #
#############################################

test_that("LLO() only accepts single numeric inputs > 0 for delta", {
  # Set up
  set.seed(47)
  n <- 100
  x <- runif(n)
  g <- 2

  # delta <= 0 - error
  d2 <- 0
  expect_error(LLO(x, d2, g))
  d3 <- -5
  expect_error(LLO(x, d3, g))

  # very large delta - no error
  d4 <- 10000
  expect_no_error(LLO(x, d4, g))
  d8 <- Inf
  expect_warning(LLO(x, d8, g))

  # character input for delta - error
  d5 <- "hello"
  expect_error(LLO(x, d5, g))
  d7 <- c("a", "b")
  expect_error(LLO(x, d7, g))

  # vector length > 1 for delta - error
  d6 <- c(1, 1)
  expect_error(LLO(x, d6, g))

  # Non-real number input for delta
  d9 <- -1+2i
  expect_error(LLO(x, d9, g))

})

test_that("LLO() only accepts single numeric inputs for gamma", {
  # Set up
  set.seed(47)
  n <- 100
  x <- runif(n)
  d8 <- 1

  # Non-real number input for gamma
  g3 <- -1+2i
  expect_error(LLO(x, d8, g3))

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

  # Single numeric value - no error
  x3 <- 0.5
  expect_no_error(p3 <- LLO(x3, d, g))

  # Numeric vector input with non [0,1] values - error
  x2 <- rnorm(n)
  expect_error(p2 <- LLO(x2, d, g))
  x5 <- 3
  expect_error(p5 <- LLO(x5, d, g))
  x6 <- -3
  expect_error(p6 <- LLO(x6, d, g))

  # Character vector input - error
  x4 <- "0.5"
  expect_error(LLO(x4, d, g))
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
  expect_vector(p <- LLO(x, d, g), ptype=numeric(), size=n)
  expect_true(check_probs(p))

  # Single numeric value - single element return
  x3 <- 0.5
  expect_vector(p <- LLO(x3, d, g), ptype=numeric(), size=1)
  expect_true(check_probs(p))

  # very large delta - vector of size n=100
  d4 <- 10000
  expect_vector(p <- LLO(x, d4, g), ptype=numeric(), size=n)
  expect_true(check_probs(p))
  d8 <- Inf
  expect_warning(p <-LLO(x, d8, g))
  expect_vector(p, ptype=numeric(), size=n)

  # very large gamma - warnings
  g2 <- Inf
  expect_warning(p <- LLO(x, d, g2))
  expect_vector(p, ptype=numeric(), size=n)
  g4 <- 10000
  expect_warning(p <- LLO(x, d, g4))
  expect_vector(p, ptype=numeric(), size=n)
  g5 <- -10000
  expect_warning(p <- LLO(x, d, g5))
  expect_vector(p, ptype=numeric(), size=n)
  g9 <- -Inf
  expect_warning(p <- LLO(x, d, g9))
  expect_vector(p, ptype=numeric(), size=n)

  # very large delta & vector of 0s and 1s - no error
  x8 <- rbinom(n, 1, prob=x)
  expect_vector(p <- LLO(x8, d4, g), ptype=numeric(), size=n)
  expect_true(check_probs(p))

})

#############################################
#  prelec() Tests                              #
#############################################

test_that("prelec() only accepts single numeric inputs > 0 for alpha", {
  # Set up
  set.seed(47)
  n <- 100
  x <- runif(n)
  b <- 2

  # alpha <= 0 - error
  a2 <- 0
  expect_error(prelec(x, a2, b))
  a3 <- -5
  expect_error(prelec(x, a3, b))

  # very large alpha - no error
  a4 <- 10000
  expect_no_condition(prelec(x, a4, b))
  a8 <- Inf
  expect_no_condition(prelec(x, a8, b))

  # character input for alpha - error
  a5 <- "hello"
  expect_error(prelec(x, a5, b))
  a7 <- c("a", "b")
  expect_error(prelec(x, a7, b))

  # vector length > 1 for alpha - error
  a6 <- c(1, 1)
  expect_error(prelec(x, a6, b))

})

test_that("Prelec() only accepts single numeric inputs > 0 for beta", {
  # Set up
  set.seed(47)
  n <- 100
  x <- runif(n)
  a <- 2

  # beta <= 0 - error
  b2 <- 0
  expect_error(prelec(x, a, b2))
  b3 <- -5
  expect_error(prelec(x, a, b3))

  # very large beta - no error
  b4 <- 10000
  expect_no_condition(prelec(x, a, b4))
  b8 <- Inf
  expect_no_condition(prelec(x, a, b8))

  # character input for beta - error
  b5 <- "hello"
  expect_error(prelec(x, a, b5))
  b7 <- c("a", "b")
  expect_error(prelec(x, a, b7))

  # vector length > 1 for beta - error
  b6 <- c(1, 1)
  expect_error(prelec(x, a, b6))

})

test_that("prelec() only accepts p in correct format", {
  # Set up
  set.seed(47)
  n <- 100
  a <- 2
  b <- 2

  # Numeric vector input - no error
  x <- runif(n)
  expect_no_error(prelec(x, a, b))
  x3 <- 0.5
  expect_no_error(prelec(x3, a, b))

  # Numeric vector input with non [0,1] values - error
  x2 <- rnorm(n)
  expect_error(prelec(x2, a, b))
  x5 <- 3
  expect_error(prelec(x5, a, b))
  x6 <- -3
  expect_error(prelec(x6, a, b))

  # Character vector input - error
  x4 <- "0.5"
  expect_error(prelec(x4, a, b))
  x7 <- c("a", "b")
  expect_error(prelec(x7, a, b))

  # Vector of 0s and 1s - no error
  x8 <- rbinom(n, 1, prob=x)
  expect_no_error(prelec(x8, a, b))


})

test_that("prelec() returns valid output", {
  # Set up
  set.seed(47)
  n <- 100
  a <- 2
  b <- 2

  # Numeric vector input - vector of size n=100
  x <- runif(n)
  expect_vector(p <- prelec(x, a, b), ptype=numeric(), size=n)
  expect_true(check_probs(p))

  # Single numeric value - single element return
  x3 <- 0.5
  expect_vector(p <- prelec(x3, a, b), ptype=numeric(), size=1)
  expect_true(check_probs(p))

  # very large alpha - vector of size n=100
  a4 <- 10000
  expect_vector(p <- prelec(x, a4, b), ptype=numeric(), size=n)
  expect_true(check_probs(p))
  a8 <- Inf
  expect_vector(p <- prelec(x, a8, b), ptype=numeric(), size=n)
  expect_true(check_probs(p))

  # very large beta - vector of size n=100
  b4 <- 10000
  expect_vector(p <- prelec(x, a, b4), ptype=numeric(), size=n)
  expect_true(check_probs(p))
  b8 <- Inf
  expect_vector(p <- prelec(x, a, b8), ptype=numeric(), size=n)
  expect_true(check_probs(p))

  # very large alpha and beta
  expect_vector(p <- prelec(x, a4, b4), ptype=numeric(), size=n)
  expect_true(check_probs(p))

  expect_warning(p <- prelec(x, a8, b8))
  expect_vector(p, ptype=numeric(), size=n)

  # vector of 0s and 1s
  x8 <- rbinom(n, 1, prob=x)
  expect_vector(prelec(x8, a, b), ptype=numeric(), size=n)

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



