
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

})

test_that("LLO() warns when NaNs returned", {
  # Set up
  set.seed(127)
  n <- 100
  d <- 2
  g <- 2
  x <- runif(n)

  # very large delta - vector of size n=100
  d4 <- 10000
  expect_vector(p <- LLO(x, d4, g), ptype=numeric(), size=n)
  expect_true(check_probs(p))
  d8 <- Inf
  expect_warning(p <- LLO(x, d8, g))
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
#  prelec() Tests                           #
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

test_that("prelec() only accepts single numeric inputs > 0 for beta", {
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
#  llo_lik() Tests                          #
#############################################


test_that("llo_lik() only accepts valid params",{
  set.seed(37)
  n <- 100
  x <- runif(n)
  y <- rbinom(n, 1, prob=x)

  # params not length 2
  params4 <- c(1)
  expect_error(llo_lik(params4, x, y))
  params5 <- c(1, 1, 2, 3)
  expect_error(llo_lik(params5, x, y))

  # delta <= 0
  params2 <- c(0, 1)
  expect_error(llo_lik(params2, x, y))
  params3 <- c(-10, 1)
  expect_error(llo_lik(params3, x, y))

  # delta non-numeric
  params6 <- c(TRUE, FALSE)
  expect_error(llo_lik(params6, x, y))
  params7 <- c("10", 1)
  expect_error(llo_lik(params7, x, y))
  params8 <- c(list(1,2), 1)
  expect_warning(expect_error(llo_lik(params8, x, y)))

  # gamma

})

test_that("llo_lik() only accepts x in correct format",{
  set.seed(37)
  n <- 100
  x <- runif(n)
  y <- rbinom(n, 1, prob=x)
  params <- c(1, 1)

  # x has values outside [0,1]
  x2 <- rnorm(n)
  expect_error(llo_lik(params, x, y2))

  # y is character vector
  y2 <- c(0, 1)
  x3 <- c("h", "f")
  expect_error(llo_lik(params, x3, y2))

  # y is logical vector
  x4 <- c(TRUE, FALSE)
  expect_error(llo_lik(params, x4, y2))

  # y is list - warning
  x5 <- list(runif(n))
  expect_warning(llo_lik(params, x5, y))
  x6 <- list(0.5, 0.2)
  expect_warning(llo_lik(params, x6, y2))

  # y is matrix - error (diff lengths) & warning (wrong type)
  x7 <- matrix(c(0.2, 1, 0.7532),ncol=1)
  expect_error(llo_lik(params, x7, y2))

  # y is matrix - warning (wrong type)
  x8 <- matrix(c(0.111, 0), ncol=1)
  expect_error(llo_lik(params, x8, y2))
})

test_that("llo_lik() only accepts y in correct format",{
  set.seed(37)
  n <- 100
  x <- runif(n)
  params <- c(1,1)

  # y has non 0 or 1s
  y2 <- rnorm(n)
  expect_error(llo_lik(params, x, y2))

  # y is character vector
  x2 <- c(0.5, 0.1)
  y3 <- c("h", "f")
  expect_error(llo_lik(params, x2, y3))

  # y is logical vector
  y4 <- c(TRUE, FALSE)
  expect_error(llo_lik(params, x2, y3))

  # y is list - warning
  y5 <- list(rbinom(n, 1, prob=x))
  expect_warning(llo_lik(params, x, y5))
  y6 <- list(1, 0)
  expect_warning(llo_lik(params, x2, y6))

  # y is matrix - error (diff lengths) & warning (wrong type)
  y7 <- matrix(c(0,1,1),ncol=1)
  expect_error(llo_lik(params, x2, y7))

  # y is matrix - warning (wrong type)
  y8 <- matrix(c(0,1),ncol=1)
  expect_error(llo_lik(params, x2, y8))
})


test_that("llo_lik() only accepts x and y of the same length",{
  set.seed(37)
  n <- 100
  x <- runif(n)
  y <- rbinom(n, 1, prob=x)
  params <- c(1,1)

  # same length - no error
  expect_no_error(llo_lik(params, x, y))

  # y shorter than x - error
  y2 <- c(0, 1,0)
  expect_error(llo_lik(params, x, y2))

  # x shorter than y - error
  x2 <-c(0.2, 0.4, 0.7)
  expect_error(llo_lik(params, x2, y))
})

test_that("llo_lik() only accepts log in correct format", {
  # set up
  set.seed(37)
  n <- 100
  x <- runif(n)
  y <- rbinom(n, 1, prob=x)
  params <- c(1,1)

  # character input - error
  expect_error(llo_lik(params, x, y, log="yes"))

  # numeric input - error
  expect_error(llo_lik(params, x, y, log=2))

  # logical input - no error
  expect_no_error(llo_lik(params, x, y, log=TRUE))
  expect_no_error(llo_lik(params, x, y, log=FALSE))
  expect_no_error(llo_lik(params, x, y, log=T))
  expect_no_error(llo_lik(params, x, y, log=F))

  # numeric but technically ok
  expect_no_error(llo_lik(params, x, y, log=0))
  expect_no_error(llo_lik(params, x, y, log=1))
})

test_that("llo_lik() only accepts neg in correct format",{
  # set up
  set.seed(37)
  n <- 100
  x <- runif(n)
  y <- rbinom(n, 1, prob=x)
  params <- c(1,1)

  # character input - error
  expect_error(llo_lik(params, x, y, neg="yes"))

  # numeric input - error
  expect_error(llo_lik(params, x, y, neg=2))

  # logical input - no error
  expect_no_error(llo_lik(params, x, y, neg=TRUE))
  expect_no_error(llo_lik(params, x, y, neg=FALSE))
  expect_no_error(llo_lik(params, x, y, neg=T))
  expect_no_error(llo_lik(params, x, y, neg=F))

  # numeric but technically ok
  expect_no_error(llo_lik(params, x, y, neg=0))
  expect_no_error(llo_lik(params, x, y, neg=1))
})

test_that("llo_lik() only accepts tau in correct format",{
  # set up
  set.seed(37)
  n <- 100
  x <- runif(n)
  y <- rbinom(n, 1, prob=x)
  params <- c(1,1)

  # character input - error
  expect_error(llo_lik(params, x, y, tau="yes"))

  # numeric input - error
  expect_error(llo_lik(params, x, y, tau=2))

  # logical input - no error
  expect_no_error(llo_lik(params, x, y, tau=TRUE))
  expect_no_error(llo_lik(params, x, y, tau=FALSE))
  expect_no_error(llo_lik(params, x, y, tau=T))
  expect_no_error(llo_lik(params, x, y, tau=F))

  # numeric but technically ok
  expect_no_error(llo_lik(params, x, y, tau=0))
  expect_no_error(llo_lik(params, x, y, tau=1))
})

#############################################
#  LLO_LRT() Tests                          #
#############################################

test_that("LLO_LRT() only accepts valid params",{
  set.seed(37)
  n <- 100
  x <- runif(n)
  y <- rbinom(n, 1, prob=x)

  # params not length 2
  params4 <- c(1)
  expect_error(LLO_LRT(x, y, params4))
  params5 <- c(1, 1, 2, 3)
  expect_error(LLO_LRT(x, y, params5))

  # delta <= 0
  params2 <- c(0, 1)
  expect_error(LLO_LRT(x, y, params2))
  params3 <- c(-10, 1)
  expect_error(LLO_LRT(x, y, params3))

  # delta non-numeric
  params6 <- c(TRUE, FALSE)
  expect_error(LLO_LRT(x, y, params6))
  params7 <- c("10", 1)
  expect_error(LLO_LRT(x, y, params7))
  params8 <- c(list(1,2), 1)
  expect_warning(expect_error(LLO_LRT(x, y, params8)))

  # gamma

})
#
# test_that("LLO_LRT() only takes valid outcomes",{
#   set.seed(37)
#   n <- 100
#   x <- runif(n)
#   params <- c(1,1)
#
#   # y has non 0 or 1s
#   y2 <- rnorm(n)
#   expect_error(llo_lik(params, x, y2))
#
#   # y is character vector
#   x2 <- c(0.5, 0.1)
#   y3 <- c("h", "f")
#   expect_error(llo_lik(params, x2, y3))
#
#   # y is logical vector
#   y4 <- c(TRUE, FALSE)
#   expect_error(llo_lik(params, x2, y3))
#
#   # y is list - warning
#   y5 <- list(rbinom(n, 1, prob=x))
#   expect_warning(llo_lik(params, x, y5))
#   y6 <- list(1, 0)
#   expect_warning(llo_lik(params, x2, y6))
#
#   # y is matrix - error (diff lengths) & warning (wrong type)
#   y7 <- matrix(c(0,1,1),ncol=1)
#   expect_error(llo_lik(params, x2, y7))
#
#   # y is matrix - warning (wrong type)
#   y8 <- matrix(c(0,1),ncol=1)
#   expect_error(llo_lik(params, x2, y8))
#
# })

test_that("LLO_LRT() gives correct p-value",{

  # number of decimal places
  dec <- 5

  # check that LLO_LRT gives correct p-value for fivethirtyeight
  lrt_538 <- LLO_LRT(hockey$x, hockey$y)
  expect_equal(round(lrt_538$pval, dec), round(0.118396594, dec))

  # check that LLO_LRT gives correct p-value for random noise
  lrt_rand <- LLO_LRT(rand_pundit$x, rand_pundit$y)
  expect_equal(round(lrt_rand$pval, dec), round(0.0000000, dec))
})

test_that("LLO_LRT() gives correct test stat",{

  # number of decimal places
  dec <- 5

  # check that LLO_LRT gives correct test_stat for fivethirtyeight
  lrt_538 <- LLO_LRT(hockey$x, hockey$y)
  expect_equal(round(lrt_538$test_stat, dec), round(4.267411, dec))

  # check that LLO_LRT gives correct test_stat for random noise
  lrt_rand <- LLO_LRT(rand_pundit$x, rand_pundit$y)
  expect_equal(round(lrt_rand$test_stat, dec), round(70.66915, dec))
})

test_that("LLO_LRT() gives correct est_params",{

  # number of decimal places
  dec <- 5

  # check that LLO_LRT gives correct est_params for fivethirtyeight
  lrt_538 <- LLO_LRT(hockey$x, hockey$y)
  expect_equal(round(lrt_538$est_params[1], dec), round(0.9453966, dec))
  expect_equal(round(lrt_538$est_params[2], dec), round(1.4005730, dec))

  # check that LLO_LRT gives correct est_params for random noise
  lrt_rand <- LLO_LRT(rand_pundit$x, rand_pundit$y)
  expect_equal(round(lrt_rand$est_params[1], dec), round(1.13946217, dec))
  expect_equal(round(lrt_rand$est_params[2], dec), round(0.07199484, dec))
})






