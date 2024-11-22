
#############################################
#  LLO() Tests                              #
#############################################

test_that("LLO() only takes valid x",{
  x <-  rnorm(10, mean=100)
  expect_error(LLO(x, 1, 1))
  
  x <-  rnorm(10, mean=-100)
  expect_error(LLO(x, 1, 1))
  
  expect_error(LLO(TRUE, 1, 1))
  expect_error(LLO("1", 1, 1))
  
  x <- runif(10)
  expect_warning(LLO(matrix(x),1,1))
  expect_error(LLO(as.data.frame(x=x),1,1))
})

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
#  llo_lrt() Tests                          #
#############################################

test_that("llo_lrt() only takes valid y",{
  x <- runif(10)
  y <- rbinom(10, 1, x)
  expect_no_condition(llo_lrt(x,y))

  expect_error(llo_lrt(x,-y))
  expect_error(llo_lrt(x, c(rep("hi", 5), rep("bye", 5))))
  expect_no_condition(llo_lrt(x, c(rep("hi", 5), rep("bye", 5)), event="hi"))
})


test_that("llo_lrt() only accepts valid input for optim_details",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_no_condition(llo_lrt(x,y,optim_details = TRUE))
  expect_no_condition(llo_lrt(x,y,optim_details = 1))
  expect_no_condition(llo_lrt(x,y,optim_details = FALSE))
  expect_no_condition(llo_lrt(x,y,optim_details = 0))
  expect_no_condition(llo_lrt(x,y,optim_details = T))
  expect_no_condition(llo_lrt(x,y,optim_details = F))
  
  expect_error(llo_lrt(x,y,optim_details = 10))
  expect_error(llo_lrt(x,y,optim_details = c(T, F)))
  expect_error(llo_lrt(x,y,optim_details = c(2, 4)))
  expect_error(llo_lrt(x,y,optim_details = "TRUE"))
  expect_error(llo_lrt(x,y,optim_details = c()))
}) 

test_that("llo_lrt() only accepts x & y of the same length",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_error(llo_lrt(x,c(y,y)))
}) 

test_that("llo_lrt() only valid input of epsilon",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_error(llo_lrt(x,y, epsilon=10))
  expect_error(llo_lrt(x,y, epsilon=-1))
  expect_error(llo_lrt(x,y, epsilon="0.3"))
  expect_error(llo_lrt(x,y, epsilon=c(0.2, 0.3)))
  expect_error(llo_lrt(x,y, epsilon=TRUE))
}) 

test_that("llo_lrt() only valid input of event",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_error(llo_lrt(x,y, event=10))
  expect_error(llo_lrt(x,y, event=-1))
  expect_error(llo_lrt(x,y, event="0.3"))
  expect_error(llo_lrt(x,y, event=c(0.2, 0.3)))
  #expect_error(llo_lrt(x,y, event=TRUE))
}) 


test_that("llo_lrt() gives correct p-value",{
  
  # number of decimal places
  dec <- 5
  
  # check that llo_lrt gives correct p-value for fivethirtyeight
  lrt_538 <- llo_lrt(hockey$x, hockey$y)
  expect_equal(round(lrt_538$pval, dec), round(0.118396594, dec))
  
  # check that llo_lrt gives correct p-value for random noise
  lrt_rand <- llo_lrt(hockey$rand, hockey$y)
  expect_equal(round(lrt_rand$pval, dec), round(0.0000000, dec))
})

test_that("llo_lrt() gives correct test stat",{
  
  # number of decimal places
  dec <- 5
  
  # check that llo_lrt gives correct test_stat for fivethirtyeight
  lrt_538 <- llo_lrt(hockey$x, hockey$y)
  expect_equal(round(lrt_538$test_stat, dec), round(4.267411, dec))
  
  # check that llo_lrt gives correct test_stat for random noise
  lrt_rand <- llo_lrt(hockey$rand, hockey$y)
  expect_equal(round(lrt_rand$test_stat, dec), round(70.66915, dec))
})

test_that("llo_lrt() gives correct MLEs",{
  
  # number of decimal places
  dec <- 5
  
  # check that llo_lrt gives correct MLEs for fivethirtyeight
  lrt_538 <- llo_lrt(hockey$x, hockey$y)
  expect_equal(round(lrt_538$MLEs[1], dec), round(0.9453966, dec))
  expect_equal(round(lrt_538$MLEs[2], dec), round(1.4005730, dec))
  
  # check that llo_lrt gives correct MLEs for random noise
  lrt_rand <- llo_lrt(hockey$rand, hockey$y)
  expect_equal(round(lrt_rand$MLEs[1], dec), round(1.13946217, dec))
  expect_equal(round(lrt_rand$MLEs[2], dec), round(0.07199484, dec))
})


#############################################
#  mle_recal() Tests                        #
#############################################

test_that("mle_recal() only accepts valid input for optim_details",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_no_condition(mle_recal(x,y,optim_details = TRUE))
  expect_no_condition(mle_recal(x,y,optim_details = 1))
  expect_no_condition(mle_recal(x,y,optim_details = FALSE))
  expect_no_condition(mle_recal(x,y,optim_details = 0))
  expect_no_condition(mle_recal(x,y,optim_details = T))
  expect_no_condition(mle_recal(x,y,optim_details = F))
  
  expect_error(mle_recal(x,y,optim_details = 10))
  expect_error(mle_recal(x,y,optim_details = c(T, F)))
  expect_error(mle_recal(x,y,optim_details = c(2, 4)))
  expect_error(mle_recal(x,y,optim_details = "TRUE"))
  expect_error(mle_recal(x,y,optim_details = c()))
}) 

test_that("mle_recal() only accepts valid input for probs_only",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_no_condition(mle_recal(x,y,probs_only = TRUE))
  expect_no_condition(mle_recal(x,y,probs_only = 1))
  expect_no_condition(mle_recal(x,y,probs_only = FALSE))
  expect_no_condition(mle_recal(x,y,probs_only = 0))
  expect_no_condition(mle_recal(x,y,probs_only = T))
  expect_no_condition(mle_recal(x,y,probs_only = F))
  
  expect_error(mle_recal(x,y,probs_only = 10))
  expect_error(mle_recal(x,y,probs_only = c(T, F)))
  expect_error(mle_recal(x,y,probs_only = c(2, 4)))
  expect_error(mle_recal(x,y,probs_only = "TRUE"))
  expect_error(mle_recal(x,y,probs_only = c()))
}) 

test_that("mle_recal() only accepts x & y of the same length",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_error(mle_recal(x,c(y,y)))
}) 


test_that("mle_recal() returns vector of correct size", {
  # Set up
  set.seed(47)
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  # Numeric vector input - vector of size n=100
  expect_no_condition(l <- mle_recal(x, y))
  expect_true(check_probs(l$probs))
  expect_vector(l$probs)
})

test_that("mle_recal() gives correct MLEs",{
  
  # number of decimal places
  dec <- 5
  
  # check that llo_lrt gives correct test_stat for fivethirtyeight
  mle_538 <- mle_recal(hockey$x, hockey$y)
  expect_equal(round(mle_538$MLEs[1], dec), round(0.9453966, dec))
  expect_equal(round(mle_538$MLEs[2], dec), round(1.4005730, dec))
})

test_that("mle_recal() accepts both probs_only and optim_details=TRUE",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_no_condition(mle_recal(x,y, optim_details = TRUE, probs_only = TRUE))
}) 

#############################################
#  logit() Tests                            #
#############################################

test_that("(INTERNAL) logit() works ok",{
  x <- runif(10)
  expect_no_condition(logit(x))
  
  expect_error(logit(x, epsilon="hi"))
})

#############################################
#  to_prob() Tests                            #
#############################################

test_that("(INTERNAL) to_prob() works ok",{
  x <- runif(10)
  expect_no_condition(to_prob(x))
  
  expect_error(to_prob(x, epsilon="hi"))
})

#############################################
#  llo_lik() Tests                          #
#############################################

test_that("(INTERNAL) llo_lik() works with log=FALSE",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_no_condition(llo_lik(params=c(1,1),x,y))
  expect_no_condition(llo_lik(params=c(1,1),x,y, log=TRUE))
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











