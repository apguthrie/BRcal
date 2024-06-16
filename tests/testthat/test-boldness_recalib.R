#############################################
#  brcal() Tests                            #
#############################################

test_that("brcal() only accepts valid input for tau",{
  set.seed(100)
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_no_condition(brcal(x,y,tau = TRUE, opts = list(maxtime=0.1)))
  expect_no_condition(brcal(x,y,tau = 1, opts = list(maxtime=0.1)))
  expect_no_condition(brcal(x,y,tau = FALSE, opts = list(maxtime=0.1)))
  expect_no_condition(brcal(x,y,tau = 0, opts = list(maxtime=0.1)))
  expect_no_condition(brcal(x,y,tau = T, opts = list(maxtime=0.1)))
  expect_no_condition(brcal(x,y,tau = F, opts = list(maxtime=0.1)))
  
  expect_error(brcal(x,y,tau = 10, opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,tau = c(T, F), opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,tau = c(2, 4), opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,tau = "TRUE", opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,tau = c(), opts = list(maxtime=0.1)))
}) 


test_that("brcal() only accepts valid input for start_at_MLEs",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_no_condition(brcal(x,y,start_at_MLEs = TRUE, opts = list(maxtime=0.1)))
  expect_no_condition(brcal(x,y,start_at_MLEs = 1, opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = FALSE, opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = 0, opts = list(maxtime=0.1)))
  expect_no_condition(brcal(x,y,start_at_MLEs = T, opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = F, opts = list(maxtime=0.1)))
  
  expect_error(brcal(x,y,start_at_MLEs = 10, opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = c(T, F), opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = c(2, 4), opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = "TRUE", opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = c(), opts = list(maxtime=0.1)))
}) 

test_that("brcal() only accepts x & y of the same length",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_error(brcal(x,c(y,y), opts = list(maxtime=0.1)))
}) 

test_that("brcal() only accepts valid starting locations",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_error(brcal(x,y,start_at_MLEs = F, x0=c(-1, 0), opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = F, x0=c(0, -2), opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = F, x0=c(1, Inf), opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = F, x0=c(1), opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = F, x0=c(1, 1, 1), opts = list(maxtime=0.1)))
  expect_error(brcal(x,y,start_at_MLEs = F, x0=c("hi"), opts = list(maxtime=0.1)))
})

test_that("brcal() doesn't need to have opts specified",{
  set.seed(47)
  x <- runif(50)
  y <- rbinom(50,1,x)
  expect_no_condition(brcal(x,y))
})

test_that("brcal() opts",{
  set.seed(47)
  x <- runif(50)
  y <- rbinom(50,1,x)
  expect_no_condition(brcal(x,y, opts=list(maxeval=1)))
  expect_no_condition(brcal(x,y, opts=list(algorithm="NLOPT_LD_MMA", maxtime=0.1)))
  expect_no_condition(brcal(x,y, opts=list(local_opts=list(eval_f = sd), maxtime=0.1)))
  
})


test_that("brcal() gives correct output",{
  
  # number of decimal places
  dec <- 5
  
  br_538 <- brcal(hockey$x, hockey$y, t=0.9, Pmc=0.25)
  
  # Check Pmc and t
  expect_equal(br_538$Pmc, 0.25)
  expect_equal(br_538$t, 0.9)
  
  # check that bayes_ms gives correct est_params for fivethirtyeight
  expect_equal(round(br_538$BR_params[1], dec), round(0.8763524, dec))
  expect_equal(round(br_538$BR_params[2], dec), round(1.9317682, dec))
  expect_equal(round(br_538$sb, dec), round(0.1634752, dec))
  # check probabilities
  # Numeric vector input - vector of size n=100
  expect_true(check_probs(br_538$probs))
  expect_vector(br_538$probs)
  expect_equal(length(br_538$probs), nrow(hockey))
  
})
