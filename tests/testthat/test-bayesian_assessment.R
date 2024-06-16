#############################################
#  bayes_ms() Tests                         #
#############################################

test_that("bayes_ms gives correct posterior_model_prob",{

  # number of decimal places
  dec <- 5

  # check that bayes_ms gives correct posterior_model_prob for fivethirtyeight
  bayes_538 <- bayes_ms(hockey$x, hockey$y)
  expect_equal(round(bayes_538$posterior_model_prob, dec), round(0.9903632, dec))

  # check that bayes_ms gives correct posterior_model_prob for random noise
  bayes_rand <- bayes_ms(rand_pundit$x, rand_pundit$y)
  expect_equal(round(bayes_rand$posterior_model_prob, dec), round(0.0000000, dec))
})

test_that("bayes_ms gives correct Bayes factor",{

  # number of decimal places
  dec <- 5

  # check that bayes_ms gives correct BF for fivethirtyeight
  bayes_538 <- bayes_ms(hockey$x, hockey$y)
  expect_equal(round(bayes_538$BF, dec), round(0.009730538, dec))

})

test_that("bayes_ms gives correct BIC",{

  # number of decimal places
  dec <- 3

  # check that bayes_ms gives correct BIC for fivethirtyeight
  bayes_538 <- bayes_ms(hockey$x, hockey$y)
  expect_equal(round(bayes_538$BIC_H0, dec), round(1148.637, dec))
  expect_equal(round(bayes_538$BIC_H1, dec), round(1157.902, dec))

  # check that bayes_ms gives correct BIC for random noise
  bayes_rand <- bayes_ms(rand_pundit$x, rand_pundit$y)
  expect_equal(round(bayes_rand$BIC_H0, dec), round(1269.664, dec))
  expect_equal(round(bayes_rand$BIC_H1, dec), round(1212.527, dec))

})

test_that("bayes_ms gives correct est_params",{

  # number of decimal places
  dec <- 5

  # check that bayes_ms gives correct est_params for fivethirtyeight
  bayes_538 <- bayes_ms(hockey$x, hockey$y)
  expect_equal(round(bayes_538$est_params[1], dec), round(0.9453966 , dec))
  expect_equal(round(bayes_538$est_params[2], dec), round(1.4005730, dec))

  # check that bayes_ms gives correct est_params for random noise
  bayes_rand <- bayes_ms(rand_pundit$x, rand_pundit$y)
  expect_equal(round(bayes_rand$est_params[1], dec), round(1.13946217, dec))
  expect_equal(round(bayes_rand$est_params[2], dec), round(0.07199484, dec))
})


test_that("bayes_ms() only accepts valid input for optim_details",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_no_condition(bayes_ms(x,y,optim_details = TRUE))
  expect_no_condition(bayes_ms(x,y,optim_details = 1))
  expect_no_condition(bayes_ms(x,y,optim_details = FALSE))
  expect_no_condition(bayes_ms(x,y,optim_details = 0))
  expect_no_condition(bayes_ms(x,y,optim_details = T))
  expect_no_condition(bayes_ms(x,y,optim_details = F))
  
  expect_error(bayes_ms(x,y,optim_details = 10))
  expect_error(bayes_ms(x,y,optim_details = c(T, F)))
  expect_error(bayes_ms(x,y,optim_details = c(2, 4)))
  expect_error(bayes_ms(x,y,optim_details = "TRUE"))
  expect_error(bayes_ms(x,y,optim_details = c()))
}) 


test_that("bayes_ms() only accepts x & y of the same length",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_error(bayes_ms(x,c(y,y)))
}) 
