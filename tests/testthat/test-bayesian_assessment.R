#############################################
#  bayes_testing() Tests                    #
#############################################

test_that("bayes_testing gives correct posterior_model_prob",{

  # number of decimal places
  dec <- 5

  # check that bayes_testing gives correct posterior_model_prob for fivethirtyeight
  bayes_538 <- bayes_testing(hockey$x, hockey$y)
  expect_equal(round(bayes_538$posterior_model_prob, dec), round(0.9903632, dec))

  # check that bayes_testing gives correct posterior_model_prob for random noise
  bayes_rand <- bayes_testing(rand_pundit$x, rand_pundit$y)
  expect_equal(round(bayes_rand$posterior_model_prob, dec), round(0.0000000, dec))
})

test_that("bayes_testing gives correct Bayes factor",{

  # number of decimal places
  dec <- 5

  # check that bayes_testing gives correct BF for fivethirtyeight
  bayes_538 <- bayes_testing(hockey$x, hockey$y)
  expect_equal(round(bayes_538$BF, dec), round(0.009730538, dec))

})

test_that("bayes_testing gives correct BIC",{

  # number of decimal places
  dec <- 3

  # check that bayes_testing gives correct BIC for fivethirtyeight
  bayes_538 <- bayes_testing(hockey$x, hockey$y)
  expect_equal(round(bayes_538$BIC_H0, dec), round(1148.637, dec))
  expect_equal(round(bayes_538$BIC_H1, dec), round(1157.902, dec))

  # check that bayes_testing gives correct BIC for random noise
  bayes_rand <- bayes_testing(rand_pundit$x, rand_pundit$y)
  expect_equal(round(bayes_rand$BIC_H0, dec), round(1269.664, dec))
  expect_equal(round(bayes_rand$BIC_H1, dec), round(1212.527, dec))

})

test_that("bayes_testing gives correct est_params",{

  # number of decimal places
  dec <- 5

  # check that bayes_testing gives correct est_params for fivethirtyeight
  bayes_538 <- bayes_testing(hockey$x, hockey$y)
  expect_equal(round(bayes_538$est_params[1], dec), round(0.9453966 , dec))
  expect_equal(round(bayes_538$est_params[2], dec), round(1.4005730, dec))

  # check that bayes_testing gives correct est_params for random noise
  bayes_rand <- bayes_testing(rand_pundit$x, rand_pundit$y)
  expect_equal(round(bayes_rand$est_params[1], dec), round(1.13946217, dec))
  expect_equal(round(bayes_rand$est_params[2], dec), round(0.07199484, dec))
})






