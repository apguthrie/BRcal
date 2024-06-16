#############################################
#  plot_params() Tests                      #
#############################################

test_that("plot_params() for hockey data doesn't break",{
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=10))
})

test_that("plot_params() only allows valid grid values",{
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=-1))
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=0))
  expect_error(plot_params(x=hockey$x, y=hockey$y, k="hi"))
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=TRUE))
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=c(5,5)))
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=Inf))
  
})

test_that("plot_params() requires either z or x&y",{
  expect_error(plot_params())
})

test_that("plot_params() only allows valid t_levels",{
  expect_error(plot_params(x=hockey$x, y=hockey$y, t_levels=c(-1)))
  expect_error(plot_params(x=hockey$x, y=hockey$y, t_levels=c(Inf)))
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, t_levels=c(0.9, 0.8, 0.4), k=10))
  expect_error(plot_params(x=hockey$x, y=hockey$y, t_levels=c(TRUE)))
  expect_error(plot_params(x=hockey$x, y=hockey$y, t_levels=c("hi")))
  expect_error(plot_params(x=hockey$x, y=hockey$y, t_levels=c(100)))

})


test_that("plot_params() only accepts valid input for return_z",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=2, return_z = TRUE))
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=2, return_z = 1))
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=2, return_z = FALSE))
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=2, return_z = 0))
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=2, return_z = T))
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=2, return_z = F))
  
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=2, return_z = 10))
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=2, return_z = "hi"))
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=2, return_z = "TRUE"))
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=2, return_z = c(TRUE, FALSE)))
}) 


test_that("plot_params() only accepts valid input for contours_only",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=2, contours_only = TRUE))
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=2, t_levels=c(0.9, 0.8, 0.4), contours_only = 1))
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=2, contours_only = FALSE))
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=2, contours_only = 0))
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=2, t_levels=c(0.9, 0.8, 0.4), contours_only = T))
  expect_no_condition(plot_params(x=hockey$x, y=hockey$y, k=2, contours_only = F))
  
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=2, contours_only = 10))
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=2, contours_only = "hi"))
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=2, contours_only = "TRUE"))
  expect_error(plot_params(x=hockey$x, y=hockey$y, k=2, contours_only = c(TRUE, FALSE)))
}) 

test_that("plot_params() only accepts x & y of the same length",{
  x <- runif(10)
  y <- rbinom(10,1,x)

  expect_error(plot_params(x,c(y,y)))
})

test_that("plot_params() options must be a list",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_error(plot_params(x,y, optim_options = c(1,2)))
  expect_error(plot_params(x,y, imgplt_options = c(1,2)))
  expect_error(plot_params(x,y, contour_options = c(1,2)))
})

test_that("plot_params() behavior at gamma = 0",{
  x <- runif(10)
  y <- rbinom(10,1,x)
  
  expect_no_condition(plot_params(x,y, k=3, glim=c(-1,1)))
})

#############################################
#  lineplot() Tests                         #
#############################################

test_that("lineplot() only accepts x & y of the same length",{
  x <- runif(10)
  y <- rbinom(10,1,x)

  expect_error(lineplot(x,c(y,y)))
})

test_that("lineplot() for hockey data doesn't break",{
  expect_no_condition(lineplot(x=hockey$x, y=hockey$y))
})

test_that("lineplot() requires either df or x&y",{
  expect_error(lineplot())
})

test_that("lineplot() returns a dataframe in list when return_df=TRUE",{
  ret <- lineplot(x=hockey$x, y=hockey$y, return_df = TRUE)
  expect_true(is.data.frame(ret$df))
  expect_true(is.list(ret))
})

test_that("lineplot() accepts t_levels",{
  set.seed(46)
  x <- runif(10)
  y <- rbinom(10,1,x)
  expect_no_condition(lineplot(x=x, y=y, t_levels=0.9, nloptr_options = list(maxtime=0.02)))
  expect_error(lineplot(x=x, y=y, t_levels=c("hi", 12, Inf)))
  expect_error(lineplot(x=x, y=y, t_levels=c(Inf)))
  expect_error(lineplot(x=x, y=y, t_levels=c(24)))
})

test_that("lineplot() thinning",{
  set.seed(46)
  x <- runif(100)
  y <- rbinom(100,1,x)
  expect_no_condition(lineplot(x=x, y=y, thin_to=10))
  expect_error(lineplot(x=x, y=y, thin_to=Inf))
  expect_error(lineplot(x=x, y=y, thin_to=-1))
  expect_error(lineplot(x=x, y=y, thin_to=1))
  expect_error(lineplot(x=x, y=y, thin_to=c(10,20)))
  
  expect_no_condition(lineplot(x=x, y=y, thin_percent=0.5))
  expect_error(lineplot(x=x, y=y, thin_percent=Inf))
  expect_error(lineplot(x=x, y=y, thin_percent=-1))
  expect_error(lineplot(x=x, y=y, thin_percent=2))
  expect_error(lineplot(x=x, y=y, thin_percent=c(0.2, 0.3)))
  
  expect_no_condition(lineplot(x=x, y=y, thin_by=2))
  expect_error(lineplot(x=x, y=y, thin_by=Inf))
  expect_error(lineplot(x=x, y=y, thin_by=-1))
  expect_error(lineplot(x=x, y=y, thin_by=c(2,3)))
})

test_that("lineplot() thinning with df",{
  set.seed(46)
  x <- runif(100)
  y <- rbinom(100,1,x)
  expect_no_condition(ret <- lineplot(x=x, y=y, return_df=TRUE))
  df <- ret$df
  expect_no_condition(lineplot(df=df, thin_to=50))
  expect_error(lineplot(df=df, thin_to=Inf))
  expect_error(lineplot(df=df, thin_to="hi"))
  expect_error(lineplot(df=df, thin_to=-1))
  expect_error(lineplot(df=df, thin_to=1))
  expect_error(lineplot(df=df, thin_to=c(10,20)))
  
  expect_no_condition(lineplot(df=df, thin_percent=0.5))
  expect_error(lineplot(df=df, thin_percent=Inf))
  expect_error(lineplot(df=df, thin_percent="hi"))
  expect_error(lineplot(df=df, thin_percent=-1))
  expect_error(lineplot(df=df, thin_percent=2))
  expect_error(lineplot(df=df, thin_percent=c(0.2, 0.3)))
  
  expect_no_condition(lineplot(df=df, thin_by=2))
  expect_error(lineplot(df=df, thin_by=Inf))
  expect_error(lineplot(df=df, thin_by="Inf"))
  expect_error(lineplot(df=df, thin_by=-1))
  expect_error(lineplot(df=df, thin_by=c(2,3)))
  
  expect_no_condition(lineplot(df=df))
})


test_that("lineplot() allows optim_options",{
  set.seed(46)
  x <- runif(100)
  y <- rbinom(100,1,x)
  expect_no_error(lineplot(x,y, optim_options = list(par=c(1,1))))
  expect_error(lineplot(x,y, optim_options = c(1,1)))
})

test_that("lineplot() allows nloptr_options",{
  set.seed(46)
  x <- runif(100)
  y <- rbinom(100,1,x)
  expect_error(lineplot(x,y, nloptr_options = c(1,1)))
})

test_that("lineplot() allows ggpoint_options",{
  set.seed(46)
  x <- runif(100)
  y <- rbinom(100,1,x)
  expect_no_error(lineplot(x,y, ggpoint_options = list(color="red")))
  expect_error(lineplot(x,y, ggpoint_options = c(1,1)))
})

test_that("lineplot() allows ggline_options",{
  set.seed(46)
  x <- runif(100)
  y <- rbinom(100,1,x)
  expect_no_error(lineplot(x,y, ggline_options = list(linewidth=10)))
  expect_error(lineplot(x,y, ggline_options = c(1,1)))
})
