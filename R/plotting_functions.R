######################################################
#  External Functions                                #
######################################################

#' Draw image plot of posterior model probability surface.
#'
#' Function to visualize the posterior model probability of the given set of
#' probabilities, `x`, after LLO-adjustment via a grid of uniformly spaced set
#' of \eqn{\delta} and \eqn{\gamma} values with optional contours.
#'
#' This function leverages the \link[fields]{image.plot} function from the
#' \link[fields]{fields} package and the \link[graphics]{contour} function from
#' the \link[graphics]{graphics} package.
#'
#' The goal of this function is to visualize how the posterior model probability
#' changes under different recalibration parameters, as this is used in
#' boldness-recalibration.  To do so, a `k` by `k` grid of uniformly spaced
#' potential values for \eqn{\delta} and \eqn{\gamma} are constructed.  Then `x`
#' is LLO-adjusted under each pair of \eqn{\delta} and \eqn{\gamma} values. The
#' posterior model probability of each LLO-adjusted set is calculated and this
#' is the quantity we use to color each grid cell in the image plot to visualize
#' change in calibration.  See below for more details on setting the grid.
#'
#' By default, only the posterior model probability surface is plotted. Argument
#' `t_levels` can be used to optionally add contours at specified levels of the
#' posterior model probability of calibration.   The goal of this is to help
#' visualize different values of \eqn{t} at which they may want to
#' boldness-recalibrate. To only draw the contours without the colored posterior
#' model probability surface, users can set `contours_only=TRUE`.
#'
#' @section Setting grid for \eqn{\delta} and \eqn{\gamma}:
#'
#'   Arguments `dlim` and `glim` are used to set the bounds of the \eqn{\delta},
#'   \eqn{\gamma} grid and the size is dictated by argument `k`. Some care is
#'   required for the selection of these arguments. The goal is to determine
#'   what range of \eqn{\delta} and \eqn{\gamma} encompasses the region of
#'   non-zero posterior probabilities of calibration.  However, it is not
#'   feasible to check the entire parameter space (as it is unbounded) and even
#'   at smaller regions it can be difficult to detect the region in which
#'   non-zero posterior probabilities are produced without as very dense grid
#'   (large `k`), as the region is often quite small relative to the entire
#'   parameter space. This is problematic, as computation time increases as `k`
#'   grows.
#'
#'   We suggest the following scheme setting `k`, `dlim`, and `glim`. First, fix
#'   `k` at some small number, less than 20 for sake of computation time. Then,
#'   center a grid with small range around the MLEs for \eqn{\delta} and
#'   \eqn{\gamma} for the given `x` and `y`. Increase the size of `k` until your
#'   grid detects approximated the probability of calibration at the MLEs that
#'   you expect. Then, expand your grid until it the region with high
#'   probability of calibration is covered or contract your grid to "zoom in" on
#'   the region. Then, increase `k` to create a fine grid of values.
#'
#'   Additionally, we caution users from including \eqn{\gamma = 0} in the grid.
#'   This setting recalibrates all values in `x` to a single value which is not
#'   desirable in practice.  Unless the single value is near the base rate, the
#'   set will be poorly calibrated and minimally bold, which does not align with
#'   the goal of boldness-recalibration.
#'
#' @section Reusing matrix `z` via `return_z`:
#'
#'   The time bottleneck for this function occurs when calculating the posterior
#'   model probabilities across the grid of parameter values. Thus it can be
#'   useful to save the resulting matrix of values to be re-used to save time
#'   when making minor cosmetic changes to your plot. If these adjustments do
#'   not change the grid bounds or density, users can set `return_z=TRUE` to
#'   return the underlying matrix of posterior mode probabilities for plotting.
#'   Then, instead of specifying `x` and `y` users can just pass the returned
#'   matrix as `z`.  Note this assumes you are NOT making any changes to `k`,
#'   `dlim`, or `glim`.  Also, it is not recommended that you construct your
#'   own matrix to pass via `z` as this function relies on the structure as
#'   returned by a previous call of `plot_params()`.
#'
#' @section Thinning:
#'
#'   Another approach to speed up the calculations of this function is to thin
#'   the data used. However, this is generally not recommended unless the sample
#'   size is very large as the calculations of the posterior model probability
#'   may change drastically under small sample sizes.  This can lead to
#'   misleading results. Under large sample sizes where thinning is used, note
#'   this is only an approximate visual of the posterior model probability. 
#'
#' @section Grid cells that show up white / inaccuracies warning message:
#'
#'   In some cases, grid cells in the plot may show up as white instead of one
#'   of the colors from red to blue shown on the legend.  A white grid cell
#'   indicates that there is no calculated posterior model probability at that
#'   cell. There are two common reasons for this: (1) that grid cell location is
#'   not covered by the `z` matrix used (i.e. you've adjusted the bounds without
#'   recalculating z) or (2) the values of the parameters at these locations
#'   cause the values in `x` to be LLO-adjusted such that they virtually equal 0
#'   or 1.  This invokes the use of `epsilon` to push them away from these
#'   boundaries for stability.   This typically happens when |gamma| is very large.
#'   However, in these extreme cases this can cause inaccuracies in this plot. 
#'   For this reason, we either throw the warning message: "Probs too close to 0 
#'   or 1 under very large |gamma|" and allow the cell to be plotted as white 
#'   to notify the user and avoid plotting artifacts.  
#'   
#'   Additionally, when gamma is very close to 0, we cannot directly calculate 
#'   the MLEs for the grid shifted prediction and thus must use `optim()`
#'   to approximate them.  In this case, we throw a `warning` to notify users there
#'   may be inaccuracies. 
#'
#' @inheritParams bayes_ms
#' @param z Matrix returned by previous call to `plot_params()` containing
#'   posterior model probabilities across k\eqn{\times}k grid of \eqn{\delta}
#'   and \eqn{\gamma}. Assumes `z` was constructed using the same `k`, `dlim`,
#'   and `glim` as the current function call.
#' @param t_levels Vector of desired level(s) of calibration at which to plot
#'   contours.
#' @param k The number of uniformly spaced \eqn{\delta} and \eqn{\gamma} values
#'   used to construct the k\eqn{\times}k grid.
#' @param dlim Vector with bounds for \eqn{\delta}, must be finite.
#' @param glim Vector with bounds for \eqn{\gamma}, must be finite.
#' @param zlim Vector with bounds for posterior probability of calibration, must
#'   be in \[0,1\].
#' @param return_z Logical.  If `TRUE`, the matrix of posterior model
#'   probabilities across the specified k\eqn{\times}k grid of \eqn{\delta} and
#'   \eqn{\gamma} will be returned.
#' @param contours_only Logical.  If `TRUE`, only the contours at the specified
#'   `t_levels` will be plotted with no color for the posterior model
#'   probability across the k\eqn{\times}k grid of \eqn{\delta} and
#'   \eqn{\gamma}.
#' @param main Plot title.
#' @param xlab Label for x-axis.
#' @param ylab Label for x-axis.
#' @param optim_options List of additional arguments to be passed to \link[stats]{optim}().
#' @param imgplt_options List of additional arguments to be passed to \link[fields]{image.plot}().
#' @param contour_options List of additional arguments to be passed to \link[graphics]{contour}().
#'
#' @return If `return_z = TRUE`, a list with the following attributes is
#'   returned: \item{\code{z}}{Matrix containing posterior model probabilities
#'   across k\eqn{\times}k grid of uniformly spaced values of \eqn{\delta} and \eqn{\gamma}
#'   in the specified ranges `dlim` and `glim`, respectively.}
#'   \item{\code{dlim}}{Vector with bounds for
#'   \eqn{\delta} used to construct z.}
#'   \item{\code{glim}}{Vector with bounds for \eqn{\gamma} used to construct
#'   z.}
#'   \item{\code{k}}{The number of uniformly spaced \eqn{\delta} and \eqn{\gamma}
#'   values used to construct z}
#' @export
#'
#' @importFrom graphics contour
#' @importFrom fields image.plot
#'
#' @references Guthrie, A. P., and Franck, C. T. (2024) Boldness-Recalibration
#'   for Binary Event Predictions, \emph{The American Statistician} 1-17.
#'
#'   Nychka, D., Furrer, R., Paige, J., Sain, S. (2021). fields: Tools for
#'   spatial data. R package version 15.2,
#'   <https://github.com/dnychka/fieldsRPackage>.
#'
#' @examples
#' 
#' # Simulate 50 predicted probabilities
#' set.seed(49)
#' x <- runif(50)
#' # Simulated 50 binary event outcomes using x
#' y <- rbinom(50, 1, x)  # By construction, x is well calibrated.
#'
#' #' # Set grid density k=20
#' plot_params(x, y, k=20)
#'
#' # Adjust bounds on delta and gamma
#' plot_params(x, y, k=20, dlim=c(0.001, 3), glim=c(0.01,2))
#'
#' # Increase grid density via k & save z matrix for faster plotting
#' zmat_list <- plot_params(x, y, k=100, dlim=c(0.001, 3), glim=c(0.01,2), return_z=TRUE)
#'
#' # Reuse z matrix
#' plot_params(z=zmat_list$z, k=100, dlim=c(0.001, 3), glim=c(0.01,2))
#'
#' # Add contours at t=0.95, 0.9, and 0.8
#' plot_params(z=zmat_list$z, k=100, dlim=c(0.001, 3), glim=c(0.01,2), t_levels=c(0.95, 0.9, 0.8))
#'
#' # Add points for 95% boldness-recalibration parameters
#' br95 <- brcal(x, y, t=0.95, print_level=0)
#' plot_params(z=zmat_list$z, k=100, dlim=c(0.001, 3), glim=c(0.01,2), t_levels=c(0.95, 0.9, 0.8))
#' points(br95$BR_params[1], br95$BR_params[2], pch=19, col="white")
#'
#' # Change color and size of contours
#' plot_params(z=zmat_list$z, k=100, dlim=c(0.001, 3), glim=c(0.01,2), t_levels = c(0.99, 0.1), 
#' contour_options=list(col="orchid", lwd=2))
#' 
#' # Plot contours only
#' plot_params(z=zmat_list$z, k=100, dlim=c(0.001, 3), glim=c(0.01,2), t_levels=c(0.95, 0.9, 0.8),
#' contours_only=TRUE)
#'
#' # Pass arguments to image.plot()
#' plot_params(z=zmat_list$z, k=100, dlim=c(0.001, 3), glim=c(0.01,2),
#'             imgplt_options=list(horizontal = TRUE, nlevel=10, 
#'             legend.lab="Posterior Model Prob"))
#'
#' # See vignette for more examples
#' 
plot_params <- function(x=NULL, y=NULL, z=NULL, t_levels = NULL,
                        Pmc = 0.5, event=1,
                        k = 100,
                        dlim = c(0.0001,5),
                        glim = c(0.0001,5),
                        zlim = c(0,1),
                        return_z = FALSE,
                        epsilon=.Machine$double.eps,
                        contours_only = FALSE,
                        main="Posterior Model Probability of Calibration",
                        xlab = "delta",
                        ylab = "gamma",
                        optim_options=NULL,
                        imgplt_options=list(legend.lab = ""),
                        contour_options=list(drawlabels=TRUE, labcex=0.6, lwd=1, 
                                             col=ifelse(contours_only, "black", "white"))){
  
  ##################
  #  Input Checks  #
  ##################
  
  # check either x and y or z are specified
  if(is.null(z) & (is.null(x) | is.null(y))) stop("must specify either x and y or z")
  
  # check t_levels are valid calibration probs
  if(!is.null(t_levels)) t_levels <- check_input_probs(t_levels, name="t_levels")
  
  # check upper and lower bounds
  check_input_delta(dlim[1], name="dlim[1]")
  check_input_delta(dlim[2], name="dlim[2]")
  check_input_gamma(glim[1], name="glim[1]")
  check_input_gamma(glim[2], name="glim[2]")
  check_value01(zlim[1], name="zlim[1]")
  check_value01(zlim[2], name="zlim[2]")
  
  # check return_z is logical
  if(!is.logical(return_z) & !(return_z %in% c(0,1))){
    stop("argument return_z must be logical")
  }
  
  # check contours_only is logical
  if(!is.logical(contours_only) & !(contours_only %in% c(0,1))){
    stop("argument contours_only must be logical")
  }
  
  if(is.null(z)) {
    # check x is vector, values in [0,1]
    x <- check_input_probs(x, name="x")
    
    # check y is vector, values are 0s or 1s, or event is specified properly
    y <- check_input_outcomes(y, name="y", event=event)
    
    # check x and y are the same length
    if(length(x) != length(y)) stop("x and y length differ")
    
    # check Pmc is valid prior model prob
    Pmc <- check_value01(Pmc, name="Pmc")
    
    # check epsilon
    epsilon <- check_value01(epsilon, name="epsilon")    
    
    # check k
    if(!is.numeric(k)) stop("k must be numeric")
    if(k < 2) stop("k must be greater than 1")
    if(is.infinite(k)) stop("k must be finite")
    
  }
  
  # check that additional options are in the form of a list
  if(!is.null(optim_options) & !is.list(optim_options)) stop("optim_options must be a list")
  if(!is.null(imgplt_options) & !is.list(imgplt_options)) stop("imgplt_options must be a list")
  if(!is.null(contour_options) & !is.list(contour_options)) stop("contour_options must be a list")
  
  # check that contour levels are specified if contours_only=TRUE
  if(is.null(t_levels) & contours_only)  stop("must provide contour levels when contours_only = TRUE") 
  
  ###################
  #  Function Code  #
  ###################
  
  # Calculate z matrix
  if(is.null(z)) {
    rows <- 1:length(x)
    x <- x[rows]
    y <- y[rows]
    z <- get_zmat(x=x, y=y, Pmc=Pmc, len.out=k, lower=c(dlim[1], glim[1]),
                  upper=c(dlim[2], glim[2]), epsilon=epsilon, optim_options=optim_options)
  }
  
  # get max z value
  max_z <- max(z[!is.na(z)])
  
  # set up delta and gamma vectors
  g <- as.numeric(colnames(z))
  d <- as.numeric(rownames(z))
  
  if(anyNA(dlim)){
    dlim <- c(min(d), max(d))
  }
  if(anyNA(glim)){
    glim <- c(min(g), max(g))
  }
  
  
  if(!contours_only){   
    # Plot color surface
    do.call(fields::image.plot, c(list(x = d, y = g, z = z,
                                       xlim = dlim, ylim = glim, zlim = zlim,
                                       main = main,
                                       xlab = xlab,
                                       ylab = ylab),
                                  imgplt_options))
    # plus contours if specified
    if(!is.null(t_levels)){ 
      do.call(contour, c(list(x = d, y = g, z = z,  
                              add = TRUE, levels = t_levels), 
                         contour_options))
    }
  } else { 
    # just plot contours
    do.call(contour, c(list(x = d, y = g, z = z,
                            xlim = dlim, ylim = glim, zlim = zlim,
                            levels = t_levels, 
                            main = main,
                            xlab = xlab,
                            ylab = ylab),
                       contour_options))
    
  } 
  
  # Return z if specified
  if(return_z){
    return(list(z = z,
                dlim = dlim,
                glim = glim,
                k = k))
  }
}


# see lineplot_dev in plotting_functions_old
# add option to control # decimal places printed in pmp
# add option for thin_unique i.e. only plot unique observations. but need to consider that outcomes could be different

#' Lineplot for LLO-adjusted Probability Predictions
#'
#' Function to visualize how predicted probabilities change under
#' MLE-recalibration and boldness-recalibration.
#'
#' This function leverages `ggplot()` and related functions from the `ggplot2`
#' package (REF).
#'
#' The goal of this function is to visualize how predicted probabilities change
#' under different recalibration parameters. By default this function only shows
#' how the original probabilities change after MLE recalibration.  Argument
#' `t_levels` can be used to specify a vector of levels of
#' boldness-recalibration to visualize in addition to MLE recalibration.
#'
#' While the x-axis shows the posterior model probabilities of each set of
#' probabilities, note the posterior model probabilities are not in ascending or
#' descending order.  Instead, they simply follow the ordering of how one might
#' use the `BRcal` package: first looking at the original predictions, then
#' maximizing calibration, then examining how far they can spread out
#' predictions while maintaining calibration with boldness-recalibration.
#'
#' @section Reusing underlying dataframe via `return_df`:
#'
#'   While this function does not typically come with a large burden on time
#'   under moderate sample sizes, there is still a call to `optim()` under the
#'   hood for MLE recalibration and a call to `nloptr()` for each level of
#'   boldness-recalibration that could cause a bottleneck on time.  With this in
#'   mind, users can specify `return_df=TRUE` to return the underlying dataframe
#'   used to build the resulting lineplot.  Then, users can pass this dataframe
#'   to `df` in subsequent calls of `lineplot` to circumvent these calls to
#'   `optim` and `nloptr` and make cosmetic changes to the plot.
#'
#'   When `return_df=TRUE`, both the plot and the dataframe are returned in a
#'   list. The dataframe contains 6 columns:
#'   * `probs`: the values of each predicted probability under each set
#'   * `outcome`: the corresponding outcome for each predicted probability
#'   * `post`: the posterior model probability of the set as a whole
#'   * `id`: the id of each individual probability used for mapping observations between sets
#'   * `set`: the set with which the probability belongs to
#'   * `label`: the label used for the x-axis in the lineplot
#'
#'   Essentially, each set of probabilities (original, MLE-, and each level of
#'   boldness-recalibration) and outcomes are "stacked" on top of each other.
#'   The `id` tells the plotting function how to connect (with line) the same
#'   observation as is changes from the original set to MLE- or
#'   boldness-recalibration.
#'
#' @section Thinning:
#'
#'   Another strategy to save time when plotting is to thin the amount of data
#'   plotted.  When sample sizes are large, the plot can become overcrowded and
#'   slow to plot.  We provide three options for thinning: `thin_to`,
#'   `thin_prop`, and `thin_by`.  By default, all three of these settings are
#'   set to `NULL`, meaning no thinning is performed.  Users can only specify
#'   one thinning strategy at a time. Care should be taken in selecting a
#'   thinning approach based on the nature of your data and problem.  Note that
#'   MLE recalibration and boldness-recalibration will be done using the full
#'   set. 
#'   
#'   Also note that if a thinning strategy is used with `return_df=TRUE`, the 
#'   returned data frame will **only contain the reduced set** (i.e. the data 
#'   *after* thinning).   
#'
#' @section Passing additional arguments to `geom_point()` and `geom_line()`:
#'
#'   To make cosmetic changes to the points and lines plotted, users can pass a
#'   list of any desired arguments of `geom_point()` and `geom_line()` to
#'   `ggpoint_options` and `ggline_options`, respectively.  These will overwrite
#'   everything passed to `geom_point()` or `geom_line()` except any aesthetic
#'   arguments in `aes()`.
#'
#'
#' @inheritParams plot_params
#' @param plot_original Logical.  If `TRUE`, the original probabilities passed 
#'   in `x` are plotted.
#' @param plot_MLE Logical.  If `TRUE`, the MLE-recalibrated probabilities are 
#'   plotted.
#' @param df Dataframe returned by previous call to lineplot() specially
#'   formatted for use in this function. Only used for faster plotting when
#'   making minor cosmetic changes to a previous call.
#' @param return_df Logical.  If `TRUE`, the dataframe used to build this plot
#'   will be returned.
#' @param title Plot title.
#' @param ylim Vector with bounds for y-axis, must be in \[0,1\].
#' @param breaks Locations along y-axis at which to draw horizontal guidelines,
#'   passed to `scale_y_continous()`.
#' @param thin_to When non-null, the observations in (x,y) are randomly sampled
#'   without replacement to form a set of size `thin_to`.
#' @param thin_prop When non-null, the observations in (x,y) are randomly
#'   sampled without replacement to form a set that is `thin_prop` * 100% of
#'   the original size of (x,y).
#' @param thin_percent This argument is deprecated, use `thin_prop` instead.
#' @param thin_by When non-null, the observations in (x,y) are thinned by
#'   selecting every `thin_by` observation.
#' @param seed Seed for random thinning.  Set to NULL for no seed.
#' @param nloptr_options List with options to be passed to `nloptr()`.
#' @param ggpoint_options List with options to be passed to `geom_point()`.
#' @param ggline_options List with options to be passed to `geom_line()`.
#'
#' @return If `return_df = TRUE`, a list with the following attributes is
#'   returned: \item{\code{plot}}{A `ggplot` object showing how the predicted
#'   probabilities under MLE recalibration and specified levels of
#'   boldness-recalibration.}
#'   \item{\code{df}}{Dataframe used to create `plot`, specially
#'   formatted for use in `lineplot()`.}
#'   Otherwise just the `ggplot` object of the plot is returned.
#' @export
#'
#' @import ggplot2
#' @importFrom lifecycle is_present
#' @importFrom lifecycle deprecated
#'
#' @references Guthrie, A. P., and Franck, C. T. (2024) Boldness-Recalibration
#'   for Binary Event Predictions, \emph{The American Statistician} 1-17.
#'
#'   Wickham, H. (2016) ggplot2: Elegant Graphics for Data Analysis.
#'   Springer-Verlag New York.
#'
#' @examples
#' 
#' set.seed(28)
#' # Simulate 100 predicted probabilities
#' x <- runif(100)
#' # Simulated 100 binary event outcomes using x
#' y <- rbinom(100, 1, x)  # By construction, x is well calibrated.
#'
#' # Lineplot show change in probabilities from original to MLE-recalibration to 
#' # specified Levels of Boldness-Recalibration via t_levels
#' # Return a list with dataframe used to construct plot with return_df=TRUE
#' lp1 <- lineplot(x, y, t_levels=c(0.98, 0.95), return_df=TRUE)
#' lp1$plot
#'
#' # Reusing the previous dataframe to save calculation time
#' lineplot(df=lp1$df)
#'
#' # Adjust geom_point cosmetics via ggpoint
#' # Increase point size and change to open circles
#' lineplot(df=lp1$df, ggpoint_options=list(size=3, shape=4))
#'
#' # Adjust geom_line cosmetics via ggline
#' # Increase line size and change transparencys
#' lineplot(df=lp1$df, ggline_options=list(linewidth=2, alpha=0.1))
#'
#' # Thinning down to 75 randomly selected observation
#' lineplot(df=lp1$df, thin_to=75)
#'
#' # Thinning down to 53% of the data
#' lineplot(df=lp1$df, thin_prop=0.53)
#'
#' # Thinning down to every 3rd observation
#' lineplot(df=lp1$df, thin_by=3)
#'
#' # Setting a different seed for thinning
#' lineplot(df=lp1$df, thin_prop=0.53, seed=47)
#'
#' # Setting NO seed for thinning (plot will be different every time)
#' lineplot(df=lp1$df, thin_to=75, seed=NULL)
#' 
lineplot <- function(x=NULL, y=NULL, t_levels=NULL, plot_original=TRUE,
                     plot_MLE=TRUE, df=NULL,
                     Pmc = 0.5, event=1, return_df=FALSE, 
                     epsilon=.Machine$double.eps,
                     title="Line Plot", ylab="Probability",
                     xlab = "Posterior Model Probability",
                     ylim=c(0,1), breaks=seq(0,1,by=0.2), 
                     thin_to=NULL,
                     thin_prop=NULL, 
                     thin_by=NULL,
                     thin_percent=deprecated(),
                     seed=0,
                     optim_options=NULL,
                     nloptr_options=NULL,
                     ggpoint_options=list(alpha=0.35, size=1.5,
                                          show.legend = FALSE),
                     ggline_options=list(alpha=0.25, linewidth=0.5,
                                         show.legend = FALSE)){
  
  ##################
  #  Input Checks  #
  ##################
  
  # check either x and y or df are specified
  if(is.null(df) & (is.null(x) | is.null(y))) stop("must specify either x and y or df")
  
  # check original and MLE toggles are logical
  if(!is.logical(plot_original) & !(plot_original %in% c(0,1))){
    stop("argument plot_original must be logical")
  }
  if(!is.logical(plot_MLE) & !(plot_MLE %in% c(0,1))){
    stop("argument plot_MLE must be logical")
  }
  
  if(is.null(df)){
    
    # check x is vector, values in [0,1]
    x <- check_input_probs(x, name="x")
    
    # check y is vector, values are 0s or 1s
    y <- check_input_outcomes(y, name="y", event=event)
    
    # check x and y are the same length
    if(length(x) != length(y)) stop("x and y length differ")
    
    # check t is valid calibration prob
    if(!is.null(t_levels)) t_levels <- check_input_probs(t_levels, name="t_levels")
    
    # check # sets specified
    if(plot_original + plot_MLE + length(t_levels) < 1) stop("Must specify at least one set to plot via plot_original, plot_MLE, and/or t_levels")
    
    # check Pmc is valid prior model prob
    Pmc <- check_input_probs(Pmc, name="Pmc")
    
    # check epsilon
    epsilon <- check_value01(epsilon, name="epsilon") 
    
    # check that additional options are in the form of a list
    if(!is.null(optim_options) & !is.list(optim_options)) stop("optim_options must be a list")
    if(!is.null(nloptr_options) & !is.list(nloptr_options)) stop("nloptr_options must be a list")
  } else{  
    if(!is.null(t_levels) & any(!(paste0(round(t_levels*100,0), "% B-R") %in% unique(df$set)))) warning("Not all t_levels found in df, will only plot those that are in df.")
  }
  
  
  
  # check return_df is logical
  if(!is.logical(return_df) & !(return_df %in% c(0,1))){
    stop("argument return_df must be logical")
  }
  
  # check thin_to
  if(!is.null(thin_to)){
    if(!is.numeric(thin_to)) stop("thin_to must be numeric")
    if(thin_to < 2) stop("thin_to must be greater than 1")
    if(is.infinite(thin_to)) stop("thin_to must be finite")
  }
  
  # Check if user using thin_percent instead of thin_prop
  if (is_present(thin_percent)) {
    lifecycle::deprecate_warn("1.0.0", "lineplot(thin_percent)", "lineplot(thin_prop)")
    thin_prop <- thin_percent
  }
  
  # check thin_prop
  if(!is.null(thin_prop)) thin_prop <- check_value01(thin_prop, name="thin_prop")
  
  # check thin_by
  if(!is.null(thin_by)){
    if(!is.numeric(thin_by)) stop("thin_by must be numeric")
    if(thin_by < 1) stop("thin_by must be greater than 0")
    if(is.infinite(thin_by)) stop("thin_by must be finite")
  }
  
  # check only one thinning strategy used
  if((!is.null(thin_to) + !is.null(thin_prop) + !is.null(thin_by)) > 1) stop("only specify one thinning strategy")
  
  # Check seed
  if(!is.null(seed)){
    if(!is.numeric(seed)) stop("seed must be numeric")
    if(is.infinite(seed)) stop("seed must be finite")
  }
  
  # Check additional options
  if(!is.null(ggpoint_options) & !is.list(ggpoint_options)) stop("ggpoint_options must be a list")
  if(!is.null(ggline_options) & !is.list(ggline_options)) stop("ggline_options must be a list")
  
  ###################
  #  Function Code  #
  ###################
  
  if(is.null(df)){
    
    rows <- 1:length(x)
    
    # Thinning if specified
    if(!is.null(thin_to)){
      if(!is.null(seed)) set.seed(seed)
      rows <- sample(1:length(x), size=thin_to)
    } else if (!is.null(thin_prop)){
      if(!is.null(seed)) set.seed(seed)
      rows <- sample(1:length(x), size=length(x)*thin_prop)
    } else if (!is.null(thin_by)){
      rows <- seq(1,length(x),thin_by)
    }  else{
      rows <- 1:length(x)
    }
    
    # extract rows to plot 
    x_plot <- x[rows]
    y_plot <- y[rows]
    
    nplot <- length(x_plot)
    ids <- 1:length(x_plot)
    
    # create empty DF
    df <- data.frame(matrix(nrow=nplot, ncol=6))
    colnames(df) <- c("probs", "outcome", "post", "id", "set", "label")
    
    # Original Set
    df$probs <- x_plot
    df$outcome <- y_plot
    # use full set to get MLEs & posterior model prob, but only plot thinned set
    bt <- do.call(bayes_ms_internal, c(list(x, y, Pmc=Pmc, epsilon=epsilon), optim_options))
    df$post <- bt$posterior_model_prob
    df$id <- ids
    df$set <- "Original"
    df$label <- paste0("Original \n(",  round(bt$posterior_model_prob,5), ")")
    
    
    temp  <- data.frame(matrix(nrow=nplot, ncol=6))
    colnames(temp) <- c("probs", "outcome", "post", "id", "set", "label")
    
    # MLE recalibrate
    temp$probs <- LLO(x_plot, bt$MLEs[1], bt$MLEs[2])
    temp$outcome <- y_plot
    bt_mle <- bayes_ms(LLO(x, bt$MLEs[1], bt$MLEs[2]), y, epsilon=epsilon)
    temp$post <- round(bt_mle$posterior_model_prob,5)
    temp$id <- ids
    temp$set <- "MLE Recal"
    temp$label <- paste0("MLE Recal. \n(",  round(bt_mle$posterior_model_prob, 5), ")")
    df <- rbind(df, temp)
    
    
    # Boldness-recalibrate at given levels
    # loop over t values
    if(!is.null(t_levels)){
      for(i in 1:length(t_levels)){
        br <- brcal(x, y, t_levels[i], Pmc=Pmc, x0 = c(bt$MLEs[1], bt$MLEs[2]),
                    start_at_MLEs=FALSE, print_level=0, epsilon=epsilon, 
                    opts=nloptr_options, optim_options=optim_options)
        temp$probs <- LLO(x=x_plot, delta=br$BR_params[1], gamma=br$BR_params[2])
        temp$outcome <- y_plot
        bt_br <- do.call(bayes_ms_internal, c(list(LLO(x, br$BR_params[1], br$BR_params[2]), y, Pmc=Pmc, epsilon=epsilon), optim_options))
        temp$post <- bt_br$posterior_model_prob
        temp$id <- ids
        temp$set <- paste0(round(t_levels[i]*100,0), "% B-R")
        temp$label <- paste0(round(t_levels[i]*100,0), "% B-R\n(",  round(bt_br$posterior_model_prob, 5), ")")
        df <- rbind(df, temp)
      }
    }
  } else{
    
    n <- max(df$id)
    rows <- 1:n
    
    # Thinning if specified
    if(!is.null(thin_to)){
      if(!is.null(seed)) set.seed(seed)
      rows <- sample(1:n, size=thin_to)
    } else if (!is.null(thin_prop)){
      if(!is.null(seed)) set.seed(seed)
      rows <- sample(1:n, size=n*thin_prop)
    } else if (!is.null(thin_by)){
      rows <- seq(1,n,thin_by)
    } else{
      rows <- 1:n
    }
    
    # Extract rows to plot 
    df <- df[df$id %in% rows,]
  }
  
  
  
  # Make sure outcome and label are factors
  df$outcome <- factor(df$outcome)
  df$set <- factor(df$set)
  df$label <- factor(df$label, levels=c(unique(df$label)))
  
  # Extract sets to plot
  cols <- c()
  if(plot_original){cols <- c(cols, "Original")}
  if(plot_MLE){cols <- c(cols, "MLE Recal")}
  if(!is.null(t_levels)){for(i in 1:length(t_levels)){cols <- c(cols,paste0(round(t_levels[i]*100,0), "% B-R"))}}
  df_plot <- df[which(df$set %in% cols),] 
  
  # Create lineplot
  lines <- ggplot2::ggplot(data = df_plot, mapping = aes(x = .data[["label"]], y = .data[["probs"]])) +
    do.call(geom_point, c(list(aes(color = as.factor(.data[["outcome"]]))), ggpoint_options)) +
    do.call(geom_line, c(list(aes(group = .data[["id"]], color = as.factor(.data[["outcome"]]))), 
                         ggline_options)) +
    labs(x = xlab,
         y = ylab) +
    ggtitle(title) +
    theme_bw() +
    scale_y_continuous(breaks = breaks,
                       limits = ylim,
                       expand = c(0, 0))+
    scale_color_manual(values = c("red", "blue")) +
    scale_x_discrete(expand = c(0, 0.075)) +
    theme(axis.text.x = element_text(hjust=0.75))
  
  # Return df if specified
  if(return_df){
    return(list(plot=lines,
                df=df))
  }
  
  # Otherwise just return the plot
  return(lines)
}



######################################################
#  Internal Functions                                #
######################################################

# Function to get matrix of posterior model probabilities across delta/gamma grid
get_zmat <- function(x, y, Pmc=0.5, len.out = 100, lower = c(0.0001,-2), upper = c(5,2), 
                     epsilon=.Machine$double.eps, optim_options=NULL){
  
  # Set up grid of Delta (d) and Gamma (g)
  d <- seq(lower[1], upper[1], length.out = len.out)
  g <- seq(lower[2], upper[2], length.out = len.out)
  grd <- expand.grid(d,g)
  
  # starting xs
  x0 <- x
  n <- length(x0)
  
  # MLE recalibrate
  xM <- do.call(mle_recal_internal, c(list(x=x0, y=y, optim_details = FALSE, probs_only=FALSE), optim_options))
  d_mle <- xM$MLEs[1]
  g_mle <- xM$MLEs[2]
  xM <- xM$probs
  
  # Set up storage
  grd.loglik <- c()
  optim.loglik <- c()
  BIC_1 <- c()
  grd.BIC_2 <- c()
  optim.BIC_2 <- c()
  warn <- NULL
  
  # Loop over grid of delta/gamma vals
  for(i in 1:nrow(grd)){
    
    if(round(grd[i,2],3) == 0){ 
      temp <- do.call(bayes_ms_internal, c(list(LLO_internal(x=x0, delta = grd[i,1], gamma = grd[i,2]), y, Pmc=Pmc), optim_options))
      BIC_1[i] <- temp$BIC_Mc
      grd.BIC_2[i] <- temp$BIC_Mu
      warn <- "There may be innaccuracies near gamma=0"
    }else{
      
      # Grid shifted xs
      xg <- LLO_internal(x=x0, delta = grd[i,1], gamma = grd[i,2])
      
      # MLEs for xg
      dg_mle <- d_mle / (grd[i,1]^(g_mle / grd[i,2]))
      gg_mle <- g_mle / grd[i,2]
      
      # Get BIC for calibrated model
      BIC_1[i] <- (-2)*llo_lik(params=c(1,1), x=xg, y=y, log=TRUE)
      
      # Get BIC for uncalibrated model
      grd.loglik[i] <- llo_lik(params=c(dg_mle, gg_mle), x=xg, y=y, log=TRUE)
      grd.BIC_2[i] <- 2*log(n) - 2*grd.loglik[i]
    }
  }
  
  # Get Bayes Factor and posterior model prob of calibration
  grd.BF <- exp(-(1/2) * (grd.BIC_2 - BIC_1)) 
  Pmu <- 1 - Pmc
  posts <- 1/(1+(grd.BF*(Pmu/Pmc))) 
  
  # Reshape vector of posterior model probs into matrix for plotting
  z_mat <- matrix(posts, nrow = length(d), ncol = length(g))
  colnames(z_mat) <- g
  rownames(z_mat) <- d
  
  if(!is.null(warn)){warning(warn)}
  if(anyNA(z_mat)){warning("Probs too close to 0 or 1 under very large |gamma|")}
  
  return(z_mat)
}

