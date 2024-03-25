######################################################
#  External Functions (Version 1)                    #
######################################################

# mention different functions being used under the hood (image, image.plot,
# contour) gamma = 0 will be approx???, option for contours only, addtnal
# options, suggestions for setting k and bounds, how to use z ad return_z
# efficiently, citations, add capability for points at B-R params???? Note
# calculations will change with thinning and is generally not recommended

#' Draw image plot of posterior model probability surface.
#'
#' Function to plot the posterior model probability of the given set of
#' probabilities, `x`, after LLO-adjustment via a grid of uniformly spaced set
#' of \eqn{\delta} and \eqn{\gamma} values with optional contours.
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
#' @param thin_to When non-null, the observations in (x,y) are randomly sampled
#'   without replacement to form a set of size `thin_to`.
#' @param thin_percent When non-null, the observations in (x,y) are randomly
#'   sampled without replacement to form a set that is `thin_percent`\% of the
#'   original size of (x,y).
#' @param thin_by When non-null, the observations in (x,y) are thinned by
#'   selecting every `thin_by` observation.
#' @param contours_only Logical.  If `TRUE`, only the contours at the specified
#'   `t_levels` will be plotted with no color for the posterior model
#'   probability across the k\eqn{\times}k grid of \eqn{\delta} and
#'   \eqn{\gamma}.
#' @param main Plot title.
#' @param xlab Label for x-axis.
#' @param ylab Label for x-axis.
#' @param legend.lab Label for legend.
#' @param drawlabels Logical.  If `TRUE`, the contour labels will be plotted
#'   with the corresponding contours passed to `contour()`.
#' @param contour_color Color of the contours and their labels passed to
#'   `contour()`.
#' @param labcex Size of contour labels passed to `contour()`.
#' @param lwd Linewidth of contours passed to `contour()`.
#' @param ... Additional arguments to be passed to `image`, `image.plot`, and
#'   `contour`.
#'
#' @return If `return_z = TRUE`, a list with the following attributes is
#'   returned: \item{\code{z}}{Matrix containing posterior model probabilities
#'   across k\eqn{\times}k grid of uniformly spaced values of \eqn{\delta} and \eqn{\gamma}
#'   in the specified ranges `dlim` and `glim`, respectively.}
#'   \item{\code{dlim}}{Vector with bounds for
#'   \eqn{\delta} used to construct z.}
#'   \item{\code{glim}}{Vector with bounds for \eqn{\gamma} used to construct
#'   z.}
#'   \item{\code{k}}{The number of uniformly spaced \eqn{\delta} and \eqn{\gamma} values used to construct
#'   z}
#' @export
#'
#' @importFrom graphics contour
#' @importFrom fields image.plot
#'
#' @examples
plot_params <- function(x, y, z=NULL, t_levels = c(0.8, 0.9),
                        Pmc = 0.5, event=1,
                        k = 100,
                        dlim = c(0.0001,5),
                        glim = c(0.0001,5),
                        zlim = c(0,1),
                        return_z = TRUE,
                        thin_to=NULL,
                        thin_percent=NULL,
                        thin_by=NULL,
                        contours_only = FALSE,
                        main="Posterior Model Probability of Calibration",
                        xlab = "delta",
                        ylab = "gamma",
                        legend.lab = "",
                        drawlabels = TRUE,
                        contour_color = "white",
                        labcex=0.6,
                        lwd=1,
                        ...){

  ##################
  #  Input Checks  #
  ##################

  # check either x and y or df are specified
  if(is.null(z) & (is.null(x) | is.null(y))) stop("must specify either x and y or z")

  # CHECK z
  # need to check contents, row and col names, match with specified grid size


  # check t_levels are valid calibration probs
  t_levels <- check_input_probs(t_levels, name="t_levels")


  # check upper and lower bounds
  check_input_delta(dlim[1], name="dlim[1]")
  check_input_delta(dlim[2], name="dlim[2]")
  check_input_gamma(glim[1], name="glim[1]")
  check_input_gamma(glim[2], name="glim[2]")


  # check return_z is logical
  if(!is.logical(return_z) & !(return_z %in% c(0,1))){
    stop("argument return_z must be logical")
  }

  # ADD MORE AFTER DECIDE ON PARAMS

  if(is.null(z)) {
    # check x is vector, values in [0,1]
    x <- check_input_probs(x, name="x")

    # check y is vector, values are 0s or 1s
    y <- check_input_outcomes(y, name="y", event=event)

    # check x and y are the same length
    if(length(x) != length(y)) stop("x and y length differ")

    # check Pmc is valid prior model prob
    Pmc <- check_input_probs(Pmc, name="Pmc")

    # check k
    if(!is.numeric(k)) stop("k must be numeric")
    if(k < 2) stop("k must be greater than 1")
    if(is.infinite(k)) stop("k must be finite")
  }


  ###################
  #  Function Code  #
  ###################


  if(is.null(z)) {
    if(!is.null(thin_to)){
      set.seed(0)
      rows <- sample(1:length(x), size=thin_to)
    } else if (!is.null(thin_percent)){
      set.seed(0)
      rows <- sample(1:length(x), size=length(x)*thin_percent)
    } else if (!is.null(thin_by)){
      rows <- seq(1,length(x),thin_by)
    }  else{
      rows <- 1:length(x)
    }

    x <- x[rows]
    y <- y[rows]

    z <- get_zmat(x=x, y=y, Pmc=Pmc, len.out=k, lower=c(dlim[1], glim[1]), upper=c(dlim[2], glim[2]))

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


    fields::image.plot(x = d, y = g, z = z,
                       xlim = dlim, ylim = glim, zlim = zlim,
                       main = main,
                       xlab = xlab,
                       ylab = ylab,
                       legend.lab =  legend.lab,
                       ...)
    # plus contours if specified
    if(!anyNA(t_levels)){
      contour(x = d, y = g, z = z,  add = TRUE, levels = t_levels, col = contour_color,
              drawlabels = drawlabels, lwd=lwd, labcex=labcex, ...)
    }


    # just plot contours
  }else{
    if(!anyNA(t_levels)){
      contour(x = d, y = g, z = z,
              xlim = dlim, ylim = glim, zlim = zlim,
              levels = t_levels, col = contour_color,
              main = main,
              xlab = xlab,
              ylab = ylab,
              drawlabels = drawlabels, lwd=lwd, labcex=labcex, ...)
    } else {
      stop("must provide contour levels when contours_only = TRUE")  # move this check up!
    }
  }
  #print("plot_params2 end")

  if(return_z){
    return(list(z = z,
                dlim = dlim,
                glim = glim,
                k = k))
  }

}



# start with simple version that you pass x and y that automatically returns OG, MLE, and desired B-R levels
# see lineplot_dev in plotting_functions_old
# add option to control # decimal places printed in pmp
# add option to thin after the fact?
# add checks for graphing params

#' Lineplot for LLO-adjusted Probability Predictions
#'
#' @inheritParams plot_params
#' @param df Dataframe returned by previous call to lineplot() specially
#'   formatted for use in this function. Only used for faster plotting when
#'   making minor cosmetic changes to a previous call.
#' @param return_df Logical.  If `TRUE`, the dataframe used to build this plot
#'   will be returned.
#' @param title Plot title.
#' @param pt_size Size of plotted points passed to `geom_point()`.
#' @param ln_size Linewidth of connecting lines passed to `geom_line()`.
#' @param pt_alpha Transparency of plotted point passed to `geom_point()`.
#' @param ln_alpha Transparency of connecting lines passed to `geom_line()`.
#' @param font_base Base font size for `ggplot()`.
#' @param ylim Vector with bounds for y-axis, must be in \[0,1\].
#' @param breaks Locations along y-axis at which to draw horizontal guidelines,
#'   passed to `scale_y_continous()`.
#'
#' @return If `return_df = TRUE`, a list with the folling attributes is
#'   returned: \item{\code{plot}}{A `ggplot` object showing how the predicted
#'   probabilities under MLE recalibration and specified levels of
#'   boldness-recalibration.}
#'   \item{\code{df}}{Dataframe used to create `plot`, specially
#'   formatted for use in `lineplot()`.}
#'   Otherwise just the `ggplot` object of the plot is returned.
#' @export
#'
#' @import ggplot2
#'
#' @examples
lineplot <- function(x, y, t_levels=NULL, df=NULL,
                     Pmc = 0.5, event=1, return_df=TRUE,
                     title="Line Plot", ylab="Probability",
                     xlab = "Posterior Model Probability",
                     pt_size = 1.5, ln_size = 0.5,
                     pt_alpha = 0.35, ln_alpha = 0.25, font_base = 10,
                     ylim=c(0,1), breaks=seq(0,1,by=0.2), thin_to=NULL,
                     thin_percent=NULL, thin_by=NULL){

  ##################
  #  Input Checks  #
  ##################

  # check either x and y or df are specified
  if(is.null(df) & (is.null(x) | is.null(y))) stop("must specify either x and y or df")

  if(is.null(df)){

    # check x is vector, values in [0,1]
    x <- check_input_probs(x, name="x")

    # check y is vector, values are 0s or 1s
    y <- check_input_outcomes(y, name="y", event=event)

    # check x and y are the same length
    if(length(x) != length(y)) stop("x and y length differ")

    # check t is valid calibration prob
    if(!is.null(t_levels)){  t_levels <- check_input_probs(t_levels, name="t_levels")}

    # check Pmc is valid prior model prob
    Pmc <- check_input_probs(Pmc, name="Pmc")

    # CHECK GRAPHING PARAMS


    ###################
    #  Function Code  #
    ###################

    rows <- 1:length(x)

    if(!is.null(thin_to)){
      set.seed(0)
      rows <- sample(1:length(x), size=thin_to)
    } else if (!is.null(thin_percent)){
      set.seed(0)
      rows <- sample(1:length(x), size=length(x)*thin_percent)
    } else if (!is.null(thin_by)){
      rows <- seq(1,length(x),thin_by)
    }  else{
      rows <- 1:length(x)
    }

    x_plot <- x[rows]
    y_plot <- y[rows]

    nplot <- length(x_plot)

    # create empty DF
    df <- data.frame(matrix(nrow=nplot, ncol=5))
    colnames(df) <- c("probs", "outcome", "post", "pairing", "label")

    pairs <- 1:length(x_plot)

    # Original Set
    df$probs <- x_plot
    df$outcome <- y_plot
    # use full set to get MLEs & posterior model prob, but only plot thinned set
    bt <- bayes_ms(x, y)
    df$post <- bt$posterior_model_prob
    df$pairing <- pairs
    df$label <- paste0("Original \n(",  round(bt$posterior_model_prob,5), ")")


    temp  <- data.frame(matrix(nrow=nplot, ncol=5))
    colnames(temp) <- c("probs", "outcome", "post", "pairing", "label")

    # MLE recalibrate
    temp$probs <- LLO(x_plot, bt$MLEs[1], bt$MLEs[2])
    temp$outcome <- y_plot
    bt_mle <- bayes_ms(LLO(x, bt$MLEs[1], bt$MLEs[2]), y)
    temp$post <- round(bt_mle$posterior_model_prob,5)
    temp$pairing <- pairs
    temp$label <- paste0("MLE Recal. \n(",  round(bt_mle$posterior_model_prob, 5), ")")
    df <- rbind(df, temp)

    # Boldness-recalibrate at given levels
    # loop over t values
    for(i in 1:length(t_levels)){
      br <- brcal(x, y, t_levels[i], Pmc=Pmc, x0 = c(bt$MLEs[1], bt$MLEs[2]),
                  start_at_MLEs=FALSE, print_level=0)
      temp$probs <- LLO(x=x_plot, delta=br$BR_params[1], gamma=br$BR_params[2])
      temp$outcome <- y_plot
      bt_br <- bayes_ms(LLO(x, br$BR_params[1], br$BR_params[2]), y)
      temp$post <- bt_br$posterior_model_prob
      temp$pairing <- pairs
      temp$label <- paste0(round(t_levels[i]*100,0), "% B-R\n(",  round(bt_br$posterior_model_prob, 5), ")")
      df <- rbind(df, temp)
    }
  }

  df$outcome <- factor(df$outcome)
  ulabs <- unique(df$label)
  df$label <- factor(df$label, levels=c(unique(df$label)))

  lines <- ggplot2::ggplot(data = df, mapping = aes_string(x = "label", y = "probs")) +
    geom_point(aes_string(color = "outcome"),
               alpha = pt_alpha, size = pt_size,
               show.legend = FALSE) +
    geom_line(aes_string(group = "pairing", color = "outcome"),
              size = ln_size, alpha = ln_alpha,
              show.legend = FALSE) +
    labs(x = xlab,
         y = ylab) +
    ggtitle(title) +
    theme_bw(base_size = font_base) +
    scale_y_continuous(breaks = breaks,
                       limits = ylim,
                       expand = c(0, 0))+
    scale_color_manual(values = c("red", "blue")) +
    scale_x_discrete(expand = c(0, 0.075)) +
    theme(axis.text.x = element_text(hjust=0.75))

  if(return_df){
    return(list(plot=lines,
                df=df))
  }

  return(lines)
}



######################################################
#  Internal Functions                                #
######################################################

# Function to get matrix of posterior model probabilities across delta/gamma grid
get_zmat <- function(x, y, Pmc=0.5, len.out = 100, lower = c(0.0001,-2), upper = c(5,2)){
  # print("get_zmat start")


  # Set up grid of Delta (d) and Gamma (g)
  d <- seq(lower[1], upper[1], length.out = len.out)
  g <- seq(lower[2], upper[2], length.out = len.out)
  grd <- expand.grid(d,g)

  #print(d)
  #print(g)

  # starting xs
  x0 <- x
  n <- length(x0)

  # MLE recalibrate
  xM <- mle_recal_internal(x0, y, optim_details = FALSE, probs_only = TRUE)

  # paramsM <- bayes_ms_internal(x0, y, optim_details=FALSE)$MLEs
  # xM <- LLO(x=x0, delta=paramsM[1], gamma=paramsM[2])

  # Set up storage
  grd.loglik <- c()
  optim.loglik <- c()
  BIC_1 <- c()
  grd.BIC_2 <- c()
  optim.BIC_2 <- c()

  # Loop over grid of delta/gamma vals
  for(i in 1:nrow(grd)){

    if(grd[i,2] == 0){  # WHAT IS THIS???
      temp <- bayes_ms_internal(LLO_internal(x=x0, delta = grd[i,1], gamma = grd[i,2]), y, Pmc=Pmc)
      BIC_1[i] <- temp$BIC_Mc
      grd.BIC_2[i] <- temp$BIC_Mu
    }else{
      # LLO-adjust based on current grid params
      xg <- LLO_internal(x=x0, delta = grd[i,1], gamma = grd[i,2])

      # grab two unique xs
      xu <- unique(xM)[1:2]

      # find their indices (make sure only grab one index for each)
      uniq_inds <- c(which(xM == xu[1])[1], which(xM == xu[2])[1])

      # Use point slope formula
      start <- logit(xg[uniq_inds])
      goal <- logit(xM[uniq_inds])
      b <- (goal[2] - goal[1]) / (start[2] - start[1])
      a <- goal[2] - b*start[2]
      # print(paste0("a=", a, ", b=", b, ", exp(a)=", exp(a)))
      # print(paste0("delta=", grd[i,1], ", gamma=", grd[i,2]))

      # Get BIC for calibrated model
      BIC_1[i] <- (-2)*llo_lik(params=c(1,1), x=xg, y=y, log=TRUE)

      # Get BIC for uncalibrated model
      grd.loglik[i] <- llo_lik(params=c(exp(a),b), x=xg, y=y, log=TRUE)
      grd.BIC_2[i] <- 2*log(n) - 2*grd.loglik[i]
    }
  }

  # Get Bayes Factor and posterior model prob of calibration
  grd.BF <- exp(-(1/2) * (grd.BIC_2 - BIC_1)) #bayes_factor(BIC1 = grd.BIC_2, BIC2 = BIC_1)
  Pmu <- 1 - Pmc
  posts <- 1/(1+(grd.BF*(Pmu/Pmc))) #post_mod_prob(grd.BF)

  # Reshape vector of posterior model probs into matrix for plotting
  z_mat <- matrix(posts, nrow = length(d), ncol = length(g))
  colnames(z_mat) <- g
  rownames(z_mat) <- d

  # print("get_zmat end")

  return(z_mat)
}

