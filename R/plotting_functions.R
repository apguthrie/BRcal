######################################################
#  External Functions (Version 1)                    #
######################################################

#' Draw image plot of posterior model probability surface.
#'
#' Function to plot the posterior model probability of the given set of
#' probabilities, `x`, after LLO-adjustment via a grid of uniformly spaced set
#' of \eqn{\delta} and \eqn{\gamma} values with optional contours.
#'
#' mention different functions being used under the hood (image, image.plot,
#' contour) gamma = 0 will be approx???, option for contours only, addtnal
#' options, suggestions for setting k and bounds, how to use z ad return_z
#' efficiently, citations, add capability for points at B-R params????
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
#' @examples
plot_params <- function(x, y, z=NULL, t_levels = c(0.8, 0.9),
                        Pmc = 0.5, event=1,
                        k = 100,
                        # lb = c(0.0001,-2),
                        # ub = c(5,2),
                        dlim = c(0.0001,5),
                        glim = c(0.0001,5),
                        zlim = c(0,1),
                        return_z = TRUE,
                        thin_to=NULL,
                        thin_percent=NULL,
                        thin_by=NULL,
                        contours_only = FALSE,
                        # no_legend=FALSE,
                        main="Posterior Model Probability of Calibration",
                        xlab = "delta",
                        ylab = "gamma",
                        drawlabels = TRUE,
                        contour_color = "white",
                        labcex=0.6,
                        lwd=1,
                        # legend.args = list(las=180),
                        # legend.mar = 9,
                        # legend.lab = "",
                        ...){

  ##################
  #  Input Checks  #
  ##################

  # check x is vector, values in [0,1]
  x <- check_input_probs(x, name="x")

  # check y is vector, values are 0s or 1s
  y <- check_input_outcomes(y, name="y", event=event)

  # check x and y are the same length
  if(length(x) != length(y)) stop("x and y length differ")

  # CHECK z
  # need to check contents, row and col names, match with specified grid size


  # check k
  if(!is.numeric(k)) stop("k must be numeric")
  if(k < 2) stop("k must be greater than 1")
  if(is.infinite(k)) stop("k must be finite")

  # check t_levels are valid calibration probs
  t_levels <- check_input_probs(t_levels, name="t_levels")

  # check Pmc is valid prior model prob
  Pmc <- check_input_probs(Pmc, name="Pmc")

  # check upper and lower bounds
  # lb <- check_input_params(lb, name="lb")
  # ub <- check_input_params(ub, name="ub")
  check_input_delta(dlim[1], name="dlim[1]")
  check_input_delta(dlim[2], name="dlim[2]")
  check_input_gamma(glim[1], name="glim[1]")
  check_input_gamma(glim[2], name="glim[2]")


  # check return_z is logical
  if(!is.logical(return_z) & !(return_z %in% c(0,1))){
    stop("argument return_z must be logical")
  }

  # ADD MORE AFTER DECIDE ON PARAMS

  ###################
  #  Function Code  #
  ###################

  #print("plot_params2 start")
  # lower <- lb
  # upper <- ub

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

    # z <- get_zmat(x=x, y=y, Pmc=Pmc, len.out=len.out, lower=lower, upper=upper)
    z <- get_zmat(x=x, y=y, Pmc=Pmc, len.out=k, lower=c(dlim[1], glim[1]), upper=c(dlim[2], glim[2]))

  }

  # get max z value
  max_z <- max(z[!is.na(z)])

  # set up delta and gamma vectors
  g <- as.numeric(colnames(z))
  d <- as.numeric(rownames(z))

  # if(anyNA(lower)){
  #   lower <- c(min(d), min(g))
  # }
  # if(anyNA(upper)){
  #   upper <- c(max(d), max(g))
  # }

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
                       ...)
    # plus contours if specified
    if(!anyNA(t_levels)){
      contour(x = d, y = g, z = z,  add = TRUE, levels = t_levels, col = contour_color,
              drawlabels = drawlabels, lwd=lwd, ...)
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
              drawlabels = drawlabels, lwd=lwd,  ...)
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


# Function to make contour plot

plot_params_OLD <- function(z, len.out = 100,
                            lower = c(0.0001,-2), upper = c(5,2),
                            cont_levels = c(0.8, 0.9),
                            sub = "",
                            zlim = c(0,1),
                            ttle_extra = "",
                            ttle = "Posterior Model Probability of Calibration",
                            contours_only = FALSE,
                            add = FALSE,
                            contour_color = "white",
                            legend.lab = "",
                            drawlabels = TRUE,
                            xlab = "delta",
                            ylab = "gamma",
                            lwd=1,
                            labcex=0.6,
                            legend.args = list(las=180),
                            legend.mar = 9, no_legend=FALSE,
                            ...){
  #library(fields)
  max_z <- max(z[!is.na(z)])

  if(anyNA(lower)){
    lower <- c(min(d), min(g))
  }
  if(anyNA(upper)){
    upper <- c(max(d), max(g))
  }

  g <- as.numeric(colnames(z))
  d <- as.numeric(rownames(z))

  if(!contours_only){

    if(no_legend){
      image(d, g, z, zlim = zlim, xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2]),
            main = paste0(ttle, ttle_extra),
            xlab = xlab,
            ylab = ylab,
            sub = sub, ...)
      if(!anyNA(cont_levels)){
        contour(d, g, z, add = TRUE, levels = cont_levels, col = contour_color,
                drawlabels = drawlabels, lwd=lwd, labcex=labcex)
      }
    }else{

      fields::image.plot(d, g, z, zlim = zlim, xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2]),
                         main = paste0(ttle, ttle_extra),
                         xlab = xlab,
                         ylab = ylab,
                         sub = sub,
                         legend.mar = legend.mar,
                         legend.lab = legend.lab,
                         legend.args = legend.args, ...)
      if(!anyNA(cont_levels)){
        contour(d, g, z, add = TRUE, levels = cont_levels, col = contour_color,
                drawlabels = drawlabels, lwd=lwd, labcex=labcex)
      }
    }
  }else{
    if(!anyNA(cont_levels)){
      contour(d, g, z, add = add, levels = cont_levels, col = contour_color,
              zlim = zlim, xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2]),
              main = paste0(ttle, ttle_extra),
              xlab = xlab,
              ylab = ylab,
              sub = sub, drawlabels = drawlabels, lwd=lwd, labcex=labcex, ...)
    } else {
      stop("must provide contour levels")
    }
  }
}


lineplot <- function(df, ttle="Line Plot", ylab="Probability", xlab = "Posterior Model Probability",
                     font_size_lp=5,
                     outside_only = FALSE, pt_size = 0.75, ln_size = 0.25,
                     pt_alpha = 0.35, ln_alpha = 0.25, font_base = 10,
                     ylim=c(0,1), breaks=seq(0,1,by=0.2)){
  font_size_lp <- font_base

  #requireNamespace(ggplot2)

  if(!outside_only){

    lines <- ggplot(data = df, mapping = aes(x = post, y = probs)) +
      geom_point(aes(color = outcome), alpha = pt_alpha, size = pt_size,
                 show.legend = FALSE) +
      geom_line(aes(group=pairing, color = outcome), size = ln_size, alpha = ln_alpha,
                show.legend = FALSE) +
      labs(x = xlab,
           y = ylab) +
      ggtitle(ttle) +
      theme_bw(base_size = font_size_lp) +
      scale_y_continuous(breaks = breaks,
                         limits = ylim,
                         expand = c(0, 0))+
      scale_color_manual(values = c("blue", "red")) +
      scale_x_discrete(expand = c(0, 0.075))
    # theme(axis.title.x = element_blank(),
    #       axis.title.y = element_blank())
    # theme(axis.title.x = element_text(size = font_size_lp),
    #       axis.title.y = element_text(size = font_size_lp),
    #       axis.text = element_text(size = font_size_lp),
    #       title = element_text(size = font_size_lp))
  } else {

    lines <- ggplot(data = df, mapping = aes(x = post, y = probs)) +
      geom_point(aes(color = pundit), alpha = pt_alpha, size = pt_size,
                 show.legend = FALSE) +
      geom_line(aes(group=c(pairing,pundit), color = pundit), size = ln_size, alpha = ln_alpha,
                show.legend = FALSE) +
      labs(x = xlab,
           y = ylab) +
      ggtitle(ttle) +
      theme_bw() +
      scale_y_continuous(breaks = breaks,
                         limits = ylim,
                         expand = c(0, 0))+
      scale_color_manual(values = c("black", "red", "blue")) +
      scale_x_discrete(expand = c(0, 0.075)) +
      theme(axis.title.x = element_text(size = font_size_lp),
            axis.title.y = element_text(size = font_size_lp),
            axis.text = element_text(size = font_size_lp),
            title = element_text(size = font_size_lp))
  }
  return(lines)

}

lineplot_dev <- function(x, y, t=NULL, delta=NULL, gamma=NULL, ttle="Line Plot", ylab="Probability",
                         xlab = "Posterior Model Probability",
                         font_size_lp=5,
                         outside_only = FALSE, pt_size = 1.5, ln_size = 0.5,
                         pt_alpha = 0.35, ln_alpha = 0.25, font_base = 10,
                         ylim=c(0,1), breaks=seq(0,1,by=0.2), thin_to=NULL,
                         thin_percent=NULL, thin_by=NULL, pmp_label=FALSE, deciles=FALSE,event=1){
  # check y only has two values
  y <- ifelse(y == event, 1, 0)

  # check validity of x,y inputs

  # numeric
  if(any(class(x) != "numeric")) stop("x must be non-empty numeric vector")
  if(any(!(class(y) %in% c("integer","numeric")))) stop("y must be non-empty numeric vector")

  # same length
  if(length(x) != length(y)) stop("x and y must be the same length")

  # y binary
  if(any(!(y %in% c(0,1)))) stop('y must be numeric vector of 1s (event) and 0s (nonevent)')

  # x between 0,1
  if(any((x < 0 | x > 1))) stop('x must be numeric vector in [0,1]')


  # check that EITHER t or delta and gamma are specified

  if(is.null(delta) != is.null(gamma)) stop('delta and gamma must either both be specified or NULL')

  if(is.null(delta) & is.null(t)) stop('EITHER t OR delta and gamma must be specified, both are NULL')

  if(!is.null(t) & !is.null(delta)) stop('EITHER t OR delta and gamma must be specified, both are specified')


  # When delta & gamma are specified....
  if(!is.null(delta)){
    # check validity of delta, gamma
    # numeric
    if(any(class(delta) != "numeric")) stop("delta must be non-empty numeric vector")
    if(any(class(gamma) != "numeric")) stop("gamma must be non-empty numeric vector")

    # same length
    if(length(delta) != length(gamma)) stop("delta and gamma must be the same length")

    # correct range for delta
    # Gamma is fine without checking because we've already checked it's numeric vector which can only hold real numbers
    if(any(delta <= 0)) stop("delta must be greater than 0")

    # ADD THINNING
    if(!is.null(thin_to) & !is.null(thin_percent)) warning('both thin_to and thin_by specified, thin_by ignored')

    if(length(x) > 5000 & is.null(thin_to) & is.null(thin_percent)) warning('plotting may be slow due to large x, considering using thin_to or thin_by')

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

    if (deciles){
      x_thin <- quantile(x,seq(0,1,.1))
      y_thin <- rep(1, length(x_thin))
      nplot <- length(x_thin)
    } else{
      nplot <- length(rows)
      x_thin <- x[rows]
      y_thin <- y[rows]
    }

    df <- data.frame(matrix(nrow=nplot, ncol=7))
    colnames(df) <- c("probs", "outcome", "post", "pairing", "delta", "gamma", "label")


    # for original
    bt <- bayes_ms(x, y)
    df$probs <- x_thin
    df$outcome <- factor(y_thin, levels = c("1", "0"))
    df$post <- rep(bt$posterior_model_prob, nplot)
    df$pairing <- factor(seq(1, nplot))
    df$delta <- df$gamma  <-  rep(1, nplot)

    if(pmp_label){
      if(round(df$post[1],4)==0){
        post = "0.0000"
      } else {
        post <- as.character(round(df$post, 4))
      }
      df$label <- paste0("Original \n(", post ,")")

    } else {
      df$label <- paste0(as.character(round(df$post, 5)), "\n (Original)")

    }

    #df$label[1:nplot] <- rep(paste0("Original \n (",
    #                                as.character(post), ")"), nplot)  # NEED TO FINISH THIS

    # NEED TO MAKE LAST PART A FACTOR AND DECIDE ON LABELS - later

    # for all sets of passed delta, gamma
    # recalibrate, then pass to bayes_ms
    for(i in 1:length(delta)){
      temp <- data.frame(matrix(nrow=nplot, ncol=7))
      colnames(temp) <- c("probs", "outcome", "post", "pairing", "delta", "gamma", "label")

      llo_probs <-  LLO(x=x, delta[i], gamma[i])


      temp$probs <- LLO(x=x_thin, delta[i], gamma[i])

      bt <- bayes_ms(llo_probs, y)

      temp$outcome <- factor(y_thin, levels = c("1", "0"))
      temp$post <- rep(bt$posterior_model_prob, nplot)
      temp$pairing <- factor(seq(1, nplot))
      temp$delta <- rep(delta[i], nplot)
      temp$gamma <- rep(gamma[i], nplot)

      if(pmp_label){
        # NEED TO GENERALIZE! usee t instead
        if(i==1){
          lab <- "MLE Recalib."
          post <- as.character(round(temp$post, 4))
        }else if(i==2){
          lab <- "95% B-R"
          post <- "0.9500"
        }

        temp$label <- paste0(lab, " \n(", post,")")

      } else {
        temp$label <- paste0(as.character(round(temp$post, 5)), "\n delta=", round(delta[i],3), "\n gamma=", round(gamma[i],3))

      }

      df <- rbind(df, temp)
    }

    df$label <- factor(df$label,
                       levels = unique(df$label))

  }


  # When t is specified...
  # check validity of t
  if(!is.null(t)){
    # numeric
    if(any(class(t) != "numeric")) stop("t must be non-empty numeric vector")

    # correct range for t
    if(any(!(t < 0 | t > 1))) stop('x must be numeric vector in [0,1]')
  }

  # CHECK GRAPHING PARAMS





  font_size_lp <- font_base

  library(ggplot2)

  if(deciles){
    lines <- ggplot(data = df, mapping = aes(x = label, y = probs)) +
      geom_point(color = "black", alpha = pt_alpha, size = pt_size,
                 show.legend = FALSE) +
      geom_line(aes(group=pairing), size = ln_size, alpha = ln_alpha,
                show.legend = FALSE, color="black") +
      labs(x = xlab,
           y = ylab) +
      ggtitle(ttle) +
      theme_bw(base_size = font_size_lp) +
      scale_y_continuous(breaks = breaks,
                         limits = ylim,
                         expand = c(0, 0))+
      #scale_color_manual(values = c("blue", "red")) +
      scale_x_discrete(expand = c(0, 0.075))
  } else if(!outside_only){

    lines <- ggplot(data = df, mapping = aes(x = label, y = probs)) +
      geom_point(aes(color = outcome), alpha = pt_alpha, size = pt_size,
                 show.legend = FALSE) +
      geom_line(aes(group=pairing, color = outcome), size = ln_size, alpha = ln_alpha,
                show.legend = FALSE) +
      labs(x = xlab,
           y = ylab) +
      ggtitle(ttle) +
      theme_bw(base_size = font_size_lp) +
      scale_y_continuous(breaks = breaks,
                         limits = ylim,
                         expand = c(0, 0))+
      scale_color_manual(values = c("blue", "red")) +
      scale_x_discrete(expand = c(0, 0.075)) +
      theme(axis.text.x = element_text(hjust=0.75))
    # theme(axis.title.x = element_blank(),
    #       axis.title.y = element_blank())
    # theme(axis.title.x = element_text(size = font_size_lp),
    #       axis.title.y = element_text(size = font_size_lp),
    #       axis.text = element_text(size = font_size_lp),
    #       title = element_text(size = font_size_lp))
  } else {

    lines <- ggplot(data = df, mapping = aes(x = post, y = probs)) +
      geom_point(aes(color = pundit), alpha = pt_alpha, size = pt_size,
                 show.legend = FALSE) +
      geom_line(aes(group=c(pairing,pundit), color = pundit), size = ln_size, alpha = ln_alpha,
                show.legend = FALSE) +
      labs(x = xlab,
           y = ylab) +
      ggtitle(ttle) +
      theme_bw() +
      scale_y_continuous(breaks = breaks,
                         limits = ylim,
                         expand = c(0, 0))+
      scale_color_manual(values = c("black", "red", "blue")) +
      scale_x_discrete(expand = c(0, 0.075)) +
      theme(axis.title.x = element_text(size = font_size_lp),
            axis.title.y = element_text(size = font_size_lp),
            axis.text = element_text(size = font_size_lp),
            title = element_text(size = font_size_lp))
  }
  #return(list(plot=lines, df=df))
  return(lines)
}

# Function to get matrix of posterior model probabilities across delta/gamma grid
get_zmat_OLD <- function(x, y, len.out = 100, lower = c(0.0001,-2), upper = c(5,2), event=1){
  # print("get_zmat start")
  # check y only has two values
  y <- ifelse(y == event, 1, 0)

  # Set up grid of Delta (d) and Gamma (g)
  d <- seq(lower[1], upper[1], length.out = len.out)
  g <- seq(lower[2], upper[2], length.out = len.out)
  grd <- expand.grid(d,g)

  #print(d)
  #print(g)
  # Loop through grid points to get posterior model probability
  #temp <- c()
  # for(i in 1:nrow(grd)){
  #   x_new <- LLO(x, delta = grd[i,1], gamma = grd[i,2])   # LLO adjust probs FIRST based on grid point
  #   pmp <- bayes_ms(x_new, y)$posterior_model_prob   # Get posterior model prob
  #   temp <- c(temp, pmp)
  # }

  x0 <- x
  paramsM <- LLO_LRT(x0, y, optim_details=FALSE)$est_params
  xM <- LLO(x=x0, delta=paramsM[1], gamma=paramsM[2])

  grd.loglik <- c()
  optim.loglik <- c()
  BIC_1 <- c()
  grd.BIC_2 <- c()
  optim.BIC_2 <- c()

  n <- length(x0)

  for(i in 1:nrow(grd)){

    if(grd[i,2] == 0){
      temp <- bayes_ms(LLO(x=x0, delta = grd[i,1], gamma = grd[i,2]), y)
      BIC_1[i] <- temp$BIC_H0
      grd.BIC_2[i] <- temp$BIC_H1
    }else{
      xg <- LLO(x=x0, delta = grd[i,1], gamma = grd[i,2])
      xu <- unique(xM)[1:2]  # grab two unique xs
      uniq_inds <- c(which(xM == xu[1])[1], which(xM == xu[2])[1]) # find their indices (make sure only grab one index for each)
      start <- logit(xg[uniq_inds])
      goal <- logit(xM[uniq_inds])
      b <- (goal[2] - goal[1]) / (start[2] - start[1])
      a <- goal[2] - b*start[2]
      # print(paste0("a=", a, ", b=", b, ", exp(a)=", exp(a)))
      # print(paste0("delta=", grd[i,1], ", gamma=", grd[i,2]))
      grd.loglik[i] <- llo_lik(params=c(exp(a),b), x=xg, y=y, log=TRUE)
      BIC_1[i] <- (-2)*llo_lik(params=c(1,1), x=xg, y=y, log=TRUE)
      grd.BIC_2[i] <- 2*log(n) - 2*grd.loglik[i]
    }
  }


  grd.BF <- bayes_factor(BIC1 = grd.BIC_2, BIC2 = BIC_1)
  temp <- post_mod_prob(grd.BF)


  # Reshape vector of posterior model probs into matrix for plotting
  z_mat <- matrix(temp, nrow = length(d), ncol = length(g))
  colnames(z_mat) <- g
  rownames(z_mat) <- d

  # print("get_zmat end")

  return(z_mat)
}



plot_params2 <- function(x, y, len.out = 100,
                         lower = c(0.0001,-2), upper = c(5,2),
                         cont_levels = c(0.8, 0.9),
                         sub = "",
                         zlim = c(0,1),
                         ttle_extra = "",
                         ttle = "Posterior Model Probability of Calibration",
                         contours_only = FALSE,
                         add = FALSE,
                         contour_color = "white",
                         legend.lab = "",
                         drawlabels = TRUE,
                         xlab = "delta",
                         ylab = "gamma",
                         lwd=1,
                         labcex=0.6,
                         legend.args = list(las=180),
                         legend.mar = 9, no_legend=FALSE,
                         thin_to=NULL,
                         thin_by=NULL,
                         thin_percent=NULL,
                         event=1,
                         ...){
  #print("plot_params2 start")

  # check y only has two values
  y <- ifelse(y == event, 1, 0)

  #rows <- 1:length(x)
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

  z <- get_zmat(x=x, y=y, len.out=len.out, lower=lower, upper=upper)

  #library(fields)
  max_z <- max(z[!is.na(z)])

  if(anyNA(lower)){
    lower <- c(min(d), min(g))
  }
  if(anyNA(upper)){
    upper <- c(max(d), max(g))
  }

  g <- as.numeric(colnames(z))
  d <- as.numeric(rownames(z))

  if(!contours_only){

    if(no_legend){
      image(d, g, z, zlim = zlim, xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2]),
            main = paste0(ttle, ttle_extra),
            xlab = xlab,
            ylab = ylab,
            sub = sub, ...)
      if(!anyNA(cont_levels)){
        contour(d, g, z, add = TRUE, levels = cont_levels, col = contour_color,
                drawlabels = drawlabels, lwd=lwd, labcex=labcex)
      }
    }else{

      fields::image.plot(d, g, z, zlim = zlim, xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2]),
                         main = paste0(ttle, ttle_extra),
                         xlab = xlab,
                         ylab = ylab,
                         sub = sub,
                         legend.mar = legend.mar,
                         legend.lab = legend.lab,
                         legend.args = legend.args, ...)
      if(!anyNA(cont_levels)){
        contour(d, g, z, add = TRUE, levels = cont_levels, col = contour_color,
                drawlabels = drawlabels, lwd=lwd, labcex=labcex)
      }
    }
  }else{
    if(!anyNA(cont_levels)){
      contour(d, g, z, add = add, levels = cont_levels, col = contour_color,
              zlim = zlim, xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2]),
              main = paste0(ttle, ttle_extra),
              xlab = xlab,
              ylab = ylab,
              sub = sub, drawlabels = drawlabels, lwd=lwd, labcex=labcex, ...)
    } else {
      stop("must provide contour levels")
    }
  }
  #print("plot_params2 end")

}
