
# Function to make contour plot

plot_params <- function(z, len.out = 100,
                        lower = c(0.0001,-2), upper = c(5,2),
                        cont_levels = c(0.8, 0.9),
                        sub = "",
                        zlim = c(0,1),
                        ttle_extra = "",
                        ttle = "Posterior Model Probability of Calibration",
                        contours_only = FALSE,
                        add = FALSE,
                        contour_color = "white",
                        leg_lab = "",
                        drawlabels = TRUE,
                        xlab = "delta",
                        ylab = "gamma",
                        lwd=1,
                        labcex=0.6,
                        legend.args = list(las=180),
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
    fields::image.plot(d, g, z, zlim = zlim, xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2]),
                       main = paste0(ttle, ttle_extra),
                       xlab = xlab,
                       ylab = ylab,
                       sub = sub,
                       legend.mar = 9,
                       legend.lab = leg_lab,
                       legend.args = legend.args, ...)
    if(!anyNA(cont_levels)){
      contour(d, g, z, add = TRUE, levels = cont_levels, col = contour_color,
              drawlabels = drawlabels, lwd=lwd, labcex=labcex)
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

  requireNamespace(ggplot2)

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
                         thin_percent=NULL, thin_by=NULL){


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
    } else{
      rows <- 1:length(x)
    }
    nplot <- length(rows)
    x <- x[rows]
    y <- y[rows]

    df <- data.frame(matrix(nrow=nplot, ncol=7))
    colnames(df) <- c("probs", "outcome", "post", "pairing", "delta", "gamma", "label")


    # for original
    bt <- bayes_testing(x, y)
    df$probs <- x
    df$outcome <- factor(y, levels = c("1", "0"))
    df$post <- rep(bt$posterior_model_prob, nplot)
    df$pairing <- factor(seq(1, nplot))
    df$delta <- df$gamma  <-  rep(1, nplot)
    df$label <- paste0(as.character(round(df$post, 5)), "\n (Original)")
    #df$label[1:nplot] <- rep(paste0("Original \n (",
    #                                as.character(post), ")"), nplot)  # NEED TO FINISH THIS

    # NEED TO MAKE LAST PART A FACTOR AND DECIDE ON LABELS - later

    # for all sets of passed delta, gamma
    # recalibrate, then pass to bayes_testing
    for(i in 1:length(delta)){
      temp <- data.frame(matrix(nrow=nplot, ncol=7))
      colnames(temp) <- c("probs", "outcome", "post", "pairing", "delta", "gamma", "label")

      temp$probs <- LLO(x, delta[i], gamma[i])

      bt <- bayes_testing(temp$probs, y)

      temp$outcome <- factor(y, levels = c("1", "0"))
      temp$post <- rep(bt$posterior_model_prob, nplot)
      temp$pairing <- factor(seq(1, nplot))
      temp$delta <- rep(delta[i], nplot)
      temp$gamma <- rep(gamma[i], nplot)
      temp$label <- paste0(as.character(round(temp$post, 5)), "\n (delta=", round(delta[i],3), ", gamma=", round(gamma[i],3), ")")
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

  if(!outside_only){

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
  return(list(plot=lines, df=df))

}

# Function to get matrix of posterior model probabilities across delta/gamma grid
get_zmat <- function(x, y, len.out = 100, lower = c(0.0001,-2), upper = c(5,2)){

  # Set up grid of Delta (d) and Gamma (g)
  d <- seq(lower[1], upper[1], length.out = len.out)
  g <- seq(lower[2], upper[2], length.out = len.out)
  grd <- expand.grid(d,g)

  # Loop through grid points to get posterior model probability
  #temp <- c()
  # for(i in 1:nrow(grd)){
  #   x_new <- LLO(x, delta = grd[i,1], gamma = grd[i,2])   # LLO adjust probs FIRST based on grid point
  #   pmp <- bayes_testing(x_new, y)$posterior_model_prob   # Get posterior model prob
  #   temp <- c(temp, pmp)
  # }

  x0 <- x
  paramsM <- LLO_LRT(x0, y, optim_details=FALSE)$est_params
  xM <- LLO(x0, delta=paramsM[1], gamma=paramsM[2])

  grd.loglik <- c()
  optim.loglik <- c()
  BIC_1 <- c()
  grd.BIC_2 <- c()
  optim.BIC_2 <- c()

  n <- length(x0)

  for(i in 1:nrow(grd)){
    xg <- LLO(x0, delta = grd[i,1], gamma = grd[i,2])
    xu <- unique(xM)[1:2]  # grab two unique xs
    uniq_inds <- c(which(xM == xu[1])[1], which(xM == xu[2])[1]) # find their indices (make sure only grab one index for each)
    start <- logit(xg[uniq_inds])
    goal <- logit(xM[uniq_inds])
    b <- (goal[2] - goal[1]) / (start[2] - start[1])
    a <- goal[2] - b*start[2]
    grd.loglik[i] <- llo_lik(params=c(exp(a),b), x=xg, y=y, log=TRUE)
    BIC_1[i] <- (-2)*llo_lik(params=c(1,1), x=xg, y=y, log=TRUE)
    grd.BIC_2[i] <- 2*log(n) - 2*grd.loglik[i]
  }

  grd.BF <- bayes_factor(BIC1 = grd.BIC_2, BIC2 = BIC_1)
  temp <- post_mod_prob(grd.BF)


  # Reshape vector of posterior model probs into matrix for plotting
  z_mat <- matrix(temp, nrow = length(d), ncol = length(g))
  colnames(z_mat) <- g
  rownames(z_mat) <- d
  return(z_mat)
}

