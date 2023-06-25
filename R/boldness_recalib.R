
# Function to get matrix of posterior model probabilities across delta/gamma grid
get_zmat <- function(x,y, len.out = 100, lower = c(0.0001,-2), upper = c(5,2)){

  # Set up grid of Delta (d) and Gamma (g)
  d <- seq(lower[1], upper[1], length.out = len.out)
  g <- seq(lower[2], upper[2], length.out = len.out)
  grd <- expand.grid(d,g)

  # Loop through grid points to get posterior model probability
  temp <- c()
  for(i in 1:nrow(grd)){
    x_new <- LLO(x, delta = grd[i,1], gamma = grd[i,2])   # LLO adjust probs FIRST based on grid point
    pmp <- bayes_testing(x_new, y)$posterior_model_prob   # Get posterior model prob
    temp <- c(temp, pmp)
  }

  # Reshape vector of posterior model probs into matrix for plotting
  z_mat <- matrix(temp, nrow = length(d), ncol = length(g))
  colnames(z_mat) <- g
  rownames(z_mat) <- d
  return(z_mat)
}


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
  library(fields)
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
    image.plot(d, g, z, zlim = zlim, xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2]),
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



bold_recalib <- function(vec, zmat, calib = 0.95, filen){
  # inds1 <- which(zmat >= calib & zmat <= (calib + 0.005), arr.ind = TRUE) # looks at a tight range around 0.9
  inds1 <- which(zmat >= calib, arr.ind = TRUE) # looks at anything greater or equal to 0.9
  adj_res <- data.frame(delta = vector(length = nrow(inds1)),
                        gamma = vector(length = nrow(inds1)),
                        post = vector(length = nrow(inds1)),
                        sd = vector(length = nrow(inds1)))

  for(i in 1:nrow(inds1)){
    adj_res$post[i] <- zmat[inds1[i,1], inds1[i,2]]
    adj_res$delta[i] <- rownames(zmat)[inds1[i,1]]
    adj_res$gamma[i] <- colnames(zmat)[inds1[i,2]]
  }

  adj_res$delta <- as.numeric(adj_res$delta)
  adj_res$gamma <- as.numeric(adj_res$gamma)

  # adj_mat <- matrix(nrow = length(vec), ncol = nrow(adj_res))

  for(i in 1:nrow(adj_res)){
    # adj_mat[,i] <- LLO(vec, delta = adj_res$delta[i], gamma = adj_res$gamma[i])
    temp <- LLO(vec, delta = adj_res$delta[i], gamma = adj_res$gamma[i])
    adj_res$sd[i] <- sd(temp)
  }


  write.csv(adj_res, paste0("simstudy_results/adjres_", filen, "_", calib, ".csv"))

  # adj_df <- data.frame(adj_mat)
  # colnames(adj_df) <- expand.grid(paste0(filen, "_"), as.character(1:ncol(adj_df))) %>%
  #   transmute(label = paste0(Var1, Var2)) %>%
  #   .$label
  # return(adj_df)
}

lineplot <- function(df, ttle="Line Plot", ylab="Probability", xlab = "Posterior Model Probability",
                     font_size_lp=5,
                     outside_only = FALSE, pt_size = 0.75, ln_size = 0.25,
                     pt_alpha = 0.35, ln_alpha = 0.25, font_base = 10,
                     ylim=c(0,1), breaks=seq(0,1,by=0.2)){
  font_size_lp <- font_base

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
