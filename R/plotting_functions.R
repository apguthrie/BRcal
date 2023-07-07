
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
