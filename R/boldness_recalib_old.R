


# Function to get matrix of posterior model probabilities across delta/gamma grid
get_zmat_old <- function(x, y, len.out = 100, lower = c(0.0001,-2), upper = c(5,2)){

  # Set up grid of Delta (d) and Gamma (g)
  d <- seq(lower[1], upper[1], length.out = len.out)
  g <- seq(lower[2], upper[2], length.out = len.out)
  grd <- expand.grid(d,g)

  # Loop through grid points to get posterior model probability
  temp <- c()
  for(i in 1:nrow(grd)){
    x_new <- LLO(x, delta = grd[i,1], gamma = grd[i,2])   # LLO adjust probs FIRST based on grid point
    pmp <- bayes_ms(x_new, y)$posterior_model_prob   # Get posterior model prob
    temp <- c(temp, pmp)
  }

  # Reshape vector of posterior model probs into matrix for plotting
  z_mat <- matrix(temp, nrow = length(d), ncol = length(g))
  colnames(z_mat) <- g
  rownames(z_mat) <- d
  return(z_mat)
}

# version where maxxed log likelihood is returned instead of pmp
get_zmat_ml <- function(x,y, len.out = 100, lower = c(0.0001,-2), upper = c(5,2)){

  # Set up grid of Delta (d) and Gamma (g)
  d <- seq(lower[1], upper[1], length.out = len.out)
  g <- seq(lower[2], upper[2], length.out = len.out)
  grd <- expand.grid(d,g)


  # # BIC under alternative
  # temp <- BIC_llo(x = x, y = y, k = k, params = params, lower = lower, upper = upper)
  # BIC2 <- temp$BIC
  #
  # optBayes <- stats::optim(c(0.5, 0.5), llo_lik, x=x, y=y, log = TRUE, neg = TRUE, method = "L-BFGS-B",
  #                          lower = lower, upper = upper)
  # max_lik <- -optBayes$value


  # Loop through grid points to get posterior model probability
  temp <- c()
  for(i in 1:nrow(grd)){
    x_new <- LLO(x, delta = grd[i,1], gamma = grd[i,2])   # LLO adjust probs FIRST based on grid point

    optBayes <- stats::optim(c(0.5, 0.5), llo_lik, x=x_new, y=y, log = TRUE, neg = TRUE, method = "L-BFGS-B",
                             lower = lower, upper = upper)
    max_lik <- -optBayes$value
    temp <- c(temp, max_lik)
  }

  # Reshape vector of posterior model probs into matrix for plotting
  z_mat <- matrix(temp, nrow = length(d), ncol = length(g))
  colnames(z_mat) <- g
  rownames(z_mat) <- d
  return(z_mat)
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
    adj_res$sd[i] <- stats::sd(temp)
  }


  utils::write.csv(adj_res, paste0("simstudy_results/adjres_", filen, "_", calib, ".csv"))

  # adj_df <- data.frame(adj_mat)
  # colnames(adj_df) <- expand.grid(paste0(filen, "_"), as.character(1:ncol(adj_df))) %>%
  #   transmute(label = paste0(Var1, Var2)) %>%
  #   .$label
  # return(adj_df)
}
