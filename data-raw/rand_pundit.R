## code to prepare `rand_pundit` dataset goes here

rand_pundit <- read.csv("../ppc/hockey/data/compiled_NHL_pundit_data20_21_FULL.csv",
                        row.names=1, stringsAsFactors=TRUE) %>%
  dplyr::transmute(y = Winner01,
                   x = with_seed(8333, x_rand <- runif(n = nrow(.),
                                                       min = min(.$HomeProb538),
                                                       max = max(.$HomeProb538))))

usethis::use_data(rand_pundit, overwrite = TRUE)

