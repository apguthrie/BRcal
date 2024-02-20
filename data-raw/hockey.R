## code to prepare `hockey` dataset goes here

hockey <- read.csv("../ppc/hockey/data/compiled_NHL_pundit_data20_21_FULL.csv",
                   row.names=1, stringsAsFactors=TRUE) %>%
  dplyr::transmute(y = Winner01,
                   x = HomeProb538,
                   rand = withr::with_seed(8333, x_rand <- runif(n = nrow(.),
                                                                 min = min(.$HomeProb538),
                                                                 max = max(.$HomeProb538))),
                   winner = Winner)

usethis::use_data(hockey, overwrite = FALSE)
