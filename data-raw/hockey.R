## code to prepare `hockey` dataset goes here

original538 <- read.csv("../ppc/hockey/data/compiled_NHL_pundit_data20_21_FULL.csv",
                        row.names=1, stringsAsFactors=TRUE)
hockey <- dplyr::transmute(original538, y = Winner01,
            x = HomeProb538)

usethis::use_data(hockey, overwrite = TRUE)
