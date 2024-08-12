## code to prepare `foreclosure` dataset goes here

foreclosure <- read.csv("../ppc/llo_method/llo_paper/case_studies/data/foreclosure_ex_oos.csv", row.names=1)
colnames(foreclosure) <- c("y", "x", "year")

# subset by just 5,000 obs from 2010
foreclosure <- foreclosure[foreclosure$year==2010,]
set.seed(47)
inds <- sample(1:nrow(foreclosure), 5000)
foreclosure <- foreclosure[inds,]
row.names(foreclosure) <- 1:nrow(foreclosure)
usethis::use_data(foreclosure, internal=TRUE, overwrite = FALSE)
