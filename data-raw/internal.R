# Foreclosure Data from Keefe et al 2017
foreclosure <- read.csv("../ppc/llo_method/llo_paper/case_studies/data/foreclosure_ex_oos.csv", row.names=1)
colnames(foreclosure) <- c("y", "x", "year")

# subset by just 5,000 obs from 2010
foreclosure <- foreclosure[foreclosure$year==2010,]
set.seed(47)
inds <- sample(1:nrow(foreclosure), 5000)
foreclosure <- foreclosure[inds,]


# CIFAR-10 predictions for Vehicle / not vehicle
cifar10 <- read.csv("../ppc/llo_method/llo_paper/case_studies/data/CIFAR-10_veh_preds.csv", row.names=1)


# Head CT scanning data from Levine et al 2013
headCT <- read.csv("../ppc/llo_method/llo_paper/case_studies/data/head_preds.csv")
colnames(headCT) <- c("y", "x")
headCT <- headCT[-which(is.na(headCT$x)),]

# Chest CT scanning data from Levine et al 2013
chestCT <- read.csv("../ppc/llo_method/llo_paper/case_studies/data/chest_preds.csv")
colnames(chestCT) <- c("y", "x")
chestCT <- chestCT[-which(is.na(chestCT$x)),]


usethis::use_data(foreclosure, cifar10, headCT, chestCT, internal=TRUE, overwrite = TRUE)
