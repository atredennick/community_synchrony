

rm(list=ls())

old_bogrs <- read.csv("../data/Montana/BOGR_unedited/growDnoNA.csv")
new_bogrs <- read.csv("../data/Montana/BOGR/growDnoNA.csv")
pose <- read.csv("../data/Montana/POSE/growDnoNA.csv")

pose_quadyrs <- with(pose, paste0(quad,year))
unique(pose_quadyrs)

bogr_quadyrs <- with(old_bogrs, paste0(quad,year))
unique(bogr_quadyrs)

new_bogr_quadyrs <- with(new_bogrs, paste0(quad,year))
unique(new_bogr_quadyrs)


new_bogrs$quadyrs <- with(new_bogrs, paste0(quad,year))
old_bogrs$quadyrs <- with(old_bogrs, paste0(quad,year))

tmp <- which(new_bogrs$quadyrs %in% old_bogrs$quadyrs)

test <- new_bogrs[tmp,]
unique(test$quadyrs)
