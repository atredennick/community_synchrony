

rm(list=ls())

# Get "old" BOGR data
old_bogrs <- read.csv("../data/Montana/BOGR_unedited/growDnoNA.csv")

# Get "new" BOGR data (which includes all the quadrats, not Chengjin's subset)
new_bogrs <- read.csv("../data/Montana/BOGR/growDnoNA.csv")

# Make new quad-year column in both data frames
new_bogrs$quadyrs <- with(new_bogrs, paste0(quad,year))
old_bogrs$quadyrs <- with(old_bogrs, paste0(quad,year))

# Only use the newdata$quadyears that are in the olddata$quadyears
ROWS_TO_USE <- which(new_bogrs$quadyrs %in% old_bogrs$quadyrs)

# Makes a subset of data that does not include "bad" BOGR quad-years (those quad-years)
#   are absent from the new data), and excludes the quad-years without enough species
#   overlap (e.g., quadyears that may be in the new data, but not in the old data)
DATA_TO_USE <- new_bogrs[ROWS_TO_USE,]

