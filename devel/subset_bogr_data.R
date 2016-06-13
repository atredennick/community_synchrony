####################################################################
##  subset_bogr_data.R: subsets the newly-edited BOGR data frames ##
##  to discard quadrat-year combinations not present in the Chu   ##
##  and Adler 2015 analysis.                                      ##
####################################################################



####
####  INDIVIDUAL GENET DATA
####
rm(list=ls()) # clear workspace each time

# Get "old" BOGR data
old_bogrs <- read.csv("../data/Montana/BOGR_unedited/BOGR_genet_xy.csv")
# Get "new" BOGR data (which includes all the quadrats, not Chengjin's subset)
new_bogrs <- read.csv("../data/Montana/BOGR_edited/BOGR_genet_xy_edited.csv")

# Make new quad-year column in both data frames
new_bogrs$quadyrs <- with(new_bogrs, paste0(quad,year))
old_bogrs$quadyrs <- with(old_bogrs, paste0(quad,year))

# Only use the newdata$quadyears that are in the olddata$quadyears
ROWS_TO_USE <- which(new_bogrs$quadyrs %in% old_bogrs$quadyrs)

# Makes a subset of data that does not include "bad" BOGR quad-years (those quad-years)
#   are absent from the new data), and excludes the quad-years without enough species
#   overlap (e.g., quadyears that may be in the new data, but not in the old data)
new_bogr_grow <- new_bogrs[ROWS_TO_USE,which(colnames(new_bogrs)!="quadyrs")]
write.csv(new_bogr_grow, "../data/Montana/BOGR/BOGR_genet_xy.csv", row.names = FALSE)



####
####  GROWTH DATA
####
rm(list=ls()) # clear workspace each time

# Get "old" BOGR data
old_bogrs <- read.csv("../data/Montana/BOGR_unedited/growDnoNA.csv")
# Get "new" BOGR data (which includes all the quadrats, not Chengjin's subset)
new_bogrs <- read.csv("../data/Montana/BOGR_edited/growDnoNA.csv")

# Make new quad-year column in both data frames
new_bogrs$quadyrs <- with(new_bogrs, paste0(quad,year))
old_bogrs$quadyrs <- with(old_bogrs, paste0(quad,year))

# Only use the newdata$quadyears that are in the olddata$quadyears
ROWS_TO_USE <- which(new_bogrs$quadyrs %in% old_bogrs$quadyrs)

# Makes a subset of data that does not include "bad" BOGR quad-years (those quad-years)
#   are absent from the new data), and excludes the quad-years without enough species
#   overlap (e.g., quadyears that may be in the new data, but not in the old data)
new_bogr_grow <- new_bogrs[ROWS_TO_USE,which(colnames(new_bogrs)!="quadyrs")]
write.csv(new_bogr_grow, "../data/Montana/BOGR/growDnoNA.csv", row.names = FALSE)



####
####  SURVIVAL DATA
####
rm(list=ls()) # clear workspace each time

# Get "old" BOGR data
old_bogrs <- read.csv("../data/Montana/BOGR_unedited/survD.csv")
# Get "new" BOGR data (which includes all the quadrats, not Chengjin's subset)
new_bogrs <- read.csv("../data/Montana/BOGR_edited/survD.csv")

# Make new quad-year column in both data frames
new_bogrs$quadyrs <- with(new_bogrs, paste0(quad,year))
old_bogrs$quadyrs <- with(old_bogrs, paste0(quad,year))

# Only use the newdata$quadyears that are in the olddata$quadyears
ROWS_TO_USE <- which(new_bogrs$quadyrs %in% old_bogrs$quadyrs)

# Makes a subset of data that does not include "bad" BOGR quad-years (those quad-years)
#   are absent from the new data), and excludes the quad-years without enough species
#   overlap (e.g., quadyears that may be in the new data, but not in the old data)
new_bogr_grow <- new_bogrs[ROWS_TO_USE,which(colnames(new_bogrs)!="quadyrs")]
write.csv(new_bogr_grow, "../data/Montana/BOGR/survD.csv", row.names = FALSE)



####
####  RECRUITMENT DATA: AREA
####
rm(list=ls()) # clear workspace each time

# Get "old" BOGR data
old_bogrs <- read.csv("../data/Montana/BOGR_unedited/recArea.csv")
# Get "new" BOGR data (which includes all the quadrats, not Chengjin's subset)
new_bogrs <- read.csv("../data/Montana/BOGR_edited/recArea.csv")

# Make new quad-year column in both data frames
new_bogrs$quadyrs <- with(new_bogrs, paste0(quad,year))
old_bogrs$quadyrs <- with(old_bogrs, paste0(quad,year))

# Only use the newdata$quadyears that are in the olddata$quadyears
ROWS_TO_USE <- which(new_bogrs$quadyrs %in% old_bogrs$quadyrs)

# Makes a subset of data that does not include "bad" BOGR quad-years (those quad-years)
#   are absent from the new data), and excludes the quad-years without enough species
#   overlap (e.g., quadyears that may be in the new data, but not in the old data)
new_bogr_grow <- new_bogrs[ROWS_TO_USE,which(colnames(new_bogrs)!="quadyrs")]
write.csv(new_bogr_grow, "../data/Montana/BOGR/recArea.csv", row.names = FALSE)



####
####  RECRUITMENT DATA: SIZE
####
rm(list=ls()) # clear workspace each time

# Get "old" BOGR data
old_bogrs <- read.csv("../data/Montana/BOGR_unedited/recSize.csv")
# Get "new" BOGR data (which includes all the quadrats, not Chengjin's subset)
new_bogrs <- read.csv("../data/Montana/BOGR_edited/recSize.csv")

# Make new quad-year column in both data frames
new_bogrs$quadyrs <- with(new_bogrs, paste0(quad,year))
old_bogrs$quadyrs <- with(old_bogrs, paste0(quad,year))

# Only use the newdata$quadyears that are in the olddata$quadyears
ROWS_TO_USE <- which(new_bogrs$quadyrs %in% old_bogrs$quadyrs)

# Makes a subset of data that does not include "bad" BOGR quad-years (those quad-years)
#   are absent from the new data), and excludes the quad-years without enough species
#   overlap (e.g., quadyears that may be in the new data, but not in the old data)
new_bogr_grow <- new_bogrs[ROWS_TO_USE,which(colnames(new_bogrs)!="quadyrs")]
write.csv(new_bogr_grow, "../data/Montana/BOGR/recSize.csv", row.names = FALSE)



####
####  QUADRAT COVER
####
rm(list=ls()) # clear workspace each time

# Get "old" BOGR data
old_bogrs <- read.csv("../data/Montana/BOGR_unedited/quadratCover.csv")
# Get "new" BOGR data (which includes all the quadrats, not Chengjin's subset)
new_bogrs <- read.csv("../data/Montana/BOGR_edited/quadratCover.csv")

# Make new quad-year column in both data frames
new_bogrs$quadyrs <- with(new_bogrs, paste0(quad,year))
old_bogrs$quadyrs <- with(old_bogrs, paste0(quad,year))

# Only use the newdata$quadyears that are in the olddata$quadyears
ROWS_TO_USE <- which(new_bogrs$quadyrs %in% old_bogrs$quadyrs)

# Makes a subset of data that does not include "bad" BOGR quad-years (those quad-years)
#   are absent from the new data), and excludes the quad-years without enough species
#   overlap (e.g., quadyears that may be in the new data, but not in the old data)
new_bogr_grow <- new_bogrs[ROWS_TO_USE,which(colnames(new_bogrs)!="quadyrs")]
write.csv(new_bogr_grow, "../data/Montana/BOGR/quadratCover.csv", row.names = FALSE)

