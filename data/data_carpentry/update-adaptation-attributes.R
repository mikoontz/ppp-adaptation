# Title: Update adaptation attributes
#
# Author: Michael Koontz
# Email: mikoontz@gmail.com
#
# Date Created: 20150414
# Last Updated: 20150414
#
# Purpose: Updates the 'adaptation attributes.csv' file with new inbreeding values in the 'observed inbreeding.csv' file (or, rather, expected loss of heterozygosity) if they have been recalculated in the 'observed bottlenecks.R' analysis script.

# I will keep the original H8 and H9 calculations because they are qualitatively the same and were used to randomly assign populations into the "mixed" group.

rm(list=ls())
setwd("/Users/mikoontz/Documents/Research/Tribolium/Demography/Clean Data")

# Read the original attributes data
attributes <- read.csv("adaptation attributes.csv")
# Read the observed inbreeding data, which may have been updated to reflect a different configuration of the inbreeding coefficient calculation
inbred <- read.csv("observed inbreeding.csv")

# Assign more relevant names to the inbred data frame (since each value represents the expected proportion of remaining heterozygosity, we call it H and then the generation number.) H0 represents the expected amount of heterozygosity lost in the initial introduction of individuals. I also change the name of the ID column to population so it matches the attributes data frame.

names(inbred)[c(1, 11:ncol(inbred))] <- c("population", paste0("H", 0:9))

# Trim the inbred data frame to include just the 3 important columns (the key and 2 value columns)
inbred.trim <- subset(inbred, select=c(population, H8, H9))

# Merge the original attributes with the trimmed inbred data frame using the 'population' column. Keep all rows of the x dataframe (attributes) but don't keep any extra rows from the y dataframe (inbred.trim). Thus, we toss the populations from blocks 1 and 2, since they weren't used for the adaptation experiment.

df <- merge(attributes, inbred.trim, by="population", all.x=TRUE, all.y=FALSE)

# Write the file. Analyses should use this adaptation attributes2.csv file
write.csv(df, "adaptation attributes2.csv", row.names=FALSE)
