# Title: generate tidy adaptation data
# 
# Author: Michael Koontz
# Email: mikoontz@gmail.com
#
# Date Created: 20150411
# Last Updated: 20150414

# This function takes the entered Tribolium flour beetle data from the "Eco-evolutionary consequences of multiple introductions" adaptation experiment (which is in long form), and puts it in a 2 dimensional form to more easily represent the time series. Each population is defined by a unique ID, which is included in each of the produced dataframes here. This makes for easy merging (using ID as a key) with the 'adaptation attributes.csv' file (which includes block number, treatment types, the degree of expected heterozygosity lost in generation 8 and 9, etc.)

# Requires the tidyr package for reshaping the data.
# Input is the long-form entered data in a dataframe object type.
# Returns a list of dataframes representing the different values that are unique to each ID/Generation combination.
# Further manipulations on these data frames are possible using other functions. 

# Load tidyr library
if (!require("tidyr"))
{install.packages("tidyr"); library(tidyr)}

tidy.adapt <- function(adapt)
{
#---------- 
# N[t+1] dataframe
# ---------
head(adapt)
# Subset to relevant parts
a <- subset(adapt, select=c(ID, Generation, Census))
# Spread it
Ntp1 <- spread(a, Generation, Census)
# Check it
head(Ntp1)
tail(Ntp1)

#---------- 
# N[t] dataframe
# ---------

a <- subset(adapt, select=c(ID, Generation, N0))
Nt <- spread(a, Generation, N0)

names(Nt) <- c("ID", paste(9:10))
head(Ntp1)
head(Nt)
tail(Nt)

#---------- 
# Census taker dataframe
# ---------

a <- subset(adapt, select=c(ID, Generation, Person))
person <- spread(a, Generation, Person)
head(person)
tail(person)


return(list(Nt=Nt, Ntp1=Ntp1, person=person))
}
