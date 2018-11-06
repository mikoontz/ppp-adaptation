# Title: observed bottlenecks
#
# Author: Michael Koontz
# Email: mikoontz@gmail.com
#
# Date Created: 20150401
# Last Updated: 20150401

# Intention: Takes actual experimental data and calculates the expected amount of heterozygosity remaining using the inbreeding function based on equation in Crow and Kimura (1970) [also referenced in McCauley and Wade (1981)]

rm(list=ls())
setwd("/Users/mikoontz/Documents/Research/Tribolium/Demography/")
source("Clean-Data/generate-tidy-data.R")
source("Clean-Analyses/inbreed-function.R")

b <- read.csv("Clean-Data/Tribolium-propagule-pressure-data.csv")
attributes <- read.csv("Clean-Data/attributes.csv")  

#-----------
# Set up objects
#-----------

# We need the number of individuals in each generation, the number of migrants, and the number of individuals in the next generation (to get the final census value)

tidy.b <- tidy.beetles(beetles=b)

Nt <- tidy.b$Nt
Ntp1 <- tidy.b$Ntp1
migrants <- tidy.b$migrants

# Correction for population 280, which had a missed survey during the first generation. This value should really be imputed, but the heterozygosity measures are very robust to a single missing data point. 

Nt[280, 3] <- 20
Ntp1[280, 2] <- 20

# This is what the population size would have been if it were used to start the next generation
Nt$'9' <- Ntp1$'9'
Nt$'9'[is.na(Nt$'9')] <- 0

# Migration rate is the migrants dataframe divided by the Nt dataframe. Ignore the ID column for the division (i.e. use the 2nd column to the last column), then add it back in using cbind()
m <- cbind(migrants$ID, migrants[, 2:ncol(migrants)]/Nt[, 2:ncol(Nt)])

#-----------
# Run inbreeding function
#-----------

inbreed.mat <- matrix(Nt$ID, nrow=nrow(Nt), ncol=ncol(Nt))
inbreed.mat <- as.data.frame(inbreed.mat)
names(inbreed.mat) <- names(Nt)

# Go through each row of the 
for (i in 1:nrow(Nt))
{
  inbreed.mat[i, 2:ncol(inbreed.mat)] <- 1 - inbreed(N=Nt[i, 2:ncol(Nt)], m=migrants[i, 2:ncol(migrants)])
}

#-----------
# Combine with attribute data
#-----------

data <- merge(attributes, inbreed.mat, by="ID")
head(data)

#-----------
# Write to file
#-----------
# setwd("/Users/mikoontz/Documents/Research/Tribolium/Demography/Clean Data")
# write.csv(data, "observed inbreeding.csv", row.names=FALSE)

#----------
# Aggregate as desired
#----------

inbred.mean <- aggregate(
                  list(gen0=data$'0', 
                       gen1=data$'1', 
                       gen2=data$'2', 
                       gen3=data$'3', 
                       gen4=data$'4', 
                       gen5=data$'5', 
                       gen6=data$'6', 
                       gen7=data$'7', 
                       gen8=data$'8', 
                       gen9=data$'9'), 
                  by=list(propagule.number=data$'number'), 
                  FUN=function(x) mean(x, na.rm=TRUE) )

inbred.x.after <- function(x, inbred.mean)
{
  plural <- ifelse(x==1, yes = '', no = 's')
  par(mar=c(5,4,5,2)+0.1)
#   layout(matrix(1:2, byrow=TRUE, ncol=2), widths=c(1,0.25))
matplot(x=0:9, t(inbred.mean[, 2:ncol(inbred.mean)]), type="l", lty=1, ylab="Heterozygosity remaining", xlab="Generation", lwd=3, main=paste0("Observed bottlenecks expressed as\nmean expected loss of heterozygosity through time\n(Dots represent heterozygosity remaining ", x, " generation", plural,"\nafter final introduction)"), las=1)
abline(h=inbred.mean[1, x+2], col=1, lty="dashed")
abline(h=inbred.mean[2, x+3], col=2, lty="dashed")
abline(h=inbred.mean[3, x+5], col=3, lty="dashed")
abline(h=inbred.mean[4, x+6], col=4, lty="dashed")
legend("topright", legend=c("1 introduction", "2 introductions", "4 introductions", "5 introductions"), bty="n", lwd=3, col=1:4)
plot.x <- rep(1,4)
y <- c(inbred.mean[1, x+2], inbred.mean[2, x+3], inbred.mean[3, x+5], inbred.mean[4, x+6])
points(x=x+c(0,1,3,4), y, pch=19, cex=1.5, col=1:4)
# plot(plot.x, y, pch=19, cex=1.5, col=1:4, ylab=NA, xlab=NA, xaxt="n", yaxt="n", ylim=range(inbred.mean[,2:ncol(inbred.mean)]))
}

for (i in 1:5)
{
  plural <- ifelse(i==1, '', 's')
  pdf(paste0('Clean Plots/observed bottlenecks ', i, ' generation', plural, ' after last intro.pdf'))
  inbred.x.after(i, inbred.mean)
  dev.off()
}

inbred.x.after(x=5, inbred.mean)
