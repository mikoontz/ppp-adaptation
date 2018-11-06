# Title: gene swamping
#
# Author: Michael Koontz
# Email: mikoontz@gmail.com
#
# Date Created: 20150404
# Last Updated: 20150626
#
# Test 'Eco-evolutionary consequences of multiple introductions for colonizing individuals' main experiment data for signatures of gene swamping
#
# Perhaps 10x2 did worse because of a large input of maladaptive alleles into the population. If so, we would expect that a greater migration rate should yield a lower population growth rate, accounting for population size.
#
# Calculate population growth rates, migration rates, and population sizes.

rm(list=ls())
setwd("/Users/mikoontz/Documents/Research/Tribolium/Demography/Clean-Data")
source("generate-tidy-data.R")

b <- read.csv("Tribolium-propagule-pressure-data.csv")
attributes <- read.csv("attributes.csv")  

#-----------
# Load libraries
#-----------

if (!require("lme4"))
{install.packages("lme4"); library(lme4)}

#-----------
# Set up objects
#-----------

# We need the number of individuals in each generation, the number of migrants, and the number of individuals in the next generation (to get the final census value)

tidy.b <- tidy.beetles(beetles=b)

Nt <- tidy.b$Nt
Ntp1 <- tidy.b$Ntp1
colnames(Ntp1) <- colnames(Nt)
migrants <- tidy.b$migrants
head(migrants)
head(Ntp1)
head(Nt)

# Lambda matrix is the population growth rate per generation
lambda <- cbind(Ntp1$ID, Ntp1[,2:ncol(Ntp1)]/Nt[,2:ncol(Nt)])
tail(lambda)
colnames(lambda)[1] <- "ID"

# Migration rate is the migrants dataframe divided by the Nt dataframe. Ignore the ID column for the division (i.e. use the 2nd column to the last column), then add it back in using cbind()

m <- cbind(migrants$ID, migrants[, 2:ncol(migrants)]/Nt[, 2:ncol(Nt)])
head(m)

# Rename columns so the migration rate, population size, and growth rates all match up to the same generation
colnames(m) <- colnames(lambda)
colnames(migrants) <- colnames(lambda)
# colnames(Nt) <- colnames(lambda)

# Put all covariates and responses into long form
m.long <- gather(m, gen, "m", -ID)
nmig.long <- gather(migrants, gen, "nmig", -ID)
l.long <- gather(lambda, gen, "lambda", -ID)
nt.long <- gather(Nt, gen, "Nt", -ID)
ntp1.long <- gather(Ntp1, gen, "Ntp1", -ID)

# Merge covariates and responses together
swamp <- merge(m.long, nmig.long)
swamp <- merge(swamp, nt.long)
swamp <- merge(swamp, l.long)
swamp <- merge(swamp, ntp1.long)

# Reorder for ease of viewing
swamp <- swamp[order(swamp$ID, swamp$gen), ]
# swamp <- swamp[complete.cases(swamp), ]
swamp$gen <- as.numeric(swamp$gen)
swamp <- merge(attributes, swamp)

swamp$scalem <- scale(swamp$m, center = TRUE, scale = TRUE)
swamp$scaleNt <- scale(swamp$Nt, center=TRUE, scale=TRUE)
swamp$scalegen <- scale(swamp$gen, center=TRUE, scale=TRUE)

swamp$ID <- as.factor(swamp$ID)
swamp$block <- as.factor(swamp$block)

swamp$migrants <- swamp$m * swamp$Nt
swamp$log.lambda <- log(swamp$lambda)
# Exclude populations with an introduction gap
# swamp <- subset(swamp, subset=!gap)

df <- subset(swamp, subset=(lambda > 0))

# Including the linearized Ricker model (with Allee effect)
names(df)


fit1 <- lmer(log.lambda ~ m*gen + Nt + I(log(Nt)) + (1 | block), data=df)
summary(fit1)
CI <- confint(fit1)
CI

fit1b <- lmer(log.lambda ~ m + Nt + I(log(Nt)) + (1 | ID) + (1 | block), data=df)
summary(fit1b)
CIb <- confint(fit1b)


fit2 <- lmer(log.lambda ~ as.factor(number)*environment + Nt + I(log(Nt)) + (1 | block) + (1 | ID), data=df)
summary(fit2)
CI2 <- confint(fit2)
# CI2boot <- confint(fit2, method="boot", nsim=1000)

# Visualize the interaction by picking 3 migration rates (0, 0.2 which is approximately the mean, and 1)
# Pick some reasonable population sizes

newdata.full <- expand.grid(m=c(0, mean(df$m, na.rm=TRUE), 1), gen=1:10, Nt=c(1:4, seq(5,70, by=5)))

preds.full <- predict(fit1, newdata=newdata.full, re.form=NA, type="response")
colors <- match(newdata.full$m, c(0, mean(df$m, na.rm=TRUE), 1))


for (i in 1:max(df$gen))
{
  plot(newdata.full$Nt, exp(preds.full), type="n", xlab="Nt", ylab="R prediction", main=paste("Predicted population growth rate\nfor different Nt and different migration rates\nGeneration", i), ylim=c(0.5, 2.0))
  
lines(newdata.full$Nt[newdata.full$m == 0 & newdata.full$gen == i], exp(preds.full[newdata.full$m == 0 & newdata.full$gen == i]), col=1, type="l", lwd=3)

lines(newdata.full$Nt[newdata.full$m == mean(df$m, na.rm=TRUE) & newdata.full$gen == i], exp(preds.full[newdata.full$m == mean(df$m, na.rm=TRUE) & newdata.full$gen == i]), col=2, type="l", lwd=3)

lines(newdata.full$Nt[newdata.full$m == 1 & newdata.full$gen == i], exp(preds.full[newdata.full$m == 1 & newdata.full$gen == i]), col=3, type="l", lwd=3)

legend("topright", col=1:3, legend=c("No migration", "Mean migration", "Only migrants"),  lwd=3)
}

names(df)

# dd <- subset(df, subset=gen%in%1)
# plot(dd$Nt, dd$Ntp1, pch=20)
# lines(newdata.full$Nt[newdata.full$m == 1 & newdata.full$gen == 1], exp(preds.full[newdata.full$m == 1 & newdata.full$gen == 1])*newdata.full$Nt[newdata.full$m == 1 & newdata.full$gen == 1])
# 
# dd <- subset(df, subset=gen%in%2)
# plot(dd$Nt, dd$Ntp1, pch=20)
# lines(newdata.full$Nt[newdata.full$m == 0 & newdata.full$gen == 2], exp(preds.full[newdata.full$m == 0 & newdata.full$gen == 2])*newdata.full$Nt[newdata.full$m == 0 & newdata.full$gen == 2])


