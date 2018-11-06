# Title: adaptation analysis
#
# Author: Michael Koontz
# Email: mikoontz@gmail.com
#
# Date Created: 20150414
# Last Updated: 20150414
#
# Test 'Eco-evolutionary consequences of multiple introductions for colonizing individuals' main experiment data for adaptation
#
# In the generation being analyzed here, all populations were started with 16 individuals and the number of individuals in the next generation was calculated.

rm(list=ls())
setwd("/Users/mikoontz/Documents/Research/Tribolium/Demography")
source("Clean-Data/generate-tidy-adaptation-data.R")

adapt <- read.csv("Clean-Data/Tribolium-adaptation-data-long.csv")
attributes <- read.csv("Clean-Data/adaptation-attributes2.csv")

#-----------
# Load libraries
#-----------

if (!require("lme4"))
{install.packages("lme4"); library(lme4)}
if (!require("dplyr"))
{install.packages("dplyr"); library(dplyr)}
if (!require("tidyr"))
{install.packages("tidyr"); library(tidyr)}
if (!require("lsmeans"))
{install.packages("lsmeans"); library(lsmeans)}

#-----------
# Set up objects
#-----------

# We need the number of individuals in each generation, the number of migrants, and the number of individuals in the next generation (to get the final census value)

tidy.adapt <- tidy.adapt(adapt)

Nt <- tidy.adapt$Nt
Ntp1 <- tidy.adapt$Ntp1

head(Nt)
head(Ntp1)

# One way to measure fitness is as pop size at generation 11 dividided by 16

W <- Ntp1$'11' / Nt$'10'

W.attr <- cbind(attributes, data.frame(Nt=Nt$'10', Ntp1=Ntp1$'11'), W)

# Remove NAs
df <- W.attr[!is.na(W.attr$W), ]
head(df)
df$number <- as.factor(df$number)

# Correlation between old calculation of H9/H0 and new. Old one was too pessimistic
plot(df$old.H9, df$H9, pch=19)
abline(b=1, a=0)
# Discrepancy tended to be biggest with more introductions (generally the additional introductions did more to reduce the loss in heterozygosity with the new calculation than with the old calculation)
# The 2-period shape here represents the 2 temporal blocks and data arranged as 20x1, 10x2, 5x4, and 4x5 within each block
plot(df$H9-df$old.H9, pch=19)

# All populations started at the same number of individuals (16) so we just use the number of individuals in the next generation as a fitness measurement and use a Poisson regression. This avoids issues with taking the log of 0 if we were to use finite rates of growth

# Just look at treatment effect (propagule number and environmental stability) using only the mixed populations. These populations should not be expressing genetic load, but would still be adapted. Question asked: Was there fundamental differences in fitness amongst experimental treatments?
m1 <- glmer(Ntp1 ~ number*env + (1 | block) + (1 | population), family="poisson", data=subset(df, subset=group == "mixed"), control=glmerControl(optCtrl=list(maxfun=1e6), optimizer="bobyqa"))
summary(m1)

# All treatments appeared to have similar maximum fitness ()

newdata <- expand.grid(number=factor(c(1,2,4,5)), env=factor(c("stable", "fluctuating")))
preds <- predict(m1, type="response", newdata=newdata, re.form=NA)

CI <- confint(m1, method="Wald")
# Really similar CI for bootstrap method
# CI2 <- confint(m3, method="boot")

mat <- matrix(model.matrix(~ number*env, data=newdata), nrow=8)

lwr <- exp(mat %*% CI[,1])
upr <- exp(mat %*% CI[,2])

stab <- which(newdata$env=="stable")
fluc <- which(newdata$env=="fluctuating")

xvals <- rep(c(1,2,4,5), times=2)
plot(xvals, preds, type="n", main="Fitness by experimental treatment", xlab="Propagule number", xaxt="n", ylim=c(0, 100), ylab=expression(bar(W)))
abline(h=16, col="blue", lty="dashed", lwd=3)
axis(side=1, at=c(1,2,4,5))
x.offset <- c(rep(-0.1, 4), rep(0.1, 4))
arrows(x0=xvals + x.offset, y0=lwr, y1=upr, code=0, lwd=3, col=c(rep(1, 4), rep(2, 4)))
legend("top", legend=c("stable", "fluctuating", "R=1"), col=c(1:2, 4), lwd=4)

# Interesting. Higher standard deviation of number of individuals for the stable environment treatment.
mu.sd.table <- aggregate(subset(df, subset=group == "mixed")$Ntp1, by=list(subset(df, subset=group == "mixed")$number, subset(df, subset=group == "mixed")$env), FUN=function(x) c(mean(x), sd(x)))
colnames(mu.sd.table)[3] <- ""
colnames(mu.sd.table[,3]) <- c("meanW", "sdW")
mu.sd.table

#------------------------

# Simple group comparison of fitness. Groups in this case are the low heterozygosity/high inbreeding, high heterozygosity/low inbreeding, mixed (no inbreeding, but adapted), and control (no inbreeding, not adapted). Groups are named based on amount of expected inbreeding (i.e. bottlenecking through time) so "low" means low inbreeding/high heterozygosity
m2 <- glmer(Ntp1 ~ group + (1 | block) + (1 | population), family="poisson", data=df, control=glmerControl(optCtrl=list(maxfun=1e6), optimizer="bobyqa"))
summary(m2)

# Clearly important differences

lsmeans::lsmeans(m2, pairwise ~ group, adjust="none")

newdata2 <- data.frame(group=c("control", "high", "low", "mixed"))

preds2 <- exp(predict(m2, newdata2, re.form=NA))
CI2 <- confint(m2, method="Wald")[-(1:2), ]

mat2 <- model.matrix(~ group, data=newdata2)

lwr2 <- exp(mat2 %*% CI2[,1])
upr2 <- exp(mat2 %*% CI2[,2])

# pdf("Clean-Plots/fitness-by-treatment.pdf")
xvals2 <- 1:4
plot(xvals2, preds2, type="p", main="Fitness by experimental treatment", xlab=NA, xaxt="n", ylim=c(12, 41), ylab=NA, pch=19, las=1)
abline(h=16, col="blue", lty="dashed", lwd=3)
axis(side=1, at=1:4, labels=c("no inbreeding\nnot adapted", "more inbred (low H9/H0)\nmaybe adapted", "less inbred (high H9/H0)\nmaybe adapted", "mixed\nmaybe adapted"), padj=0.5, cex.axis=0.75)
mtext(side=2, text=expression(bar(W)), las=1, line=3)
arrows(x0=xvals2, y0=lwr2, y1=upr2, code=0, lwd=3)
text(x=1:4, y=rep(40, 4), labels=c("a", "ab", "b", "c"))
# dev.off()

#-----------------------------
# Okay, we see that adaptation likely occurred. We see that adaptation didn't occur differently depending on environmental and propagule number treatment. But we also know that there are some differences in expected remaining heterozygosity by treatment. Let's explore that.

# More introductions had much less loss in heterozygosity by the 9th generation. Doesn't include mixed or control groups.
aggregate(H9 ~ number, FUN=mean, data=df)

# Visualize correlation between introduction scenarios and H9/H0. 
# Get propagule number as an integer
num <- c(0:2, 4:5)[as.numeric(df$number)]
# Only plot when propagule number > 0 (i.e. not control groups). This plot excludes mixed groups also (because their H9 value is NA)

plot(num[num>0], df$H9[num>0], pch=20, xlab="Propagule number", ylab=NA, las=1, yaxt="n", main="Expected heterozygosity remaining vs. introduction scenario")
axis(side=2, at=seq(0.6, 0.9, by=0.1), las=1)
mtext(side=2, text=expression(frac(H[9], H[0])), line=2.75, las=1)

plot(df$H9, df$Ntp1, pch=20, xlab=NA, ylab=NA, yaxt="n", main="Fitness vs. expected heterozygosity remaining")
mtext(side=1, text=expression(frac(H[9], H[0])), line=4)
axis(side=2, at=seq(0, 50, by=10), las=1)
mtext(side=2, text=expression(N[t+1]), las=1, line=2.5)

# First test the effect of propagule number, environmental stability, and expected heterozygosity on fitness. Subset to exclude control populations and populations that were mixed.
m3 <- glmer(Ntp1 ~ H9*number*env + (1 | block) + (1 | population), family="poisson", data=subset(df, subset=!(group %in% c("mixed", "control"))), control=glmerControl(optCtrl=list(maxfun=1e6), optimizer="bobyqa"), na.action=na.exclude)
summary(m3)
# Only expected heterozygosity is signficant given other covariates

# pdf("Clean Plots/fitness vs heterozygosity remaining.pdf")
preds3 <- predict(m3, type="response")
plot(df$H9[!is.na(df$H9)], df$Ntp1[!is.na(df$H9)], pch=20, xlab=NA, ylab=NA, yaxt="n", main="Fitness vs. expected heterozygosity remaining")
mtext(side=1, text=expression(frac(H[9], H[0])), line=4)
axis(side=2, at=seq(0, 50, by=10), las=1)
mtext(side=2, text=expression(N[t+1]), las=1, line=2.5)
points(df$H9[!is.na(df$H9)], preds3, col="red", pch=20)
legend("topleft", legend=c("Observed", "Modeled"), col=1:2, pch=20)

# Just use expected heterozygosity remaining for plotting purposes
pdf("Clean-Plots/fitness-vs-heterozygosity-remaining-simple.pdf")
m4 <- glmer(Ntp1 ~ H9 + (1 | block) + (1 | population), family="poisson", data=subset(df, subset=!(group %in% c("mixed", "control"))), control=glmerControl(optCtrl=list(maxfun=1e6), optimizer="bobyqa"))
coefs <- summary(m4)$coefficients[, "Estimate"]
x <- c(0.6, 0.95)
y <- exp(coefs[1] + coefs[2]*x)

plot(df$H9[!is.na(df$H9)], df$Ntp1[!is.na(df$H9)], pch=20, xlab=NA, ylab=NA, yaxt="n", main="Fitness vs. expected heterozygosity remaining")
mtext(side=1, text=expression(frac(H[9], H[0])), line=4)
axis(side=2, at=seq(0, 50, by=10), las=1)
axis(side=2, at=16, labels=16, col.axis="blue", las=1)
abline(h=16, col="blue", lty="dashed")
mtext(side=2, text=expression(N[t+1]), las=1, line=2.5)

lines(x,y, lwd=3, col="red", lty="solid")
dev.off()

#-----------------------