# Script for model fitting

library(poweRlaw)
library(ggplot2) # pretty plotting
library(quantmod) # finding peaks function
library(poweRlaw) # power law fitting
library(cowplot) # for grids
library(pracma) 
library(dplyr)
library(stringr)




# St1

steps.ninv1 <- read.csv("all_steps_non_invaders1_newold.csv")
steps.ninv1 <- steps.ninv1$x
pl.ninv1 <- conpl$new(steps.ninv1[steps.ninv1>0])
xmin.ninv1 <- estimate_xmin(pl.ninv1)
pl.ninv1$setXmin(xmin.ninv1)
pars.ninv1 <- estimate_pars(pl.ninv1, pl.ninv1$pars)
str_c("The power law exponent for invadors in movement 1 is ",pl.ninv1$pars-1)

pl.ll <- dist_ll(pl.ninv1)
str_c("The negative likelihood for power law is ", pl.ll)

pl.AIC <- 2*length(pl.ninv1$pars) - 2*pl.ll
str_c("The Akaike information for power law is ", pl.AIC)

exp.ninv1 <- conexp$new(steps.ninv1[steps.ninv1>0])

xmin.ninv1 <- estimate_xmin(exp.ninv1)
exp.ninv1$setXmin(xmin.ninv1)
pars.ninv1 <- estimate_pars(exp.ninv1, exp.ninv1$pars)
str_c("The exponential for invadors in movement 1 is ",exp.ninv1$pars)

exp.ll <- dist_ll(exp.ninv1)
str_c("The negative likelihood for exponential is ", exp.ll)

exp.AIC <- 2*length(exp.ninv1$pars) - 2*exp.ll
str_c("The Akaike information for exponential is ", exp.AIC)

exp.AIC <- 2*length(exp.ninv1$pars) - 2*exp.ll
AICs <- c(exp.AIC, pl.AIC)
AICscores <- tibble(AIC_exp = exp.AIC, AIC_PL = pl.AIC)
AICscores

pl.delta.AIC <- pl.AIC - min(AICs)
exp.delta.AIC <- exp.AIC - min(AICs)
pl.weight.AIC <- exp(-pl.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
exp.weight.AIC <- exp(-exp.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
str_c('AIC Weight for power law = ', pl.weight.AIC)
str_c('AIC Weight for exponential = ', exp.weight.AIC)

pdf(file="ninv1_comp_newold.pdf", height = 3.6, width = 3.6)
plot(pl.ninv1, pch = 20, font.main = 1,col = "firebrick2",
     xaxt = "n", yaxt = "n",
     ylab = "", xlab = "")
lines(exp.ninv1, col = "darkslategray3", lw = 3, lty = 4)
lines(pl.ninv1, col = "darkslategray", lw = 3, lty = 5)
legend(0.0007, 0.0025, legend = c(paste("Exp, μ = ", num2str(exp.ninv1$pars)),#, ", lmin = ", num2str(exp.ninv1$xmin), "μm"), 
                                paste("PL, μ = ", num2str(pl.ninv1$pars-1))),#, ", lmin = ", num2str(pl.ninv1$xmin), "μm")), 
                   col = c("darkslategray3", "darkslategrey"), lty = c(4,5), lw = 3, box.lty = 0)
legend(0.0007, 0.0065, legend = paste("Total steps =", length(steps.ninv1)), box.lty = 0)
axis(1, at = c(0.01, 0.1, 1, 10, 100, 300),
     labels = c("0.01", "0.1", "1", "10", "100", "300"))
axis(2, at = c(0.001, 0.01, 0.1, 1), labels = c("0.001", "0.01", "0.1", "1"))
dev.off()




# St2

steps.ninv2 <- read.csv("all_steps_non_invaders2_newold.csv")
steps.ninv2 <- steps.ninv2$x
pl.ninv2 <- conpl$new(steps.ninv2[steps.ninv2>0])

xmin.ninv2 <- estimate_xmin(pl.ninv2)
pl.ninv2$setXmin(xmin.ninv2)
str_c("The power law exponent for non-invadors in movement 2 is ",pl.ninv2$pars-1)

pl.ll <- dist_ll(pl.ninv2)
str_c("The negative likelihood for power law is ", pl.ll)

pl.AIC <- 2*length(pl.ninv2$pars) - 2*pl.ll
str_c("The Akaike information for power law is ", pl.AIC)

exp.ninv2 <- conexp$new(steps.ninv2[steps.ninv2>0])

xmin.ninv2 <- estimate_xmin(exp.ninv2)
exp.ninv2$setXmin(xmin.ninv2)
pars.ninv2 <- estimate_pars(exp.ninv2, exp.ninv2$pars)
str_c("The exponential for non-invadors in movement 2 is ",exp.ninv2$pars)

exp.ll <- dist_ll(exp.ninv2)
str_c("The negative likelihood for exponential is ", exp.ll)

exp.AIC <- 2*length(exp.ninv2$pars) - 2*exp.ll
str_c("The Akaike information for exponential is ", exp.AIC)

exp.AIC <- 2*length(exp.ninv2$pars) - 2*exp.ll
AICs <- c(exp.AIC, pl.AIC)
AICscores <- tibble(AIC_exp = exp.AIC, AIC_PL = pl.AIC)
AICscores

pl.delta.AIC <- pl.AIC - min(AICs)
exp.delta.AIC <- exp.AIC - min(AICs)
pl.weight.AIC <- exp(-pl.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
exp.weight.AIC <- exp(-exp.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
str_c('AIC Weight for power law = ', pl.weight.AIC)
str_c('AIC Weight for exponential = ', exp.weight.AIC)

pdf(file="ninv2_comp_newold.pdf", height = 5.3, width = 5.3)
plot(pl.ninv2, pch = 20, font.main = 1,col = "firebrick2", main = "NINV St2",
     xaxt = "n", yaxt = "n",
     ylab = "Probability", xlab = "Step Length (μm)")
lines(exp.ninv2, col = "darkslategray3", lw = 3, lty = 4)
lines(pl.ninv2, col = "darkslategray", lw = 3, lty = 5)
legend(0.0015, 0.001, legend = c(paste("Exp, μ = ", num2str(exp.ninv2$pars), ", lmin = ", num2str(exp.ninv2$xmin), "μm"), 
                                paste("PL, μ = ", num2str(pl.ninv2$pars-1), ", lmin = ", num2str(pl.ninv2$xmin), "μm")), 
                   col = c("darkslategray3", "darkslategrey"), lty = c(4,5), lw = 3, box.lty = 0)
legend(0.0015, 0.0019, legend = paste("Total steps =", length(steps.ninv2)), box.lty = 0)
axis(1, at = c(0.01, 0.1, 1, 10, 100, 250),
     labels = c("0.01", "0.1", "1", "10", "100", "250"))
axis(2, at = c(0.001, 0.01, 0.1, 0.5), c("0.001", "0.01", "0.1", "0.5"))
dev.off()




# create the data of step lenghts st1+st2 with mixed states allowed

data_non_invader <- read.csv("ninv_newold_states.csv")
all.steps.non.invaders <- c()

for(parasite in unique(data_non_invader$ID))
{
peaks.x <- findPeaks(data_non_invader$x[data_non_invader$ID == parasite]) - 1
peaks.y <- findPeaks(data_non_invader$y[data_non_invader$ID == parasite]) - 1
valleys.x <- findValleys(data_non_invader$x[data_non_invader$ID == parasite]) - 1
valleys.y <- findValleys(data_non_invader$y[data_non_invader$ID == parasite]) - 1
peaks.x <- sort(c(peaks.x, valleys.x), decreasing = FALSE)
peaks.y <- sort(c(peaks.y, valleys.y), decreasing = FALSE)
step.data_non_invader.x <- replicate(length(peaks.x) - 1, 0)
step.data_non_invader.y <- replicate(length(peaks.y) - 1, 0)
StepX.ind <- data_non_invader$StepX[data_non_invader$ID == parasite]
StepY.ind <- data_non_invader$StepY[data_non_invader$ID == parasite]
Cluster.ind <- as.numeric(data_non_invader$State[data_non_invader$ID == parasite])
for(ind in c(1:(length(peaks.x) - 1)))
{
step.data_non_invader.x[ind] = sum(StepX.ind[peaks.x[ind]:peaks.x[ind + 1]], na.rm=TRUE)
}
for(ind in c(1:(length(peaks.y) - 1)))
{
step.data_non_invader.y[ind] = sum(StepY.ind[peaks.y[ind]:peaks.y[ind + 1]], na.rm=TRUE)
}

combined_steps <- c(step.data_non_invader.x,step.data_non_invader.y)
combined_steps <- combined_steps[combined_steps > 0]

all.steps.non.invaders <- append(all.steps.non.invaders, combined_steps)
}

write.csv(all.steps.non.invaders, "all_steps_non_invaders_newold.csv")

# Plot the power law for ninv st1+st2

steps.ninv <- read.csv("all_steps_non_invaders_newold.csv")
steps.ninv <- steps.ninv$x
pl.ninv <- conpl$new(steps.ninv[steps.ninv>0.00001])

xmin.ninv <- estimate_xmin(pl.ninv)
pl.ninv$setXmin(xmin.ninv)
pars.ninv <- estimate_pars(pl.ninv, pl.ninv$pars)
str_c("The power law exponent for invadors in movement 1 is ",pl.ninv$pars-1)

pl.ll <- dist_ll(pl.ninv)
str_c("The negative likelihood for power law is ", pl.ll)

pl.AIC <- 2*length(pl.ninv$pars) - 2*pl.ll
str_c("The Akaike information for power law is ", pl.AIC)

exp.ninv <- conexp$new(steps.ninv[steps.ninv>0])

xmin.ninv <- estimate_xmin(exp.ninv)
exp.ninv$setXmin(xmin.ninv)
pars.ninv <- estimate_pars(exp.ninv, exp.ninv$pars)
str_c("The exponential for invadors in movement 1 is ",exp.ninv$pars)

exp.ll <- dist_ll(exp.ninv)
str_c("The negative likelihood for exponential is ", exp.ll)

exp.AIC <- 2*length(exp.ninv$pars) - 2*exp.ll
str_c("The Akaike information for exponential is ", exp.AIC)

exp.AIC <- 2*length(exp.ninv$pars) - 2*exp.ll
AICs <- c(exp.AIC, pl.AIC)
AICscores <- tibble(AIC_exp = exp.AIC, AIC_PL = pl.AIC)
AICscores

pl.delta.AIC <- pl.AIC - min(AICs)
exp.delta.AIC <- exp.AIC - min(AICs)
pl.weight.AIC <- exp(-pl.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
exp.weight.AIC <- exp(-exp.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
str_c('AIC Weight for power law = ', pl.weight.AIC)
str_c('AIC Weight for exponential = ', exp.weight.AIC)

pdf(file="ninv_comp_newold_v2.pdf", height = 5.3, width = 5.3)
plot(pl.ninv, pch = 20, font.main = 1,col = "firebrick2", main = "NINV St1+St2",
     xaxt = "n", yaxt = "n",
     ylab = "Probability", xlab = "Step Length (μm)")
lines(exp.ninv, col = "darkslategray3", lw = 3, lty = 4)
lines(pl.ninv, col = "darkslategray", lw = 3, lty = 5)
legend(0.001, 0.0008, legend = c(paste("Exp, μ = ", num2str(exp.ninv$pars), ", lmin = ", num2str(exp.ninv$xmin), "μm"), 
                                paste("PL, μ = ", num2str(pl.ninv$pars-1), ", lmin = ", num2str(pl.ninv$xmin), "μm")), 
                               col = c("darkslategray3", "darkslategrey"), lty = c(4,5), lw = 3, box.lty = 0)
legend(0.001, 0.0015, legend = paste("Total steps =", length(steps.ninv)), box.lty = 0)
axis(1, at = c(0.01, 0.1, 1, 10, 100, 300),
     labels = c("0.01", "0.1", "1", "10", "100", "300"))
axis(2, at = c(0.001, 0.01, 0.1, 1), labels = c("0.001", "0.01", "0.1", "1"))
dev.off()