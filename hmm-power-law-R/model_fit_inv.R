library(poweRlaw)
library(ggplot2) # pretty plotting
library(quantmod) # finding peaks function
library(poweRlaw) # power law fitting
library(cowplot) # for grids
library(pracma) 
library(dplyr)
library(stringr)

steps.inv1 <- read.csv("all_steps_invaders1_old.csv")
steps.inv1 <- steps.inv1$x
pl.inv1 <- conpl$new(steps.inv1[steps.inv1>0])

xmin.inv1 <- estimate_xmin(pl.inv1)
pl.inv1$setXmin(xmin.inv1)
pars.inv1 <- estimate_pars(pl.inv1, pl.inv1$pars)
str_c("The power law exponent for invadors in movement 1 is ",pl.inv1$pars-1)

pl.ll <- dist_ll(pl.inv1)
str_c("The negative likelihood for power law is ", pl.ll)

pl.AIC <- 2*length(pl.inv1$pars) - 2*pl.ll
str_c("The Akaike information for power law is ", pl.AIC)

exp.inv1 <- conexp$new(steps.inv1[steps.inv1>0])

xmin.inv1 <- estimate_xmin(exp.inv1)
exp.inv1$setXmin(xmin.inv1)
pars.inv1 <- estimate_pars(exp.inv1, exp.inv1$pars)
str_c("The exponential for invadors in movement 1 is ",exp.inv1$pars)

exp.ll <- dist_ll(exp.inv1)
str_c("The negative likelihood for exponential is ", exp.ll)

exp.AIC <- 2*length(exp.inv1$pars) - 2*exp.ll
str_c("The Akaike information for exponential is ", exp.AIC)

exp.AIC <- 2*length(exp.inv1$pars) - 2*exp.ll
AICs <- c(exp.AIC, pl.AIC)
AICscores <- tibble(AIC_exp = exp.AIC, AIC_PL = pl.AIC)
AICscores

pl.delta.AIC <- pl.AIC - min(AICs)
exp.delta.AIC <- exp.AIC - min(AICs)
pl.weight.AIC <- exp(-pl.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
exp.weight.AIC <- exp(-exp.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
str_c('AIC Weight for power law = ', pl.weight.AIC)
str_c('AIC Weight for exponential = ', exp.weight.AIC)

pdf(file="inv1_comp_old.pdf", height = 3.6, width = 3.6)
plot(pl.inv1, pch = 20, font.main = 1,col = "blue3",
     xaxt = "n",
     ylab = "", xlab = "")
lines(exp.inv1, col = "darkslategray3", lw = 3, lty = 4)
lines(pl.inv1, col = "darkslategray", lw = 3, lty = 5)
legend(0.0007, 0.012, legend = c(paste("Exp, μ = ", num2str(exp.inv1$pars)),#", lmin = ", num2str(exp.inv1$xmin), "μm"), 
                                paste("PL, μ = ", num2str(pl.inv1$pars-1))),#, ", lmin = ", num2str(pl.inv1$xmin), "μm")), 
                   col = c("darkslategray3", "darkslategrey"), lty = c(4,5), lw = 3, box.lty = 0)
legend(0.0007, 0.024, legend = paste("Total steps =", length(steps.inv1)), box.lty = 0)
axis(1, at = c(0.01, 0.1, 1, 10, 100),
     labels = c("0.01", "0.1", "1", "10", "100"))
dev.off()

steps.inv2 <- read.csv("all_steps_invaders2_old.csv")
steps.inv2 <- steps.inv2$x
pl.inv2 <- conpl$new(steps.inv2[steps.inv2>0])

xmin.inv2 <- estimate_xmin(pl.inv2)
pl.inv2$setXmin(xmin.inv2)
pars.inv2 <- estimate_pars(pl.inv2, pl.inv2$pars)
str_c("The power law exponent for invadors in movement 2 is ",pl.inv2$pars-1)

pl.ll <- dist_ll(pl.inv2)
str_c("The negative likelihood for power law is ", pl.ll)

pl.AIC <- 2*length(pl.inv2$pars) - 2*pl.ll
str_c("The Akaike information for power law is ", pl.AIC)

exp.inv2 <- conexp$new(steps.inv2[steps.inv2>0])

xmin.inv2 <- estimate_xmin(exp.inv2)
exp.inv2$setXmin(xmin.inv2)
pars.inv2 <- estimate_pars(exp.inv2, exp.inv2$pars)
str_c("The exponential for invadors in movement 2 is ",exp.inv2$pars)

exp.ll <- dist_ll(exp.inv2)
str_c("The negative likelihood for exponential is ", exp.ll)

exp.AIC <- 2*length(exp.inv2$pars) - 2*exp.ll
str_c("The Akaike information for exponential is ", exp.AIC)

pl.delta.AIC <- pl.AIC - min(AICs)
exp.delta.AIC <- exp.AIC - min(AICs)
pl.weight.AIC <- exp(-pl.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
exp.weight.AIC <- exp(-exp.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
str_c('AIC Weight for power law = ', pl.weight.AIC)
str_c('AIC Weight for exponential = ', exp.weight.AIC)

pdf(file="inv2_comp_old.pdf", height = 5.3, width = 5.3)
plot(pl.inv2, pch = 20, font.main = 1,col = "blue3", main = "INV St2",
     xaxt = "n",
     ylab = "Probability", xlab = "Step Length (μm)")
lines(exp.inv2, col = "darkslategray3", lw = 3, lty = 4)
lines(pl.inv2, col = "darkslategray", lw = 3, lty = 5)
legend(0.02, 0.01, legend = c(paste("Exp, μ = ", num2str(exp.inv2$pars),", lmin = ", num2str(exp.inv2$xmin), "μm"), 
                                paste("PL, μ = ", num2str(pl.inv2$pars-1), ", lmin = ", num2str(pl.inv2$xmin), "μm")), 
                   col = c("darkslategray3", "darkslategrey"), lty = c(4,5), lw = 3, box.lty = 0)
legend(0.02, 0.015, legend = paste("Total steps =", length(steps.inv2)), box.lty = 0)
axis(1, at = c(0.01, 0.1, 1, 10, 100),
     labels = c("0.01", "0.1", "1", "10", "100"))
dev.off()

data_invader <- read.csv("inv_old_states.csv")
all.steps.invaders <- c()

for (parasite in unique(data_invader$ID))
{
peaks.x <- findPeaks(data_invader$x[data_invader$ID == parasite]) - 1
peaks.y <- findPeaks(data_invader$y[data_invader$ID == parasite]) - 1
valleys.x <- findValleys(data_invader$x[data_invader$ID == parasite]) - 1
valleys.y <- findValleys(data_invader$y[data_invader$ID == parasite]) - 1
peaks.x <- sort(c(peaks.x, valleys.x), decreasing = FALSE)
peaks.y <- sort(c(peaks.y, valleys.y), decreasing = FALSE)
step.data_invader.x <- replicate(length(peaks.x) - 1, 0)
step.data_invader.y <- replicate(length(peaks.y) - 1, 0)
StepX.ind <- data_invader$StepX[data_invader$ID == parasite]
StepY.ind <- data_invader$StepY[data_invader$ID == parasite]
Cluster.ind <- as.numeric(data_invader$State[data_invader$ID == parasite])
for(ind in c(1:(length(peaks.x) - 1)))
{
step.data_invader.x[ind] = sum(StepX.ind[peaks.x[ind]:peaks.x[ind + 1]], na.rm=TRUE)
}
for(ind in c(1:(length(peaks.y) - 1)))
{
step.data_invader.y[ind] = sum(StepY.ind[peaks.y[ind]:peaks.y[ind + 1]], na.rm=TRUE)
}

combined_steps <- c(step.data_invader.x, step.data_invader.y)
combined_steps <- combined_steps[combined_steps > 0]

all.steps.invaders <- append(all.steps.invaders, combined_steps)
}

write.csv(all.steps.invaders, "all_steps_invaders_newold.csv")

steps.inv <- read.csv("all_steps_invaders_newold.csv")
steps.inv <- steps.inv$x

pl.inv <- conpl$new(steps.inv[steps.inv>0])

xmin.inv <- estimate_xmin(pl.inv)
pl.inv$setXmin(xmin.inv)
pars.inv <- estimate_pars(pl.inv, pl.inv$pars)
str_c("The power law exponent for invadors in movement 2 is ",pl.inv$pars)

pl.ll <- dist_ll(pl.inv)
str_c("The negative likelihood for power law is ", pl.ll)

pl.AIC <- 2*length(pl.inv$pars) - 2*pl.ll
str_c("The Akaike information for power law is ", pl.AIC)

exp.inv <- conexp$new(steps.inv[steps.inv>0])

xmin.inv <- estimate_xmin(exp.inv)
exp.inv$setXmin(xmin.inv)
pars.inv <- estimate_pars(exp.inv, exp.inv$pars)
str_c("The exponential for invadors in movement 2 is ",exp.inv$pars)

exp.ll <- dist_ll(exp.inv)
str_c("The negative likelihood for exponential is ", exp.ll)

exp.AIC <- 2*length(exp.inv$pars) - 2*exp.ll
str_c("The Akaike information for power law is ", exp.AIC)

exp.AIC <- 2*length(exp.inv$pars) - 2*exp.ll
AICs <- c(exp.AIC, pl.AIC)
AICscores <- tibble(AIC_exp = exp.AIC, AIC_PL = pl.AIC)
AICscores

pl.delta.AIC <- pl.AIC - min(AICs)
exp.delta.AIC <- exp.AIC - min(AICs)
pl.weight.AIC <- exp(-pl.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
exp.weight.AIC <- exp(-exp.delta.AIC/2)/(exp(-pl.delta.AIC/2)+exp(-exp.delta.AIC/2))
str_c('AIC Weight for power law = ', pl.weight.AIC)
str_c('AIC Weight for exponential = ', exp.weight.AIC)

pdf(file="inv_comp_old_v2.pdf", height = 5.3, width = 5.3)
plot(pl.inv, pch = 20, font.main = 1,col = "blue3", main = "INV St1+St2",
     xaxt = "n", yaxt = "n",
     ylab = "Probability", xlab = "Step Length (μm)")
lines(exp.inv, col = "darkslategray3", lw = 3, lty = 4)
lines(pl.inv, col = "darkslategray", lw = 3, lty = 5)
legend(0.001, 0.003, legend = c(paste("Exp, μ = ", num2str(exp.inv$pars), ", lmin = ", num2str(exp.inv$xmin), "μm"), 
                                paste("PL, μ = ", num2str(pl.inv$pars-1), ", lmin = ", num2str(pl.inv$xmin), "μm")), 
                   col = c("darkslategray3", "darkslategrey"), lty = c(4,5), lw = 3, box.lty = 0)
legend(0.001, 0.005, legend = paste("Total steps =", length(steps.inv)), box.lty = 0)
axis(1, at = c(0.01, 0.1, 1, 10, 100, 250),
     labels = c("0.01", "0.1", "1", "10", "100", "250"))
axis(2, at = c(0.001, 0.01, 0.1, 0.5), c("0.001", "0.01", "0.1", "0.5"))
dev.off()