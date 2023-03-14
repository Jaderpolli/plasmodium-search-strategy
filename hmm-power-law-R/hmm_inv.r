# This is the final notebook of the analysis of sporozoite data
# with objective of applying the hmm model into the data

# Main Code and idealization
# Marina E Wosniack
# Max Planck Instutute for Brain Research
# 2019-2020

# Structuring and Combining codes
# Jaderson G. Polli
# Universidade Federal do Paraná
# 2022

# loading libraries
library(moveHMM) # has the hmm package
library(ggplot2) # pretty plotting
library(quantmod) # finding peaks function
library(poweRlaw) # power law fitting
library(cowplot) # for grids
library(pracma) 
library(dplyr)

# reading the data files and correcting

data_invader <- read.csv("raw-data/invader_data.csv")
data_invader <- subset(data_invader, select = c("PARASITE", "x..micron.", "y..micron.", "t..sec."))
names(data_invader) <- c("ID", "x", "y", "t")

data_invader$ID <- as.factor(data_invader$ID)


data_invader$ID <- as.factor(data_invader$ID)
data_invader$x <- as.numeric(data_invader$x)
data_invader$y <- as.numeric(data_invader$y)
data_invader$t <- as.numeric(data_invader$t)

data_invader$StepX <- NaN
data_invader$StepY <- NaN

for (parasite in unique(data_invader$ID))
{
    x <- data_invader$x[data_invader$ID == parasite]
    y <- data_invader$y[data_invader$ID == parasite]
    t <- data_invader$t[data_invader$ID == parasite]
    data_invader$StepX[data_invader$ID == parasite][2:length(data_invader$StepX[data_invader$ID == parasite])] <- x[2:length(x)]-x[1:length(x)-1]
    data_invader$StepY[data_invader$ID == parasite][2:length(data_invader$StepY[data_invader$ID == parasite])] <- y[2:length(y)]-y[1:length(y)-1]
}

data_invader$v <- NaN
vel <- ((data_invader$StepX[2:nrow(data_invader)]/(data_invader$t[2]-data_invader$t[1]))^2+(data_invader$StepY[2:nrow(data_invader)]/(data_invader$t[2]-data_invader$t[1]))^2)^(1/2)
data_invader$v[2:nrow(data_invader)] <- vel

data_invader <- na.omit(data_invader)


threshold_speed <- 0.0
#maxthreshold_speed <- 20

# Now filtering the data entries whose speed was lower than a given threshold
v_min_inv = min(na.omit(data_invader$v))
#print(v_min_inv)
data_invader <- data_invader[data_invader$v > threshold_speed, ]
#data_invader <- data_invader[data_invader$v < maxthreshold_speed, ]

data_invader <- data_invader[!is.na(data_invader$ID), ]

track.invaders <- prepData(data_invader, type="UTM", coordNames = c("x", "y"))

# reference parameters for each state
mu0 <- c(1, 10) # step mean (two parameters: one for each state)
sigma0 <- c(1, 5) # step SD
zeromass0 <- c(0, 0) # step zero-mass
stepPar0 <- c(mu0, sigma0)
angleMean0 <- c(pi ,0) # angle mean
kappa0 <- c(0.3, 3) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

m.invaders <- fitHMM(data=track.invaders, nbStates=2, stepPar0=stepPar0, anglePar0=anglePar0)
#m.invaders

#m.invaders
m.invaders.states <- viterbi(m.invaders)
data_invader$State <- as.factor(m.invaders.states)
#plot(m.invaders, plotCI = FALSE, log = "y")
#data_invader = track.invaders$angle
write.csv(data_invader, "inv_old_states.csv")

all.steps.invaders1 <- c()
all.steps.invaders2 <- c()

for (parasite in unique(data_invader$ID)){
peaks.x <- findPeaks(data_invader$x[data_invader$ID == parasite]) - 1
peaks.y <- findPeaks(data_invader$y[data_invader$ID == parasite]) - 1
valleys.x <- findValleys(data_invader$x[data_invader$ID == parasite]) - 1
valleys.y <- findValleys(data_invader$y[data_invader$ID == parasite]) - 1
peaks.x <- sort(c(peaks.x, valleys.x), decreasing = FALSE)
peaks.y <- sort(c(peaks.y, valleys.y), decreasing = FALSE)
step.data_invader.x1 <- replicate(length(peaks.x) - 1, 0)
step.data_invader.y1 <- replicate(length(peaks.y) - 1, 0)
step.data_invader.x2 <- replicate(length(peaks.x) - 1, 0)
step.data_invader.y2 <- replicate(length(peaks.y) - 1, 0)
StepX.ind <- data_invader$StepX[data_invader$ID == parasite]
StepY.ind <- data_invader$StepY[data_invader$ID == parasite]
Cluster.ind <- as.numeric(data_invader$State[data_invader$ID == parasite])
for(ind in c(1:(length(peaks.x) - 1)))
{
# If the next if passes, it means that the step we are computing is either pure state 1 or
# pure state 2
if (length(unique(Cluster.ind[peaks.x[ind]:peaks.x[ind + 1]])) == 1)
{
# If we are talking about state 1...
if (Cluster.ind[peaks.x[ind]] == 1)
{
# Now computing the step size in case state 1...
step.data_invader.x1[ind] = sum(StepX.ind[peaks.x[ind]:peaks.x[ind + 1]], na.rm=TRUE)

}
# Talking about state 2...
if (Cluster.ind[peaks.x[ind]] == 2)
{
step.data_invader.x2[ind] = sum(StepX.ind[peaks.x[ind]:peaks.x[ind + 1]], na.rm=TRUE)

}
}
}
for(ind in c(1:(length(peaks.y) - 1)))
{
if (length(unique(Cluster.ind[peaks.y[ind]:peaks.y[ind + 1]])) == 1)
{
if (Cluster.ind[peaks.y[ind]] == 1)
{
step.data_invader.y1[ind] = sum(StepY.ind[peaks.y[ind]:peaks.y[ind + 1]], na.rm=TRUE)

}
if (Cluster.ind[peaks.y[ind]] == 2)
{
step.data_invader.y2[ind] = sum(StepY.ind[peaks.y[ind]:peaks.y[ind + 1]], na.rm=TRUE)

}
}
}
combined_steps1 <- c(step.data_invader.x1,step.data_invader.y1)
combined_steps2 <- c(step.data_invader.x2,step.data_invader.y2)
combined_steps1 <- combined_steps1[combined_steps1 > 0]
combined_steps2 <- combined_steps2[combined_steps2 > 0]

all.steps.invaders1 <- append(all.steps.invaders1, combined_steps1)
all.steps.invaders2 <- append(all.steps.invaders2, combined_steps2)
}

plot(all.steps.invaders1)
plot(all.steps.invaders2)

write.csv(all.steps.invaders1, "all_steps_invaders1_old.csv")
write.csv(all.steps.invaders2, "all_steps_invaders2_old.csv")

tt.pl <- conpl$new(all.steps.invaders1[all.steps.invaders1>0])
tt.xmin <- estimate_xmin(tt.pl)
tt.pl$setXmin(tt.xmin)
data.pl <- data.frame(l = tt.pl$dat, P = tt.pl$internal$cum_n/length(tt.pl$dat))
write.csv(data.pl, "invaders1_power_law_old.csv")
results.pl <- data.frame(mu2 = num2str(tt.pl$pars-1), xmin2 = num2str(tt.pl$xmin), total_steps2 = num2str(length(all.steps.invaders1)))
write.csv(results.pl, "invaders1_power_law_results_old.csv")

tt.pl <- conpl$new(all.steps.invaders2[all.steps.invaders2>0])
tt.xmin <- estimate_xmin(tt.pl)
tt.pl$setXmin(tt.xmin)
data.pl <- data.frame(l = tt.pl$dat, P = tt.pl$internal$cum_n/length(tt.pl$dat))
write.csv(data.pl, "invaders2_power_law_old.csv")
results.pl <- data.frame(mu2 = num2str(tt.pl$pars-1), xmin2 = num2str(tt.pl$xmin), total_steps2 = num2str(length(all.steps.invaders2)))
write.csv(results.pl, "invaders2_power_law_results_old.csv")
dev.off()

all.steps.invaders <- c(all.steps.invaders1, all.steps.invaders2)
length(all.steps.invaders)

tt.pl <- conpl$new(all.steps.invaders)
tt.xmin <- estimate_xmin(tt.pl)
tt.pl$setXmin(tt.xmin)
data.pl <- data.frame(l = tt.pl$dat, P = tt.pl$internal$cum_n/length(tt.pl$dat))
write.csv(data.pl, "invaders_power_law_old.csv")
results.pl <- data.frame(mu2 = num2str(tt.pl$pars-1), xmin2 = num2str(tt.pl$xmin), total_steps2 = num2str(length(all.steps.invaders)))
write.csv(results.pl, "invaders_power_law_results_old.csv")