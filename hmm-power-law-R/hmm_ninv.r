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

data_non_invader <- read.csv("raw-data/non_invader_data.csv", sep = ";")
data_non_invader <- subset(data_non_invader, select = c("PARASITE", "x..micron.", "y..micron.", "t..sec."))
names(data_non_invader) <- c("ID", "x", "y", "t")

data_non_invader2 <- read.csv("ninv_20221220.csv")
data_non_invader2 <- subset(data_non_invader2, select = c("ID", "x", "y", "t"))
names(data_non_invader2) <- c("ID", "x", "y", "t")
data_non_invader2$ID <- paste("20221220",data_non_invader2$ID, sep = "")

data_non_invader3 <- read.csv("ninv_20230116.csv")
data_non_invader3 <- subset(data_non_invader3, select = c("ID", "x", "y", "t"))
names(data_non_invader3) <- c("ID", "x", "y", "t")
data_non_invader3$ID <- paste("20230116",data_non_invader3$ID, sep = "")

data_non_invader1 <- read.csv("ninv_20230111.csv")
data_non_invader1 <- subset(data_non_invader1, select = c("TID", "x", "y", "t"))
names(data_non_invader1) <- c("ID", "x", "y", "t")
data_non_invader1$ID <- paste("20230111",data_non_invader1$ID, sep = "")

data_non_invader <- rbind(data_non_invader, data_non_invader2, data_non_invader3, data_non_invader1)

IDs <- c(1:length(unique(data_non_invader$ID)))


#fixing IDs 

j = 0
for (parasite in unique(data_non_invader$ID))
{
    j <- j+1
    data_non_invader$ID[data_non_invader$ID == parasite] <- rep(paste("NINV",IDs[j], sep = ""), each = length(data_non_invader$ID[data_non_invader$ID == parasite]))
}
data_non_invader$ID <- as.factor(data_non_invader$ID)

data_non_invader$ID <- as.factor(data_non_invader$ID)
data_non_invader$x <- as.numeric(data_non_invader$x)
data_non_invader$y <- as.numeric(data_non_invader$y)
data_non_invader$t <- as.numeric(data_non_invader$t)

data_non_invader$StepX <- NaN
data_non_invader$StepY <- NaN

for (parasite in unique(data_non_invader$ID))
{
    x <- data_non_invader$x[data_non_invader$ID == parasite]
    y <- data_non_invader$y[data_non_invader$ID == parasite]
    t <- data_non_invader$t[data_non_invader$ID == parasite]
    data_non_invader$StepX[data_non_invader$ID == parasite][2:length(data_non_invader$StepX[data_non_invader$ID == parasite])] <- x[2:length(x)]-x[1:length(x)-1]
    data_non_invader$StepY[data_non_invader$ID == parasite][2:length(data_non_invader$StepY[data_non_invader$ID == parasite])] <- y[2:length(y)]-y[1:length(y)-1]
}

data_non_invader$v <- NaN
vel <- ((data_non_invader$StepX[2:nrow(data_non_invader)]/(data_non_invader$t[2]-data_non_invader$t[1]))^2+(data_non_invader$StepY[2:nrow(data_non_invader)]/(data_non_invader$t[2]-data_non_invader$t[1]))^2)^(1/2)
data_non_invader$v[2:nrow(data_non_invader)] <- vel

data_non_invader <- na.omit(data_non_invader)

threshold_speed <- 0

# Now filtering the data entries whose speed was lower than a given threshold
v_min_non_inv = min(na.omit(data_non_invader$v))
data_non_invader <- data_non_invader[data_non_invader$v > threshold_speed, ]

data_non_invader <- data_non_invader[!is.na(data_non_invader$ID), ]

track.non.invaders <- prepData(data_non_invader, type="UTM", coordNames = c("x", "y"))

# reference parameters for each state
mu0 <- c(1, 10) # step mean (two parameters: one for each state)
sigma0 <- c(1, 5) # step SD
zeromass0 <- c(0.1, 0.1) # step zero-mass
stepPar0 <- c(mu0, sigma0)
angleMean0 <- c(pi ,0) # angle mean
kappa0 <- c(0.3, 3) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

m.non.invaders <- fitHMM(data=track.non.invaders, nbStates=2, stepPar0=stepPar0, anglePar0=anglePar0)
#m.non.invaders

#m.non.invaders
m.non.invaders.states <- viterbi(m.non.invaders)
data_non_invader$State <- as.factor(m.non.invaders.states)
#plot(m.non.invaders, plotCI = FALSE, log = "y")
#data_non_invader$angle <- track.non.invaders$angle

# Saving the csv file with the states appended to the data:

write.csv(data_non_invader, "ninv_newold_states.csv")

#creating the steps distribution

all.steps.non.invaders1 <- c()
all.steps.non.invaders2 <- c()



for (parasite in unique(data_non_invader$ID)){
peaks.x <- findPeaks(data_non_invader$x[data_non_invader$ID == parasite]) - 1
peaks.y <- findPeaks(data_non_invader$y[data_non_invader$ID == parasite]) - 1
valleys.x <- findValleys(data_non_invader$x[data_non_invader$ID == parasite]) - 1
valleys.y <- findValleys(data_non_invader$y[data_non_invader$ID == parasite]) - 1
peaks.x <- sort(c(peaks.x, valleys.x), decreasing = FALSE)
peaks.y <- sort(c(peaks.y, valleys.y), decreasing = FALSE)
step.data_non_invader.x1 <- replicate(length(peaks.x) - 1, 0)
step.data_non_invader.y1 <- replicate(length(peaks.y) - 1, 0)
step.data_non_invader.x2 <- replicate(length(peaks.x) - 1, 0)
step.data_non_invader.y2 <- replicate(length(peaks.y) - 1, 0)
StepX.ind <- data_non_invader$StepX[data_non_invader$ID == parasite]
StepY.ind <- data_non_invader$StepY[data_non_invader$ID == parasite]
Cluster.ind <- as.numeric(data_non_invader$State[data_non_invader$ID == parasite])
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
step.data_non_invader.x1[ind] = sum(StepX.ind[peaks.x[ind]:peaks.x[ind + 1]], na.rm=TRUE)

}
# Talking about state 2...
if (Cluster.ind[peaks.x[ind]] == 2)
{
step.data_non_invader.x2[ind] = sum(StepX.ind[peaks.x[ind]:peaks.x[ind + 1]], na.rm=TRUE)

}
}
}
for(ind in c(1:(length(peaks.y) - 1)))
{
if (length(unique(Cluster.ind[peaks.y[ind]:peaks.y[ind + 1]])) == 1)
{
if (Cluster.ind[peaks.y[ind]] == 1)
{
step.data_non_invader.y1[ind] = sum(StepY.ind[peaks.y[ind]:peaks.y[ind + 1]], na.rm=TRUE)

}
if (Cluster.ind[peaks.y[ind]] == 2)
{
step.data_non_invader.y2[ind] = sum(StepY.ind[peaks.y[ind]:peaks.y[ind + 1]], na.rm=TRUE)

}
}
}
combined_steps1 <- c(step.data_non_invader.x1,step.data_non_invader.y1)
combined_steps2 <- c(step.data_non_invader.x2,step.data_non_invader.y2)
combined_steps1 <- combined_steps1[combined_steps1 > 0]
combined_steps2 <- combined_steps2[combined_steps2 > 0]

all.steps.non.invaders1 <- append(all.steps.non.invaders1, combined_steps1)
all.steps.non.invaders2 <- append(all.steps.non.invaders2, combined_steps2)
}

write.csv(all.steps.non.invaders1, "all_steps_non_invaders1.csv")
write.csv(all.steps.non.invaders2, "all_steps_non_invaders2.csv")

tt.pl <- conpl$new(all.steps.non.invaders1[all.steps.non.invaders1>0])
tt.xmin <- estimate_xmin(tt.pl)
tt.pl$setXmin(tt.xmin)
data.pl <- data.frame(l = tt.pl$dat, P = tt.pl$internal$cum_n/length(tt.pl$dat))
write.csv(data.pl, "power_law_non-invaders_st1.csv")

tt.pl <- conpl$new(all.steps.non.invaders2[all.steps.non.invaders2>0])
tt.xmin <- estimate_xmin(tt.pl)
tt.pl$setXmin(tt.xmin)
data.pl <- data.frame(l = tt.pl$dat, P = tt.pl$internal$cum_n/length(tt.pl$dat))
write.csv(data.pl, "power_law_non-invaders_st2.csv")

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

pl.ninv <- conpl$new(all.steps.non.invaders[all.steps.non.invaders>0.00001])
tt.xmin <- estimate_xmin(pl.ninv)
pl.ninv$setXmin(pl.ninv)
data.pl <- data.frame(l = pl.ninv$dat, P = pl.ninv$internal$cum_n/length(pl.ninv$dat))
write.csv(data.pl, "power_law_non_invaders_st1_st2.csv")