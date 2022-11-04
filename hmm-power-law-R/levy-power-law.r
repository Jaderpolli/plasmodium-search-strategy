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

# reading the data files
data_invader <- read.csv("../raw-data/invader_data.csv")
data_non_invader <- read.csv("../raw-data/non_invader_data.csv", sep = ";")

# now the power law analysis for invaders and non-invaders

# fixing IDs
all.steps.invaders1 <- c()
all.steps.invaders2 <- c()

# calculating steps as described in suplemental material

for (parasite in unique(data_invader$ID)){
    peaks.x <- findPeaks(data_invader$X[data_invader$ID == parasite]) - 1
    peaks.y <- findPeaks(data_invader$Y[data_invader$ID == parasite]) - 1
    valleys.x <- findValleys(data_invader$X[data_invader$ID == parasite]) - 1
    valleys.y <- findValleys(data_invader$Y[data_invader$ID == parasite]) - 1
    peaks.x <- sort(c(peaks.x, valleys.x), decreasing = FALSE)
    peaks.y <- sort(c(peaks.y, valleys.y), decreasing = FALSE)
    step.data_invader.x1 <- replicate(length(peaks.x) - 1, 0)
    step.data_invader.y1 <- replicate(length(peaks.y) - 1, 0)
    step.data_invader.x2 <- replicate(length(peaks.x) - 1, 0)
    step.data_invader.y2 <- replicate(length(peaks.y) - 1, 0)
    StepX.ind <- data_invader$StepX[data_invader$ID == parasite]
    StepY.ind <- data_invader$StepY[data_invader$ID == parasite]
    Cluster.ind <- as.numeric(data_invader$ClusterID[data_invader$ID == parasite])
    
    for(ind in c(1:(length(peaks.x) - 1)))
    {
        # If the next if passes, it means that the step we are computing is either pure state 1 or
        # pure state 2
        if (length(unique(Cluster.ind[peaks.x[ind]:peaks.x[ind + 1]])) == 1)
        {
            # If we are talking about state 1...
            if (Cluster.ind[peaks.x[1]] == 1)
            {
                # Now computing the step size in case state 1...
                step.data_invader.x1[ind] = sum(StepX.ind[peaks.x[ind]:peaks.x[ind+ 1]], na.rm=TRUE)
            }
            # Talking about state 2...
            if (Cluster.ind[peaks.x[1]] == 2)
            {
                step.data_invader.x2[ind] = sum(StepX.ind[peaks.x[ind]:peaks.x[ind+ 1]], na.rm=TRUE)
            }
        }
    }
    for(ind in c(1:(length(peaks.y) - 1)))
    {
        if (length(unique(Cluster.ind[peaks.y[ind]:peaks.y[ind + 1]])) == 1)
        {
            if (Cluster.ind[peaks.y[1]] == 1)
        {
            step.data_invader.y1[ind] = sum(StepY.ind[peaks.y[ind]:peaks.y[ind+ 1]], na.rm=TRUE)
        }
        if (Cluster.ind[peaks.y[1]] == 2)
        {
            step.data_invader.y2[ind] = sum(StepY.ind[peaks.y[ind]:peaks.y[ind+ 1]], na.rm=TRUE)
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

# calculating the power law via Continuous Power Law fit and plotting:

tt.pl <- conpl$new(all.steps.invaders1[all.steps.invaders1>0])
tt.xmin <- estimate_xmin(tt.pl)
tt.pl$setXmin(tt.xmin)
text.plot <- paste('Invaders state 1', 'mu1 = ', num2str(tt.pl$pars), ', xmin1= ', num2str(tt.pl$xmin), ', total_steps1 = ', num2str(length(all.steps.invaders1)))
plot(tt.pl, ylab=text.plot,
xlab='Step length (1D)')
lines(tt.pl, col=4, lwd=2)

for (parasite in unique(data_non_invader$ID)){
    peaks.x <- findPeaks(data_non_invader$X[data_non_invader$ID == parasite]) - 1
    peaks.y <- findPeaks(data_non_invader$Y[data_non_invader$ID == parasite]) - 1
    valleys.x <- findValleys(data_non_invader$X[data_non_invader$ID == parasite]) - 1
    valleys.y <- findValleys(data_non_invader$Y[data_non_invader$ID == parasite]) - 1
    peaks.x <- sort(c(peaks.x, valleys.x), decreasing = FALSE)
    peaks.y <- sort(c(peaks.y, valleys.y), decreasing = FALSE)
    step.data_non_invader.x1 <- replicate(length(peaks.x) - 1, 0)
    step.data_non_invader.y1 <- replicate(length(peaks.y) - 1, 0)
    step.data_non_invader.x2 <- replicate(length(peaks.x) - 1, 0)
    step.data_non_invader.y2 <- replicate(length(peaks.y) - 1, 0)
    StepX.ind <- data_non_invader$StepX[data_non_invader$ID == parasite]
    StepY.ind <- data_non_invader$StepY[data_non_invader$ID == parasite]
    Cluster.ind <- as.numeric(data_non_invader$ClusterID[data_non_invader$ID ==parasite])
        for(ind in c(1:(length(peaks.x) - 1)))
    {
        # If the next if passes, it means that the step we are computing is either pure state 1 or
        # pure state 2
        if (length(unique(Cluster.ind[peaks.x[ind]:peaks.x[ind + 1]])) == 1)
        {
            # If we are talking about state 1...
            if (Cluster.ind[peaks.x[1]] == 1)
            {
                # Now computing the step size in case state 1...
                step.data_non_invader.x1[ind] = sum(StepX.ind[peaks.x[ind]:peaks.x[ind + 1]], na.rm=TRUE)

            }
            # Talking about state 2...
            if (Cluster.ind[peaks.x[1]] == 2)
            {
                step.data_non_invader.x2[ind] = sum(StepX.ind[peaks.x[ind]:peaks.x[ind + 1]], na.rm=TRUE)

            }
        }
    }
    for(ind in c(1:(length(peaks.y) - 1)))
    {
        if (length(unique(Cluster.ind[peaks.y[ind]:peaks.y[ind + 1]])) == 1)
            {
            if (Cluster.ind[peaks.y[1]] == 1)
            {
                step.data_non_invader.y1[ind] = sum(StepY.ind[peaks.y[ind]:peaks.y[ind + 1]], na.rm=TRUE)
            }
            if (Cluster.ind[peaks.y[1]] == 2)
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

# Combined stats state 1 FOR NON-INVADERSS
tt.pl <- conpl$new(all.steps.non.invaders1[all.steps.non.invaders1>0])
tt.xmin <- estimate_xmin(tt.pl)
tt.pl$setXmin(tt.xmin)
text.plot <- paste('Non-Invaders state 1', 'mu1 = ', num2str(tt.pl$pars), ',xmin1 = ', num2str(tt.pl$xmin), ', total_steps1 = ', num2str(length(all.steps.non.invaders1)))
plot(tt.pl, ylab=text.plot,
xlab='Step length (1D)')
lines(tt.pl, col=4, lwd=2)


# Combined stats state 2 FOR NON-INVADERSS
tt.pl <- conpl$new(all.steps.non.invaders2[all.steps.non.invaders2>0])
tt.xmin <- estimate_xmin(tt.pl)
tt.pl$setXmin(tt.xmin)
text.plot <- paste('Non-Invaders state 2', 'mu2 = ', num2str(tt.pl$pars), ',xmin2 = ', num2str(tt.pl$xmin), ', total_steps2 = ', num2str(length(all.steps.non.invaders2)))
plot(tt.pl, ylab=text.plot,
xlab='Step length (1D)')
lines(tt.pl, col=5, lwd=2)