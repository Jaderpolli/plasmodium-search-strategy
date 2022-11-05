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
data_invader <- read.csv("raw-data/invader_data.csv")
data_non_invader <- read.csv("raw-data/non_invader_data.csv", sep = ";")


# corrections
data_invader <- subset(data_invader, select = c("CID", "TID", "x..micron.", 
                    "y..micron.", "v..micron.sec.", "movie", "type"))
data_non_invader <- subset(data_non_invader, select = c("CID", "TID", "x..micron.", 
                    "y..micron.", "v..micron.sec.", "movie", "type"))
#
names(data_invader) <- c("ClusterID", "TID", "X", "Y", "V", "Movie", "Type")
names(data_non_invader) <- c("ClusterID", "TID", "X", "Y", "V", "Movie", "Type")

#
# fixing factors
#
data_invader$ClusterID <- as.factor(data_invader$ClusterID)
data_invader$TID <- as.factor(data_invader$TID)
data_invader$X <- as.numeric(data_invader$X)
data_invader$Y <- as.numeric(data_invader$Y)
data_invader$Movie <- as.factor(data_invader$Movie)
data_invader$Type <- as.factor(data_invader$Type)
#
#
data_non_invader$ClusterID <- as.factor(data_non_invader$ClusterID)
data_non_invader$TID <- as.factor(data_non_invader$TID)
data_non_invader$X <- as.numeric(data_non_invader$X)
data_non_invader$Y <- as.numeric(data_non_invader$Y)
data_non_invader$Movie <- as.factor(data_non_invader$Movie)
data_non_invader$Type <- as.factor(data_non_invader$Type)

# here ID should be as factors... and we will add the cluster number in 
# front of ID number just to have unique iDs..
data_non_invader$ID <- as.factor(paste(data_non_invader$ClusterID, data_non_invader$ID, sep = ""))
data_non_invader$ID <- as.factor(paste(data_non_invader$Movie, '-', data_non_invader$ClusterID,'-',data_non_invader$TID, '-', data_non_invader$Type, sep = ""))
data_invader$ID <- as.factor(paste(data_invader$ClusterID, data_invader$ID, sep = ""))
data_invader$ID <- as.factor(paste(data_invader$Movie, '-', data_invader$ClusterID,'-',data_invader$TID, '-', data_invader$Type, sep = ""))

# extra quantities
data_invader$StepX <- NaN
vec_dx <- data_invader$X[2:nrow(data_invader)] - data_invader$X[1:nrow(data_invader) - 1]
data_invader$StepX[2:nrow(data_invader)] <- abs(vec_dx)
data_invader$StepY <- NaN
vec_dy <- data_invader$Y[2:nrow(data_invader)] - data_invader$Y[1:nrow(data_invader) - 1]
data_invader$StepY[2:nrow(data_invader)] <- abs(vec_dy)
#
data_non_invader$StepX <- NaN
vec_dx <- data_non_invader$X[2:nrow(data_non_invader)] - data_non_invader$X[1:nrow(data_non_invader) - 1]
data_non_invader$StepX[2:nrow(data_non_invader)] <- abs(vec_dx)
data_non_invader$StepY <- NaN
vec_dy <- data_non_invader$Y[2:nrow(data_non_invader)] - data_non_invader$Y[1:nrow(data_non_invader) - 1]
data_non_invader$StepY[2:nrow(data_non_invader)] <- abs(vec_dy)

# Now filtering the data entries whose speed was lower than a given threshold
#v_min_non_inv = min(na.omit(data_non_invader$V))
#print(v_min_non_inv)
#threshold_speed <- 0.1
#data_non_invader <- data_non_invader[data_non_invader$V > threshold_speed, ]
#v_min_inv = min(na.omit(data_invader$V))
#print(v_min_inv)
#threshold_speed <- 0.1
#data_invader <- data_invader[data_invader$V > threshold_speed, ]
#
# fixing ID issues
#
data_non_invader <- data_non_invader[!is.na(data_non_invader$ID), ]
data_invader <- data_invader[!is.na(data_invader$ID), ]

# Step max threshold

threshold_step <- 20

data_invader$StepX[data_invader$StepX > threshold_step ] <- NA
data_invader$StepY[data_invader$StepY > threshold_step ] <- NA
data_non_invader$StepX[data_non_invader$StepX > threshold_step ] <- NA
data_non_invader$StepY[data_non_invader$StepY > threshold_step ] <- NA

data_non_invader <- data_non_invader[-c(18),]

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

tt.pl <- conpl$new(all.steps.invaders1[all.steps.invaders1>0])
tt.xmin <- estimate_xmin(tt.pl)
tt.pl$setXmin(tt.xmin)
pdf(file="invaders1_power_law.pdf")
text.plot <- paste('Invaders state 1', 'mu1 = ', num2str(tt.pl$pars), ', xmin1= ', num2str(tt.pl$xmin), ', total_steps1 = ', num2str(length(all.steps.invaders1)))
plot(tt.pl, ylab=text.plot,
xlab='Step length (1D)')
lines(tt.pl, col=4, lwd=2)
dev.off()

all.steps.non.invaders1 <- c()
all.steps.non.invaders2 <- c()

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

tt.pl <- conpl$new(all.steps.non.invaders1[all.steps.non.invaders1>0])
tt.xmin <- estimate_xmin(tt.pl)
tt.pl$setXmin(tt.xmin)
pdf(file="non-invaders1_power_law.pdf")
text.plot <- paste('Non-Invaders state 1', 'mu1 = ', num2str(tt.pl$pars), ',xmin1 = ', num2str(tt.pl$xmin), ', total_steps1 = ', num2str(length(all.steps.non.invaders1)))
plot(tt.pl, ylab=text.plot,
xlab='Step length (1D)')
lines(tt.pl, col=4, lwd=2)
dev.off()

tt.pl <- conpl$new(all.steps.non.invaders2[all.steps.non.invaders2>0])
tt.xmin <- estimate_xmin(tt.pl)
tt.pl$setXmin(tt.xmin)
pdf(file="non-invaders2_power_law.pdf")
text.plot <- paste('Non-Invaders state 2', 'mu2 = ', num2str(tt.pl$pars), ',xmin2 = ', num2str(tt.pl$xmin), ', total_steps2 = ', num2str(length(all.steps.non.invaders2)))
plot(tt.pl, ylab=text.plot,
xlab='Step length (1D)')
lines(tt.pl, col=5, lwd=2)
dev.off()