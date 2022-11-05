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

threshold_speed <- 0.1

# Now filtering the data entries whose speed was lower than a given threshold
v_min_non_inv = min(na.omit(data_non_invader$V))
print(v_min_non_inv)
v_min_inv = min(na.omit(data_invader$V))
print(v_min_inv)
data_non_invader <- data_non_invader[data_non_invader$V > threshold_speed, ]
data_invader <- data_invader[data_invader$V > threshold_speed, ]
#
# fixing ID issues
#
data_non_invader <- data_non_invader[!is.na(data_non_invader$ID), ]
data_invader <- data_invader[!is.na(data_invader$ID), ]

# reference parameters for each state
mu0 <- c(0.1, 0.5) # step mean (two parameters: one for each state)
sigma0 <- c(0.1, 0.5) # step SD
zeromass0 <- c(0.1, 0.05) # step zero-mass
stepPar0 <- c(mu0, sigma0)
angleMean0 <- c(pi ,0) # angle mean
kappa0 <- c(1, 1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

# analyzing the invaders file first
track.invaders <- prepData(data_invader, type="UTM", coordNames = c("X", "Y"))
m.invaders <- fitHMM(data=track.invaders, nbStates=2, stepPar0=stepPar0, anglePar0=anglePar0)
# adding the info about the states to the dataframe now
m.invaders.states <- viterbi(m.invaders)
data_invader$State <- as.factor(m.invaders.states)
track.invaders$State <- as.factor(m.invaders.states)

# saving the results
write.csv(data_invader, "invader-hmm-states-results.csv")

# analyzing non-invaders...

data_non_invader <- data_non_invader[-c(18),]

track.non.invaders <- prepData(data_non_invader, type="UTM", coordNames = c("X", "Y"))
m.non.invaders <- fitHMM(data=track.non.invaders, nbStates=2, stepPar0=stepPar0, anglePar0=anglePar0)
# adding the info about the states to the dataframe now
m.non.invaders.states <- viterbi(m.non.invaders)
data_non_invader$State <- as.factor(m.non.invaders.states)
track.non.invaders$State <- as.factor(m.non.invaders.states)

# saving the results
write.csv(data_non_invader, "non-invader-hmm-states-results.csv")