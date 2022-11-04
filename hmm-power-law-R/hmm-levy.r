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

threshold_speed <- 0.1

# Now filtering the data entries whose speed was lower than a given threshold
data_non_invader <- data_non_invader[data_non_invader$V > threshold_speed, ]
data_invader <- data_invader[data_invader$V > threshold_speed, ]

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

track.non.invaders <- prepData(data_non_invader, type="UTM", coordNames = c("X", "Y"))
m.non.invaders <- fitHMM(data=track.non.invaders, nbStates=2, stepPar0=stepPar0, anglePar0=anglePar0)
# adding the info about the states to the dataframe now
m.non.invaders.states <- viterbi(m.non.invaders)
data_non_invader$State <- as.factor(m.non.invaders.states)
track.non.invaders$State <- as.factor(m.non.invaders.states)

# saving the results
write.csv(data_non_invader, "invader-hmm-states-results.csv")