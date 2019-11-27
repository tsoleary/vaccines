# Analysis of data for EvoComp Project
# P. Alexander Burnham
# November 25, 2019



# read in our data (using partial data right now (only three reps)) - Alex
RealData <- read.csv("realNWdata.csv", 
                              header = FALSE,
                              stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)

names(RealData) <- c('Fitness', 'Func.Calls', "Transcendence", 'Network.Size')

RealData$Transcendence <- as.character(RealData$Transcendence)
RealData$Network.Size <- as.character(RealData$Network.Size)

ggplot(RealData, aes(y = Fitness, x = Network.Size, fill = Transcendence)) +
  geom_boxplot(alpha=0.7) + theme_bw()




# read in our data (using partial data right now (only three reps)) - Alex
toyData <- read.csv("toyNWdata.csv", 
                     header = FALSE,
                     stringsAsFactors = FALSE)

names(toyData) <- c('Fitness', 'Func.Calls', 'Network' , "Transcendence")

toyData$Network <- as.factor(toyData$Network)
toyData$Transcendence <- as.character(toyData$Transcendence)

toyData$Network <- factor(toyData$Network, levels = c('1', '2', '3','4'), 
       labels = c('lattice', 'star', 'chain', 'Erdos-Renyi'))

ggplot(toyData, aes(y = Fitness, x = Network, fill = Transcendence)) +
  geom_boxplot(alpha=0.7) + theme_bw()


