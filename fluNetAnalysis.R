# Analysis of data for EvoComp Project
# P. Alexander Burnham
# November 25, 2019


library(dplyr)
library(ggplot2)


# read in our data for 20 reps of 3 trans values for 9 different nws
RealData <- read.csv("realNWdata1.csv",
                     header = FALSE,
                     stringsAsFactors = FALSE)

# add names
names(RealData) <- c('Fitness', 'Func.Calls', "Transcendence", 'Network.Size')

# convert to character values
RealData$Transcendence <- as.character(RealData$Transcendence)
RealData$Network.Size <- as.factor(RealData$Network.Size)

# plot figures for fittness
fr <- ggplot(RealData, aes(y = Fitness, x = Network.Size, fill = Transcendence)) + 
  geom_boxplot() +
  theme_minimal(base_size = 18) +
  labs(x="Network Size", y = "Fitness (prop. supercritical)", title = 'Fitness by Real Networks') +
  theme(legend.position=c(.2, .8)) +
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "trans.")




# read in our data 20 reps for four toy nworks for 2 trans values
toyData <- read.csv("toyNWdata.csv",
                    header = FALSE,
                    stringsAsFactors = FALSE)

# add names
names(toyData) <- c('Fitness', 'Func.Calls', 'Network' , "Transcendence")

# make numeric values character strings
toyData$Network <- as.factor(toyData$Network)
toyData$Transcendence <- as.character(toyData$Transcendence)

# change factor names to real names of nws
toyData$Network <- factor(toyData$Network, levels = c('1', '2', '3','4'),
                          labels = c('lattice', 'star', 'chain', 'Erdos-Renyi'))
# plot figures for fittness
ft <- ggplot(toyData, aes(y = Fitness, x = Network, fill = Transcendence)) +
  geom_boxplot() +
  theme_minimal(base_size = 18) +
  labs(x="Network Type", y = "Fitness (prop. supercritical)", title = 'Fitness by Toy Networks') +
  theme(legend.position=c(.8, .8)) + 
  scale_fill_manual(values = c('grey', 'steelblue'), name = "trans.")


# plot figures for calls
ct <- ggplot(toyData, aes(y = Func.Calls, x = Network, fill = Transcendence)) +
  geom_boxplot() +
  theme_minimal(base_size = 18) +
  labs(x="Network Type", y = "# Function Calls", title = 'Function Calls by Toy Networks') +
  theme(legend.position=c(.8, .8)) + 
  scale_fill_manual(values = c('grey', 'steelblue'), name = "trans.")


# plot figures for calls
cr <- ggplot(RealData, aes(y = Func.Calls, x = Network.Size, fill = Transcendence)) +
  geom_boxplot() +
  theme_minimal(base_size = 18) +
  labs(x="Network Type", y = "# Function Calls", title = 'Function Calls by Real Networks') +
  theme(legend.position=c(.8, .8)) + 
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "trans.")

library(cowplot)

cowplot::plot_grid(ft, ct, fr, cr, labels = c("A", "B", "C", "D"), align = "v")


