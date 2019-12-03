# Analysis of data for EvoComp Project
# P. Alexander Burnham
# November 25, 2019


library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)


# read in our data for 20 reps of 3 trans values for 9 different nws
RealData <- read.csv("realNWdata2.csv",
                     header = FALSE,
                     stringsAsFactors = FALSE)

# add names
names(RealData) <- c('Fitness', 'Func.Calls', "Transcendence", 'Network.Size')

# convert to character values
RealData$Transcendence <- as.character(RealData$Transcendence)
RealData$Network.Size <- as.factor(RealData$Network.Size)

RealData$Log_fitness <- -log10(RealData$Fitness + 10^-5) # log10(RealData$Fitness+1)



# plot figures for fittness
fr <- ggplot(RealData, aes(y = Fitness, x = Network.Size, fill = Transcendence)) + 
  geom_boxplot() +
  theme_minimal(base_size = 18) +
  labs(x="Network Size", y = "prop. supercritical") +
  theme(legend.position = 'top',axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "Trans:") + 
  coord_cartesian(ylim = c(0,.015)) #+ scale_y_log10()




# read in our data 20 reps for four toy nworks for 2 trans values
toyData <- read.csv("toyNWdataFinal.csv",
                    header = FALSE,
                    stringsAsFactors = FALSE)

# add names
names(toyData) <- c('Fitness', 'Func.Calls', 'Network' , "Transcendence")

# make numeric values character strings
toyData$Network <- as.factor(toyData$Network)
toyData$Transcendence <- as.character(toyData$Transcendence)

# change factor names to real names of nws
toyData$Network <- factor(toyData$Network, levels = c('1', '2', '3','4'),
                          labels = c('lattice', 'star', 'chain', 'E-R'))
# plot figures for fittness
ft <- ggplot(toyData, aes(y = log10(Fitness), x = Network, fill = Transcendence)) +
  geom_boxplot() +
  theme_minimal(base_size = 18) +
  labs(x="Network Type", y = "prop. supercritical") +
  theme(legend.position = 'top',axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "Trans:") + coord_cartesian(ylim = c(0,.4))


# plot figures for calls
ct <- ggplot(toyData, aes(y = log10(Func.Calls), x = Network, fill = Transcendence)) +
  geom_boxplot() +
  theme_minimal(base_size = 18) +
  labs(x="Network Type", y = "log10(func. calls)") +
  theme(legend.position = 'top',axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "Trans:") + coord_cartesian(ylim = c(2,4.5))

# plot figures for calls
cr <- ggplot(RealData, aes(y = log10(Func.Calls), x = Network.Size, fill = Transcendence)) +
  geom_boxplot() +
  theme_minimal(base_size = 18) +
  labs(x="Network Size", y = "log10(func. calls)") +
  theme(legend.position = 'top',axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "Trans:") + coord_cartesian(ylim = c(2,4.5))

# plot the four plots together
cowplot::plot_grid(ft,fr,ct,cr, labels = c("A", "B", "C", "D"), align = "v")



# vector between ranges of real net sizes
N <- 20:1430

# vector of function calls required to brut force networks with 1 to 4 vaccines
kn <- c(choose(k=1, n=N), choose(k=2, n=N), choose(k=3, n=N), choose(k=4, n=N))

# rep the number of N vals for plotting
Nvec <- rep(N, 4)

# rep K values
Kvec <- as.character(c(rep(1, length(N)),rep(2, length(N)),rep(3, length(N)),rep(4, length(N))))

# create a df
chooseDF <- data.frame(kn, Nvec, as.character(Kvec))

library(plyr)
# summarise values for actual function calls
DF3 <- ddply(RealData, c("Network.Size"), summarise, 
             n = length(Func.Calls),
             mean = mean(Func.Calls, na.rm=TRUE),
             sd = sd(Func.Calls, na.rm=TRUE),
             se = sd / sqrt(n),
             upper = mean+sd,
             lower = mean-sd)


# work on function call plot
ggplot(chooseDF, aes(x=Nvec, y = kn, group = Kvec, color=Kvec)) + 
  geom_line(size=1.5) + theme_classic(base_size = 18) +
  labs(y='log10(function calls)', x='network size', color='# Vaccines:') +
  theme(legend.position = 'top') +
  annotate("point", x = as.numeric(as.character(DF3$Network.Size)), 
           y = as.numeric(DF3$mean), 
           colour = "blue") +
  annotate('segment', x = as.numeric(as.character(DF3$Network.Size)), 
           y = DF3$lower, xend = as.numeric(as.character(DF3$Network.Size)), 
                   yend = DF3$upper) + scale_y_log10()


x=split(RealData, RealData$Transcendence)
temp <- x$`1`[x$`1`$Network.Size=='1430', ]
temp$Fitness
