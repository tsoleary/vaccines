# Analysis of data for EvoComp Project
# P. Alexander Burnham
# November 25, 2019


library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)
library(latex2exp)
library(plyr)
library(ggforce)
library(matrixStats)
library(tidyr)




##############################################################################################
# read in random data for toy networks
randToyMat <- read.csv("randToyMat.csv",
                       header = FALSE,
                       stringsAsFactors = FALSE)
# make a matrix
randToyMat <- as.matrix(randToyMat)

# summary as dataframe
randTDF <- data.frame(randToyMat[,c(1:2)],
                      rowMeans(randToyMat[,c(3:1002)]),
                      rowSds(randToyMat[,c(3:1002)]))

# add names to dataframe
names(randTDF) <- c('nwTypes', 'Trans', "mean", 'sd')
##############################################################################################






##############################################################################################
# read in random data for toy networks
randRealMat3 <- read.csv("randRealMat3.csv",
                         header = FALSE,
                         stringsAsFactors = FALSE)
# make a matrix
randRealMat3 <- as.matrix(randRealMat3)

# summary as dataframe
randRDF3 <- data.frame(randRealMat3[,c(1:2)],
                       rowMeans(randRealMat3[,c(3:1002)]),
                       rowSds(randRealMat3[,c(3:1002)]))

# add names to dataframe
names(randRDF3) <- c('nwTypes', 'Trans', "mean", 'sd')

randRDF3<-randRDF3[!(randRDF3$nwTypes %in% c(7, 8, 9)),]
randRDF3$nwTypes <- rep(6:1,3)
##############################################################################################






##############################################################################################
# read in random data for toy networks
randRealMat4 <- read.csv("randRealMat4.csv",
                         header = FALSE,
                         stringsAsFactors = FALSE)
# make a matrix
randRealMat4 <- as.matrix(randRealMat4)

# summary as dataframe
randRDF4 <- data.frame(randRealMat4[,c(1:2)],
                       rowMeans(randRealMat4[,c(3:1002)]),
                       rowSds(randRealMat4[,c(3:1002)]))

# add names to dataframe
names(randRDF4) <- c('nwTypes', 'Trans', "mean", 'sd')

randRDF4<-randRDF4[!(randRDF4$nwTypes %in% c(7, 8, 9)),]
randRDF4$nwTypes <- rep(6:1,3)
##############################################################################################





# read in our data for 20 reps of 3 trans values for 9 different nws
RealData <- read.csv("realNWdata4vac.csv",
                     header = FALSE,
                     stringsAsFactors = FALSE)

# add names
names(RealData) <- c('Fitness', 'Func.Calls', "Transcendence", 'Network.Size')

# convert to character values
RealData$Transcendence <- as.character(RealData$Transcendence*2)
RealData$Network.Size <- as.factor(RealData$Network.Size)

RealData$Log_fitness <- -log10(RealData$Fitness + 10^-5) # log10(RealData$Fitness+1)
RealDataLarge<-RealData[!(RealData$Network.Size %in% c(20, 43, 52)),]

split4<-split(randRDF4, randRDF4$Trans)

# plot figures for fittness
fr <- ggplot(RealDataLarge, aes(y = Fitness, x = Network.Size, fill = Transcendence)) +
  geom_boxplot() +
  theme_light(base_size = 20) +
  labs(x=NULL, y = "prop. supercritical") +
  theme(legend.position = 'none',axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "Trans:") +
  annotate("point", x = split4$`2`$nwTypes, y =split4$`2`$mean, colour = "steelblue", size=3, shape=18) +
  annotate('segment', x = split4$`2`$nwTypes, y = split4$`2`$mean-split4$`2`$sd,
           xend = split4$`2`$nwTypes, yend = split4$`2`$mean+split4$`2`$sd, size=.5, colour = "steelblue") +
  annotate("point", x = split4$`1`$nwTypes-.25, y =split4$`1`$mean, colour = "grey", size=3, shape=18) +
  annotate('segment', x = split4$`1`$nwTypes-.25, y = split4$`1`$mean-split4$`1`$sd,
           xend = split4$`1`$nwTypes-.25, yend = split4$`1`$mean+split4$`1`$sd, size=.5, colour = "grey") +
  annotate("point", x = split4$`3`$nwTypes+.25, y =split4$`3`$mean, colour = "black", size=3, shape=18) +
  annotate('segment', x = split4$`3`$nwTypes+.25, y = split4$`3`$mean-split4$`3`$sd,
           xend = split4$`3`$nwTypes+.25, yend = split4$`3`$mean+split4$`3`$sd, size=.5, colour = "black") +facet_zoom(ylim = c(0, .015), zoom.size=2, show.area=F)

fr


# read in our data 20 reps for four toy nworks for 2 trans values
toyData <- read.csv("toyNWdataFinal2.csv",
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


toysplit <- split(randTDF, randTDF$Trans)


# plot figures for fittness
ft <- ggplot(toyData, aes(y = Fitness, x = Network, fill = Transcendence)) +
  geom_boxplot() +
  theme_light(base_size = 20) +
  labs(x=NULL, y = "prop. supercritical") +
  theme(legend.position = 'none',axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "Trans:") +
  annotate("point", x = toysplit$`2.5`$nwTypes, y =toysplit$`2.5`$mean, colour = "steelblue", size=3, shape=18) +
  annotate('segment', x = toysplit$`2.5`$nwTypes, y = toysplit$`2.5`$mean-toysplit$`2.5`$sd,
           xend = toysplit$`2.5`$nwTypes, yend = toysplit$`2.5`$mean+toysplit$`2.5`$sd, size=.5, colour = "steelblue") +
  annotate("point", x = toysplit$`2`$nwTypes-.25, y =toysplit$`2`$mean, colour = "grey", size=3, shape=18) +
  annotate('segment', x = toysplit$`2`$nwTypes-.25, y = toysplit$`2`$mean-toysplit$`2`$sd,
           xend = toysplit$`2`$nwTypes-.25, yend = toysplit$`2`$mean+toysplit$`2`$sd, size=.5, colour = "grey") +
  annotate("point", x = toysplit$`3`$nwTypes+.25, y =toysplit$`3`$mean, colour = "black", size=3, shape=18) +
  annotate('segment', x = toysplit$`3`$nwTypes+.25, y = toysplit$`3`$mean-toysplit$`3`$sd,
           xend = toysplit$`3`$nwTypes+.25, yend = toysplit$`3`$mean+toysplit$`3`$sd, size=.5, colour = "black")+
  facet_zoom(ylim = c(0, .075), zoom.size=2, show.area=F)

ft




















# read in our data for 20 reps of 3 trans values for 9 different nws
RealData3Vac <- read.csv("realNWdata3vac.csv",
                         header = FALSE,
                         stringsAsFactors = FALSE)

# add names
names(RealData3Vac) <- c('Fitness', 'Func.Calls', "Transcendence", 'Network.Size')

# convert to character values
RealData3Vac$Transcendence <- as.character(RealData3Vac$Transcendence)
RealData3Vac$Network.Size <- as.factor(RealData3Vac$Network.Size)



split3<-split(randRDF3, randRDF3$Trans)

# plot figures for fittness
fr3 <- ggplot(RealData3Vac, aes(y = Fitness, x = Network.Size, fill = Transcendence)) +
  geom_boxplot() +
  theme_light(base_size = 20) +
  labs(x=NULL, y = "prop. supercritical") +
  theme(legend.position = 'none',axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "Trans:") +
  annotate("point", x = split3$`2`$nwTypes, y =split3$`2`$mean, colour = "steelblue", size=3, shape=20) +
  annotate('segment', x = split3$`2`$nwTypes, y = split3$`2`$mean-split3$`2`$sd,
           xend = split3$`2`$nwTypes, yend = split3$`2`$mean+split3$`2`$sd, size=.5, colour = "steelblue") +
  annotate("point", x = split3$`1`$nwTypes-.25, y =split3$`1`$mean, colour = "grey", size=3, shape=18) +
  annotate('segment', x = split3$`1`$nwTypes-.25, y = split3$`1`$mean-split3$`1`$sd,
           xend = split3$`1`$nwTypes-.25, yend = split3$`1`$mean+split3$`1`$sd, size=.5, colour = "grey") +
  annotate("point", x = split3$`3`$nwTypes+.25, y =split3$`3`$mean, colour = "black", size=3, shape=18) +
  annotate('segment', x = split3$`3`$nwTypes+.25, y = split3$`3`$mean-split3$`3`$sd,
           xend = split3$`3`$nwTypes+.25, yend = split3$`3`$mean+split3$`3`$sd, size=.5, colour = "black") +facet_zoom(ylim = c(0, .0083), zoom.size=2, show.area=F)

fr3










# plot figures for calls
ct <- ggplot(toyData, aes(y = log10(Func.Calls), x = Network, fill = Transcendence)) +
  geom_boxplot() +
  theme_light(base_size = 20) +
  labs(x=NULL, y = "log10(func. calls)") +
  theme(legend.position = c(.74,.74),axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "Trans:") +
  coord_cartesian(ylim = c(2.25,6))
ct

# plot figures for calls
cr <- ggplot(RealDataLarge, aes(y = log10(Func.Calls), x = Network.Size, fill = Transcendence)) +
  geom_boxplot() +
  theme_light(base_size = 20) +
  labs(x=NULL, y = "log10(func. calls)") +
  theme(legend.position = c(.74,.74),axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "Trans:") +
  coord_cartesian(ylim = c(2.25,6))
cr



# plot figures for calls
cr3 <- ggplot(RealData3Vac, aes(y = log10(Func.Calls), x = Network.Size, fill = Transcendence)) +
  geom_boxplot() +
  theme_light(base_size = 20) +
  labs(x=NULL, y = "log10(func. calls)") +
  theme(legend.position = c(.74,.74),axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c('grey', 'steelblue', 'black'), name = "Trans:") +
  coord_cartesian(ylim = c(2.25,6))
cr3






# plot the four plots together
cowplot::plot_grid(ct,ft, cr3, fr3, cr,fr, labels = c("A", "B", "C","D","E","F"),
                   align = "hv", ncol=2, rel_widths = c(1,2))





# vector between ranges of real net sizes
N <- 20:1500

# vector of function calls required to brut force networks with 1 to 4 vaccines
kn <- c(choose(k=1, n=N), choose(k=2, n=N), choose(k=3, n=N), choose(k=4, n=N))

# rep the number of N vals for plotting
Nvec <- rep(N, 4)

# rep K values
Kvec <- as.character(c(rep(1, length(N)),rep(2, length(N)),rep(3, length(N)),rep(4, length(N))))

# create a df
chooseDF <- data.frame(kn, Nvec, as.character(Kvec))

# summarise values for actual function calls
DF3 <- ddply(RealData, c("Network.Size"), summarise,
             n = length(Func.Calls),
             mean = log10(mean(Func.Calls, na.rm=TRUE)),
             sd = log10(sd(Func.Calls, na.rm=TRUE)),
             se = sd / sqrt(n),
             upper = mean+se,
             lower = mean-se)

# change -INF values to 0
DF3[1,6:7] <- 0

# work on function call plot
ggplot(chooseDF, aes(x=Nvec, y = log10(kn), group = Kvec, color=Kvec)) +
  geom_line(size=1.8) + theme_classic(base_size = 20) +
  labs(y='func. calls (log10(N choose k))', x='network size (N)', color='# Vaccines (k):') +
  theme(legend.position = 'top') +
  annotate("point", x = as.numeric(as.character(DF3$Network.Size)),
           y = as.numeric(DF3$mean),
           colour = "black", size=3) +
  annotate('segment', x = as.numeric(as.character(DF3$Network.Size)),
           y = DF3$lower, xend = as.numeric(as.character(DF3$Network.Size)),
           yend = DF3$upper, size=1) + coord_cartesian(xlim = c(0, 1500), ylim = c(0,12))


# growing network analysis, ga v. rand -----------------------------------------

growData <- read.csv("growMat.csv",
                     header = FALSE,
                     stringsAsFactors = FALSE)

colnames(growData) <- c("days", rep.int("Random", 1800), rep.int("GA", 10), "Max")

growData <- pivot_longer(growData,
                         -contains("days"),
                         names_to = "group",
                         values_to = "fitness") %>%
  filter(days != 270 & days != 450)

growData$days <- as.factor(growData$days)

MaxData <- growData[growData$group=="Max",]
NoMaxData <- growData[!growData$group=="Max",]

gr <- ggplot(NoMaxData, aes(y = fitness, x = days, fill = group)) +
  geom_boxplot() +
  theme_minimal(base_size = 18) +
  labs(x = "Days after vaccine selection",
       y = "Fitness") +
  theme(legend.position = 'top',axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c('grey', 'steelblue'), name = "") +
  annotate("point", x = MaxData$days, y =MaxData$fitness, colour = "red", size=3, shape=18)
gr

# stats for time mod
timeMod <- aov(data=NoMaxData, fitness~group*days)
summary(timeMod)

# growing network analysis, ga v. rand -----------------------------------------

brutData <- read.csv("brutForceGAfig.csv",
                     header = FALSE,
                     stringsAsFactors = FALSE)

colnames(brutData) <- c("fitness", "func_calls", "trans", "net")

brutData$net <- as.factor(brutData$net)
brutData$trans <- as.factor(brutData$trans)

br <- ggplot(brutData, aes(y = fitness, x = net, fill = trans)) +
  geom_boxplot() +
  theme_minimal(base_size = 18) +
  labs(x = "Days after vaccine selection",
       y = "Fitness") +
  theme(legend.position = 'top',axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c('grey', 'black', 'steelblue'), name = "")
br
