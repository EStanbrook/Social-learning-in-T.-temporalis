##### Code for zero-inflated models, GLMM models, and figures 1 and 3 #####

##### Packages #####
# Load necessary library packages

library(MASS)
library(lme4)
library(pscl)
library(growthcurver)
library(ggplot2)
library(ggpubr)
library(lmerTest)
library(car)

##### zero-inflated model: factors affecting learning #####
# model only applied to groups 1-8 i.e. groups
  # that had undergone the full 25 trials

feed <- read.csv("ZeroInflated.csv")


zeromod <- zeroinfl(Boxes.Opened ~ 
                   Cumulative.Feeds+Cumulative.Observed.Feeds
                 + Feeds.Observed+Opens.observed|Cumulative.Observed.Opens+
                   Cumulative.Observed.Feeds+
                   Opens.observed, data=feed, dist = "poisson")
summary(zeromod)


##### GLMMs: the effect of dominance on learning #####

all1=read.csv("NodeMetrics.csv")

# Ag/As/Af = agonistic/ associative/ affiliative
# Dom = dominance (David's score)
# Qav = average number of quartiles entered by fish (measure of movement)

##### Latency to first feed #####
latfeed_net<- glmer(log10(LatFeedS+10)~ 
                      AsCloseness+AfCloseness+
                      Dom+ (1|Group), data =  all1)
summary(latfeed_net)
Anova(latfeed_net)

##### Number of boxes fed from #####
feed1 <- glmer(sqrt(NoFeeds)~AsCloseness+
                 AfCloseness+Dom+(1|Group),
               data = all1)
summary(feed1)
Anova(feed1)

##### Latency to open first box #####
open1 <- lmer(log10(LatOpenS+10)~
                AsCloseness+AfCloseness+
                Dom+(1|Group),
              data = all1)
summary(open1)
Anova(open1)

##### Number of boxes opened #####

open2 <- lmer(log10(NoOpens+10)~
                AsCloseness+AfCloseness+
                Dom+(1|Group),
              data = all1)
summary(open2)
Anova(open2)


##### Figure 1 #####

fig1 <- read.csv("Fig1.csv")

ggplot(fig1, aes(x=Fish_Day, y = Proportion_1_3)) +  
  ylab("Proportion of fish that opened a box") +
  xlab("Fish Days")+
  geom_point(size=2,colour="#344d90",
             position = position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=Proportion_1_3-
                      Proportion_1_3_SE, 
                    ymax=Proportion_1_3+
                      Proportion_1_3_SE),width=0.25,
                colour="#344d90") +
  xlim(0,10)+ylim(0,1)+
  scale_x_continuous(breaks=1:10)+
  geom_smooth(se = FALSE, method = "gam", formula = y ~ x + I(x^2),
              colour="#344d90",size=0.3)+
  theme_bw()

ggplot(fig1, aes(x=Fish_Day, 
                 y = Cumulative_Proportion_1_3)) +  
  ylab("Proportion of fish that opened a box") +
  xlab("Fish Days")+
  geom_point(size=2,colour="#5cc5ef",
             position = position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=Cumulative_Proportion_1_3-
                      Cumulative_Proportion_1_3_SE, 
                    ymax=Cumulative_Proportion_1_3+
                      Cumulative_Proportion_1_3_SE),width=0.25,
                colour="#5cc5ef") +
  xlim(0,10)+ylim(0,1)+
  scale_x_continuous(breaks=1:10)+
  geom_smooth(se = FALSE, method = "gam", 
              formula = y ~ x + I(x^2),
              colour="#5cc5ef",size=0.3)+
  theme_bw()


ggplot(fig1, aes(x=Fish_Day, y = Proportion_4_10)) +  
  ylab("Proportion of fish that opened a box") +
  xlab("Fish Days")+
  geom_point(size=2,colour="#e7552c",
             position = position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=Proportion_4_10-
                      Proportion_4_10_SE, 
                    ymax=Proportion_4_10+
                      Proportion_4_10_SE),width=0.25,
                colour="#e7552c") +
  xlim(0,10)+ylim(0,1)+
  scale_x_continuous(breaks=1:10)+
  geom_smooth(se = FALSE, method = "gam", 
              formula = y ~ s(log(x)),
              colour="#e7552c",size=0.3)+
  theme_bw()

ggplot(fig1, aes(x=Fish_Day, 
                 y = Cumulative_Proportion_4_10)) +  
  ylab("Proportion of fish that opened a box") +
  xlab("Fish Days")+
  geom_point(size=2,colour="#ffb745",
             position = position_dodge(), stat="identity") + 
  geom_errorbar(aes(ymin=Cumulative_Proportion_4_10-
                      Cumulative_Proportion_4_10_SE, 
                    ymax=Cumulative_Proportion_4_10+
                      Cumulative_Proportion_4_10_SE),width=0.25,
                colour="#ffb745") +
  xlim(0,10)+ylim(0,1)+
  scale_x_continuous(breaks=1:10)+
  geom_smooth(se = FALSE, method = "gam", 
              formula = y ~ s(log(x)),
              colour="#ffb745",size=0.3)+
  theme_bw()



##### Figure 3 #####

fig3=read.csv("Fig3.csv")

ggplot(fig3,aes(x=DominanceScore,y=LatencyToFirstOpenS))+
  geom_point(size=2,colour="#344d90")+
  xlab('Dominance score')+
  ylab('Latency to open a box (s)')+
  geom_smooth(method='lm',se=FALSE,colour="#344d90")+
  theme_bw(base_size=13)


ggplot(fig3,aes(x=DominanceScore,y=LogNumberOfBoxesOpened))+
  geom_point(size=2,colour="#344d90",shape=17)+
  xlab('Dominance score')+
  ylab('Log number of boxes opened/ fed from')+
  ylim(0,1.6)+
  geom_smooth(method='lm',se=FALSE,colour="#344d90")+
  theme_bw(base_size=13)

ggplot(fig3,aes(x=DominanceScore,y=LogNumberOfBoxesFedFrom))+
  geom_point(size=2,shape=15,colour="#5cc5ef")+
  xlab('Dominance score')+
  ylab('Log number of boxes fed from')+
  ylim(0,1.6)+
  geom_smooth(method='lm',se=FALSE,colour="#5cc5ef")+
  theme_bw(base_size=13)















