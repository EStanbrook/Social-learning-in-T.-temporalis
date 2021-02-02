##### Code for network-based diffusion analysis
    # (NBDA) #####
# 168 NBDAs: one NBDA including and one excluding
# asocial variables (dominance and propensity 
# to feed/open) for each of the 6 social networks
# (symmetrical and directed for associative, 
# affiliative and agonistic) for each of the
# 7 groups (groups 4-10)
      # AICc value differences >2 determine 
      # which of the social or asocial transmission
      # models best explains the order of 
      # acquisition (smaller AIC=better model)
    
#####NB: BEFORE RUNNING NBDAs #####
    # Before running the NBDAs you must run the 
    # code provided (Hoppitt NBDA Code 1.2.13)
    # available at 
    # https://lalandlab.st-andrews.ac.uk/freeware/

##### Packages #####
# Load necessary library packages, 

library(statnet)
library(igraph)
library(asnipe)
library(tnet)
library(assortnet)
library(remotes)
library(NBDA)




##### Input orders of acquisition #####
    # of feeding and opening behaviour as a 
    # vector for subsequent NBDAs, where numbers 
    # indicate the row/column representing that 
    # individual in the network matrix 

OOAFeed4<-c(1,3,4,2,5)
OOAOpen4<-c(1,3,5,4,2)

OOAFeed5<-c(1,2,3,4,5)
OOAOpen5<-c(1,4,5,3,2)

OOAFeed6<-c(2,5,4,1,3)
OOAOpen6<-c(2,3,4,5,1)

OOAFeed7<-c(1,2,5,4,3)
OOAOpen7<-c(1,2,3,5)

OOAFeed8<-c(1,2,3,4)
OOAOpen8<-c(1,2,4)

OOAFeed9<-c(3,5,2,1,6,7,4)
OOAOpen9<-c(3,1,6,7)

OOAFeed10<-c(2,3,8,1,6)
OOAOpen10<-c(2,3,6)




# Generate matrix of individual-level [asocial] variables for each group, including propensity to feed [mean-centred, normalised feeds], to open, and dominance rank




PropFeed4<-c(10.25,1,7.5,1,1.25)
PropOpen4<-c(36,3,14,1,1)
dom4<-c(0.299864007,0,1,0.546237536,0.805757029)

PropFeed5<-c(1.85,1.25,1,1.4,1.3)
PropOpen5<-c(7.5,1.5,1,3.25,2.75)
dom5<-c(0.332092499,0,1,0.548770012,0.860254241)

PropFeed6<-c(1.125,3.5,1,4.875,5)
PropOpen6<-c(1,10.5,1.5,11,13.5)
dom6<-c(0,0.641118,0.450927,0.641118,1)

PropFeed7<-c(4.5,9,5,1,4.75)
PropOpen7<-c(1.75,3.875,1.125,0,1)
dom7<-c(0,0.387507378,0.883806974,0.464031974,1)

PropFeed8<-c(4,5.75,1,2.75)
PropOpen8<-c(3.33,5.83,0,1)
dom8<-c(0,0.931507,0.534247,1)

PropFeed9<-c(5.33,2,6,1,2.33,5.33,4)
PropOpen9<-c(1,0,3,0,0,1.67,1)
dom9<-c(0.49084997,0.180601,0.482329,0.288331,0,0.696928,1)

PropFeed10<-c(1,16,7,0,0,5,0,4)
PropOpen10<-c(0,9,4,0,0,1,0,0)
dom10<-c(0,0.71954181,0.65411764,0.263529423,0.24423529,1,0.375764715,0.803046458)





asoc4<-cbind(PropFeed4,PropOpen4,dom4)
asoc5<-cbind(PropFeed5,PropOpen5,dom5)
asoc6<-cbind(PropFeed6,PropOpen6,dom6)
asoc7<-cbind(PropFeed7,PropOpen7,dom7)
asoc8<-cbind(PropFeed8,PropOpen8,dom8)
asoc9<-cbind(PropFeed9,PropOpen9,dom9)
asoc10<-cbind(PropFeed10,PropOpen10,dom10)




##### Load associative matrices #####
    # [both directed and symmetrical] and format 
    # for NBDA

group4AsD <- read.csv("Group4AssociativeDir.csv")
rownames(group4AsD)<-group4AsD$X
group4AsD$X<-NULL
m4AsD<-as.matrix(group4AsD, header=T, row.names=1)
a4AsD <- as.array(m4AsD)
group4AsD_network<-as.network(m4AsD, data_format = "SP", association_index = "HWI")
t4AsD<- as.tnet(m4AsD)
d4AsD<- dichotomise_w(t4AsD, GT = 90)

group4AsS <- read.csv("Group4AssociativeSym.csv")
rownames(group4AsS)<-group4AsS$X
group4AsS$X<-NULL
m4AsS<-as.matrix(group4AsS, header=T, row.names=1)
a4AsS <- as.array(m4AsS)
group4AsS_network<-as.network(m4AsS, data_format = "SP", association_index = "HWI")
t4AsS<- as.tnet(m4AsS)
d4AsS<- dichotomise_w(t4AsS, GT = 90)

group5AsD <- read.csv("Group5AssociativeDir.csv")
rownames(group5AsD)<-group5AsD$X
group5AsD$X<-NULL
m5AsD<-as.matrix(group5AsD, header=T, row.names=1)
a5AsD <- as.array(m5AsD)
group5AsD_network<-as.network(m5AsD, data_format = "SP", association_index = "HWI")
t5AsD<- as.tnet(m5AsD)
d5AsD<- dichotomise_w(t5AsD, GT = 90)

group5AsS <- read.csv("Group5AssociativeSym.csv")
rownames(group5AsS)<-group5AsS$X
group5AsS$X<-NULL
m5AsS<-as.matrix(group5AsS, header=T, row.names=1)
a5AsS <- as.array(m5AsS)
group5AsS_network<-as.network(m5AsS, data_format = "SP", association_index = "HWI")
t5AsS<- as.tnet(m5AsS)
d5AsS<- dichotomise_w(t5AsS, GT = 90)

group6AsD <- read.csv("Group6AssociativeDir.csv")
rownames(group6AsD)<-group6AsD$X
group6AsD$X<-NULL
m6AsD<-as.matrix(group6AsD, header=T, row.names=1)
a6AsD <- as.array(m6AsD)
group6AsD_network<-as.network(m6AsD, data_format = "SP", association_index = "HWI")
t6AsD<- as.tnet(m6AsD)
d6AsD<- dichotomise_w(t6AsD, GT = 90)

group6AsS <- read.csv("Group6AssociativeSym.csv")
rownames(group6AsS)<-group6AsS$X
group6AsS$X<-NULL
m6AsS<-as.matrix(group6AsS, header=T, row.names=1)
a6AsS <- as.array(m6AsS)
group6AsS_network<-as.network(m6AsS, data_format = "SP", association_index = "HWI")
t6AsS<- as.tnet(m6AsS)
d6AsS<- dichotomise_w(t6AsS, GT = 90)

group7AsD <- read.csv("Group7AssociativeDir.csv")
rownames(group7AsD)<-group7AsD$X
group7AsD$X<-NULL
m7AsD<-as.matrix(group7AsD, header=T, row.names=1)
a7AsD <- as.array(m7AsD)
group7AsD_network<-as.network(m7AsD, data_format = "SP", association_index = "HWI")
t7AsD<- as.tnet(m7AsD)
d7AsD<- dichotomise_w(t7AsD, GT = 90)

group7AsS <- read.csv("Group7AssociativeSym.csv")
rownames(group7AsS)<-group7AsS$X
group7AsS$X<-NULL
m7AsS<-as.matrix(group7AsS, header=T, row.names=1)
a7AsS <- as.array(m7AsS)
group7AsS_network<-as.network(m7AsS, data_format = "SP", association_index = "HWI")
t7AsS<- as.tnet(m7AsS)
d7AsS<- dichotomise_w(t7AsS, GT = 90)

group8AsD <- read.csv("Group8AssociativeDir.csv")
rownames(group8AsD)<-group8AsD$X
group8AsD$X<-NULL
m8AsD<-as.matrix(group8AsD, header=T, row.names=1)
a8AsD <- as.array(m8AsD)
group8AsD_network<-as.network(m8AsD, data_format = "SP", association_index = "HWI")
t8AsD<- as.tnet(m8AsD)
d8AsD<- dichotomise_w(t8AsD, GT = 90)

group8AsS <- read.csv("Group8AssociativeSym.csv")
rownames(group8AsS)<-group8AsS$X
group8AsS$X<-NULL
m8AsS<-as.matrix(group8AsS, header=T, row.names=1)
a8AsS <- as.array(m8AsS)
group8AsS_network<-as.network(m8AsS, data_format = "SP", association_index = "HWI")
t8AsS<- as.tnet(m8AsS)
d8AsS<- dichotomise_w(t8AsS, GT = 90)

group9AsD <- read.csv("Group9AssociativeDir.csv")
rownames(group9AsD)<-group9AsD$X
group9AsD$X<-NULL
m9AsD<-as.matrix(group9AsD, header=T, row.names=1)
a9AsD <- as.array(m9AsD)
group9AsD_network<-as.network(m9AsD, data_format = "SP", association_index = "HWI")
t9AsD<- as.tnet(m9AsD)
d9AsD<- dichotomise_w(t9AsD, GT = 90)

group9AsS <- read.csv("Group9AssociativeSym.csv")
rownames(group9AsS)<-group9AsS$X
group9AsS$X<-NULL
m9AsS<-as.matrix(group9AsS, header=T, row.names=1)
a9AsS <- as.array(m9AsS)
group9AsS_network<-as.network(m9AsS, data_format = "SP", association_index = "HWI")
t9AsS<- as.tnet(m9AsS)
d9AsS<- dichotomise_w(t9AsS, GT = 90)

group10AsD <- read.csv("Group10AssociativeDir.csv")
rownames(group10AsD)<-group10AsD$X
group10AsD$X<-NULL
m10AsD<-as.matrix(group10AsD, header=T, row.names=1)
a10AsD <- as.array(m10AsD)
group10AsD_network<-as.network(m10AsD, data_format = "SP", association_index = "HWI")
t10AsD<- as.tnet(m10AsD)
d10AsD<- dichotomise_w(t10AsD, GT = 90)

group10AsS <- read.csv("Group10AssociativeSym.csv")
rownames(group10AsS)<-group10AsS$X
group10AsS$X<-NULL
m10AsS<-as.matrix(group10AsS, header=T, row.names=1)
a10AsS <- as.array(m10AsS)
group10AsS_network<-as.network(m10AsS, data_format = "SP", association_index = "HWI")
t10AsS<- as.tnet(m10AsS)
d10AsS<- dichotomise_w(t10AsS, GT = 90)




##### Associative NBDAs (without asocial variables) ##### 
  # for the acquisition of feeding and opening 
  # behaviour, based on the above matrices, 
  # without the inclusion of individual-level 
  # [asocial] variables


OADA4AsDFeed<-oaData(assMatrix=m4AsD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AsDFeed<-addFit(oadata=OADA4AsDFeed)
summary(ModelOADA4AsDFeed)

OADA4AsSFeed<-oaData(assMatrix=m4AsS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AsSFeed<-addFit(oadata=OADA4AsSFeed)
summary(ModelOADA4AsSFeed)

OADA4AsDOpen<-oaData(assMatrix=m4AsD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AsDOpen<-addFit(oadata=OADA4AsDOpen)
summary(ModelOADA4AsDOpen)

OADA4AsSOpen<-oaData(assMatrix=m4AsS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AsSOpen<-addFit(oadata=OADA4AsSOpen)
summary(ModelOADA4AsSOpen)



OADA5AsDFeed<-oaData(assMatrix=m5AsD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AsDFeed<-addFit(oadata=OADA5AsDFeed)
summary(ModelOADA5AsDFeed)

OADA5AsSFeed<-oaData(assMatrix=m5AsS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AsSFeed<-addFit(oadata=OADA5AsSFeed)
summary(ModelOADA5AsSFeed)

OADA5AsDOpen<-oaData(assMatrix=m5AsD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AsDOpen<-addFit(oadata=OADA5AsDOpen)
summary(ModelOADA5AsDOpen)

OADA5AsSOpen<-oaData(assMatrix=m5AsS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AsSOpen<-addFit(oadata=OADA5AsSOpen)
summary(ModelOADA5AsSOpen)



OADA6AsDFeed<-oaData(assMatrix=m6AsD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AsDFeed<-addFit(oadata=OADA6AsDFeed)
summary(ModelOADA6AsDFeed)

OADA6AsSFeed<-oaData(assMatrix=m6AsS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AsSFeed<-addFit(oadata=OADA6AsSFeed)
summary(ModelOADA6AsSFeed)

OADA6AsDOpen<-oaData(assMatrix=m6AsD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AsDOpen<-addFit(oadata=OADA6AsDOpen)
summary(ModelOADA6AsDOpen)

OADA6AsSOpen<-oaData(assMatrix=m6AsS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AsSOpen<-addFit(oadata=OADA6AsSOpen)
summary(ModelOADA6AsSOpen)



OADA7AsDFeed<-oaData(assMatrix=m7AsD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AsDFeed<-addFit(oadata=OADA7AsDFeed)
summary(ModelOADA7AsDFeed)

OADA7AsSFeed<-oaData(assMatrix=m7AsS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AsSFeed<-addFit(oadata=OADA7AsSFeed)
summary(ModelOADA7AsSFeed)

OADA7AsDOpen<-oaData(assMatrix=m7AsD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AsDOpen<-addFit(oadata=OADA7AsDOpen)
summary(ModelOADA7AsDOpen)

OADA7AsSOpen<-oaData(assMatrix=m7AsS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AsSOpen<-addFit(oadata=OADA7AsSOpen)
summary(ModelOADA7AsSOpen)



OADA8AsDFeed<-oaData(assMatrix=m8AsD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AsDFeed<-addFit(oadata=OADA8AsDFeed)
summary(ModelOADA8AsDFeed)

OADA8AsSFeed<-oaData(assMatrix=m8AsS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AsSFeed<-addFit(oadata=OADA8AsSFeed)
summary(ModelOADA8AsSFeed)

OADA8AsDOpen<-oaData(assMatrix=m8AsD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AsDOpen<-addFit(oadata=OADA8AsDOpen)
summary(ModelOADA8AsDOpen)

OADA8AsSOpen<-oaData(assMatrix=m8AsS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AsSOpen<-addFit(oadata=OADA8AsSOpen)
summary(ModelOADA8AsSOpen)



OADA9AsDFeed<-oaData(assMatrix=m9AsD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AsDFeed<-addFit(oadata=OADA9AsDFeed)
summary(ModelOADA9AsDFeed)

OADA9AsSFeed<-oaData(assMatrix=m9AsS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AsSFeed<-addFit(oadata=OADA9AsSFeed)
summary(ModelOADA9AsSFeed)

OADA9AsDOpen<-oaData(assMatrix=m9AsD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AsDOpen<-addFit(oadata=OADA9AsDOpen)
summary(ModelOADA9AsDOpen)

OADA9AsSOpen<-oaData(assMatrix=m9AsS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AsSOpen<-addFit(oadata=OADA9AsSOpen)
summary(ModelOADA9AsSOpen)



OADA10AsDFeed<-oaData(assMatrix=m10AsD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AsDFeed<-addFit(oadata=OADA10AsDFeed)
summary(ModelOADA10AsDFeed)

OADA10AsSFeed<-oaData(assMatrix=m10AsS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AsSFeed<-addFit(oadata=OADA10AsSFeed)
summary(ModelOADA10AsSFeed)

OADA10AsDOpen<-oaData(assMatrix=m10AsD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AsDOpen<-addFit(oadata=OADA10AsDOpen)
summary(ModelOADA10AsDOpen)

OADA10AsSOpen<-oaData(assMatrix=m10AsS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AsSOpen<-addFit(oadata=OADA10AsSOpen)
summary(ModelOADA10AsSOpen)




##### Associative NBDAs (with asocial variables) #####
  # by attributing a numerical vector to the 
  # argument "asocialVar", specifying the columns 
  # in the "asoc" matrix containing the desired 
  # variables. In this case, this includes 
  # propensity to open, to feed, and dominance rank.




OADA4AsDFeed<-oaData(assMatrix=m4AsD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AsDFeed<-addFit(oadata=OADA4AsDFeed, asocialVar=c(1,2,3))
summary(ModelOADA4AsDFeed)

OADA4AsSFeed<-oaData(assMatrix=m4AsS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AsSFeed<-addFit(oadata=OADA4AsSFeed, asocialVar=c(1,2,3))
summary(ModelOADA4AsSFeed)

OADA4AsDOpen<-oaData(assMatrix=m4AsD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AsDOpen<-addFit(oadata=OADA4AsDOpen, asocialVar=c(1,2,3))
summary(ModelOADA4AsDOpen)

OADA4AsSOpen<-oaData(assMatrix=m4AsS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AsSOpen<-addFit(oadata=OADA4AsSOpen, asocialVar=c(1,2,3))
summary(ModelOADA4AsSOpen)



OADA5AsDFeed<-oaData(assMatrix=m5AsD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AsDFeed<-addFit(oadata=OADA5AsDFeed, asocialVar=c(1,2,3))
summary(ModelOADA5AsDFeed)

OADA5AsSFeed<-oaData(assMatrix=m5AsS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AsSFeed<-addFit(oadata=OADA5AsSFeed, asocialVar=c(1,2,3))
summary(ModelOADA5AsSFeed)

OADA5AsDOpen<-oaData(assMatrix=m5AsD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AsDOpen<-addFit(oadata=OADA5AsDOpen, asocialVar=c(1,2,3))
summary(ModelOADA5AsDOpen)

OADA5AsSOpen<-oaData(assMatrix=m5AsS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AsSOpen<-addFit(oadata=OADA5AsSOpen, asocialVar=c(1,2,3))
summary(ModelOADA5AsSOpen)



OADA6AsDFeed<-oaData(assMatrix=m6AsD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AsDFeed<-addFit(oadata=OADA6AsDFeed, asocialVar=c(1,2,3))
summary(ModelOADA6AsDFeed)

OADA6AsSFeed<-oaData(assMatrix=m6AsS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AsSFeed<-addFit(oadata=OADA6AsSFeed, asocialVar=c(1,2,3))
summary(ModelOADA6AsSFeed)

OADA6AsDOpen<-oaData(assMatrix=m6AsD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AsDOpen<-addFit(oadata=OADA6AsDOpen, asocialVar=c(1,2,3))
summary(ModelOADA6AsDOpen)

OADA6AsSOpen<-oaData(assMatrix=m6AsS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AsSOpen<-addFit(oadata=OADA6AsSOpen, asocialVar=c(1,2,3))
summary(ModelOADA6AsSOpen)



OADA7AsDFeed<-oaData(assMatrix=m7AsD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AsDFeed<-addFit(oadata=OADA7AsDFeed, asocialVar=c(1,2,3))
summary(ModelOADA7AsDFeed)

OADA7AsSFeed<-oaData(assMatrix=m7AsS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AsSFeed<-addFit(oadata=OADA7AsSFeed, asocialVar=c(1,2,3))
summary(ModelOADA7AsSFeed)

OADA7AsDOpen<-oaData(assMatrix=m7AsD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AsDOpen<-addFit(oadata=OADA7AsDOpen, asocialVar=c(1,2,3))
summary(ModelOADA7AsDOpen)

OADA7AsSOpen<-oaData(assMatrix=m7AsS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AsSOpen<-addFit(oadata=OADA7AsSOpen, asocialVar=c(1,2,3))
summary(ModelOADA7AsSOpen)



OADA8AsDFeed<-oaData(assMatrix=m8AsD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AsDFeed<-addFit(oadata=OADA8AsDFeed, asocialVar=c(1,2,3))
summary(ModelOADA8AsDFeed)

OADA8AsSFeed<-oaData(assMatrix=m8AsS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AsSFeed<-addFit(oadata=OADA8AsSFeed, asocialVar=c(1,2,3))
summary(ModelOADA8AsSFeed)

OADA8AsDOpen<-oaData(assMatrix=m8AsD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AsDOpen<-addFit(oadata=OADA8AsDOpen, asocialVar=c(1,2,3))
summary(ModelOADA8AsDOpen)

OADA8AsSOpen<-oaData(assMatrix=m8AsS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AsSOpen<-addFit(oadata=OADA8AsSOpen, asocialVar=c(1,2,3))
summary(ModelOADA8AsSOpen)



OADA9AsDFeed<-oaData(assMatrix=m9AsD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AsDFeed<-addFit(oadata=OADA9AsDFeed, asocialVar=c(1,2,3))
summary(ModelOADA9AsDFeed)

OADA9AsSFeed<-oaData(assMatrix=m9AsS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AsSFeed<-addFit(oadata=OADA9AsSFeed, asocialVar=c(1,2,3))
summary(ModelOADA9AsSFeed)

OADA9AsDOpen<-oaData(assMatrix=m9AsD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AsDOpen<-addFit(oadata=OADA9AsDOpen, asocialVar=c(1,2,3))
summary(ModelOADA9AsDOpen)

OADA9AsSOpen<-oaData(assMatrix=m9AsS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AsSOpen<-addFit(oadata=OADA9AsSOpen, asocialVar=c(1,2,3))
summary(ModelOADA9AsSOpen)



OADA10AsDFeed<-oaData(assMatrix=m10AsD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AsDFeed<-addFit(oadata=OADA10AsDFeed, asocialVar=c(1,2,3))
summary(ModelOADA10AsDFeed)

OADA10AsSFeed<-oaData(assMatrix=m10AsS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AsSFeed<-addFit(oadata=OADA10AsSFeed, asocialVar=c(1,2,3))
summary(ModelOADA10AsSFeed)

OADA10AsDOpen<-oaData(assMatrix=m10AsD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AsDOpen<-addFit(oadata=OADA10AsDOpen, asocialVar=c(1,2,3))
summary(ModelOADA10AsDOpen)

OADA10AsSOpen<-oaData(assMatrix=m10AsS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AsSOpen<-addFit(oadata=OADA10AsSOpen, asocialVar=c(1,2,3))
summary(ModelOADA10AsSOpen)




##### REPEAT for agonistic interactions #####


group4AgD <- read.csv("Group4AgonisticDir.csv")
rownames(group4AgD)<-group4AgD$X
group4AgD$X<-NULL
m4AgD<-as.matrix(group4AgD, header=T, row.names=1)
a4AgD <- as.array(m4AgD)
group4AgD_network<-as.network(m4AgD, data_format = "SP", association_index = "HWI")
t4AgD<- as.tnet(m4AgD)
d4AgD<- dichotomise_w(t4AgD, GT = 90)

group4AgS <- read.csv("Group4AgonisticSym.csv")
rownames(group4AgS)<-group4AgS$X
group4AgS$X<-NULL
m4AgS<-as.matrix(group4AgS, header=T, row.names=1)
a4AgS <- as.array(m4AgS)
group4AgS_network<-as.network(m4AgS, data_format = "SP", association_index = "HWI")
t4AgS<- as.tnet(m4AgS)
d4AgS<- dichotomise_w(t4AgS, GT = 90)

group5AgD <- read.csv("Group5AgonisticDir.csv")
rownames(group5AgD)<-group5AgD$X
group5AgD$X<-NULL
m5AgD<-as.matrix(group5AgD, header=T, row.names=1)
a5AgD <- as.array(m5AgD)
group5AgD_network<-as.network(m5AgD, data_format = "SP", association_index = "HWI")
t5AgD<- as.tnet(m5AgD)
d5AgD<- dichotomise_w(t5AgD, GT = 90)

group5AgS <- read.csv("Group5AgonisticSym.csv")
rownames(group5AgS)<-group5AgS$X
group5AgS$X<-NULL
m5AgS<-as.matrix(group5AgS, header=T, row.names=1)
a5AgS <- as.array(m5AgS)
group5AgS_network<-as.network(m5AgS, data_format = "SP", association_index = "HWI")
t5AgS<- as.tnet(m5AgS)
d5AgS<- dichotomise_w(t5AgS, GT = 90)

group6AgD <- read.csv("Group6AgonisticDir.csv")
rownames(group6AgD)<-group6AgD$X
group6AgD$X<-NULL
m6AgD<-as.matrix(group6AgD, header=T, row.names=1)
a6AgD <- as.array(m6AgD)
group6AgD_network<-as.network(m6AgD, data_format = "SP", association_index = "HWI")
t6AgD<- as.tnet(m6AgD)
d6AgD<- dichotomise_w(t6AgD, GT = 90)

group6AgS <- read.csv("Group6AgonisticSym.csv")
rownames(group6AgS)<-group6AgS$X
group6AgS$X<-NULL
m6AgS<-as.matrix(group6AgS, header=T, row.names=1)
a6AgS <- as.array(m6AgS)
group6AgS_network<-as.network(m6AgS, data_format = "SP", association_index = "HWI")
t6AgS<- as.tnet(m6AgS)
d6AgS<- dichotomise_w(t6AgS, GT = 90)

group7AgD <- read.csv("Group7AgonisticDir.csv")
rownames(group7AgD)<-group7AgD$X
group7AgD$X<-NULL
m7AgD<-as.matrix(group7AgD, header=T, row.names=1)
a7AgD <- as.array(m7AgD)
group7AgD_network<-as.network(m7AgD, data_format = "SP", association_index = "HWI")
t7AgD<- as.tnet(m7AgD)
d7AgD<- dichotomise_w(t7AgD, GT = 90)

group7AgS <- read.csv("Group7AgonisticSym.csv")
rownames(group7AgS)<-group7AgS$X
group7AgS$X<-NULL
m7AgS<-as.matrix(group7AgS, header=T, row.names=1)
a7AgS <- as.array(m7AgS)
group7AgS_network<-as.network(m7AgS, data_format = "SP", association_index = "HWI")
t7AgS<- as.tnet(m7AgS)
d7AgS<- dichotomise_w(t7AgS, GT = 90)

group8AgD <- read.csv("Group8AgonisticDir.csv")
rownames(group8AgD)<-group8AgD$X
group8AgD$X<-NULL
m8AgD<-as.matrix(group8AgD, header=T, row.names=1)
a8AgD <- as.array(m8AgD)
group8AgD_network<-as.network(m8AgD, data_format = "SP", association_index = "HWI")
t8AgD<- as.tnet(m8AgD)
d8AgD<- dichotomise_w(t8AgD, GT = 90)

group8AgS <- read.csv("Group8AgonisticSym.csv")
rownames(group8AgS)<-group8AgS$X
group8AgS$X<-NULL
m8AgS<-as.matrix(group8AgS, header=T, row.names=1)
a8AgS <- as.array(m8AgS)
group8AgS_network<-as.network(m8AgS, data_format = "SP", association_index = "HWI")
t8AgS<- as.tnet(m8AgS)
d8AgS<- dichotomise_w(t8AgS, GT = 90)

group9AgD <- read.csv("Group9AgonisticDir.csv")
rownames(group9AgD)<-group9AgD$X
group9AgD$X<-NULL
m9AgD<-as.matrix(group9AgD, header=T, row.names=1)
a9AgD <- as.array(m9AgD)
group9AgD_network<-as.network(m9AgD, data_format = "SP", association_index = "HWI")
t9AgD<- as.tnet(m9AgD)
d9AgD<- dichotomise_w(t9AgD, GT = 90)

group9AgS <- read.csv("Group9AgonisticSym.csv")
rownames(group9AgS)<-group9AgS$X
group9AgS$X<-NULL
m9AgS<-as.matrix(group9AgS, header=T, row.names=1)
a9AgS <- as.array(m9AgS)
group9AgS_network<-as.network(m9AgS, data_format = "SP", association_index = "HWI")
t9AgS<- as.tnet(m9AgS)
d9AgS<- dichotomise_w(t9AgS, GT = 90)

group10AgD <- read.csv("Group10AgonisticDir.csv")
rownames(group10AgD)<-group10AgD$X
group10AgD$X<-NULL
m10AgD<-as.matrix(group10AgD, header=T, row.names=1)
a10AgD <- as.array(m10AgD)
group10AgD_network<-as.network(m10AgD, data_format = "SP", association_index = "HWI")
t10AgD<- as.tnet(m10AgD)
d10AgD<- dichotomise_w(t10AgD, GT = 90)

group10AgS <- read.csv("Group10AgonisticSym.csv")
rownames(group10AgS)<-group10AgS$X
group10AgS$X<-NULL
m10AgS<-as.matrix(group10AgS, header=T, row.names=1)
a10AgS <- as.array(m10AgS)
group10AgS_network<-as.network(m10AgS, data_format = "SP", association_index = "HWI")
t10AgS<- as.tnet(m10AgS)
d10AgS<- dichotomise_w(t10AgS, GT = 90)




###### Agonistic NBDAs (without asocial variables #####

OADA4AgDFeed<-oaData(assMatrix=m4AgD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AgDFeed<-addFit(oadata=OADA4AgDFeed)
summary(ModelOADA4AgDFeed)

OADA4AgSFeed<-oaData(assMatrix=m4AgS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AgSFeed<-addFit(oadata=OADA4AgSFeed)
summary(ModelOADA4AgSFeed)

OADA4AgDOpen<-oaData(assMatrix=m4AgD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AgDOpen<-addFit(oadata=OADA4AgDOpen)
summary(ModelOADA4AgDOpen)

OADA4AgSOpen<-oaData(assMatrix=m4AgS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AgSOpen<-addFit(oadata=OADA4AgSOpen)
summary(ModelOADA4AgSOpen)



OADA5AgDFeed<-oaData(assMatrix=m5AgD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AgDFeed<-addFit(oadata=OADA5AgDFeed)
summary(ModelOADA5AgDFeed)

OADA5AgSFeed<-oaData(assMatrix=m5AgS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AgSFeed<-addFit(oadata=OADA5AgSFeed)
summary(ModelOADA5AgSFeed)

OADA5AgDOpen<-oaData(assMatrix=m5AgD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AgDOpen<-addFit(oadata=OADA5AgDOpen)
summary(ModelOADA5AgDOpen)

OADA5AgSOpen<-oaData(assMatrix=m5AgS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AgSOpen<-addFit(oadata=OADA5AgSOpen)
summary(ModelOADA5AgSOpen)



OADA6AgDFeed<-oaData(assMatrix=m6AgD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AgDFeed<-addFit(oadata=OADA6AgDFeed)
summary(ModelOADA6AgDFeed)

OADA6AgSFeed<-oaData(assMatrix=m6AgS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AgSFeed<-addFit(oadata=OADA6AgSFeed)
summary(ModelOADA6AgSFeed)

OADA6AgDOpen<-oaData(assMatrix=m6AgD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AgDOpen<-addFit(oadata=OADA6AgDOpen)
summary(ModelOADA6AgDOpen)

OADA6AgSOpen<-oaData(assMatrix=m6AgS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AgSOpen<-addFit(oadata=OADA6AgSOpen)
summary(ModelOADA6AgSOpen)



OADA7AgDFeed<-oaData(assMatrix=m7AgD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AgDFeed<-addFit(oadata=OADA7AgDFeed)
summary(ModelOADA7AgDFeed)

OADA7AgSFeed<-oaData(assMatrix=m7AgS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AgSFeed<-addFit(oadata=OADA7AgSFeed)
summary(ModelOADA7AgSFeed)

OADA7AgDOpen<-oaData(assMatrix=m7AgD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AgDOpen<-addFit(oadata=OADA7AgDOpen)
summary(ModelOADA7AgDOpen)

OADA7AgSOpen<-oaData(assMatrix=m7AgS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AgSOpen<-addFit(oadata=OADA7AgSOpen)
summary(ModelOADA7AgSOpen)



OADA8AgDFeed<-oaData(assMatrix=m8AgD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AgDFeed<-addFit(oadata=OADA8AgDFeed)
summary(ModelOADA8AgDFeed)

OADA8AgSFeed<-oaData(assMatrix=m8AgS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AgSFeed<-addFit(oadata=OADA8AgSFeed)
summary(ModelOADA8AgSFeed)

OADA8AgDOpen<-oaData(assMatrix=m8AgD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AgDOpen<-addFit(oadata=OADA8AgDOpen)
summary(ModelOADA8AgDOpen)

OADA8AgSOpen<-oaData(assMatrix=m8AgS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AgSOpen<-addFit(oadata=OADA8AgSOpen)
summary(ModelOADA8AgSOpen)



OADA9AgDFeed<-oaData(assMatrix=m9AgD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AgDFeed<-addFit(oadata=OADA9AgDFeed)
summary(ModelOADA9AgDFeed)

OADA9AgSFeed<-oaData(assMatrix=m9AgS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AgSFeed<-addFit(oadata=OADA9AgSFeed)
summary(ModelOADA9AgSFeed)

OADA9AgDOpen<-oaData(assMatrix=m9AgD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AgDOpen<-addFit(oadata=OADA9AgDOpen)
summary(ModelOADA9AgDOpen)

OADA9AgSOpen<-oaData(assMatrix=m9AgS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AgSOpen<-addFit(oadata=OADA9AgSOpen)
summary(ModelOADA9AgSOpen)



OADA10AgDFeed<-oaData(assMatrix=m10AgD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AgDFeed<-addFit(oadata=OADA10AgDFeed)
summary(ModelOADA10AgDFeed)

OADA10AgSFeed<-oaData(assMatrix=m10AgS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AgSFeed<-addFit(oadata=OADA10AgSFeed)
summary(ModelOADA10AgSFeed)

OADA10AgDOpen<-oaData(assMatrix=m10AgD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AgDOpen<-addFit(oadata=OADA10AgDOpen)
summary(ModelOADA10AgDOpen)

OADA10AgSOpen<-oaData(assMatrix=m10AgS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AgSOpen<-addFit(oadata=OADA10AgSOpen)
summary(ModelOADA10AgSOpen)




##### Agonistic NBDAs (with asocial variables) #####


OADA4AgDFeed<-oaData(assMatrix=m4AgD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AgDFeed<-addFit(oadata=OADA4AgDFeed, asocialVar=c(1,2,3))
summary(ModelOADA4AgDFeed)

OADA4AgSFeed<-oaData(assMatrix=m4AgS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AgSFeed<-addFit(oadata=OADA4AgSFeed, asocialVar=c(1,2,3))
summary(ModelOADA4AgSFeed)

OADA4AgDOpen<-oaData(assMatrix=m4AgD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AgDOpen<-addFit(oadata=OADA4AgDOpen, asocialVar=c(1,2,3))
summary(ModelOADA4AgDOpen)

OADA4AgSOpen<-oaData(assMatrix=m4AgS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AgSOpen<-addFit(oadata=OADA4AgSOpen, asocialVar=c(1,2,3))
summary(ModelOADA4AgSOpen)



OADA5AgDFeed<-oaData(assMatrix=m5AgD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AgDFeed<-addFit(oadata=OADA5AgDFeed, asocialVar=c(1,2,3))
summary(ModelOADA5AgDFeed)

OADA5AgSFeed<-oaData(assMatrix=m5AgS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AgSFeed<-addFit(oadata=OADA5AgSFeed, asocialVar=c(1,2,3))
summary(ModelOADA5AgSFeed)

OADA5AgDOpen<-oaData(assMatrix=m5AgD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AgDOpen<-addFit(oadata=OADA5AgDOpen, asocialVar=c(1,2,3))
summary(ModelOADA5AgDOpen)

OADA5AgSOpen<-oaData(assMatrix=m5AgS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AgSOpen<-addFit(oadata=OADA5AgSOpen, asocialVar=c(1,2,3))
summary(ModelOADA5AgSOpen)



OADA6AgDFeed<-oaData(assMatrix=m6AgD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AgDFeed<-addFit(oadata=OADA6AgDFeed, asocialVar=c(1,2,3))
summary(ModelOADA6AgDFeed)

OADA6AgSFeed<-oaData(assMatrix=m6AgS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AgSFeed<-addFit(oadata=OADA6AgSFeed, asocialVar=c(1,2,3))
summary(ModelOADA6AgSFeed)

OADA6AgDOpen<-oaData(assMatrix=m6AgD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AgDOpen<-addFit(oadata=OADA6AgDOpen, asocialVar=c(1,2,3))
summary(ModelOADA6AgDOpen)

OADA6AgSOpen<-oaData(assMatrix=m6AgS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AgSOpen<-addFit(oadata=OADA6AgSOpen, asocialVar=c(1,2,3))
summary(ModelOADA6AgSOpen)



OADA7AgDFeed<-oaData(assMatrix=m7AgD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AgDFeed<-addFit(oadata=OADA7AgDFeed, asocialVar=c(1,2,3))
summary(ModelOADA7AgDFeed)

OADA7AgSFeed<-oaData(assMatrix=m7AgS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AgSFeed<-addFit(oadata=OADA7AgSFeed, asocialVar=c(1,2,3))
summary(ModelOADA7AgSFeed)

OADA7AgDOpen<-oaData(assMatrix=m7AgD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AgDOpen<-addFit(oadata=OADA7AgDOpen, asocialVar=c(1,2,3))
summary(ModelOADA7AgDOpen)

OADA7AgSOpen<-oaData(assMatrix=m7AgS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AgSOpen<-addFit(oadata=OADA7AgSOpen, asocialVar=c(1,2,3))
summary(ModelOADA7AgSOpen)



OADA8AgDFeed<-oaData(assMatrix=m8AgD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AgDFeed<-addFit(oadata=OADA8AgDFeed, asocialVar=c(1,2,3))
summary(ModelOADA8AgDFeed)

OADA8AgSFeed<-oaData(assMatrix=m8AgS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AgSFeed<-addFit(oadata=OADA8AgSFeed, asocialVar=c(1,2,3))
summary(ModelOADA8AgSFeed)

OADA8AgDOpen<-oaData(assMatrix=m8AgD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AgDOpen<-addFit(oadata=OADA8AgDOpen, asocialVar=c(1,2,3))
summary(ModelOADA8AgDOpen)

OADA8AgSOpen<-oaData(assMatrix=m8AgS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AgSOpen<-addFit(oadata=OADA8AgSOpen, asocialVar=c(1,2,3))
summary(ModelOADA8AgSOpen)



OADA9AgDFeed<-oaData(assMatrix=m9AgD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AgDFeed<-addFit(oadata=OADA9AgDFeed, asocialVar=c(1,2,3))
summary(ModelOADA9AgDFeed)

OADA9AgSFeed<-oaData(assMatrix=m9AgS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AgSFeed<-addFit(oadata=OADA9AgSFeed, asocialVar=c(1,2,3))
summary(ModelOADA9AgSFeed)

OADA9AgDOpen<-oaData(assMatrix=m9AgD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AgDOpen<-addFit(oadata=OADA9AgDOpen, asocialVar=c(1,2,3))
summary(ModelOADA9AgDOpen)

OADA9AgSOpen<-oaData(assMatrix=m9AgS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AgSOpen<-addFit(oadata=OADA9AgSOpen, asocialVar=c(1,2,3))
summary(ModelOADA9AgSOpen)



OADA10AgDFeed<-oaData(assMatrix=m10AgD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AgDFeed<-addFit(oadata=OADA10AgDFeed, asocialVar=c(1,2,3))
summary(ModelOADA10AgDFeed)

OADA10AgSFeed<-oaData(assMatrix=m10AgS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AgSFeed<-addFit(oadata=OADA10AgSFeed, asocialVar=c(1,2,3))
summary(ModelOADA10AgSFeed)

OADA10AgDOpen<-oaData(assMatrix=m10AgD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AgDOpen<-addFit(oadata=OADA10AgDOpen, asocialVar=c(1,2,3))
summary(ModelOADA10AgDOpen)

OADA10AgSOpen<-oaData(assMatrix=m10AgS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AgSOpen<-addFit(oadata=OADA10AgSOpen, asocialVar=c(1,2,3))
summary(ModelOADA10AgSOpen)




##### REPEAT for the affiliative interactions #####


group4AfD <- read.csv("Group4AffiliativeDir.csv")
rownames(group4AfD)<-group4AfD$X
group4AfD$X<-NULL
m4AfD<-as.matrix(group4AfD, header=T, row.names=1)
a4AfD <- as.array(m4AfD)
group4AfD_network<-as.network(m4AfD, data_format = "SP", association_index = "HWI")
t4AfD<- as.tnet(m4AfD)
d4AfD<- dichotomise_w(t4AfD, GT = 90)

group4AfS <- read.csv("Group4AffiliativeSym.csv")
rownames(group4AfS)<-group4AfS$X
group4AfS$X<-NULL
m4AfS<-as.matrix(group4AfS, header=T, row.names=1)
a4AfS <- as.array(m4AfS)
group4AfS_network<-as.network(m4AfS, data_format = "SP", association_index = "HWI")
t4AfS<- as.tnet(m4AfS)
d4AfS<- dichotomise_w(t4AfS, GT = 90)

group5AfD <- read.csv("Group5AffiliativeDir.csv")
rownames(group5AfD)<-group5AfD$X
group5AfD$X<-NULL
m5AfD<-as.matrix(group5AfD, header=T, row.names=1)
a5AfD <- as.array(m5AfD)
group5AfD_network<-as.network(m5AfD, data_format = "SP", association_index = "HWI")
t5AfD<- as.tnet(m5AfD)
d5AfD<- dichotomise_w(t5AfD, GT = 90)

group5AfS <- read.csv("Group5AffiliativeSym.csv")
rownames(group5AfS)<-group5AfS$X
group5AfS$X<-NULL
m5AfS<-as.matrix(group5AfS, header=T, row.names=1)
a5AfS <- as.array(m5AfS)
group5AfS_network<-as.network(m5AfS, data_format = "SP", association_index = "HWI")
t5AfS<- as.tnet(m5AfS)
d5AfS<- dichotomise_w(t5AfS, GT = 90)

group6AfD <- read.csv("Group6AffiliativeDir.csv")
rownames(group6AfD)<-group6AfD$X
group6AfD$X<-NULL
m6AfD<-as.matrix(group6AfD, header=T, row.names=1)
a6AfD <- as.array(m6AfD)
group6AfD_network<-as.network(m6AfD, data_format = "SP", association_index = "HWI")
t6AfD<- as.tnet(m6AfD)
d6AfD<- dichotomise_w(t6AfD, GT = 90)

group6AfS <- read.csv("Group6AffiliativeSym.csv")
rownames(group6AfS)<-group6AfS$X
group6AfS$X<-NULL
m6AfS<-as.matrix(group6AfS, header=T, row.names=1)
a6AfS <- as.array(m6AfS)
group6AfS_network<-as.network(m6AfS, data_format = "SP", association_index = "HWI")
t6AfS<- as.tnet(m6AfS)
d6AfS<- dichotomise_w(t6AfS, GT = 90)

group7AfD <- read.csv("Group7AffiliativeDir.csv")
rownames(group7AfD)<-group7AfD$X
group7AfD$X<-NULL
m7AfD<-as.matrix(group7AfD, header=T, row.names=1)
a7AfD <- as.array(m7AfD)
group7AfD_network<-as.network(m7AfD, data_format = "SP", association_index = "HWI")
t7AfD<- as.tnet(m7AfD)
d7AfD<- dichotomise_w(t7AfD, GT = 90)

group7AfS <- read.csv("Group7AffiliativeSym.csv")
rownames(group7AfS)<-group7AfS$X
group7AfS$X<-NULL
m7AfS<-as.matrix(group7AfS, header=T, row.names=1)
a7AfS <- as.array(m7AfS)
group7AfS_network<-as.network(m7AfS, data_format = "SP", association_index = "HWI")
t7AfS<- as.tnet(m7AfS)
d7AfS<- dichotomise_w(t7AfS, GT = 90)

group8AfD <- read.csv("Group8AffiliativeDir.csv")
rownames(group8AfD)<-group8AfD$X
group8AfD$X<-NULL
m8AfD<-as.matrix(group8AfD, header=T, row.names=1)
a8AfD <- as.array(m8AfD)
group8AfD_network<-as.network(m8AfD, data_format = "SP", association_index = "HWI")
t8AfD<- as.tnet(m8AfD)
d8AfD<- dichotomise_w(t8AfD, GT = 90)

group8AfS <- read.csv("Group8AffiliativeSym.csv")
rownames(group8AfS)<-group8AfS$X
group8AfS$X<-NULL
m8AfS<-as.matrix(group8AfS, header=T, row.names=1)
a8AfS <- as.array(m8AfS)
group8AfS_network<-as.network(m8AfS, data_format = "SP", association_index = "HWI")
t8AfS<- as.tnet(m8AfS)
d8AfS<- dichotomise_w(t8AfS, GT = 90)

group9AfD <- read.csv("Group9AffiliativeDir.csv")
rownames(group9AfD)<-group9AfD$X
group9AfD$X<-NULL
m9AfD<-as.matrix(group9AfD, header=T, row.names=1)
a9AfD <- as.array(m9AfD)
group9AfD_network<-as.network(m9AfD, data_format = "SP", association_index = "HWI")
t9AfD<- as.tnet(m9AfD)
d9AfD<- dichotomise_w(t9AfD, GT = 90)

group9AfS <- read.csv("Group9AffiliativeSym.csv")
rownames(group9AfS)<-group9AfS$X
group9AfS$X<-NULL
m9AfS<-as.matrix(group9AfS, header=T, row.names=1)
a9AfS <- as.array(m9AfS)
group9AfS_network<-as.network(m9AfS, data_format = "SP", association_index = "HWI")
t9AfS<- as.tnet(m9AfS)
d9AfS<- dichotomise_w(t9AfS, GT = 90)

group10AfD <- read.csv("Group10AffiliativeDir.csv")
rownames(group10AfD)<-group10AfD$X
group10AfD$X<-NULL
m10AfD<-as.matrix(group10AfD, header=T, row.names=1)
a10AfD <- as.array(m10AfD)
group10AfD_network<-as.network(m10AfD, data_format = "SP", association_index = "HWI")
t10AfD<- as.tnet(m10AfD)
d10AfD<- dichotomise_w(t10AfD, GT = 90)

group10AfS <- read.csv("Group10AffiliativeSym.csv")
rownames(group10AfS)<-group10AfS$X
group10AfS$X<-NULL
m10AfS<-as.matrix(group10AfS, header=T, row.names=1)
a10AfS <- as.array(m10AfS)
group10AfS_network<-as.network(m10AfS, data_format = "SP", association_index = "HWI")
t10AfS<- as.tnet(m10AfS)
d10AfS<- dichotomise_w(t10AfS, GT = 90)




##### Affiliative NBDAs (without asocial variables) #####




OADA4AfDFeed<-oaData(assMatrix=m4AfD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AfDFeed<-addFit(oadata=OADA4AfDFeed)
summary(ModelOADA4AfDFeed)

OADA4AfSFeed<-oaData(assMatrix=m4AfS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AfSFeed<-addFit(oadata=OADA4AfSFeed)
summary(ModelOADA4AfSFeed)

OADA4AfDOpen<-oaData(assMatrix=m4AfD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AfDOpen<-addFit(oadata=OADA4AfDOpen)
summary(ModelOADA4AfDOpen)

OADA4AfSOpen<-oaData(assMatrix=m4AfS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AfSOpen<-addFit(oadata=OADA4AfSOpen)
summary(ModelOADA4AfSOpen)



OADA5AfDFeed<-oaData(assMatrix=m5AfD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AfDFeed<-addFit(oadata=OADA5AfDFeed)
summary(ModelOADA5AfDFeed)

OADA5AfSFeed<-oaData(assMatrix=m5AfS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AfSFeed<-addFit(oadata=OADA5AfSFeed)
summary(ModelOADA5AfSFeed)

OADA5AfDOpen<-oaData(assMatrix=m5AfD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AfDOpen<-addFit(oadata=OADA5AfDOpen)
summary(ModelOADA5AfDOpen)

OADA5AfSOpen<-oaData(assMatrix=m5AfS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AfSOpen<-addFit(oadata=OADA5AfSOpen)
summary(ModelOADA5AfSOpen)



OADA6AfDFeed<-oaData(assMatrix=m6AfD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AfDFeed<-addFit(oadata=OADA6AfDFeed)
summary(ModelOADA6AfDFeed)

OADA6AfSFeed<-oaData(assMatrix=m6AfS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AfSFeed<-addFit(oadata=OADA6AfSFeed)
summary(ModelOADA6AfSFeed)

OADA6AfDOpen<-oaData(assMatrix=m6AfD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AfDOpen<-addFit(oadata=OADA6AfDOpen)
summary(ModelOADA6AfDOpen)

OADA6AfSOpen<-oaData(assMatrix=m6AfS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AfSOpen<-addFit(oadata=OADA6AfSOpen)
summary(ModelOADA6AfSOpen)



OADA7AfDFeed<-oaData(assMatrix=m7AfD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AfDFeed<-addFit(oadata=OADA7AfDFeed)
summary(ModelOADA7AfDFeed)

OADA7AfSFeed<-oaData(assMatrix=m7AfS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AfSFeed<-addFit(oadata=OADA7AfSFeed)
summary(ModelOADA7AfSFeed)

OADA7AfDOpen<-oaData(assMatrix=m7AfD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AfDOpen<-addFit(oadata=OADA7AfDOpen)
summary(ModelOADA7AfDOpen)

OADA7AfSOpen<-oaData(assMatrix=m7AfS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AfSOpen<-addFit(oadata=OADA7AfSOpen)
summary(ModelOADA7AfSOpen)



OADA8AfDFeed<-oaData(assMatrix=m8AfD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AfDFeed<-addFit(oadata=OADA8AfDFeed)
summary(ModelOADA8AfDFeed)

OADA8AfSFeed<-oaData(assMatrix=m8AfS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AfSFeed<-addFit(oadata=OADA8AfSFeed)
summary(ModelOADA8AfSFeed)

OADA8AfDOpen<-oaData(assMatrix=m8AfD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AfDOpen<-addFit(oadata=OADA8AfDOpen)
summary(ModelOADA8AfDOpen)

OADA8AfSOpen<-oaData(assMatrix=m8AfS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AfSOpen<-addFit(oadata=OADA8AfSOpen)
summary(ModelOADA8AfSOpen)



OADA9AfDFeed<-oaData(assMatrix=m9AfD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AfDFeed<-addFit(oadata=OADA9AfDFeed)
summary(ModelOADA9AfDFeed)

OADA9AfSFeed<-oaData(assMatrix=m9AfS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AfSFeed<-addFit(oadata=OADA9AfSFeed)
summary(ModelOADA9AfSFeed)

OADA9AfDOpen<-oaData(assMatrix=m9AfD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AfDOpen<-addFit(oadata=OADA9AfDOpen)
summary(ModelOADA9AfDOpen)

OADA9AfSOpen<-oaData(assMatrix=m9AfS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AfSOpen<-addFit(oadata=OADA9AfSOpen)
summary(ModelOADA9AfSOpen)



OADA10AfDFeed<-oaData(assMatrix=m10AfD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AfDFeed<-addFit(oadata=OADA10AfDFeed)
summary(ModelOADA10AfDFeed)

OADA10AfSFeed<-oaData(assMatrix=m10AfS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AfSFeed<-addFit(oadata=OADA10AfSFeed)
summary(ModelOADA10AfSFeed)

OADA10AfDOpen<-oaData(assMatrix=m10AfD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AfDOpen<-addFit(oadata=OADA10AfDOpen)
summary(ModelOADA10AfDOpen)

OADA10AfSOpen<-oaData(assMatrix=m10AfS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AfSOpen<-addFit(oadata=OADA10AfSOpen)
summary(ModelOADA10AfSOpen)




###### Affiliative NBDAs (with asocial variables) #####



OADA4AfDFeed<-oaData(assMatrix=m4AfD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AfDFeed<-addFit(oadata=OADA4AfDFeed, asocialVar=c(1,2,3))
summary(ModelOADA4AfDFeed)

OADA4AfSFeed<-oaData(assMatrix=m4AfS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed4, groupid="4", taskid="F", weights=PropFeed4)
ModelOADA4AfSFeed<-addFit(oadata=OADA4AfSFeed, asocialVar=c(1,2,3))
summary(ModelOADA4AfSFeed)

OADA4AfDOpen<-oaData(assMatrix=m4AfD, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AfDOpen<-addFit(oadata=OADA4AfDOpen, asocialVar=c(1,2,3))
summary(ModelOADA4AfDOpen)

OADA4AfSOpen<-oaData(assMatrix=m4AfS, asoc=asoc4, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen4, groupid="4", taskid="O", weights=PropOpen4)
ModelOADA4AfSOpen<-addFit(oadata=OADA4AfSOpen, asocialVar=c(1,2,3))
summary(ModelOADA4AfSOpen)



OADA5AfDFeed<-oaData(assMatrix=m5AfD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AfDFeed<-addFit(oadata=OADA5AfDFeed, asocialVar=c(1,2,3))
summary(ModelOADA5AfDFeed)

OADA5AfSFeed<-oaData(assMatrix=m5AfS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed5, groupid="5", taskid="F", weights=PropFeed5)
ModelOADA5AfSFeed<-addFit(oadata=OADA5AfSFeed, asocialVar=c(1,2,3))
summary(ModelOADA5AfSFeed)

OADA5AfDOpen<-oaData(assMatrix=m5AfD, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AfDOpen<-addFit(oadata=OADA5AfDOpen, asocialVar=c(1,2,3))
summary(ModelOADA5AfDOpen)

OADA5AfSOpen<-oaData(assMatrix=m5AfS, asoc=asoc5, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen5, groupid="5", taskid="O", weights=PropOpen5)
ModelOADA5AfSOpen<-addFit(oadata=OADA5AfSOpen, asocialVar=c(1,2,3))
summary(ModelOADA5AfSOpen)



OADA6AfDFeed<-oaData(assMatrix=m6AfD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AfDFeed<-addFit(oadata=OADA6AfDFeed, asocialVar=c(1,2,3))
summary(ModelOADA6AfDFeed)

OADA6AfSFeed<-oaData(assMatrix=m6AfS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed6, groupid="6", taskid="F", weights=PropFeed6)
ModelOADA6AfSFeed<-addFit(oadata=OADA6AfSFeed, asocialVar=c(1,2,3))
summary(ModelOADA6AfSFeed)

OADA6AfDOpen<-oaData(assMatrix=m6AfD, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AfDOpen<-addFit(oadata=OADA6AfDOpen, asocialVar=c(1,2,3))
summary(ModelOADA6AfDOpen)

OADA6AfSOpen<-oaData(assMatrix=m6AfS, asoc=asoc6, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen6, groupid="6", taskid="O", weights=PropOpen6)
ModelOADA6AfSOpen<-addFit(oadata=OADA6AfSOpen, asocialVar=c(1,2,3))
summary(ModelOADA6AfSOpen)



OADA7AfDFeed<-oaData(assMatrix=m7AfD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AfDFeed<-addFit(oadata=OADA7AfDFeed, asocialVar=c(1,2,3))
summary(ModelOADA7AfDFeed)

OADA7AfSFeed<-oaData(assMatrix=m7AfS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed7, groupid="7", taskid="F", weights=PropFeed7)
ModelOADA7AfSFeed<-addFit(oadata=OADA7AfSFeed, asocialVar=c(1,2,3))
summary(ModelOADA7AfSFeed)

OADA7AfDOpen<-oaData(assMatrix=m7AfD, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AfDOpen<-addFit(oadata=OADA7AfDOpen, asocialVar=c(1,2,3))
summary(ModelOADA7AfDOpen)

OADA7AfSOpen<-oaData(assMatrix=m7AfS, asoc=asoc7, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen7, groupid="7", taskid="O", weights=PropOpen7)
ModelOADA7AfSOpen<-addFit(oadata=OADA7AfSOpen, asocialVar=c(1,2,3))
summary(ModelOADA7AfSOpen)



OADA8AfDFeed<-oaData(assMatrix=m8AfD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AfDFeed<-addFit(oadata=OADA8AfDFeed, asocialVar=c(1,2,3))
summary(ModelOADA8AfDFeed)

OADA8AfSFeed<-oaData(assMatrix=m8AfS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed8, groupid="8", taskid="F", weights=PropFeed8)
ModelOADA8AfSFeed<-addFit(oadata=OADA8AfSFeed, asocialVar=c(1,2,3))
summary(ModelOADA8AfSFeed)

OADA8AfDOpen<-oaData(assMatrix=m8AfD, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AfDOpen<-addFit(oadata=OADA8AfDOpen, asocialVar=c(1,2,3))
summary(ModelOADA8AfDOpen)

OADA8AfSOpen<-oaData(assMatrix=m8AfS, asoc=asoc8, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen8, groupid="8", taskid="O", weights=PropOpen8)
ModelOADA8AfSOpen<-addFit(oadata=OADA8AfSOpen, asocialVar=c(1,2,3))
summary(ModelOADA8AfSOpen)



OADA9AfDFeed<-oaData(assMatrix=m9AfD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AfDFeed<-addFit(oadata=OADA9AfDFeed, asocialVar=c(1,2,3))
summary(ModelOADA9AfDFeed)

OADA9AfSFeed<-oaData(assMatrix=m9AfS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed9, groupid="9", taskid="F", weights=PropFeed9)
ModelOADA9AfSFeed<-addFit(oadata=OADA9AfSFeed, asocialVar=c(1,2,3))
summary(ModelOADA9AfSFeed)

OADA9AfDOpen<-oaData(assMatrix=m9AfD, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AfDOpen<-addFit(oadata=OADA9AfDOpen, asocialVar=c(1,2,3))
summary(ModelOADA9AfDOpen)

OADA9AfSOpen<-oaData(assMatrix=m9AfS, asoc=asoc9, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen9, groupid="9", taskid="O", weights=PropOpen9)
ModelOADA9AfSOpen<-addFit(oadata=OADA9AfSOpen, asocialVar=c(1,2,3))
summary(ModelOADA9AfSOpen)



OADA10AfDFeed<-oaData(assMatrix=m10AfD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AfDFeed<-addFit(oadata=OADA10AfDFeed, asocialVar=c(1,2,3))
summary(ModelOADA10AfDFeed)

OADA10AfSFeed<-oaData(assMatrix=m10AfS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAFeed10, groupid="10", taskid="F", weights=PropFeed10)
ModelOADA10AfSFeed<-addFit(oadata=OADA10AfSFeed, asocialVar=c(1,2,3))
summary(ModelOADA10AfSFeed)

OADA10AfDOpen<-oaData(assMatrix=m10AfD, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AfDOpen<-addFit(oadata=OADA10AfDOpen, asocialVar=c(1,2,3))
summary(ModelOADA10AfDOpen)

OADA10AfSOpen<-oaData(assMatrix=m10AfS, asoc=asoc10, ties=c(0,0,0,0,0,0,0), orderAcq=OOAOpen10, groupid="10", taskid="O", weights=PropOpen10)
ModelOADA10AfSOpen<-addFit(oadata=OADA10AfSOpen, asocialVar=c(1,2,3))
summary(ModelOADA10AfSOpen)



