##### Code for creation of figure 2: #####
# social networks based on spatial association,
# affiliative interactions, and agonistic
# interactions of each of the stable groups (4-10)
# of T. temporalis

##### Packages #####
# Load necessary library packages

library(statnet)
library(igraph)
library(asnipe)
library(tnet)
library(assortnet)
library(remotes)
library(NBDA)


##### Produce social network matrices from 
          # spatial association data #####

group4As <- read.csv("Group4As.csv")
rownames(group4As)<-group4As$X
group4As$X<-NULL
m4As<-as.matrix(group4As, header=T, row.names=1)
group4As_network<-as.network(m4As, data_format = "SP", association_index = "HWI")
graph4As<- graph_from_adjacency_matrix(m4As, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group5As <- read.csv("group5As.csv")
rownames(group5As)<-group5As$X
group5As$X<-NULL
m5As<-as.matrix(group5As, header=T, row.names=1)
group5As_network<-as.network(m5As, data_format = "SP", association_index = "HWI")
graph5As<- graph_from_adjacency_matrix(m5As, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group6As <- read.csv("group6As.csv")
rownames(group6As)<-group6As$X
group6As$X<-NULL
m6As<-as.matrix(group6As, header=T, row.names=1)
group6As_network<-as.network(m6As, data_format = "SP", association_index = "HWI")
graph6As<- graph_from_adjacency_matrix(m6As, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group7As <- read.csv("group7As.csv")
rownames(group7As)<-group7As$X
group7As$X<-NULL
m7As<-as.matrix(group7As, header=T, row.names=1)
group7As_network<-as.network(m7As, data_format = "SP", association_index = "HWI")
graph7As<- graph_from_adjacency_matrix(m7As, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group8As <- read.csv("group8As.csv")
rownames(group8As)<-group8As$X
group8As$X<-NULL
m8As<-as.matrix(group8As, header=T, row.names=1)
group8As_network<-as.network(m8As, data_format = "SP", association_index = "HWI")
graph8As<- graph_from_adjacency_matrix(m8As, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group9As <- read.csv("group9As.csv")
rownames(group9As)<-group9As$X
group9As$X<-NULL
m9As<-as.matrix(group9As, header=T, row.names=1)
group9As_network<-as.network(m9As, data_format = "SP", association_index = "HWI")
graph9As<- graph_from_adjacency_matrix(m9As, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group10As <- read.csv("group10As.csv")
rownames(group10As)<-group10As$X
group10As$X<-NULL
m10As<-as.matrix(group10As, header=T, row.names=1)
group10As_network<-as.network(m10As, data_format = "SP", association_index = "HWI")
graph10As<- graph_from_adjacency_matrix(m10As, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)




##### Calculate associative network metrics #####

graph4As_degree <- degree(graph4As, mode="all", loops=F)
graph4As_weighted_degree <- graph.strength(graph4As)
graph4As_weighted_degree_standardised <- graph4As_weighted_degree/max(graph4As_weighted_degree)
graph4As_trans<- transitivity(graph4As)
graph4As_closeness <- closeness(graph4As)
graph4As_betweenness <- betweenness(graph4As)
graph4As_density<-edge_density(graph4As)
graph4As_assortment<-assortment.continuous(m4As, weighted=T, vertex_values = E(graph4As))
t4As<- as.tnet(m4As)
graph4As_diameter<-max(distance_w(t4As), na.rm=T)
graph4As_clust_global<- clustering_w(t4As, measure = "am")

graph5As_degree <- degree(graph5As, mode="all", loops=F)
graph5As_weighted_degree <- graph.strength(graph5As)
graph5As_weighted_degree_standardised <- graph5As_weighted_degree/max(graph5As_weighted_degree)
graph5As_trans<- transitivity(graph5As)
graph5As_closeness <- closeness(graph5As)
graph5As_betweenness <- betweenness(graph5As)
graph5As_density<-edge_density(graph5As)
graph5As_assortment<-assortment.continuous(m5As, weighted=T, vertex_values = E(graph5As))
t5As<- as.tnet(m5As)
graph5As_diameter<-max(distance_w(t5As), na.rm=T)
graph5As_clust_global<- clustering_w(t5As, measure = "am")

graph6As_degree <- degree(graph6As, mode="all", loops=F)
graph6As_weighted_degree <- graph.strength(graph6As)
graph6As_weighted_degree_standardised <- graph6As_weighted_degree/max(graph6As_weighted_degree)
graph6As_trans<- transitivity(graph6As)
graph6As_closeness <- closeness(graph6As)
graph6As_betweenness <- betweenness(graph6As)
graph6As_density<-edge_density(graph6As)
graph6As_assortment<-assortment.continuous(m6As, weighted=T, vertex_values = E(graph6As))
t6As<- as.tnet(m6As)
graph6As_diameter<-max(distance_w(t6As), na.rm=T)
graph6As_clust_global<- clustering_w(t6As, measure = "am")

graph7As_degree <- degree(graph7As, mode="all", loops=F)
graph7As_weighted_degree <- graph.strength(graph7As)
graph7As_weighted_degree_standardised <- graph7As_weighted_degree/max(graph7As_weighted_degree)
graph7As_trans<- transitivity(graph7As)
graph7As_closeness <- closeness(graph7As)
graph7As_betweenness <- betweenness(graph7As)
graph7As_density<-edge_density(graph7As)
graph7As_assortment<-assortment.continuous(m7As, weighted=T, vertex_values = E(graph7As))
t7As<- as.tnet(m7As)
graph7As_diameter<-max(distance_w(t7As), na.rm=T)
graph7As_clust_global<- clustering_w(t7As, measure = "am")

graph8As_degree <- degree(graph8As, mode="all", loops=F)
graph8As_weighted_degree <- graph.strength(graph8As)
graph8As_weighted_degree_standardised <- graph8As_weighted_degree/max(graph8As_weighted_degree)
graph8As_trans<- transitivity(graph8As)
graph8As_closeness <- closeness(graph8As)
graph8As_betweenness <- betweenness(graph8As)
graph8As_density<-edge_density(graph8As)
graph8As_assortment<-assortment.continuous(m8As, weighted=T, vertex_values = E(graph8As))
t8As<- as.tnet(m8As)
graph8As_diameter<-max(distance_w(t8As), na.rm=T)
graph8As_clust_global<- clustering_w(t8As, measure = "am")

graph9As_degree <- degree(graph9As, mode="all", loops=F)
graph9As_weighted_degree <- graph.strength(graph9As)
graph9As_weighted_degree_standardised <- graph9As_weighted_degree/max(graph9As_weighted_degree)
graph9As_trans<- transitivity(graph9As)
graph9As_closeness <- closeness(graph9As)
graph9As_betweenness <- betweenness(graph9As)
graph9As_density<-edge_density(graph9As)
graph9As_assortment<-assortment.continuous(m9As, weighted=T, vertex_values = E(graph9As))
t9As<- as.tnet(m9As)
graph9As_diameter<-max(distance_w(t9As), na.rm=T)
graph9As_clust_global<- clustering_w(t9As, measure = "am")

graph10As_degree <- degree(graph10As, mode="all", loops=F)
graph10As_weighted_degree <- graph.strength(graph10As)
graph10As_weighted_degree_standardised <- graph10As_weighted_degree/max(graph10As_weighted_degree)
graph10As_trans<- transitivity(graph10As)
graph10As_closeness <- closeness(graph10As)
graph10As_betweenness <- betweenness(graph10As)
graph10As_density<-edge_density(graph10As)
graph10As_assortment<-assortment.continuous(m10As, weighted=T, vertex_values = E(graph10As))
t10As<- as.tnet(m10As)
graph10As_diameter<-max(distance_w(t10As), na.rm=T)
graph10As_clust_global<- clustering_w(t10As, measure = "am")

##### Produce social network matrices 
          # from agonistic interactions data #####

group4Ag <- read.csv("group4Ag.csv")
rownames(group4Ag)<-group4Ag$X
group4Ag$X<-NULL
m4Ag<-as.matrix(group4Ag, header=T, row.names=1)
group4Ag_network<-as.network(m4Ag, data_format = "SP", association_index = "HWI")
graph4Ag<- graph_from_adjacency_matrix(m4Ag, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group5Ag <- read.csv("group5Ag.csv")
rownames(group5Ag)<-group5Ag$X
group5Ag$X<-NULL
m5Ag<-as.matrix(group5Ag, header=T, row.names=1)
group5Ag_network<-as.network(m5Ag, data_format = "SP", association_index = "HWI")
graph5Ag<- graph_from_adjacency_matrix(m5Ag, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group6Ag <- read.csv("group6Ag.csv")
rownames(group6Ag)<-group6Ag$X
group6Ag$X<-NULL
m6Ag<-as.matrix(group6Ag, header=T, row.names=1)
group6Ag_network<-as.network(m6Ag, data_format = "SP", association_index = "HWI")
graph6Ag<- graph_from_adjacency_matrix(m6Ag, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group7Ag <- read.csv("group7Ag.csv")
rownames(group7Ag)<-group7Ag$X
group7Ag$X<-NULL
m7Ag<-as.matrix(group7Ag, header=T, row.names=1)
group7Ag_network<-as.network(m7Ag, data_format = "SP", association_index = "HWI")
graph7Ag<- graph_from_adjacency_matrix(m7Ag, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group8Ag <- read.csv("group8Ag.csv")
rownames(group8Ag)<-group8Ag$X
group8Ag$X<-NULL
m8Ag<-as.matrix(group8Ag, header=T, row.names=1)
group8Ag_network<-as.network(m8Ag, data_format = "SP", association_index = "HWI")
graph8Ag<- graph_from_adjacency_matrix(m8Ag, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group9Ag <- read.csv("group9Ag.csv")
rownames(group9Ag)<-group9Ag$X
group9Ag$X<-NULL
m9Ag<-as.matrix(group9Ag, header=T, row.names=1)
group9Ag_network<-as.network(m9Ag, data_format = "SP", association_index = "HWI")
graph9Ag<- graph_from_adjacency_matrix(m9Ag, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group10Ag <- read.csv("group10Ag.csv")
rownames(group10Ag)<-group10Ag$X
group10Ag$X<-NULL
m10Ag<-as.matrix(group10Ag, header=T, row.names=1)
group10Ag_network<-as.network(m10Ag, data_format = "SP", association_index = "HWI")
graph10Ag<- graph_from_adjacency_matrix(m10Ag, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)


##### Calculate agonistic network metrics #####

graph4Ag_degree <- degree(graph4Ag, mode="all", loops=F)
graph4Ag_weighted_degree <- graph.strength(graph4Ag)
graph4Ag_weighted_degree_standardised <- graph4Ag_weighted_degree/max(graph4Ag_weighted_degree)
graph4Ag_trans<- transitivity(graph4Ag)
graph4Ag_closeness <- closeness(graph4Ag)
graph4Ag_betweenness <- betweenness(graph4Ag)
graph4Ag_density<-edge_density(graph4Ag)
graph4Ag_assortment<-assortment.continuous(m4Ag, weighted=T, vertex_values = E(graph4Ag))
t4Ag<- as.tnet(m4Ag)
graph4Ag_diameter<-max(distance_w(t4Ag), na.rm=T)
graph4Ag_clust_global<- clustering_w(t4Ag, measure = "am")

graph5Ag_degree <- degree(graph5Ag, mode="all", loops=F)
graph5Ag_weighted_degree <- graph.strength(graph5Ag)
graph5Ag_weighted_degree_standardised <- graph5Ag_weighted_degree/max(graph5Ag_weighted_degree)
graph5Ag_trans<- transitivity(graph5Ag)
graph5Ag_closeness <- closeness(graph5Ag)
graph5Ag_betweenness <- betweenness(graph5Ag)
graph5Ag_density<-edge_density(graph5Ag)
graph5Ag_assortment<-assortment.continuous(m5Ag, weighted=T, vertex_values = E(graph5Ag))
t5Ag<- as.tnet(m5Ag)
graph5Ag_diameter<-max(distance_w(t5Ag), na.rm=T)
graph5Ag_clust_global<- clustering_w(t5Ag, measure = "am")

graph6Ag_degree <- degree(graph6Ag, mode="all", loops=F)
graph6Ag_weighted_degree <- graph.strength(graph6Ag)
graph6Ag_weighted_degree_standardised <- graph6Ag_weighted_degree/max(graph6Ag_weighted_degree)
graph6Ag_trans<- transitivity(graph6Ag)
graph6Ag_closeness <- closeness(graph6Ag)
graph6Ag_betweenness <- betweenness(graph6Ag)
graph6Ag_density<-edge_density(graph6Ag)
graph6Ag_assortment<-assortment.continuous(m6Ag, weighted=T, vertex_values = E(graph6Ag))
t6Ag<- as.tnet(m6Ag)
graph6Ag_diameter<-max(distance_w(t6Ag), na.rm=T)
graph6Ag_clust_global<- clustering_w(t6Ag, measure = "am")

graph7Ag_degree <- degree(graph7Ag, mode="all", loops=F)
graph7Ag_weighted_degree <- graph.strength(graph7Ag)
graph7Ag_weighted_degree_standardised <- graph7Ag_weighted_degree/max(graph7Ag_weighted_degree)
graph7Ag_trans<- transitivity(graph7Ag)
graph7Ag_closeness <- closeness(graph7Ag)
graph7Ag_betweenness <- betweenness(graph7Ag)
graph7Ag_density<-edge_density(graph7Ag)
graph7Ag_assortment<-assortment.continuous(m7Ag, weighted=T, vertex_values = E(graph7Ag))
t7Ag<- as.tnet(m7Ag)
graph7Ag_diameter<-max(distance_w(t7Ag), na.rm=T)
graph7Ag_clust_global<- clustering_w(t7Ag, measure = "am")

graph8Ag_degree <- degree(graph8Ag, mode="all", loops=F)
graph8Ag_weighted_degree <- graph.strength(graph8Ag)
graph8Ag_weighted_degree_standardised <- graph8Ag_weighted_degree/max(graph8Ag_weighted_degree)
graph8Ag_trans<- transitivity(graph8Ag)
graph8Ag_closeness <- closeness(graph8Ag)
graph8Ag_betweenness <- betweenness(graph8Ag)
graph8Ag_density<-edge_density(graph8Ag)
graph8Ag_assortment<-assortment.continuous(m8Ag, weighted=T, vertex_values = E(graph8Ag))
t8Ag<- as.tnet(m8Ag)
graph8Ag_diameter<-max(distance_w(t8Ag), na.rm=T)
graph8Ag_clust_global<- clustering_w(t8Ag, measure = "am")

graph9Ag_degree <- degree(graph9Ag, mode="all", loops=F)
graph9Ag_weighted_degree <- graph.strength(graph9Ag)
graph9Ag_weighted_degree_standardised <- graph9Ag_weighted_degree/max(graph9Ag_weighted_degree)
graph9Ag_trans<- transitivity(graph9Ag)
graph9Ag_closeness <- closeness(graph9Ag)
graph9Ag_betweenness <- betweenness(graph9Ag)
graph9Ag_density<-edge_density(graph9Ag)
graph9Ag_assortment<-assortment.continuous(m9Ag, weighted=T, vertex_values = E(graph9Ag))
t9Ag<- as.tnet(m9Ag)
graph9Ag_diameter<-max(distance_w(t9Ag), na.rm=T)
graph9Ag_clust_global<- clustering_w(t9Ag, measure = "am")

graph10Ag_degree <- degree(graph10Ag, mode="all", loops=F)
graph10Ag_weighted_degree <- graph.strength(graph10Ag)
graph10Ag_weighted_degree_standardised <- graph10Ag_weighted_degree/max(graph10Ag_weighted_degree)
graph10Ag_trans<- transitivity(graph10Ag)
graph10Ag_closeness <- closeness(graph10Ag)
graph10Ag_betweenness <- betweenness(graph10Ag)
graph10Ag_density<-edge_density(graph10Ag)
graph10Ag_assortment<-assortment.continuous(m10Ag, weighted=T, vertex_values = E(graph10Ag))
t10Ag<- as.tnet(m10Ag)
graph10Ag_diameter<-max(distance_w(t10Ag), na.rm=T)
graph10Ag_clust_global<- clustering_w(t10Ag, measure = "am")


##### Produce social network matrices from 
          # affiliative interactions data #####

group4Af <- read.csv("group4Af.csv")
rownames(group4Af)<-group4Af$X
group4Af$X<-NULL
m4Af<-as.matrix(group4Af, header=T, row.names=1)
group4Af_network<-as.network(m4Af, data_format = "SP", association_index = "HWI")
graph4Af<- graph_from_adjacency_matrix(m4Af, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group5Af <- read.csv("group5Af.csv")
rownames(group5Af)<-group5Af$X
group5Af$X<-NULL
m5Af<-as.matrix(group5Af, header=T, row.names=1)
group5Af_network<-as.network(m5Af, data_format = "SP", association_index = "HWI")
graph5Af<- graph_from_adjacency_matrix(m5Af, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group6Af <- read.csv("group6Af.csv")
rownames(group6Af)<-group6Af$X
group6Af$X<-NULL
m6Af<-as.matrix(group6Af, header=T, row.names=1)
group6Af_network<-as.network(m6Af, data_format = "SP", association_index = "HWI")
graph6Af<- graph_from_adjacency_matrix(m6Af, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group7Af <- read.csv("group7Af.csv")
rownames(group7Af)<-group7Af$X
group7Af$X<-NULL
m7Af<-as.matrix(group7Af, header=T, row.names=1)
group7Af_network<-as.network(m7Af, data_format = "SP", association_index = "HWI")
graph7Af<- graph_from_adjacency_matrix(m7Af, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group8Af <- read.csv("group8Af.csv")
rownames(group8Af)<-group8Af$X
group8Af$X<-NULL
m8Af<-as.matrix(group8Af, header=T, row.names=1)
group8Af_network<-as.network(m8Af, data_format = "SP", association_index = "HWI")
graph8Af<- graph_from_adjacency_matrix(m8Af, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group9Af <- read.csv("group9Af.csv")
rownames(group9Af)<-group9Af$X
group9Af$X<-NULL
m9Af<-as.matrix(group9Af, header=T, row.names=1)
group9Af_network<-as.network(m9Af, data_format = "SP", association_index = "HWI")
graph9Af<- graph_from_adjacency_matrix(m9Af, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)

group10Af <- read.csv("group10Af.csv")
rownames(group10Af)<-group10Af$X
group10Af$X<-NULL
m10Af<-as.matrix(group10Af, header=T, row.names=1)
group10Af_network<-as.network(m10Af, data_format = "SP", association_index = "HWI")
graph10Af<- graph_from_adjacency_matrix(m10Af, mode = "directed", weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NULL)




# Calculate affiliative network metrics




graph4Af_degree <- degree(graph4Af, mode="all", loops=F)
graph4Af_weighted_degree <- graph.strength(graph4Af)
graph4Af_weighted_degree_standardised <- graph4Af_weighted_degree/max(graph4Af_weighted_degree)
graph4Af_trans<- transitivity(graph4Af)
graph4Af_closeness <- closeness(graph4Af)
graph4Af_betweenness <- betweenness(graph4Af)
graph4Af_density<-edge_density(graph4Af)
graph4Af_assortment<-assortment.continuous(m4Af, weighted=T, vertex_values = E(graph4Af))
t4Af<- as.tnet(m4Af)
graph4Af_diameter<-max(distance_w(t4Af), na.rm=T)
graph4Af_clust_global<- clustering_w(t4Af, measure = "am")

graph5Af_degree <- degree(graph5Af, mode="all", loops=F)
graph5Af_weighted_degree <- graph.strength(graph5Af)
graph5Af_weighted_degree_standardised <- graph5Af_weighted_degree/max(graph5Af_weighted_degree)
graph5Af_trans<- transitivity(graph5Af)
graph5Af_closeness <- closeness(graph5Af)
graph5Af_betweenness <- betweenness(graph5Af)
graph5Af_density<-edge_density(graph5Af)
graph5Af_assortment<-assortment.continuous(m5Af, weighted=T, vertex_values = E(graph5Af))
t5Af<- as.tnet(m5Af)
graph5Af_diameter<-max(distance_w(t5Af), na.rm=T)
graph5Af_clust_global<- clustering_w(t5Af, measure = "am")

graph6Af_degree <- degree(graph6Af, mode="all", loops=F)
graph6Af_weighted_degree <- graph.strength(graph6Af)
graph6Af_weighted_degree_standardised <- graph6Af_weighted_degree/max(graph6Af_weighted_degree)
graph6Af_trans<- transitivity(graph6Af)
graph6Af_closeness <- closeness(graph6Af)
graph6Af_betweenness <- betweenness(graph6Af)
graph6Af_density<-edge_density(graph6Af)
graph6Af_assortment<-assortment.continuous(m6Af, weighted=T, vertex_values = E(graph6Af))
t6Af<- as.tnet(m6Af)
graph6Af_diameter<-max(distance_w(t6Af), na.rm=T)
graph6Af_clust_global<- clustering_w(t6Af, measure = "am")

graph7Af_degree <- degree(graph7Af, mode="all", loops=F)
graph7Af_weighted_degree <- graph.strength(graph7Af)
graph7Af_weighted_degree_standardised <- graph7Af_weighted_degree/max(graph7Af_weighted_degree)
graph7Af_trans<- transitivity(graph7Af)
graph7Af_closeness <- closeness(graph7Af)
graph7Af_betweenness <- betweenness(graph7Af)
graph7Af_density<-edge_density(graph7Af)
graph7Af_assortment<-assortment.continuous(m7Af, weighted=T, vertex_values = E(graph7Af))
t7Af<- as.tnet(m7Af)
graph7Af_diameter<-max(distance_w(t7Af), na.rm=T)
graph7Af_clust_global<- clustering_w(t7Af, measure = "am")

graph8Af_degree <- degree(graph8Af, mode="all", loops=F)
graph8Af_weighted_degree <- graph.strength(graph8Af)
graph8Af_weighted_degree_standardised <- graph8Af_weighted_degree/max(graph8Af_weighted_degree)
graph8Af_trans<- transitivity(graph8Af)
graph8Af_closeness <- closeness(graph8Af)
graph8Af_betweenness <- betweenness(graph8Af)
graph8Af_density<-edge_density(graph8Af)
graph8Af_assortment<-assortment.continuous(m8Af, weighted=T, vertex_values = E(graph8Af))
t8Af<- as.tnet(m8Af)
graph8Af_diameter<-max(distance_w(t8Af), na.rm=T)
graph8Af_clust_global<- clustering_w(t8Af, measure = "am")

graph9Af_degree <- degree(graph9Af, mode="all", loops=F)
graph9Af_weighted_degree <- graph.strength(graph9Af)
graph9Af_weighted_degree_standardised <- graph9Af_weighted_degree/max(graph9Af_weighted_degree)
graph9Af_trans<- transitivity(graph9Af)
graph9Af_closeness <- closeness(graph9Af)
graph9Af_betweenness <- betweenness(graph9Af)
graph9Af_density<-edge_density(graph9Af)
graph9Af_assortment<-assortment.continuous(m9Af, weighted=T, vertex_values = E(graph9Af))
t9Af<- as.tnet(m9Af)
graph9Af_diameter<-max(distance_w(t9Af), na.rm=T)
graph9Af_clust_global<- clustering_w(t9Af, measure = "am")

graph10Af_degree <- degree(graph10Af, mode="all", loops=F)
graph10Af_weighted_degree <- graph.strength(graph10Af)
graph10Af_weighted_degree_standardised <- graph10Af_weighted_degree/max(graph10Af_weighted_degree)
graph10Af_trans<- transitivity(graph10Af)
graph10Af_closeness <- closeness(graph10Af)
graph10Af_betweenness <- betweenness(graph10Af)
graph10Af_density<-edge_density(graph10Af)
graph10Af_assortment<-assortment.continuous(m10Af, weighted=T, vertex_values = E(graph10Af))
t10Af<- as.tnet(m10Af)
graph10Af_diameter<-max(distance_w(t10Af), na.rm=T)
graph10Af_clust_global<- clustering_w(t10Af, measure = "am")




# Establish layout and colour palette, and assign weighted degree and dominance for the production of igraph social network figures


l <- layout_in_circle

dom4<-cbind(0.299864007, 0, 1, 0.546237536, 0.805757029)
dom5<-cbind(0.332092499, 0, 1, 0.548770012, 0.860254241)
dom6<-cbind(0, 0.641118, 0.450927, 0.641118, 1)
dom7<-cbind(0, 0.387507378, 0.883806974, 0.464031974, 1)
dom8<-cbind(0, 0.931507, 0.534247, 1)
dom9<-cbind(0.49084997, 0.180601, 0.482329, 0.288331, 0, 0.696928, 1)
dom10<-cbind(0, 0.71954181, 0.65411764, 0.263529423, 0.24423529, 1, 0.375764715, 0.803046458)

my_resolution = 100
my_paletteAs    = colorRampPalette(c('#FCB59B', '#FF8C00'))
my_max4    = max(dom4, na.rm=TRUE)
my_vector4 = dom4 / my_max4
my_colors4As = my_paletteAs(my_resolution)[as.numeric(cut(my_vector4, breaks=my_resolution))]

my_max5    = max(dom5, na.rm=TRUE)
my_vector5 = dom5 / my_max5
my_colors5As = my_paletteAs(my_resolution)[as.numeric(cut(my_vector5, breaks=my_resolution))]

my_max6    = max(dom6, na.rm=TRUE)
my_vector6 = dom6 / my_max6
my_colors6As = my_paletteAs(my_resolution)[as.numeric(cut(my_vector6, breaks=my_resolution))]

my_max7    = max(dom7, na.rm=TRUE)
my_vector7 = dom7 / my_max7
my_colors7As = my_paletteAs(my_resolution)[as.numeric(cut(my_vector7, breaks=my_resolution))]

my_max8    = max(dom8, na.rm=TRUE)
my_vector8 = dom8 / my_max8
my_colors8As = my_paletteAs(my_resolution)[as.numeric(cut(my_vector8, breaks=my_resolution))]

my_max9    = max(dom9, na.rm=TRUE)
my_vector9 = dom9 / my_max9
my_colors9As = my_paletteAs(my_resolution)[as.numeric(cut(my_vector9, breaks=my_resolution))]

my_max10    = max(dom10, na.rm=TRUE)
my_vector10 = dom10 / my_max10
my_colors10As = my_paletteAs(my_resolution)[as.numeric(cut(my_vector10, breaks=my_resolution))]


WD4As<-cbind(0.715994, 0.3692078, 0.9730942, 0.7802691, 1)
WD5As<-cbind(0.5410959, 0.5684932, 0.9060665, 0.6144814, 1)
WD6As<-cbind(0.5931198, 0.6453144, 0.5670225, 1, 0.886121)
WD7As<-cbind(0.5149834, 0.4506104, 1, 0.3651498, 0.945616)
WD8As<-cbind(0.7267176, 0.9160305, 0.9480916, 1)
WD9As<-cbind(0.9024927, 0.9464809, 1, 0.8049853, 0.9347507, 0.936217, 0.9178886)
WD10As<-cbind(0.7086331, 0.5671463, 0.9820144, 0.6282974, 0.8645084, 1, 0.7422062, 0.6438849)


my_paletteAg    = colorRampPalette(c("#85D0FF", "#0320FC"))
my_colors4Ag = my_paletteAs(my_resolution)[as.numeric(cut(my_vector4, breaks=my_resolution))]
my_colors5Ag = my_paletteAs(my_resolution)[as.numeric(cut(my_vector5, breaks=my_resolution))]
my_colors6Ag = my_paletteAs(my_resolution)[as.numeric(cut(my_vector6, breaks=my_resolution))]
my_colors7Ag = my_paletteAs(my_resolution)[as.numeric(cut(my_vector7, breaks=my_resolution))]
my_colors8Ag = my_paletteAs(my_resolution)[as.numeric(cut(my_vector8, breaks=my_resolution))]
my_colors9Ag = my_paletteAs(my_resolution)[as.numeric(cut(my_vector9, breaks=my_resolution))]
my_colors10Ag = my_paletteAs(my_resolution)[as.numeric(cut(my_vector10, breaks=my_resolution))]


WD4Ag<-cbind(0.5, 0.276316, 1, 0.407895, 0.631579)
WD5Ag<-cbind(0.835616, 0.534247, 0.835616, 1, 0.547945)
WD6Ag<-cbind(0.859649, 1, 0.473684, 0.614035, 0.877193)
WD7Ag<-cbind(0.602041, 0.602041, 1, 0.520408, 0.622449)
WD8Ag<-cbind(0.305556, 0.388889, 1, 0.805556)
WD9Ag<-cbind(0.609329, 0.664723, 0.673469, 0.548105, 0.571429, 0.580175, 1)
WD10Ag<-cbind(0.328125, 1, 0.421875, 0.28125, 0.328125, 0.828125, 0.328125, 0.703125)

my_paletteAf    = colorRampPalette(c("#FFBAFD", "#D119B8"))
                                     my_colors4Af = my_paletteAs(my_resolution)[as.numeric(cut(my_vector4, breaks=my_resolution))]
                                     my_colors5Af = my_paletteAs(my_resolution)[as.numeric(cut(my_vector5, breaks=my_resolution))]
                                     my_colors6Af = my_paletteAs(my_resolution)[as.numeric(cut(my_vector6, breaks=my_resolution))]
                                     my_colors7Af = my_paletteAs(my_resolution)[as.numeric(cut(my_vector7, breaks=my_resolution))]
                                     my_colors8Af = my_paletteAs(my_resolution)[as.numeric(cut(my_vector8, breaks=my_resolution))]
                                     my_colors9Af = my_paletteAs(my_resolution)[as.numeric(cut(my_vector9, breaks=my_resolution))]
                                     my_colors10Af = my_paletteAs(my_resolution)[as.numeric(cut(my_vector10, breaks=my_resolution))]
                                     
                                     
WD4Af<-cbind(0.416667, 0.166667, 0.25, 0.833333, 1)
WD5Af<-cbind(0.230769, 0.128205, 0.948718, 0.153846, 1)
WD6Af<-cbind(0.166667, 0.1, 0.1, 0.966667, 1)
WD7Af<-cbind(0.058824, 0.039216, 0.960784, 0.019608, 1)
WD8Af<-cbind(0.01, 1, 0.181818, 0.818182)
WD9Af<-cbind(0.888889, 0.888889, 0.888889, 0.925926, 1, 0.888889, 0.962963)
WD10Af<-cbind(0.01, 0.043478, 0.869565, 0.01, 0.086957, 1, 0.130435, 0.130435)
                                     
                                     
                                     
                                     
##### Produce associative igraph network images #####
                                     
                                     
                                     
                                     
plot.igraph(graph4As, edge.width = (log(E(graph4As)$weight)),  vertex.size=(WD4As*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors4As, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors4As)
                                     
plot.igraph(graph5As, edge.width = (log(E(graph5As)$weight)),  vertex.size=(WD5As*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors5As, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors5As)
                                     
plot.igraph(graph6As, edge.width = (log(E(graph6As)$weight)),  vertex.size=(WD6As*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors6As, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors6As)
                                     
plot.igraph(graph7As, edge.width = (log(E(graph7As)$weight)),  vertex.size=(WD7As*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors7As, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors7As)
                                     
plot.igraph(graph8As, edge.width = (log(E(graph8As)$weight)),  vertex.size=(WD8As*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors8As, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors8As)
                                     
plot.igraph(graph9As, edge.width = (log(E(graph9As)$weight)),  vertex.size=(WD9As*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors9As, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors9As)
                                     
plot.igraph(graph10As, edge.width = (log(E(graph10As)$weight)),  vertex.size=(WD10As*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors10As, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors10As)
                                     
                                     
                                     
                                     
##### Produce agonistic igraph network images #####
                                     
                                     
                                     
plot.igraph(graph4Ag, edge.width = (log(E(graph4Ag)$weight)),  vertex.size=(WD4Ag*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors4Ag, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors4Ag)
                                     
plot.igraph(graph5Ag, edge.width = (log(E(graph5Ag)$weight)),  vertex.size=(WD5Ag*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors5Ag, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors5Ag)
                                     
plot.igraph(graph6Ag, edge.width = (log(E(graph6Ag)$weight)),  vertex.size=(WD6Ag*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors6Ag, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors6Ag)
                                     
plot.igraph(graph7Ag, edge.width = (log(E(graph7Ag)$weight)),  vertex.size=(WD7Ag*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors7Ag, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors7Ag)
                                     
plot.igraph(graph8Ag, edge.width = (log(E(graph8Ag)$weight)),  vertex.size=(WD8Ag*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors8Ag, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors8Ag)
                                     
plot.igraph(graph9Ag, edge.width = (log(E(graph9Ag)$weight)),  vertex.size=(WD9Ag*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors9Ag, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors9Ag)
                                     
plot.igraph(graph10Ag, edge.width = (log(E(graph10Ag)$weight)),  vertex.size=(WD10Ag*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors10Ag, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors10Ag)
                                     
                                     
                                     
                                     
##### Produce affiliative igraph network images #####
                                     
                                     
plot.igraph(graph4Af, edge.width = (log(E(graph4Af)$weight)),  vertex.size=(WD4Af*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors4Af, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors4Af)
                                     
plot.igraph(graph5Af, edge.width = (log(E(graph5Af)$weight)),  vertex.size=(WD5Af*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors5Af, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors5Af)
                                     
                                     plot.igraph(graph6Af, edge.width = (log(E(graph6Af)$weight)),  vertex.size=(WD6Af*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors6Af, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors6Af)
                                     
plot.igraph(graph7Af, edge.width = (log(E(graph7Af)$weight)),  vertex.size=(WD7Af*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors7Af, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors7Af)
                                     
plot.igraph(graph8Af, edge.width = (log(E(graph8Af)$weight)),  vertex.size=(WD8Af*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors8Af, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors8Af)
                                     
plot.igraph(graph9Af, edge.width = (log(E(graph9Af)$weight)),  vertex.size=(WD9Af*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors9Af, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors9Af)
                                     
plot.igraph(graph10Af, edge.width = (log(E(graph10Af)$weight)),  vertex.size=(WD10Af*30), vertex.label.dist=4, edge.arrow.size=0.45, edge.curved=.2, layout=l, vertex.color=my_colors10Af, vertex.label.color="black", vertex.label.cex=1, vertex.frame.color=my_colors10Af)





