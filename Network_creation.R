######## Script to construct network 
library("igraph")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("reshape2")

# Read 
## 293T cellline PPI network
BioPlex_293T_Network<- read.delim("BioPlex_293T_Network.tsv")
 ###### SARS-COV-2 PPI data and extracting VIP
HEK293T_SARS_CoV_2_node <- read.csv("HEK293T_SARS-CoV-2 node.csv")
SARS_proteins<-HEK293T_SARS_CoV_2_node%>%
  filter(Bait_Boolean==0)%>%
  dplyr::select(shared.name)

##### Human proteins in network 
a<-data.frame(BioPlex_293T_Network$SymbolA)
b<-data.frame(BioPlex_293T_Network$SymbolB)
names(b)<-names(a)
PPI_protein<-rbind(a,b)%>%
  unique()
names(PPI_protein)[names(PPI_protein)=="BioPlex_293T_Network.SymbolA"] <- "id"


###### Creating the network
nodes<-PPI_protein
links<-BioPlex_293T_Network %>%
  dplyr::select(SymbolA,SymbolB)

##### renamimg columns
names(links)[names(links) == "SymbolA"] <- "From"
names(links)[names(links) == "SymbolB"] <- "To"  

### Intital network 
network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
#### Check if the graph is simple 
is_simple(network)
###### Creating the SARS-CoV-2 network 
SARSnodes <- SARS_proteins[['shared.name']]
V(network)$id<-as_ids(V(network))
selnodes <- V(network)[id %in% SARSnodes]
selegoV <- ego(network, order=1, nodes = selnodes, mode = "all", mindist = 0)
SARS_graph <- induced_subgraph(network,unlist(selegoV))
V(SARS_graph)$id = V(SARS_graph)$name
vcount(SARS_graph) ### number of nodes 
ecount(SARS_graph) ### number of edges 

#### calculating the topology 
deg_net <- degree(SARS_graph, mode="all")
betweeness_net<-betweenness(SARS_graph,normalized = TRUE)

######## Read created attribute file that's provided as supplementray table S1b 
SARS_network_nodes_attribute <- node_attributes <- read_excel("node_attributes.xlsx")
edge_list <- read_excel("edge_list.xlsx")


######## Create the network for analysis
links<-edge_list
nodes<-SARS_network_nodes_attribute  

#### renaming column 
names(nodes)[names(nodes) == "Target class"] <- "Target_class"


COVID_net <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
vcount(COVID_net)
ecount(COVID_net)
#############################Extracting giant component 
components <- igraph::clusters(COVID_net, mode="weak")
biggest_cluster_id <- which.max(components$csize)
# ids
vert_ids <- V(COVID_net)[components$membership == biggest_cluster_id]

# subgraph

gaint_component<-igraph::induced_subgraph(COVID_net, vert_ids)

########## Calculating the shortest path between NonVIP and VIP

NonVIPtoVIP_shortest<-distances(gaint_component, v=V(gaint_component)[Target_class==c("Non-VIP","NIP")], to=V(gaint_component)[which(Target_class=="VIP")], weights=NA,mode ="all")

avg_VIP_shoterst<-as.data.frame(rowMeans(NonVIPtoVIP_shortest,na.rm =TRUE))

mean(avg_VIP_shoterst[,1],na.rm = TRUE)

#  convert from matrix to data frame and preserve row names as column
NonVIPtoVIP_shortest <- data.frame(NonVIPs = row.names(NonVIPtoVIP_shortest), as.data.frame(NonVIPtoVIP_shortest), row.names = NULL)

# gather so in a tidy format for ease of use in ggplot2

NonVIPtoVIP_shortest<-gather(as.data.frame(NonVIPtoVIP_shortest),VIP_genes,value,-1)

###### Counting number of VIPS in different neigbhourhood

VIP_distance_count<-NonVIPtoVIP_shortest %>%
  group_by(value) %>%
  summarise(Total_interactions=n()) %>%
  dplyr::rename(Neighbourhood_size=value)

#################### For random Path############################################################

Random_VIP_shoterst<-distances(gaint_component,v=sample(V(gaint_component)[which(Target_class==c("Non-VIP","NIP"))],298), to=V(gaint_component)[which(Target_class=="VIP")], weights=NA,mode ="all")

Random_VIP_shoterst[!is.finite(Random_VIP_shoterst)]<-NA

avg_Random_VIP_shoterst<-as.data.frame(rowMeans(Random_VIP_shoterst,na.rm =TRUE))

mean(avg_Random_VIP_shoterst[,1],na.rm = TRUE)

#  convert from matrix to data frame and preserve row names as column
Random_VIP<- data.frame(NonVIPs = row.names(Random_VIP_shoterst), as.data.frame(Random_VIP_shoterst), row.names = NULL)

# gather so in a tidy format for ease of use in ggplot2
Random_VIP <- gather(as.data.frame(Random_VIP),VIP_genes,value,-1)

###### Counting number of VIPS in different neigbhourhood

random_distance_count<-Random_VIP %>%
  group_by(value) %>%
  summarise(Total_interactions=n()) %>%
  dplyr::rename(Neighbourhood_size=value)


####### Statistical tests for distance between VIP and random nodes

data_shortest_paths<-cbind(avg_VIP_shoterst$`rowMeans(NonVIPtoVIP_shortest, na.rm = TRUE)`,avg_Random_VIP_shoterst$`rowMeans(Random_VIP_shoterst, na.rm = TRUE)`)

data_shortest_paths<-as.data.frame(data_shortest_paths) %>% 
  dplyr::rename(SARS2_targeted=V1,
         Random_genes=V2)

data_shortest_paths<-as.data.frame(data_shortest_paths)

data_shortest_paths_melt<-melt(data_shortest_paths,variable.names="key",value.names="value")

distance_VIPrandom<-data_shortest_paths_melt %>% 
  group_by(variable) %>%
  summarise(
    count = n(),
    median = median(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE),
    mean=mean(value, na.rm = TRUE),
    sd_dist=sd(value),
    se_dist=sd(value)/sqrt(n()),
    test=wilcox.test(value~variable,data=.)$p.value)

###########

require(scales)

ggplot(VIP_distance_count,aes(x=as.factor(Neighbourhood_size),y=Total_interactions,fill=as.factor(Neighbourhood_size)))+geom_col(fill="#999999",width = 0.9)+theme(axis.title.y=element_text(face = "bold",size = 12),axis.title.x=element_text(face = "bold",size=10))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size=12,face = "bold"),plot.margin = unit(c(1,0.5,1.5,1.2),"cm"))+xlab("Neighbourhood\nsize")+ylab("Number of\ninteractions")+scale_y_continuous(expand = c(0,0),trans='log10',breaks = c(0,10,100,1000,10000,100000),labels = function(x) format(x, scientific = FALSE))

ggplot(mapping = aes(x, y))+geom_bar(data=data.frame(x=as.factor(VIP_distance_count$Neighbourhood_size),y=VIP_distance_count$Total_interactions),width=0.4,stat = 'identity',fill="grey7") +
  geom_bar(data=data.frame(x=as.factor(random_distance_count$Neighbourhood_size),y=random_distance_count$Total_interactions),width=0.4, stat = 'identity', fill = '#999999',alpha=0.8) +
  theme_classic() + scale_y_continuous(expand = c(0, 0))+theme(axis.title.y=element_text(size = 12),axis.title.x=element_text(size=12))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size=12,family="ArialMT"),axis.text = element_text(color = "black", size = 12))+xlab("Neighbourhood size")+ylab("Number of interactions")+scale_y_continuous(expand = c(0,0),trans='log10',breaks = c(0,10,100,1000,10000,100000),labels = function(x) format(x, scientific = FALSE))

##############  network Topology analysis

network_data_analysis<-nodes %>% mutate(Target_class=ifelse(Target_class=='NIP','Non-VIP',
                       ifelse(Target_class=='VIP','VIP','Non-VIP')))

deg_VIP<-network_data_analysis%>% 
  group_by(Target_class) %>%
  summarise(
    count = n(),
    median = median(Degree, na.rm = TRUE),
    IQR = IQR(Degree, na.rm = TRUE),
    mean=mean(Degree, na.rm = TRUE),
    sd_deg=sd(Degree),
    se_deg=sd(Degree)/sqrt(n()),
    test=wilcox.test(Degree~Target_class,data=.)$p.value)

deg<-ggplot(deg_VIP,aes(Target_class,mean,fill=Target_class))+ geom_bar(position=position_dodge(0.4),fill=c("#999999","#FF3300"),width = 0.4,stat="identity")+ geom_errorbar(aes(ymin = mean - se_deg, ymax = mean+ se_deg), width=0.1)+
  theme(axis.title.x = element_blank(),legend.title = element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size=12,family="ArialMT"),plot.margin = unit(c(1,0.1,1.5,1.1),"cm"),axis.text = element_text(color = "black", size = 12))+ylab("Connectivity (k)")+scale_y_continuous(expand = c(0,0) )+xlab("Node type")

deg1<-deg+geom_text(data=tibble(x=1.35,y=24),aes(x=x,y=y,label="p=0.003"),inherit.aes = FALSE,size=4)
ggsave("deg1.pdf",plot=deg1,width = 8,height = 7,units = "cm")

names(network_data_analysis)[names(network_data_analysis) == "Betweenness centrality"] <- "BetweennessCentrality"

bet_VIP<-network_data_analysis %>% 
  group_by(Target_class) %>%
  summarise(
    count = n(),
    median = median(BetweennessCentrality, na.rm = TRUE),
    IQR = IQR(BetweennessCentrality, na.rm = TRUE),
    mean=mean(BetweennessCentrality, na.rm = TRUE),
    sd_deg=sd(BetweennessCentrality),
    se_deg=sd(BetweennessCentrality)/sqrt(n()),
    test=wilcox.test(BetweennessCentrality~Target_class,data=.)$p.value)

bet<-ggplot(bet_VIP,aes(Target_class,mean,fill=Target_class))+ geom_bar(fill=c("#999999","#FF3300"),stat="identity",width = 0.4)+ geom_errorbar(aes(ymin = mean - se_deg, ymax = mean+ se_deg), width=0.1)+theme(axis.title.x = element_blank(),axis.title.y=element_text(size = 12))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size=12,family="ArialMT"),plot.margin = unit(c(1,0.1,1.5,1.1),"cm"),axis.text = element_text(color = "black", size = 12))+ylab("Betweenness\ncentrality (b)")+scale_y_continuous(expand = c(0,0))+xlab("Node type")

bet1<-bet+geom_text(data=tibble(x=1.35,y=0.00128),aes(x=x,y=y,label="p=3.09e-11"),size=4,inherit.aes = FALSE)

