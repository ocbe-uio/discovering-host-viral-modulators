if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)

node_attributes <- read_excel("~/Documents/discovering-host-viral-modulators/data/node_attributes.xlsx")
names(node_attributes)[names(node_attributes) == "Gene name"] <- "Gene_name"
names(node_attributes)[names(node_attributes) == "Target class"] <- "Target_class"
names(node_attributes)[names(node_attributes) == "Entrez ID"] <- "Entrez_ID"

NIP_nodes<-node_attributes %>% 
  filter(Target_class=='NIP') %>% 
  dplyr::select(Entrez_ID)

VIP_nodes<-node_attributes %>% 
  filter(Target_class=='VIP') %>% 
  dplyr::select(Entrez_ID)

NIP_nodes<-NIP_nodes$Entrez_ID
VIP_nodes<-VIP_nodes$Entrez_ID

######### GO Biological process #################################################
clusterContainer_VIP_NIP = list(VIP_nodes,NIP_nodes)
names(clusterContainer_VIP_NIP) <- c("VIP","NIP")

geneset_BP_VIP_NIP<-compareCluster(geneCluster = clusterContainer_VIP_NIP, fun = "enrichGO",OrgDb='org.Hs.eg.db', ont="BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",readable = TRUE)

Compre_BP<- setReadable(geneset_BP_VIP_NIP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

Compare_BP<-as.data.frame(Compre_BP@compareClusterResult)

dotplot(geneset_BP_VIP_NIP)

############## Independent analysis ##################




############################# Compare KEGG#########################################
geneset_kegg_VIP_NIP <- compareCluster(clusterContainer_VIP_NIP, fun="enrichKEGG",
                                  organism="hsa", pvalueCutoff=0.05)
dotplot(geneset_kegg_VIP_NIP)

compare_KEGG<-as.data.frame(geneset_kegg_VIP_NIP@compareClusterResult)

VIP_kegg <- enrichKEGG(gene   = VIP_nodes,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05
)
VIP_symbol<- setReadable(VIP_kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
VIP_kegg<-VIP_symbol@result

NIP_kegg <- enrichKEGG(gene   = NIP_nodes,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05
)
NIP_symbol<- setReadable(NIP_kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
NIP_kegg<-NIP_symbol@result

dotplot(NIP_kegg)

########################### For the plot as in the paper download the supplementray table S2 ########

Compare_BP_kegg <- read.csv("../data/GO_analysis.csv", sep=";")

all <- Compare_BP_kegg %>% 
  pivot_wider(
    names_from = Cluster, 
    values_from = "p.adjust", 
    id_cols = c("Description")
  )

all_new <- Compare_BP_kegg %>% 
  pivot_wider(
    names_from = Cluster, 
    values_from = "p.adjust", 
    id_cols = c("Description", "ID")
  ) %>% 
  mutate(id = case_when(
    grepl("GO", ID) ~ "Biological Process",
    grepl("hsa", ID) ~ "KEGG pathway"
  ))
common <- all_new %>% 
  dplyr::select(-ID) %>% 
  na.omit() %>% 
  pivot_longer(
    names_to = "Cluster", 
    values_to = "p.adjust", 
    cols = c(-Description, -id)
  )
################ Plot ############################
common %>%
  filter(p.adjust <= 0.05) %>%
  group_by(id) %>% 
  ggplot(aes(Cluster,
             reorder(Description, as.numeric(as.factor(id)) + p.adjust),
             fill = p.adjust)) +
  geom_tile(
    colour = "grey",
    linetype = 1,
    size = 0.1
  ) +
  theme_void(base_size = 12, base_family = "Arial") +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 12, hjust = 1),
    axis.title.x = element_text(size = 12),
    text = element_text(size = 12, family = "Arial"),
    axis.text = element_text(color = "black", size = 12)
  ) +
  scale_fill_viridis_c(option = "D") +
  #guide = guide_colorbar(barheight = unit(0.5, "npc"))) +
  scale_y_discrete(
    labels = function(x)
      stringr::str_to_sentence(x) %>% 
      stringr::str_wrap(width = 28)
  ) +
  labs(x = NULL, y = NULL, fill = "p-adjusted") +
  coord_equal() +
  annotate(
    geom = "tile",
    x = 2.75,
    y = 13.5,
    height = 4,
    width = 0.5,
    fill = "lightgrey"
  ) +
  annotate(
    geom = "tile",
    x = 2.75,
    y = 6.5,
    height = 12,
    width = 0.5,
    fill = "lightgrey"
  ) +
  annotate(
    geom = "segment",
    x = 2.5,
    xend = 3,
    y = 11.5,
    yend = 11.5,
    color = "white",
    size = 1
  ) +
  annotate(geom = "text", x = 2.75, y = 13.5, label = "KEGG pathway", size = 4,  angle = -90, hjust = 0.5, vjust = 0.5, family = "Arial") +
  annotate(geom = "text", x = 2.75, y = 6.5, label = "Biological process", size = 4,  angle = -90, hjust = 0.5, vjust = 0.5, family = "Arial") +
  expand_limits(x = c(0.5, 3))


