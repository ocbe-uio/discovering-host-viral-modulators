#### expression data 
install.packages(HPAAnalyze)
library(HPAanalyze)

node_attributes <- read_excel("~/Desktop/codes_submission/node_attributes.xlsx")
names(node_attributes)[names(node_attributes) == "Gene name"] <- "Gene_name"
names(node_attributes)[names(node_attributes) == "Target class"] <- "Target_class"

network_nodes_VIP<-as.vector(node_attributes %>% 
                  filter(Target_class=="VIP") %>% 
                  dplyr::select(Gene_name))

network_nodes_NIP<-node_attributes %>% 
  filter(Target_class=="NIP") %>% 
  dplyr::select(Gene_name)


downloadedData <- hpaDownload(downloadList='histology', version='latest')

########## VIP list ##########################

geneList_VIP<-network_nodes_VIP$Gene_name
geneList_NIP<-network_nodes_NIP$Gene_name
tissueList <- c('lung','bronchus','nasopharynx')

hpa_VIP_tissue<-hpaVis(data=downloadedData,
                             targetGene=geneList_VIP,
                             targetTissue =tissueList,
                             visType = c("Tissue"))

subset_VIP<-hpaSubset(data=downloadedData,
                      targetGene=geneList_VIP,
                      targetTissue =tissueList)

hpaExport(data=subset_VIP,
          fileName='VIP_tissue_atlas.xlsx',
          fileType='xlsx')

############ For priotrized proteins ##############

hpa_NIP<-hpaVis(data=downloadedData,
                      targetGene=geneList_NIP,
                      targetTissue =tissueList,
                      visType = c("Tissue"))

subset_NIP<-hpaSubset(data=downloadedData,
                      targetGene=geneList_NIP,
                      targetTissue =tissueList)


hpaExport(data=subset_NIP,
          fileName='NIP_tissue_selected.xlsx',
          fileType='xlsx')



##################################### ####################

### Visulization plot

library(readxl)
VIP_tissue_atlas_xlsx <- read_excel("VIP_tissue_atlas.xlsx", 
                                    sheet = "normal_tissue", col_types = c("skip", "text", "text", "text", "text", "skip"))

VIP_set$level<-factor(VIP_set$level,levels = c("High","Medium","Low","Not detected","Not avaiable"))

VIP_set<-VIP_tissue_atlas_xlsx
VIP_set$tissue_cell<-paste(VIP_set$tissue, "/", VIP_set$cell_type)

colours<-c("#810f7c","#e66101","#fdb863","#fff7bc","#f0f0f0")

VIP_plot<-ggplot(VIP_set,aes(gene,tissue_cell,fill=level))+ geom_tile(color = "grey",lwd = 0.3,linetype = 1)+ scale_fill_manual(values = colours,na.value = "white")+theme(axis.text.x = element_text(family="ArialMT",size=12, angle=90,hjust = 0.95,vjust = 0.17),axis.text.y = element_text(family="ArialMT",size = 12),axis.title.x = element_blank(),legend.title = element_blank())+labs(y="Tissue/Cell",family="ArialMT")+coord_equal()

RWR_tissue_selected_xlsx <- read_excel("NIP_tissue_selected.xlsx", 
                                       +     col_types = c("skip", "text", "text", 
                                                           +         "text", "text", "skip"))

RWR_data<-RWR_tissue_selected_xlsx

RWR_data$level<-factor(RWR_data$level,levels = c("High","Medium","Low","Not detected","Not avaiable"))


RWR_data$tissue_cell<-paste(RWR_data$tissue, "/", RWR_data$cell_type)

RWR_plot<-ggplot(RWR_data,aes(gene,tissue_cell,fill=level))+ geom_tile(color = "grey",lwd = 0.3,linetype = 1)+ scale_fill_manual(values = colours,na.value = "white")+theme(axis.text.x = element_text(face="bold",size=10, angle=90,hjust = 0.95,vjust = 0.17),axis.text.y = element_text(face="bold",size = 10),axis.title.x = element_blank(),legend.title = element_blank())+labs(y="Tissue/Cell",face="bold")+coord_equal()

library(ggpubr)

ggarrange(VIP_plot,RWR_plot,labels = c("A", "B"),
          ncol = 2,
          nrow=2, common.legend = TRUE,legend = "bottom")

