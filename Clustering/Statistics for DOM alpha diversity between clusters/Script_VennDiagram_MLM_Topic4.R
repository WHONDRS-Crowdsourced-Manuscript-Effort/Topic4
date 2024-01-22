############################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Venn Diagrams
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Why? Checking how many molecules are shared among clusters, how many are exclusive
#Create List of molecules of each cluster for surface water (sw) and sediment (sed)


#load libraries
library("ggvenn")
library(ggplot2)

#Folder - Michaela's pc (change yours)
setwd("C:/Users/micha/Google Drive/WHONDRS/Paper-Topic 4-Statistics/Statistics for DOM alpha diversity between clusters")

#Import labels
clust_water=read.csv("water_cluster_label.csv", row.names = 1)
clust_sed=read.csv("sediment_cluster_label.csv", row.names = 1)

#Import DOM metadata and richness
rich_water=read.csv("water_merged_metadataindices.csv", row.names = 1)
rich_sed=read.csv("sediment_merged_metadataindices.csv", row.names = 1)


#----> WATER SAMPLES
#Merge both tables
merge_water=merge(clust_water, rich_water, by="row.names")
#rename clusters
merge_water$cluster <-str_replace_all(merge_water$cluster, c("sw_0" = "Cluster 0", "sw_1" = "Cluster 1", "sw_2" = "Cluster 2"))

#----> SEDIMENT SAMPLES
#Merge both tables
merge_sed=merge(clust_sed, rich_sed, by="row.names")
#rename clusters
merge_sed$cluster <-str_replace_all(merge_sed$cluster, c("sed_0" = "Cluster 0", "sed_1" = "Cluster 1", "sed_2" = "Cluster 2"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Analyzes:

#-->Water

#Cluster 0
c0_sw=filter(merge_water, cluster=="Cluster 0") #147 samples
colnames(c0_sw) #options(max.print=1000000)
row.names(c0_sw)=c0_sw$Row.names
c0_sw.mol =c0_sw[,-c(1:2,4938:4979)] #remove what isnt molecules
#Clean MF with 0
c0_sw.mol2=c0_sw.mol[, colSums (c0_sw.mol) != 0] #4935
#Create a list of molecules
list_c0.sw=list(row.names(t(c0_sw.mol2)))

#Cluster 1
c1_sw=filter(merge_water, cluster=="Cluster 1") #28 samples
c1_sw.mol =c1_sw[,-c(1:2,4938:4979)] 
#Clean MF with 0
c1_sw.mol2=c1_sw.mol[, colSums (c1_sw.mol) != 0] #4043
#Create a list of molecules
list_c1.sw=list(row.names(t(c1_sw.mol2)))

#Cluster 2
c2_sw=filter(merge_water, cluster=="Cluster 2") #90 samples
c2_sw.mol =c2_sw[,-c(1:2,4938:4979)] 
#Clean MF with 0
c2_sw.mol2=c2_sw.mol[, colSums (c2_sw.mol) != 0] #4848
#Create a list of molecules
list_c2.sw=list(row.names(t(c2_sw.mol2)))


#Venn Diagrams
list_venn.sw <- c(list_c0.sw,list_c1.sw, list_c2.sw)
list_venn.sw 
names(list_venn.sw)[1] <- "C0 - sw"
names(list_venn.sw)[2] <- "C1 - sw"
names(list_venn.sw)[3] <- "C2 - sw"


#Plot
ven.sw=ggvenn(list_venn.sw, c("C0 - sw", "C1 - sw", "C2 - sw"), stroke_size = 0.5) +
  ggplot2:: scale_fill_brewer(palette="Dark2")+
  theme( plot.background = element_rect(fill = "white", colour = "white"))

ven.sw
ggsave("Venn_diagram_SW.png", dpi=300, height = 5, width = 5) #Save it!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-->Sediment

#Cluster 0
c0_sed=filter(merge_sed, cluster=="Cluster 0") #76 samples
colnames(c0_sed) #options(max.print=1000000)
colnames(c0_sed)[1:5]
row.names(c0_sed)=c0_sed$Row.names
c0_sed.mol =c0_sed[,-c(1:2,4055:4096)] #select only molecules
#Clean MF with 0
c0_sed.mol2=c0_sed.mol[, colSums (c0_sed.mol) != 0] #4051
#Create a list of molecules
list_c0.sed=list(row.names(t(c0_sed.mol2)))

#Cluster 1
c1_sed=filter(merge_sed, cluster=="Cluster 1") #17 samples
c1_sed.mol =c1_sed[,-c(1:2,4055:4096)] 
#Clean MF with 0
c1_sed.mol2=c1_sed.mol[, colSums (c1_sed.mol) != 0] #2841
#Create a list of molecules
list_c1.sed=list(row.names(t(c1_sed.mol2)))

#Cluster 2
c2_sed=filter(merge_sed, cluster=="Cluster 2") #146 samples
c2_sed.mol =c2_sed[,-c(1:2,4055:4096)] 
#Clean MF with 0
c2_sed.mol2=c2_sed.mol[, colSums (c2_sed.mol) != 0] #4052
#Create a list of molecules
list_c2.sed=list(row.names(t(c2_sed.mol2)))


#Venn Diagrams
list_venn.sed <- c(list_c0.sed,list_c1.sed, list_c2.sed)
list_venn.sed 
names(list_venn.sed)[1] <- "C0 - sed"
names(list_venn.sed)[2] <- "C1 - sed"
names(list_venn.sed)[3] <- "C2 - sed"

#plot
ven.sed=ggvenn(list_venn.sed, c("C0 - sed", "C1 - sed", "C2 - sed"), stroke_size = 0.5) +
  ggplot2:: scale_fill_brewer(palette="Dark2")+
  theme( plot.background = element_rect(fill = "white", colour = "white"))

ven.sed
ggsave("Venn_diagram_SED.png", dpi=300, height = 5, width = 5) #Save it!
#################################################################################
