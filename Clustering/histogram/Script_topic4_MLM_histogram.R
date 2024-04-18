#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Script for Topic 4: DOM properties in sample clusters

#Author: Michaela de Melo
#Date: July 2022

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)


### MOLECULAR PROPERTIES OF GROUPS HIGHLIGHTED BY THE HEATMAP

######################
## SEDIMENTS SAMPLES
######################

## Set path
setwd("") #change it!
#Upload table with molecular properties
cross_sed=read.csv("Sediment_Prevalence_10_crosstable_2021-09-29.csv", row.names = 1)
cross_wat=read.csv("Water_Prevalence_10_crosstable_2021-09-29.csv", row.names = 1)


#CFs identified in the heatmap (I used a Text extracted tool online to copy the CFs from heatmap image and past to a .csv file) )
setwd("")
heat_sed=read.csv("Heatmap_CF_sediment.csv", row.names = 1)
heat_wat=read.csv("Heatmap_CF_water.csv", row.names = 1)



#Merge both tables
merge_sed=merge(heat_sed, cross_sed, by = "row.names")
merge_water=merge(heat_wat, cross_wat, by = "row.names")
unique(merge_water$Class) #There is Carb in water but not in Sed

#Boxplots to investigate differences in molecular properties between groups of dominant CFs
library(grid)
library(gridExtra)

merge_sed$group=factor(merge_sed$group, levels=c("A", "B", "C"),
                       labels=c("Sed-2", "Sed-1", "Sed-0")) #order data

merge_water$group=factor(merge_water$group, levels=c("D", "F", "E"),
                       labels=c("Wat-2", "Wat-1", "Wat-0")) #order data


## Histogram
merge_sed$Class=factor(merge_sed$Class, levels=c("AminoSugar", "Carb", "ConHC", "Lignin", "Lipid","Protein", "Tannin", "UnsatHC", "Other" ))
gsed=ggplot(merge_sed, aes(x = group, fill = Class)) +theme_classic()+
  geom_bar(position = "fill") +labs( x=" ", y="Relative Contribution", title="Sediment clusters")+
  theme( aspect.ratio=1, strip.text.x = element_text(size = 12), axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.x = element_blank(), axis.text.x = element_text(size = 12))+
  scale_fill_brewer(palette="Set3")
gsed

## Histogram
merge_water$Class=factor(merge_water$Class, levels=c("AminoSugar", "ConHC", "Lignin", "Lipid","Protein", "Tannin", "UnsatHC", "Other", "Carb" ))
gwat=ggplot(merge_water, aes(x = group, fill = Class)) +theme_classic()+
  geom_bar(position = "fill") +labs( x=" ", y="Relative Contribution", title="Water clusters")+
  theme( aspect.ratio=1, strip.text.x = element_text(size = 12), axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.x = element_blank(), axis.text.x = element_text(size = 12))+
  scale_fill_brewer(palette="Set3")
gwat




##Arrange and save it!
g=ggarrange(gsed, gwat, ncol=2, labels = c("A", "B"))
g
ggsave("Histogram_heatmap.png", dpi=300, height = 5, width=12, g)


