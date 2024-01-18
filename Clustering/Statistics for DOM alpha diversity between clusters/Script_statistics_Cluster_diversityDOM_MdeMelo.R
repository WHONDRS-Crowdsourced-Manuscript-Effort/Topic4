
#############################################################################################
# WHONDRS - Clustering

#Goal: Run statistics to compare variability in DOm moleuclar composition between clusters
#Author: Michaela de Melo
#Date: 12/27/202
#############################################################################################
#Libraries
library(ggplot2)
library(car)
library(FSA)
library(stringr)
library(ggpubr)



#Folder - Michaela's pc (change yours)
setwd("C:/Users/micha/Google Drive/WHONDRS/Paper-Topic 3-Statistics")

#Import labels
clust_water=read.csv("water_cluster_label.csv", row.names = 1)
clust_sed=read.csv("sediment_cluster_label.csv", row.names = 1)

#Import DOM metadata and richness
rich_water=read.csv("water_merged_metadataindices.csv", row.names = 1)
rich_sed=read.csv("sediment_merged_metadataindices.csv", row.names = 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----> WATER SAMPLES

#Merge both tables
merge_water=merge(clust_water, rich_water, by="row.names")
#rename clusters
merge_water$cluster <-str_replace_all(merge_water$cluster, c("sw_0" = "Cluster 0", "sw_1" = "Cluster 1", "sw_2" = "Cluster 2"))


#Boxplot - Richness

colnames(merge_water[, 4900:4977])
g1a=ggplot()+geom_boxplot(data=merge_water, aes(x=cluster, y=Richness.observed, fill=cluster))+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+ labs(y="Observed Richness")+
  theme(legend.position="none", axis.title.y = element_text(size=12), axis.title.x = element_blank(),  
        axis.text.x = element_text( size=10), axis.text.y = element_text(  size=8))+
  annotate("text", x = 2.8, y = 3800, label = "Kruskal–Wallis test, p-value<0.01", size=3)+
annotate("text", x = 1, y = 3550, label = "a", size=4)+
annotate("text", x = 2, y = 2350, label = "b", size=4)+
annotate("text", x = 3, y = 3250, label = "c", size=4)

ggsave("Boxplot_Richness_Water.png", dpi=200, width = 6, height = 5)


#-->Statistics
#test for homogeneity of variances
leveneTest(Richness.observed ~ as.factor(cluster), data = merge_water)
#test for data normality
shapiro.test(merge_water$Richness.observed)
#Result: non-normal and non-homogeneous (p<0.05)

#-->Non-parametric test : Kruskal-Wallis
#A Kruskal-Wallis test is used to determine whether or not there is a statistically significant difference between the medians of three or more independent groups. 
#It is considered to be the non-parametric equivalent of the One-Way ANOVA.
kruskal.test(Richness.observed ~ as.factor(cluster), data = merge_water) #p-value < 2.2e-16

#If the results of a Kruskal-Wallis test are statistically significant, then it's appropriate to conduct Dunn's Test to determine exactly which groups are different.
dunnTest(Richness.observed ~ as.factor(cluster), data = merge_water,    method="holm")
#There is difference in richness among clusters!



#Boxplot - Simpson

colnames(merge_water[, 4900:4977])
g1b=ggplot()+geom_boxplot(data=merge_water, aes(x=cluster, y=Simpson.s.Diversity.Index, fill=cluster))+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+ labs(y="Simpson's Index")+
  theme(legend.position="none", axis.title.y = element_text(size=12), axis.title.x = element_blank(),  
        axis.text.x = element_text( size=10), axis.text.y = element_text(  size=8))+
 annotate("text", x = 2.8, y = 1, label = "Kruskal–Wallis test, p-value<0.01", size=3)+
  annotate("text", x = 1, y = 0.99975, label = "a", size=4)+
  annotate("text", x = 2, y = 0.99965, label = "b", size=4)+
  annotate("text", x = 3, y = 0.99975, label = "c", size=4)

water_panel=ggarrange(g1a, g1b, nrow = 1, ncol=2)
water_panel

ggsave("Boxplot_Diversity_Water.png", dpi=200, width = 10, height = 4)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----> SEDIMENT SAMPLES

#Merge both tables
merge_sed=merge(clust_sed, rich_sed, by="row.names")
#rename clusters
merge_sed$cluster <-str_replace_all(merge_sed$cluster, c("sed_0" = "Cluster 0", "sed_1" = "Cluster 1", "sed_2" = "Cluster 2"))


#Boxplot - Richness

g2a=ggplot()+geom_boxplot(data=merge_sed, aes(x=cluster, y=Richness.observed, fill=cluster))+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+ labs(y="Observed Richness")+
  theme(legend.position="none", axis.title.y = element_text(size=12), axis.title.x = element_blank(),  
        axis.text.x = element_text( size=10), axis.text.y = element_text(  size=8))+
  annotate("text", x = 1.1, y = 2800, label = "Kruskal–Wallis test, p-value<0.01", size=3)+
  annotate("text", x = 1, y = 2450, label = "a", size=4)+
  annotate("text", x = 2, y = 1550, label = "b", size=4)+
  annotate("text", x = 3, y = 2850, label = "c", size=4)
g2a

#ggsave("Boxplot_Richness_Water.png", dpi=200, width = 6, height = 5)


#-->Statistics
#test for homogeneity of variances
leveneTest(Richness.observed ~ as.factor(cluster), data = merge_sed) #OK
#test for data normality
shapiro.test(merge_sed$Richness.observed)
#Result: non-normal and homogeneous (p<0.05)

#-->Non-parametric test : Kruskal-Wallis
#A Kruskal-Wallis test is used to determine whether or not there is a statistically significant difference between the medians of three or more independent groups. 
#It is considered to be the non-parametric equivalent of the One-Way ANOVA.
kruskal.test(Richness.observed ~ as.factor(cluster), data = merge_sed) #p-value < 2.2e-16

#If the results of a Kruskal-Wallis test are statistically significant, then it's appropriate to conduct Dunn's Test to determine exactly which groups are different.
dunnTest(Richness.observed ~ as.factor(cluster), data = merge_sed,    method="holm")
#There is difference in richness among clusters!



#Boxplot - Simpson

#colnames(merge_water[, 4900:4977])
g2b=ggplot()+geom_boxplot(data=merge_sed, aes(x=cluster, y=Simpson.s.Diversity.Index, fill=cluster))+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+ labs(y="Simpson's Index")+
  theme(legend.position="none", axis.title.y = element_text(size=12), axis.title.x = element_blank(),  
        axis.text.x = element_text( size=10), axis.text.y = element_text(  size=8))+
  annotate("text", x = 2.8, y = 1, label = "Kruskal–Wallis test, p-value<0.01", size=3)+
  annotate("text", x = 1, y = 0.99975, label = "a", size=4)+
  annotate("text", x = 2, y = 0.99965, label = "b", size=4)+
  annotate("text", x = 3, y = 0.99975, label = "c", size=4)
g2b
sed_panel=ggarrange(g2a, g2b, nrow = 1)
sed_panel

ggsave("Boxplot_Diversity_Sediment.png", dpi=200, width = 10, height = 4)



#-->Statistics
#test for homogeneity of variances
leveneTest(Simpson.s.Diversity.Index ~ as.factor(cluster), data = merge_sed)
#test for data normality
shapiro.test(merge_sed$Simpson.s.Diversity.Index)
#Result: non-normal and non-homogeneous (p<0.05)

#-->Non-parametric test : Kruskal-Wallis
kruskal.test(Simpson.s.Diversity.Index~ as.factor(cluster), data = merge_sed) #p-value < 2.2e-16
dunnTest(Simpson.s.Diversity.Index~ as.factor(cluster), data = merge_sed,    method="holm")
#There is difference in richness among clusters!


############################################################################################

