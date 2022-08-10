library(sp)
library(rgdal)
library(isoWater)

setwd("..")

# Well records -----

#Read wells
combi = readOGR("us_wells")

#how many wells have no well depth?
noWD1 = sum(is.na(combi$Wll_Dpt))
noWD1

#use water depth if nothing else available
for(i in length(combi)){
  if(is.na(combi$Wll_Dpt[i])){combi$Wll_Dpt[i] = combi$Wtr_Dpt[i]}
}
#for how many did w sub in water depth?
noWD2 = sum(is.na(combi$Wll_Dpt))
noWD1 - noWD2
#answer is none.

#Screen out really shallow and deep ones, missing depths
wd = combi[!is.na(combi$Wll_Dpt),]
wd = wd[wd$Wll_Dpt > 1,]
wd = wd[wd$Wll_Dpt < 2000,]
wd = wd[!is.na(wd$Wll_Dpt),]

#Work with depths in log space
wd$lnWD = log(wd$Wll_Dpt)

#table for reporting
t = table(wd$State)
sort(t)
tdf = data.frame(t)
names(tdf) = c("State", "Count")
s = paste(wd$State, wd$Source, sep = ";")
s = unique(s)
sdf = strsplit(s, ";")
sdf = data.frame(sdf)
sdf = t(sdf)
row.names(sdf) = c()
sdf = data.frame(sdf)
names(sdf) = c("State", "Source")
ts = merge(tdf, sdf)
write.csv(ts, "wellSources.csv", row.names = "FALSE")

save(wd, file = "wellDepths.rda")