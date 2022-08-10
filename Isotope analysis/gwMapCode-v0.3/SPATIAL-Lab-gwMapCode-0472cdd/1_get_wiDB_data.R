#####
#Get GW data from wiDB and save for use in analysis

library(sp)
library(rgdal)
setwd("..")

#load wiDB wells data
library(RODBC)
ch = odbcConnect("WIDB") ##Only works for mapped data source, i.e. on Gabe's machines
wigw = sqlQuery(ch, "SELECT Sites.Site_ID, Sites.Longitude, Sites.Latitude, 
                Sites.Elevation_mabsl, Sites.State_or_Province, 
                Samples.Sample_ID, Samples.Depth_meters, 
                Samples.Collection_Date, AVG(Water_Isotope_Data.d2H), 
                AVG(Water_Isotope_Data.d18O), Samples.Project_ID FROM 
                Water_Isotope_Data INNER JOIN Samples ON 
                Water_Isotope_Data.Sample_ID = Samples.Sample_ID INNER JOIN 
                Sites ON Samples.Site_ID = Sites.Site_ID WHERE Samples.Type = 
                'Ground' AND Water_Isotope_Data.WI_Analysis_Ignore = 0 
                GROUP BY Samples.Sample_ID")

close(ch)

#remove some NA coords or depths plus really shallow wells
wigw = wigw[!is.na(wigw$Depth_meters),]
max(wigw$Depth_meters)
wigw = wigw[wigw$Depth_meters < 2000,]
min(wigw$Depth_meters)
wigw = wigw[wigw$Depth_meters > 1,]
dex = wigw$`AVG(Water_Isotope_Data.d2H)` - 8 * wigw$`AVG(Water_Isotope_Data.d18O)`
wigw = wigw[dex > -10,]
dex = wigw$`AVG(Water_Isotope_Data.d2H)` - 8 * wigw$`AVG(Water_Isotope_Data.d18O)`
wigw = wigw[dex < 25,]
wigw = wigw[!is.na(wigw$Latitude),]
wigw = wigw[!is.na(wigw$Longitude),]

#make geographic
wigw.spdf = SpatialPointsDataFrame(data.frame(wigw$Longitude, wigw$Latitude), data = wigw, 
                                   proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

#project and extract only USA
states = readOGR("states_shapefile/states.shp")
states = states[2:50,] #CONUS only
states = spTransform(states, CRS(aea))
wigw.spdf = spTransform(wigw.spdf, CRS(aea))
wigw.usa = wigw.spdf[states,]
names(wigw.usa)[9:10] = c("d2H", "d18O")

plot(states)
points(wigw.usa)
wigwStates = over(wigw.usa, states)
table(wigwStates$STATE_ABBR)

#save it
writeOGR(wigw.usa, paste0("wiDB_wells/", Sys.Date(), "_wiDB.shp"), paste0("wiDB_wells/", Sys.Date(), "_wiDB.shp"), driver = "ESRI Shapefile")
write.csv(wigw.usa@data, paste0("wiDB_wells/", Sys.Date(), "_wiDB.csv"), row.names = FALSE)
lf = list.files("wiDB_wells", pattern = "wiDB.shp")
wigw.usa = readOGR(paste0("wiDB_wells/", sort(lf, decreasing = TRUE)[1]))
wigw.usa$lnWD = log(wigw.usa$Dpth_mt)
save(wigw.usa, file = "wellIsotopes.rda")

#get and write out sources for reporting
ch = odbcConnect("WIDB")

sources = unique(wigw.usa$Project_ID)
s = as.character(sources[1])
for(i in sources[-1]){
  s = paste(s, i, sep = ", ")
}
qtext = paste0("SELECT * FROM Projects WHERE Project_ID IN (", s, ")")
wis = sqlQuery(ch, qtext)
write.csv(wis, "wiSources.csv", row.names = FALSE)

close(ch)

#Get tap water data -----

library(isoWater)
library(sp)
library(rgdal)

#get 2019 susa transect
susa = wiDB_data(minDate = "2019-03-01", maxDate = "2019-04-15", projects = "00225")
susa = susa[[1]]

#remove a couple sites with surface water source
susa = susa[susa$Site_Name != "GWC airport",]
susa = susa[susa$Site_Name != "DNL airport",]

#get 2007 water data
ousa = wiDB_data(projects = "00077")
ousa = ousa[[1]]

#get only the GW towns
gwt = read.csv("gw_tap_sites.csv")
ousa = ousa[ousa$Site_ID %in% gwt$Site_ID,]

tw = rbind(susa, ousa)

tw = SpatialPointsDataFrame(data.frame(tw$Longitude, tw$Latitude), 
                              data = tw, 
                              proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
tw = spTransform(tw, CRS(aea))

#Remove sites not on the GW layer
load("gw_iso.rda")
gw = extract(gw_iso[[1]], tw)
tw = tw[!is.na(gw),]

sum(tw$Project_ID == 77)
sum(tw$Project_ID == 225)

save(tw, file = "tapData.rda")
