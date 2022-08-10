library(gstat)
library(spacetime)
library(rgdal)
library(future)
library(future.apply)
library(raster)

setwd("..")

#Read well data
lf = list.files("wiDB_wells", pattern = "wiDB.shp")
wigw.usa = readOGR(paste0("wiDB_wells/", sort(lf, decreasing = TRUE)[1]))

#Pull out spatial points and check for duplicate coords
wi.sp = SpatialPoints(wigw.usa, proj4string = CRS(proj4string(wigw.usa)))
dup = zerodist(wi.sp)

#11,738 co-located wells; 16,000 unique locations
length(unique(dup[,1]))
nrow(wigw.usa) - length(unique(dup[,1]))

#There are a bunch. Check for equal depths
dep1 = wigw.usa$Dpth_mt[dup[,1]]
dep2 = wigw.usa$Dpth_mt[dup[,2]]
dep = cbind(dep1, dep2)
dep.id = apply(dep, 1, diff)
sum(dep.id == 0)

#9,318 resamples of same well + depth
dup = dup[dep.id == 0,]
length(unique(dup[,1]))

#Now get a list of all samples that appear in the dupes list columns
duplist1 = unique(dup[,1])
duplist2 = unique(dup[,2])

#These should be a unique list of samples w/ a single one in each dup set
keeplist = !(duplist1 %in% duplist2)
keeplist = duplist1[keeplist]

#Calculate averages for each keeper
avesH = numeric(length(keeplist))
avesO = numeric(length(keeplist))
for(i in 1:length(keeplist)){
  matches = dup[dup[,1] %in% keeplist[i],2]
  avesH[i] = mean(c(wigw.usa$d2H[keeplist[i]], wigw.usa$d2H[matches]), na.rm = TRUE)
  avesO[i] = mean(c(wigw.usa$d18O[keeplist[i]], wigw.usa$d18O[matches]), na.rm = TRUE)
}

#Replace the keepers with averages
wigw.usa$d2H[keeplist] = avesH
wigw.usa$d18O[keeplist] = avesO

#Get rid of the dupes
dupall = unique(c(duplist1, duplist2))
droplist = dupall[!(dupall %in% keeplist)]
wigw.usa = wigw.usa[-droplist,]
wi.sp = wi.sp[-droplist,]

#Pull out depth dimension
dpth = log(wigw.usa$Dpth_mt)
dpth.date = as.Date(dpth, origin = "1900-01-01")

#Pull out isotope obs
d2H = wigw.usa$d2H
d18O = wigw.usa$d18O

#Combine
d2H.st = STIDF(wi.sp, time = dpth.date, data = data.frame(d2H))
d18O.st = STIDF(wi.sp, time = dpth.date, data = data.frame(d18O))

#reduce data
indx = seq_along(d18O.st)
indx = cbind(sample(indx, 4000), sample(indx, 4000))
d18O.st.sub.1 = d18O.st[indx[,1]]
d18O.st.sub.2 = d18O.st[indx[,2]]

#st semivariogram for 1st O sub-dataset
t.start = proc.time()
o.var1 = variogramST(d18O~1, data = d18O.st.sub.1, tunit = "days", tlags = seq(0, 6, by = 1),
                    cutoff = 3.5e6, cores = 10)
proc.time()[3] - t.start[3]

#st semivariogram for 2nd O sub-dataset
t.start = proc.time()
o.var2 = variogramST(d18O~1, data = d18O.st.sub.2, tunit = "days", tlags = seq(0, 6, by = 1),
                     cutoff = 3.5e6, cores = 10)
proc.time()[3] - t.start[3]

#convert timelag values back to meters
o.var1$timelag = exp(o.var1$timelag)
o.var2$timelag = exp(o.var2$timelag)

save(o.var1, o.var2, file = "stVario.rda")

#Supplementary Figure 3
png("S3_Fig.png", width = 5.5, height = 3.5, units = "in", res = 600)
plot(o.var1, ylab = "Vertical distance (m)", xlab = "Horizontal distance (m)", 
     log = "y", col = heat.colors(16))
dev.off()

#Now 2d -----

load("depths.rda")

#Read well data
lf = list.files("wiDB_wells", pattern = "wiDB.shp")
wigw.usa = readOGR(paste0("wiDB_wells/", sort(lf, decreasing = TRUE)[1]))

wis = list()

for(i in seq_along(depths$ud)){
  wi = wigw.usa[wigw.usa$Dpth_mt >= exp(depths$ud[i]) & 
                        wigw.usa$Dpth_mt < exp(depths$ld[i]),]
  dup = zerodist(wi)
  
  #Now get a list of all samples that appear in the dupes list columns
  duplist1 = unique(dup[,1])
  duplist2 = unique(dup[,2])
  
  #These should be a unique list of samples w/ a single one in each dup set
  keeplist = !(duplist1 %in% duplist2)
  keeplist = duplist1[keeplist]
  
  #Calculate averages for each keeper
  avesH = numeric(length(keeplist))
  avesO = numeric(length(keeplist))
  for(j in 1:length(keeplist)){
    matches = dup[dup[,1] %in% keeplist[j],2]
    avesH[j] = mean(c(wi$d2H[keeplist[j]], wi$d2H[matches]), na.rm = TRUE)
    avesO[j] = mean(c(wi$d18O[keeplist[j]], wi$d18O[matches]), na.rm = TRUE)
  }
  
  #Replace the keepers with averages
  wi$d2H[keeplist] = avesH
  wi$d18O[keeplist] = avesO
  
  #Get rid of the dupes
  dupall = unique(c(duplist1, duplist2))
  droplist = dupall[!(dupall %in% keeplist)]
  wi = wi[-droplist,]
  
  wis[[i]] = wi
}

save(wis, file = "wis.rda")

#Fit variograms for oxygen ----
v.2d.samp = v.2d.mod = list()

for(i in seq_along(depths$ud)){
  pts = wis[[i]]
  #sample semivariogram
  v.2d.samp[[i]] = variogram(d18O ~ 1, pts, cutoff = 4e6, width = 1e5)
  
  if(i < 7){
    #initial parameter estimates
    v.2d.mod[[i]] = vgm(psill = 35, "Mat", range = 9e5, nugget = 3, kappa = 1, cutoff = 4e6)
    
    #fit the model
    v.2d.mod[[i]] = fit.variogram(v.2d.samp[[i]], v.2d.mod[[i]], fit.kappa = TRUE,
                                  fit.method = 7)
  } else{
    #initial parameter estimates
    v.2d.mod[[i]] = vgm(psill = 80, "Mat", range = 2e6, nugget = 3, kappa = 1, cutoff = 4e6)
    
    #fit the model
    v.2d.mod[[i]] = fit.variogram(v.2d.samp[[i]], v.2d.mod[[i]], fit.ranges = FALSE, 
                                  fit.kappa = TRUE, fit.method = 7)
  }

  print(plot(v.2d.samp[[i]], v.2d.mod[[i]], main = paste(exp(depths$ud[i]), "-", 
                                                         exp(depths$ld[i]), "m,",
                                                         length(pts), "samples")))
}

save(v.2d.samp, v.2d.mod, file = "2dVario.rda")

#Supplementary figure 2
png("S2_Fig.png", width = 5.5, height = 7, units = "in", res = 600)
layout(matrix(c(1, 2), nrow = 2))
par(mar = c(5, 5, 1, 1))

i = 2
vs = v.2d.samp[[i]]
vm = v.2d.mod[[i]]
plot(vs$dist, vs$gamma, xlab = "Distance (m)", ylab = "Semivariance",
     xlim = c(0, 3.5e6))
lines(variogramLine(vm, 4e6))
text(200, 40, "10 - 25 m", pos = 4)

i = 5
vs = v.2d.samp[[i]]
vm = v.2d.mod[[i]]
plot(vs$dist, vs$gamma, xlab = "Distance (m)", ylab = "Semivariance",
     xlim = c(0, 3.5e6))
lines(variogramLine(vm, 4e6))
text(200, 26.7, "100 - 200 m", pos = 4)

dev.off()

#save parameters
vparms = data.frame(numeric(0), numeric(0), numeric(0), numeric(0))
for(i in seq_along(depths$ud)){
  vm = v.2d.mod[[i]]
  vparms = rbind(vparms, data.frame(vm$psill[1], vm$psill[2], vm$range[2] / 1e3,
                                    vm$kappa[2]))
}
np = unlist(lapply(wis, length))
vparms = cbind(vparms, np)
names(vparms) = c("Nugget", "Partial Sill", "Range (km)", "Kappa", "Observations")
write.csv(vparms, "2dVarioParms.csv", row.names = FALSE)

#Exploratory plot
ranges = kappa = double(length(v.2d.mod))
np = unlist(lapply(wis, length))
for(i in 1:length(v.2d.mod)){
  ranges[i] = v.2d.mod[[i]]$range[2]
  kappa[i] = v.2d.mod[[i]]$kappa[2]
}
plot(depths$ud, ranges, ylim = c(0.9e6, 7e6))
text(depths$ud, ranges, np, pos = 3)
nr = lm(ranges~np)
plot(exp(depths$ud), nr$residuals)

#Fit variograms for hydrogen ----
v.2d.samp = v.2d.mod = list()

for(i in seq_along(depths$ud)){
  pts = wis[[i]]
  #sample semivariogram
  v.2d.samp[[i]] = variogram(d2H ~ 1, pts, cutoff = 4e6, width = 1e5)
  
  if(i < 7){
    #initial parameter estimates
    v.2d.mod[[i]] = vgm(psill = 2000, "Mat", range = 9e5, nugget = 50, kappa = 1, cutoff = 4e6)
    
    #fit the model
    v.2d.mod[[i]] = fit.variogram(v.2d.samp[[i]], v.2d.mod[[i]], fit.kappa = TRUE,
                                  fit.method = 7)
  } else{
    #initial parameter estimates
    v.2d.mod[[i]] = vgm(psill = 2000, "Mat", range = 2e6, nugget = 50, kappa = 1, cutoff = 4e6)
    
    #fit the model
    v.2d.mod[[i]] = fit.variogram(v.2d.samp[[i]], v.2d.mod[[i]], fit.ranges = FALSE, 
                                  fit.kappa = TRUE, fit.method = 7)
  }
  
  print(plot(v.2d.samp[[i]], v.2d.mod[[i]], main = paste(exp(depths$ud[i]), "-", 
                                                         exp(depths$ld[i]), "m,",
                                                         length(pts), "samples")))
}

save(v.2d.samp, v.2d.mod, file = "2dVario_d2H.rda")

#save parameters
vparms = data.frame(numeric(0), numeric(0), numeric(0), numeric(0))
for(i in seq_along(depths$ud)){
  vm = v.2d.mod[[i]]
  vparms = rbind(vparms, data.frame(vm$psill[1], vm$psill[2], vm$range[2] / 1e3,
                                    vm$kappa[2]))
}
np = unlist(lapply(wis, length))
vparms = cbind(vparms, np)
names(vparms) = c("Nugget", "Partial Sill", "Range (km)", "Kappa", "Observations")
write.csv(vparms, "2dVarioParms_d2H.csv", row.names = FALSE)
