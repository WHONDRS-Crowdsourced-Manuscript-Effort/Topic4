library(rgdal)
library(sp)
library(raster)
library(rgeos)

setwd("..")

# Well records -----
load("wellDepths.rda")

#dimensions and bounds of area
res = 25000
xmn = bbox(wd)[1,1]
xmx = bbox(wd)[1,2]
ymn = bbox(wd)[2,1]
ymx = bbox(wd)[2,2]
xdim = ceiling((xmx - xmn)/res)
ydim = ceiling((ymx - ymn)/res)

#multiband raster to store subsurface data
gdr = raster(nrows = ydim, ncols = xdim, xmn = xmn, xmx = xmn + (xdim) * res, 
             ymn = ymn, ymx = ymn + (ydim) * res, crs = proj4string(wd), vals=NA)
gdr2 = gdr
for(i in 2:7){
  gdr = brick(c(gdr, gdr2))
}

#corresponding sppoly layer to use in extraction
gd = GridTopology(c(bbox(wd)[1,1]+res/2, bbox(wd)[2,1]+res/2), c(res,res), c(xdim, ydim))
polys = as.SpatialPolygons.GridTopology(gd)
proj4string(polys) = proj4string(wd)

#check extents, should match!
extent(gdr[[1]])
extent(polys)

#depths for raster bands, from 1 m to 2000 m
#ud = seq(0, 7.601, length.out = 17)[1:16]
#ld = seq(0, 7.601, length.out = 17)[2:17]
#round(exp(ud),1); round(exp(ld),1)
ud = c(1, 10, 25, 50, 100, 200, 500)
ld = c(10, 25, 50, 100, 200, 500, 2000)
ud = log(ud)
ld = log(ld)

depths = data.frame(ud, ld)
save(depths, file = "depths.rda")

#find all grid depths w/ wells
for(i in 1:length(polys)){
  w = wd$lnWD[which(!is.na(over(x=wd, y=polys[i])))]
  if(length(w) > 0){
    for(j in 1:length(ud)){
      if(any(w > ud[j] & w < ld[j])){
        gdr[[j]][i] = 1
      } else{
        gdr[[j]][i] = 0
      }
    }
  } else{
    for(j in 1:length(ud)){
      gdr[[j]][i] = NA
    }
  } 
}

names(gdr) = paste0("UD_", round(exp(ud), 1))

save(gdr, file = "gdr.rda")

# Isotope wells -----

#load wiDB wells data
load("wellIsotopes.rda")

#quick plot
plot(wigw.usa)
lines(states)

#multiband raster to store subsurface data
idr = raster(nrows = ydim, ncols = xdim, xmn = xmn, xmx = xmn + (xdim) * res, 
             ymn = ymn, ymx = ymn + (ydim) * res, crs = proj4string(wd), vals=NA)
idr2 = idr
for(i in 2:7){
  idr = brick(c(idr, idr2))
}

#find all grid depths w/ isotope data
for(i in 1:length(polys)){
  w = wigw.usa$lnWD[which(!is.na(over(x=wigw.usa, y=polys[i])))]
  if(length(w) > 0){
    for(j in 1:length(ud)){
      if(any(w > ud[j] & w < ld[j])){
        idr[[j]][i] = 1
      } else {
        idr[[j]][i] = 0
      }
    }
  } else{
    for(j in 1:length(ud)){
      idr[[j]][i] = NA
    }  
  }
}

names(idr) = paste0("UD_", round(exp(ud), 1))

save(idr, file = "idr.rda")

# Merge wells -----
library(doParallel)
load("gdr.rda")
load("idr.rda")

cdr = gdr[[1]]

registerDoParallel(7)

cdrL = foreach(i = 1:7) %dopar% {
  for(j in 1:raster::ncell(cdr)){
    #No well depth data at gridcell
    if(is.na(gdr[[i]][j])){
      # + No isotope data at gridcell = NA
      if(is.na(idr[[i]][j])){
        cdr[j] = NA
        # only isotope data at depth = 2
      } else if(idr[[i]][j] == 1){
        cdr[j] = 2
        # only isotope data in gridcell = 0
      } else{
        cdr[j] = 0
      }
      #No isotope data at gridcell
    } else if(is.na(idr[[i]][j])){
      # only well depth data at depth = 1
      if(gdr[[i]][j] == 1){
        cdr[j] = 1
        # only well depth data in gridcell = 0
      } else{
        cdr[j] = 0
      }
      #Well depth data at depth
    } else if(gdr[[i]][j] == 1){
      # + isotope data at depth = 3
      if(idr[[i]][j] == 1){
        cdr[j] = 3
        # + no isotope data at depth = 1
      } else{
        cdr[j] = 1
      }
      #Isotope data at depth, well depth data at gridcell = 2
    } else if(idr[[i]][j] == 1){
      cdr[j] = 2
      #Neither at depth, one or both at gridcell = 0
    } else{
      cdr[j] = 0
    }
  }
  cdr
}
stopImplicitCluster()

cdr = stack(cdrL)

names(cdr) = paste0("UD_", round(exp(ud), 1))

save(cdr, file = "cdr.rda")

# Subsurface averages for d18O -----

#multiband raster to store subsurface data
isor = raster(nrows = ydim, ncols = xdim, xmn = xmn, xmx = xmn + (xdim) * res, 
              ymn = ymn, ymx = ymn + (ydim) * res, crs = proj4string(wd), vals=NA)
isor2 = isor
for(i in 2:7){
  isor = brick(c(isor, isor2))
}
isosdr = isoctr = isor

#find all grid depths w/ isotope data
for(i in 1:length(polys)){
  w = wigw.usa[which(!is.na(over(x=wigw.usa, y=polys[i]))),]
  if(length(w) > 0){
    for(j in 1:length(ud)){
      if(any(w$lnWD > ud[j] & w$lnWD <= ld[j])){
        isor[[j]][i] = mean(w$d18O[w$lnWD > ud[j] & w$lnWD <= ld[j]], na.rm = TRUE)
        isosdr[[j]][i] = sd(w$d18O[w$lnWD > ud[j] & w$lnWD <= ld[j]], na.rm = TRUE)
        isoctr[[j]][i] = length(w$d18O[w$lnWD > ud[j] & w$lnWD <= ld[j]])
      }
    }
  }
}

names(isor) = names(isosdr) = names(isoctr) = paste0("UD_", round(exp(ud), 1))

save(isor, isosdr, isoctr, file = "iso3d.rda")

#Stats for reporting
localSD = cellStats(isosdr, mean)
mean(localSD)
cellStats(isoctr, max)
save(localSD, file = "localSD.rda")

#Repeat with different resolutions to test number and variance within voxels
res = c(2500, 10000, 25000, 50000, 100000)
xdim = ceiling((xmx - xmn)/res)
ydim = ceiling((ymx - ymn)/res)

sts = data.frame("SD" = numeric(5), "CT" = numeric(5), "NSD" = numeric(5), "N" = numeric(5))

for(k in 1:5) {
  stime = proc.time()
  gd = GridTopology(c(bbox(wd)[1,1]+res[k]/2, bbox(wd)[2,1]+res[k]/2), 
                    c(res[k],res[k]), c(xdim[k], ydim[k]))
  polys = as.SpatialPolygons.GridTopology(gd)
  proj4string(polys) = proj4string(wd)
  
  isosdr = isoctr = isonsd = ison = 0
  
  #find all grid depths w/ isotope data
  for(i in 1:length(polys)){
    w = wigw.usa[which(!is.na(over(x=wigw.usa, y=polys[i]))),]
    if(length(w) > 0){
      for(j in 1:length(ud)){
        wl = w$d18O[w$lnWD > ud[j] & w$lnWD <= ld[j]]
        if(length(wl)  > 0){
          if(length(wl) > 1) {
            isosdr = isosdr + sd(w$d18O[w$lnWD > ud[j] & w$lnWD <= ld[j]], na.rm = TRUE)
            isonsd = isonsd + 1
          }
          isoctr = isoctr + length(w$d18O[w$lnWD > ud[j] & w$lnWD <= ld[j]])
          ison = ison + 1
        }
      }
    }
  }
  sts[k,] = c(isosdr / isonsd, isoctr / ison, isonsd, ison)
  print(paste("Grid", res[k], "done in", proc.time()[3] - stime[3], "seconds"))
}
View(sts)

# Subsurface averages for d2H -----
res = 25000
xdim = ceiling((xmx - xmn)/res)
ydim = ceiling((ymx - ymn)/res)

#multiband raster to store subsurface data
isor = raster(nrows = ydim, ncols = xdim, xmn = xmn, xmx = xmn + (xdim) * res, 
              ymn = ymn, ymx = ymn + (ydim) * res, crs = proj4string(wd), vals=NA)
isor2 = isor
for(i in 2:7){
  isor = brick(c(isor, isor2))
}
isosdr = isoctr = isor

#find all grid depths w/ isotope data
for(i in 1:length(polys)){
  w = wigw.usa[which(!is.na(over(x=wigw.usa, y=polys[i]))),]
  if(length(w) > 0){
    for(j in 1:length(ud)){
      if(any(w$lnWD > ud[j] & w$lnWD <= ld[j])){
        isor[[j]][i] = mean(w$d2H[w$lnWD > ud[j] & w$lnWD <= ld[j]], na.rm = TRUE)
        isosdr[[j]][i] = sd(w$d2H[w$lnWD > ud[j] & w$lnWD <= ld[j]], na.rm = TRUE)
        isoctr[[j]][i] = length(w$d2H[w$lnWD > ud[j] & w$lnWD <= ld[j]])
      }
    }
  }
}

names(isor) = names(isosdr) = names(isoctr) = paste0("UD_", round(exp(ud), 1))

save(isor, isosdr, isoctr, file = "iso3d_d2H.rda")

#Stats for reporting
localSD = cellStats(isosdr, mean)
mean(localSD)
cellStats(isoctr, max)
save(localSD, file = "localSD_d2H.rda")

# 3d aquifer map -----

aq = cdr

for(i in 1:nlayers(aq)){
  for(j in 1:ncell(aq[[1]])){
    if(!is.na(cdr[[i]][j])){
      if(cdr[[i]][j] > 0){
        aq[[i]][j] = 1
      } else{
        aq[[i]][j] = 0
      }
    } else(
      aq[[i]][j] = 0
    )
  }  
}

names(aq) = paste0("UD_", round(exp(ud), 1))

save(aq, file = "aquifer3d.rda")
