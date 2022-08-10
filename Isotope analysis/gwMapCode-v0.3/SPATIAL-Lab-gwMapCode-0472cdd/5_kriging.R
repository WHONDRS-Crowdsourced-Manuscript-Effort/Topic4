library(gstat)
library(raster)
library(doParallel)
library(foreach)

setwd("..")

#raster layers as input
load("aquifer3d.rda")
#load("iso3d.rda")
load("2dVario.rda")
load("wis.rda")

# kriging for d18O -----

#space for output
vs = vms = ps = r.ps = r.vs = cv = list()

registerDoParallel(7)

ko = foreach(i = 1:nlayers(aq), .packages = c("gstat", "raster")) %dopar% {
  #grid for prediction
  grd = aq[[i]]
  grdpts = rasterToPolygons(grd, fun = function(x){x==1})
  
  #prep data
  pts = spTransform(wis[[i]], proj4string(grdpts))

  #predict at points
  p = krige(d18O ~ 1, pts, grdpts, v.2d.mod[[i]], nmin = 3, 
                   nmax = 50, maxdist = 2.5e6)

  #cross validate
  c = krige.cv(d18O~1, pts, model = v.2d.mod[[i]], nfold = nrow(pts))
  
  #points to grid
  r.p = rasterize(p, grd, p$var1.pred)
  r.v = rasterize(p, grd, p$var1.var)
  
  list(c, r.p, r.v)
}
stopImplicitCluster()

load("localSD.rda")
for(i in 1:7){
  cv[[names(aq)[i]]] = ko[[i]][[1]]@data
  r.ps[[names(aq)[i]]] = ko[[i]][[2]]
  r.vs[[names(aq)[i]]] = sqrt(ko[[i]][[3]] + localSD[i]^2)
}

r.ps.b = brick(r.ps)
plot(r.ps.b)

r.vs.b = brick(r.vs)
plot(r.vs.b)

#Error statistics
for(i in 1:7){
  print(shapiro.test(cv[[i]]$residual))
}
for(i in 1:7){
  print(round(c(exp(depths$ud[i]), mean(cv[[i]]$residual), mean(abs(cv[[i]]$residual)),
              mean(cv[[i]]$residual^2), sqrt(mean(cv[[i]]$residual^2))), 2))
}

save(r.ps.b, file = "isoscape.rda")
save(r.vs.b, file = "isovar.rda")
save(cv, file = "krigCV.rda")

writeRaster(r.ps.b, "gwMapCode/out/isoscape_d18O.nc", overwrite = TRUE)
writeRaster(r.vs.b, "gwMapCode/out/isovar_d18O.nc", overwrite = TRUE)

#collapse to 2d stats -----

gw.mean = mean(r.ps.b, na.rm = TRUE)
gw.max = max(r.ps.b, na.rm = TRUE)
gw.min = min(r.ps.b, na.rm = TRUE)
gw.sd = calc(r.ps.b, fun = sd, na.rm = TRUE)
gw = brick(gw.mean, gw.max, gw.min, gw.sd)
names(gw) = c("mean", "max", "min", "sd")
save(gw, file = "gw_iso.rda")

# kriging for d2H -----

load("iso3d_d2H.rda")

#space for output
vs = vms = ps = r.ps = r.vs = cv = list()

registerDoParallel(7)

ko = foreach(i = 1:nlayers(aq), .packages = c("gstat", "raster")) %dopar% {
  #grid for prediction
  grd = aq[[i]]
  grdpts = rasterToPolygons(grd, fun = function(x){x==1})
  
  #prep data
  pts = spTransform(wis[[i]], proj4string(grdpts))
  
  #predict at points
  p = krige(d2H ~ 1, pts, grdpts, v.2d.mod[[i]], nmin = 3, 
            nmax = 50, maxdist = 2.5e6)
  
  #cross validate
  c = krige.cv(d2H~1, pts, model = v.2d.mod[[i]], nfold = nrow(pts))

  #points to grid
  r.p = rasterize(p, grd, p$var1.pred)
  r.v = rasterize(p, grd, p$var1.var)

  list(c, r.p, r.v)
}
stopImplicitCluster()

load("localSD_d2H.rda")
for(i in 1:7){
  cv[[names(aq)[i]]] = ko[[i]][[1]]@data
  r.ps[[names(aq)[i]]] = ko[[i]][[2]]
  r.vs[[names(aq)[i]]] = sqrt(ko[[i]][[3]] + localSD[i]^2)
}

r.ps.b = brick(r.ps)
plot(r.ps.b)

r.vs.b = brick(r.vs)
plot(r.vs.b)

#Error statistics
for(i in 1:7){
  print(shapiro.test(cv[[i]]$residual))
}
for(i in 1:7){
  print(round(c(exp(depths$ud[i]), mean(cv[[i]]$residual), mean(abs(cv[[i]]$residual)),
                mean(cv[[i]]$residual^2), sqrt(mean(cv[[i]]$residual^2))), 2))
}

save(r.ps.b, file = "isoscape_d2H.rda")
save(r.vs.b, file = "isovar_d2H.rda")
save(cv, file = "krigCV_d2H.rda")

writeRaster(r.ps.b, "gwMapCode/out/isoscape_d2H.nc", overwrite = TRUE)
writeRaster(r.vs.b, "gwMapCode/out/isovar_d2H.nc", overwrite = TRUE)
