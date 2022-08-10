library(RColorBrewer)
library(cubeview)
library(htmlwidgets)
library(raster)
library(rgdal)

setwd("..")

#####
#Load data
#####

load("gdr.rda")
load("idr.rda")
load("cdr.rda")
load("wellDepths.rda")
load("wellIsotopes.rda")
load("depths.rda")

#Facet plots

plot(gdr, col = c("dark grey", "blue"))
plot(idr, col = c("dark grey", "red"))
plot(cdr, col = c("grey", "blue", "red", "purple"))

#3D models

saveWidget(cubeview(gdr[[7:1]], c(-0.5, 0.5, 1.5), col.regions = c("grey", "blue"), na.color = "white"), "gdr.html")
saveWidget(cubeview(idr[[7:1]], c(-0.5, 0.5, 1.5), col.regions = c("grey", "red"), na.color = "white"), "idr.html")
saveWidget(cubeview(cdr[[7:1]], c(-0.5, 0.5, 1.5, 2.5, 3.5), col.regions = c("grey", "blue", "red", "purple"), 
                    na.color = "white"), "S1_Fig.html")

#Plots...
states = readOGR("states_shapefile/states.shp")
states = spTransform(states, proj4string(cdr))
states = states[2:50,]

#FIGURE1 wd database and wigw database depth distributions
png("Fig1.png", width = 7, height = 5, units = "in", res = 600)
layout(matrix(c(1,3,2,3), nrow = 2, byrow = TRUE), 
       widths = c(lcm(4 * 2.54), lcm(3 * 2.54)),
       heights = rep(lcm(2.5 * 2.54), 2))

par(mar = rep(0, 4))
plot(states, col = "grey90", border = "grey50")
points(wd, pch = ".", col = "steelblue2")
text(-2.4e6, 3e6, "A", cex = 1.5)

plot(states, col = "grey90", border = "grey50")
points(wigw.usa, pch = ".", col = "tomato1")
points(wigw.usa[wigw.usa$Prjc_ID == 225,], pch = ".")
text(-2.4e6, 3e6, "B", cex = 1.5)

par(mar = c(7,5,5,1))
plot(density(wd$lnWD, adjust = 4), main = "", 
     xlab = "Well depth (m)", col = "steelblue2", axes = FALSE)
axis(2)
axis(1, log(c(1, 10, 25, 100, 500, 2000)), 
     c(1, 10, 25, 100, 500, 2000))
box()
lines(density(wigw.usa$lnWD), col = "tomato1")
text(0, 0.43, "C", cex = 1.5)
dev.off()

#Explore comparisons between USGS prinicipal aquifers and grid 
aqs = readOGR("USGS_layers/aquifrp025.shp")
proj4string(aqs) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
aqs = spTransform(aqs, proj4string(cdr))
aqs = aqs[states,]

aqz = unique(aqs$AQ_NAME)
i=3
for(aq in aqz){
  plot(cdr[[i]], col = c("grey75", "steelblue2", "tomato1", "mediumpurple2"), legend = FALSE, 
       axes = FALSE, box = FALSE)
  lines(aqs[aqs$AQ_NAME == aq,])
  text(1.2e6, 2.8e6, aq, adj = 0.5)
}

#Function for figure
shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
          labels, col=bg, ... )
  }
  text(xy$x, xy$y, labels, col=col, ... )
}

#FIGURE2 showing example layers
png("Fig2.png", width = 8, height = 12, units = "in", res = 600)
layout(matrix(c(1,2,3), nrow = 3))

par(mai = c(0.1, 0.1, 0.1, 0.1))

i = 1
plot(cdr[[i]], col = c("grey75", "steelblue2", "tomato1", "mediumpurple2"), 
     legend = FALSE, axes = FALSE, box = FALSE)
lines(states, col = "grey25")
text(-2.4e6, 3e6, "A", cex = 3)
lt = paste(round(exp(depths$ud[i]), 0), "-", round(exp(depths$ld[i]), 0), "m")
text(1.2e6, 2.8e6, lt, adj = 0.5, cex = 2)
lines(aqs[aqs$AQ_CODE %in% c(107),], lw = 2)
shadowtext(-0.3e6, 2.1e6, "HP", cex = 2, col = "black", bg = "lightyellow")
lines(c(-2e6, -2e6, -1e6, -1e6), c(5.5e5, 5e5, 5e5, 5.5e5))
text(-1.5e6, 6e5, "1000 km", cex = 1.5)

i = 4
plot(cdr[[i]], col = c("grey75", "steelblue2", "tomato1", "mediumpurple2"), 
     legend = FALSE, axes = FALSE, box = FALSE)
lines(states, col = "grey25")
text(-2.4e6, 3e6, "B", cex = 3)
lt = paste(round(exp(depths$ud[i]), 0), "-", round(exp(depths$ld[i]), 0), "m")
text(1.2e6, 2.8e6, lt, adj = 0.5, cex = 2)
lines(aqs[aqs$AQ_CODE %in% c(107, 315, 115),], lw = 2)
shadowtext(-1.8e6, 2.8e6, "WV", cex = 2, col = "black", bg = "lightyellow")
shadowtext(-0.3e6, 2.1e6, "HP", cex = 2, col = "black", bg = "lightyellow")
shadowtext(-0.1e6, 2.6e6, "UC", cex = 2, col = "black", bg = "lightyellow")

i = 6
plot(cdr[[i]], col = c("grey75", "steelblue2", "tomato1", "mediumpurple2"), 
     legend = FALSE, axes = FALSE, box = FALSE)
lines(states, col = "grey25")
text(-2.4e6, 3e6, "C", cex = 3)
lt = paste(round(exp(depths$ud[i]), 0), "-", round(exp(depths$ld[i]), 0), "m")
text(1.2e6, 2.8e6, lt, adj = 0.5, cex = 2)
lines(aqs[aqs$AQ_CODE %in% c(405, 106, 302, 202, 115, 112, 606),], lw = 2)

shadowtext(-1.8e6, 2.8e6, "WPS", cex = 2, col = "black", bg = "lightyellow")
shadowtext(-2.3e6, 1.7e6, "CV", cex = 2, col = "black", bg = "lightyellow")
shadowtext(-1.45e6, 2.55e6, "SRP", cex = 2, col = "black", bg = "lightyellow")
shadowtext(-0.5e6, 1.9e6, "DB", cex = 2, col = "black", bg = "lightyellow")
shadowtext(0.1e6, 1.85e6, "OP", cex = 2, col = "black", bg = "lightyellow")
shadowtext(0.2e6, 0.8e6, "TCU", cex = 2, col = "black", bg = "lightyellow")

dev.off()

#####
#summarization
#####
wel = cdr == 1
wel.s = sum(raster::cellStats(wel, "sum"))
iso = cdr == 2
iso.s = sum(raster::cellStats(iso, "sum"))
bot = cdr == 3
bot.s = sum(raster::cellStats(bot, "sum"))
anyWell = cdr > 0
anyWell.s = sum(raster::cellStats(anyWell, "sum"))
n = !is.na(cdr)
n.s = sum(raster::cellStats(n, "sum"))

#fraction of all cells with aquifer
anyWell.s / n.s

#now do that by depth
for(i in 1:7){
  print(paste(exp(depths$ud[i]), ":", exp(depths$ld[i])))
  print(sum(getValues(cdr[[i]]) > 0, na.rm = TRUE) / sum(!is.na(getValues(cdr[[i]]))))
}

#fraction of wd cells with isotope data, ~11.5%
bot.s / wel.s

#fraction of isotope cells w/o wd, ~6.5%
iso.s / (iso.s + bot.s)


library(gstat)
load("isoscape.rda")
load("isovar.rda")
load("variograms.rda")
load("variomodels.rda")
saveWidget(cubeview(r.ps.b[[7:1]], seq(-20, 4, by=2)), "S5_Fig.html")
saveWidget(cubeview(r.vs.b[[7:1]], c(0.4, 0.6, 0.9, 1.3, 1.8, 2.4, 3.1, 4)), 
           "S6_Fig.html")

#FIGURE 3 showing example layers from isoscape
png("Fig3.png", res = 600, units = "in", width = 8, height = 12)

layout(matrix(c(1,2,3), nrow = 3))
par(mai = c(0.1, 0.1, 0.1, 0.1))

i = 1
plot(states, col = "grey")
plot(r.ps.b[[i]], col = rev(heat.colors(15)), add = TRUE, legend = FALSE)
lines(states, col = "grey25")
text(-2.4e6, 3e6, "A", cex = 3)
lt = paste(round(exp(depths$ud[i]), 0), "-", round(exp(depths$ld[i]), 0), "m")
text(1.2e6, 2.8e6, lt, adj = 0.5, cex = 2)
lines(c(-2e6, -2e6, -1e6, -1e6), c(5.5e5, 5e5, 5e5, 5.5e5))
text(-1.5e6, 6e5, "1000 km", cex = 1.5)

i = 3
plot(states, col = "grey")
plot(r.ps.b[[i]], col = rev(heat.colors(15)), add = TRUE, legend = FALSE)
lines(states, col = "grey25")
plot(r.ps.b[[i]], col = rev(heat.colors(15)), legend.only = TRUE,
     smallplot = c(0.9, 0.92, 0.2, 0.8),
     axis.args = list(cex.axis = 1.5),
     legend.args = list(text = expression("Modeled "*delta^{18}*"O"), 
                        side = 2, line = 0.5, cex = 1.2))
text(-2.4e6, 3e6, "B", cex = 3)
lt = paste(round(exp(depths$ud[i]), 0), "-", round(exp(depths$ld[i]), 0), "m")
text(1.2e6, 2.8e6, lt, adj = 0.5, cex = 2)

i = 5
plot(states, col = "grey")
plot(r.ps.b[[i]], col = rev(heat.colors(15)), add = TRUE, legend = FALSE)
lines(states, col = "grey25")
text(-2.4e6, 3e6, "C", cex = 3)
lt = paste(round(exp(depths$ud[i]), 0), "-", round(exp(depths$ld[i]), 0), "m")
text(1.2e6, 2.8e6, lt, adj = 0.5, cex = 2)

dev.off()

#####
#compare with precip isoscape
#####

pcp = raster("Oma.asc")
pcp = projectRaster(pcp, r.ps.b)

#first compare our map w/ the most recent Terzer product, downloaded from 
#https://isotopehydrologynetwork.iaea.org/
pcp.t1 = raster("O18_W120_N0_13.tif")
pcp.t2 = raster("O18_W120_N30_13.tif")
pcp.t3 = raster("O18_W90_N0_13.tif")
pcp.t4 = raster("O18_W90_N30_13.tif")
pcp.t5 = raster("O18_W150_N30_13.tif")

pcp.t = mosaic(pcp.t1, pcp.t2, pcp.t3, pcp.t4, pcp.t5, fun = "mean")
pcp.t = projectRaster(pcp.t, pcp)
plot(pcp.t - pcp)
plot(r.ps.b[[c(1, 3, 5)]] - pcp.t)

#Figure 4 showing GWI data v Precip isoscape values
pwg = extract(pcp, wigw.usa)
pal = brewer.pal(8, "YlGnBu")
pal = pal[2:8]
coli = findInterval(wigw.usa$lnWD, depths$ud)

png("Fig4.png", units = "in", res = 600, width = 6, height = 4)
layout(matrix(c(1,2), nrow = 1), widths = c(5,1))
par(mar = c(5, 5, 2, 1))
plot(pwg, wigw.usa$d18O, col = pal[coli],
     xlab = expression("Precipitation "*delta^{18}*"O (\u2030)"),
     ylab = expression("Groundwater "*delta^{18}*"O (\u2030)"))
abline(0, 1)
par(mar = rep(0, 4))
plot(0, type = "n", axes = FALSE, xlab = "")
legend("left", legend = exp(c(depths$ud)), pch = 1, col = pal, 
       bty = "n", title = "Depth\nInterval\nTop (m)")
dev.off()

#Some statistics for the layers
for(i in 1:7){
  cat(paste(exp(depths$ud[i])), mean(wigw.usa$d18O[coli == i] - 
                                       pwg[coli == i], na.rm = TRUE),
      "\n")
}
mean(wigw.usa$d18O[coli %in% c(1,2,3)] - pwg[coli %in% c(1,2,3)], 
     na.rm = TRUE)
mean(wigw.usa$d18O[coli %in% c(6,7)] - pwg[coli %in% c(6,7)], 
     na.rm = TRUE)

#now GW-PCP isoscape differences
gwdif = r.ps.b - pcp

save(gwdif, file = "gwdiff.rda")

load("gwdiff.rda")
breaks = c(-13, -4, -2, -1, 0, 1, 2, 4, 10)

cols = brewer.pal(8, "RdYlBu")
saveWidget(cubeview(gwdif[[7:1]], breaks, rev(cols)), "S7_Fig.html")

#Figure 5 showing precip-gw layers
png("Fig5.png", res = 600, units = "in", width = 8.2, height = 12)

layout(matrix(c(1,2,3), nrow = 3))
par(mai = c(0.1, 0.1, 0.1, 0.1))

i = 1
plot(states, col = "grey")
plot(gwdif[[i]], breaks = breaks, col = rev(cols), 
     add = TRUE, legend = FALSE)
lines(states, col = "grey25")
text(-2.4e6, 3e6, "A", cex = 3)
lt = paste(round(exp(depths$ud[i]), 0), "-", round(exp(depths$ld[i]), 0), "m")
text(1.2e6, 2.8e6, lt, adj = 0.5, cex = 2)
lines(c(-2e6, -2e6, -1e6, -1e6), c(5.5e5, 5e5, 5e5, 5.5e5))
text(-1.5e6, 6e5, "1000 km", cex = 1.5)

i = 3
plot(states, col = "grey")
plot(gwdif[[i]], breaks = breaks, col = rev(cols), 
     add = TRUE, legend = FALSE)
lines(states, col = "grey25")
plot(gwdif[[i]], breaks = breaks, col = rev(cols), 
     legend.only = TRUE, add = TRUE,
     smallplot = c(0.9, 0.92, 0.2, 0.8),
     axis.args = list(cex.axis = 1.5),
     legend.args = list(text = expression("Groundwater - precipitation "*delta^{18}*"O"), 
                        side = 2, line = 0.5, cex = 1.2))
text(-2.4e6, 3e6, "B", cex = 3)
lt = paste(round(exp(depths$ud[i]), 0), "-", round(exp(depths$ld[i]), 0), "m")
text(1.2e6, 2.8e6, lt, adj = 0.5, cex = 2)

i = 5
plot(states, col = "grey")
plot(gwdif[[i]], breaks = breaks, col = rev(cols), 
     add = TRUE, legend = FALSE)
lines(states, col = "grey25")
text(-2.4e6, 3e6, "C", cex = 3)
lt = paste(round(exp(depths$ud[i]), 0), "-", round(exp(depths$ld[i]), 0), "m")
text(1.2e6, 2.8e6, lt, adj = 0.5, cex = 2)

dev.off()

#Plot residuals -----
load("krigCV.rda")
load("depths.rda")
png("FigS4.png", width = 5, height = 4, units = "in", res = 600)
par(mar = c(5, 5, 1, 1))
plot(density(cv[[1]]$residual), xlab = "Cross-validation error", main = "", 
     ylim = c(0, 1))
for(i in 2:7){
  lines(density(cv[[i]]$residual), col = i)
}
legend(-7.7, 0.8, paste(exp(depths[,1]), "-", exp(depths[,2]), "m"), 
       lty = 1, col = c(1:7), bty = "n")
dev.off()
  
#Figure 6: 2d summary stats and compare with tap -----
load("tapData.rda")
load("gw_iso.rda")

png("Fig6.png", width = 8, height = 8, units = "in", res = 600)
layout(matrix(c(1, 2), nrow = 2))
par(mar = c(0, 0, 0, 5))

plot(states, col = "grey")
vrange = summary(gw$mean)[c(1, 5)]
breaks = seq(-20, 2, by = 2)
cols = rev(heat.colors(length(breaks) - 1))
plot(gw$mean, breaks = breaks, col = cols, 
     axes = FALSE, box = FALSE, add = TRUE,
     legend.args = 
       list(text = expression("Ground or tap water "*delta^{18}*"O"), 
            side = 2, line = 0.5))
lines(states, col = "grey25")
pcol = ceiling((tw$d18O - min(breaks)) / 2)
points(tw@coords, pch = 21, bg = cols[pcol])
text(-2.4e6, 3e6, "A", cex = 2.5)
lines(c(-2e6, -2e6, -1e6, -1e6), c(5.5e5, 5e5, 5e5, 5.5e5))
text(-1.5e6, 6e5, "1000 km")

plot(states, col = "grey")
plot(gw$max - gw$min, col = rev(heat.colors(15)), axes = FALSE, box = FALSE, 
     add = TRUE,
     legend.args = 
       list(text = expression("Groundwater "*delta^{18}*"O range"), 
            side = 2, line = 0.5))
lines(states, col = "grey25")
text(-2.4e6, 3e6, "B", cex = 2.5)
dev.off()

#Now xy comparisons
gwe = extract(gw, tw)

#get precip & uncertainty, from https://wateriso.utah.edu/waterisotopes/pages/data_access/ArcGrids.html
pcp = raster("Oma.asc")
pcp.ci = raster("OmaCI.asc")
proj4string(pcp.ci) = proj4string(pcp)
pcp = brick(pcp, pcp.ci / 1.96)
names(pcp) = c("pcp", "pcp_sd")
pcp = projectRaster(pcp, gw)
pwe = extract(pcp, tw)

#Stahl isoscape
gw.stahl = raster("Stahl/RF_model_GW_isotopes_O.gri")
proj4string(gw.stahl) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
gw.stahl = projectRaster(gw.stahl, gw)

gwe.stahl = extract(gw.stahl, tw)
plot(gwe.stahl, tw$d18O)

##Figure 7: validation plots
png("Fig7.png", res = 600, units = "in", width = 4, height = 12)
layout(matrix(c(1, 2, 3), nrow = 3))
par("mar" = c(5.1, 5.1, 0.5, 0.5), cex.axis = 1.5, cex.lab = 1.5)
plot(gwe[,1], tw$d18O, pch = 20, xlim = c(-20, 1), ylim = c(-20, 1), 
     xlab = expression("Groundwater "*delta^{18}*"O (\u2030)"),
     ylab = expression("Tap water "*delta^{18}*"O (\u2030)"), cex = 1.5)
arrows(gwe[,1] - gwe[,4], tw$d18O, gwe[,1] + gwe[,4], tw$d18O, 
       length = 0.05, angle = 90, code = 3, col = "grey")
points(gwe[,1], tw$d18O, pch = 20, col = "blue")
gwlm = lm(tw$d18O ~ gwe[,1])
abline(0, 1, lty = 3)
abline(gwlm, col = "blue")
text(-3, -16, paste("R^2 =", round(summary(gwlm)$adj.r.squared, 2)), cex = 1.5)
text(-3, -18, paste("MAE =", round(mean(tw$d18O - gwe[,1], na.rm = TRUE), 2)), 
     cex = 1.5)
text(-19.5, 0.6, "A", cex = 1.5)

plot(gwe.stahl, tw$d18O, pch = 20, xlim = c(-20, 1), ylim = c(-20, 1), 
     xlab = expression("Groundwater "*delta^{18}*"O (\u2030)"),
     ylab = expression("Tap water "*delta^{18}*"O (\u2030)"), 
     col = "dark green", cex = 1.5)
gwlm.stahl = lm(tw$d18O ~ gwe.stahl)
abline(0, 1, lty = 3)
abline(gwlm, col = "dark green")
text(-3, -16, paste("R^2 =", round(summary(gwlm.stahl)$adj.r.squared, 2)), cex = 1.5)
text(-3, -18, paste("MAE =", round(mean(tw$d18O - gwe.stahl, na.rm = TRUE), 2)), 
     cex = 1.5)
text(-19.5, 0.6, "B", cex = 1.5)

plot(pwe[,1], tw$d18O, pch = 20, xlim = c(-20, 1), ylim = c(-20, 1), 
     xlab = expression("Precipitation "*delta^{18}*"O (\u2030)"),
     ylab = expression("Tap water "*delta^{18}*"O (\u2030)"), cex = 1.5)
arrows(pwe[,1] - pwe[,2], tw$d18O, pwe[,1] + pwe[,2], tw$d18O, 
       length = 0.05, angle = 90, col = "grey", code = 3)
points(pwe[,1], tw$d18O, pch = 20, col = "red")
abline(0, 1, lty = 3)
plm = lm(tw$d18O ~ pwe[,1])
abline(plm, col = "red")
text(-3, -16, paste("R^2 =", round(summary(plm)$adj.r.squared, 2)), cex = 1.5)
text(-3, -18, paste("MAE =", round(mean(tw$d18O - pwe[,1], na.rm = TRUE), 2)), 
     cex = 1.5)
text(-19.5, 0.6, "C", cex = 1.5)
dev.off()

#stats for prediction
load("isoscape.rda")
load("isovar.rda")
gwe3 = extract(r.ps.b, tw)
gwe3sd = extract(r.vs.b, tw)

gws95 = gws1s = gwe3
for(i in 1:7){
  gws1s[,i] = tw$d18O > (gwe3[,i] - gwe3sd[,i]) &
    tw$d18O < (gwe3[,i] + gwe3sd[,i]) 
  gws95[,i] = tw$d18O > (gwe3[,i] - gwe3sd[,i] * 1.96) & 
    tw$d18O < (gwe3[,i] + gwe3sd[,i] * 1.96)
}

mch1s = mch95 = rep(FALSE, (nrow(gws1s)))
for(i in 1:nrow(gws1s)){
  if(any(gws1s[i, !is.na(gws1s[i,])] == 1)) mch1s[i] = TRUE
  if(any(gws95[i, !is.na(gws95[i,])] == 1)) mch95[i] = TRUE
}
sum(mch1s) / 273
sum(mch95) / 273

#Not clear what error values to use here, their reported prediction variance
#is very low and the spatial layer for this not available; Value below
#is RMSE for their validation dataset
sum(tw$d18O > gwe.stahl - 1.17 & tw$d18O < gwe.stahl + 1.17, na.rm = TRUE) / 
  sum(!is.na(gwe.stahl))
sum(tw$d18O > gwe.stahl - 1.17 * 1.96 & tw$d18O < gwe.stahl + 1.17 * 1.96, na.rm = TRUE) / 
  sum(!is.na(gwe.stahl))

#this is using their average reported prediction error
sum(tw$d18O > gwe.stahl - 0.2 & tw$d18O < gwe.stahl + 0.2, na.rm = TRUE) / 
  sum(!is.na(gwe.stahl))
sum(tw$d18O > gwe.stahl - 0.2 * 1.96 & tw$d18O < gwe.stahl + 0.2 * 1.96, na.rm = TRUE) / 
  sum(!is.na(gwe.stahl))

#for precipitation
sum(tw$d18O > pwe[,1] - pwe[,2] & tw$d18O < pwe[,1] + pwe[,2], na.rm = TRUE) / 
  sum(!is.na(pwe[,1]))
sum(tw$d18O > pwe[,1] - 1.96 * pwe[,2] & tw$d18O < pwe[,1] + 1.96 * pwe[,2], na.rm = TRUE) / 
  sum(!is.na(pwe[,1]))

#exploratory, where is our model closer than the ML SW one?
twc95.stahl = tw@data[!(tw$d18O > gwe.stahl - 2 * 1.17 & 
                          tw$d18O < gwe.stahl + 2 * 1.17),]
twc95.stahl = tw[tw$Site_ID %in% twc95.stahl$Site_ID,]
twc68.stahl = tw@data[!(tw$d18O > gwe.stahl - 1.17 & 
                          tw$d18O < gwe.stahl + 1.17),]
twc68.stahl = tw[tw$Site_ID %in% twc68.stahl$Site_ID,]
twc95 = tw[!mch95,]
twc68 = tw[!mch1s,]
plot(states, col = "light grey")
plot(tw, pch = 19, add = TRUE, cex = 0.5)
plot(twc68, pch = 21, add = TRUE, col = "light blue", bg = "white")
plot(twc95, pch = 21, add = TRUE, col = "blue", bg = "white")
plot(twc68.stahl, add = TRUE, col = "pink")
plot(twc95.stahl, add = TRUE, col = "red")

#A couple of examples
i = match("Gaullup", tw$Site_Name)
gwe3[i,]
gwe.stahl[i]
tw$d18O[i]

i = match("Monmouth", tw$Site_Name)
gwe3[i,]
gwe.stahl[i]
tw$d18O[i]

#stats for shallow GW missmatch
gwe3.stahl68 = gwe3[tw$Site_ID %in% twc68.stahl$Site_ID,]
gwe3.stahl95 = gwe3[tw$Site_ID %in% twc95.stahl$Site_ID,]
mean(apply(apply(gwe3, 1, range, na.rm = TRUE), 2, diff))
mean(apply(apply(gwe3.stahl68, 1, range, na.rm = TRUE), 2, diff))
mean(apply(apply(gwe3.stahl95, 1, range, na.rm = TRUE), 2, diff))
