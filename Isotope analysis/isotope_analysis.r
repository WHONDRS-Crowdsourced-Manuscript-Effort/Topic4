# attempt to infer sources of streamwater
# code written by AJ Tanentzap (ajt65@cam.ac.uk) on 3 Aug 2022
#############################################################
# 1) PREPARE WORKSPACE
# read in data
mdw <- read.csv('metadata_whondrs.csv')

# load packages
library(rstan)
library(raster)
library(ncdf4)
library(sf)

#############################################################
# 2) CALCULATE D-EXCESS
# d-excess says somethin about precip sources
# more negative values suggest old recharge and greater evaporation
# beneath 10‰ may indicate secondary evapration and more humid air masses
mdw$dexc <- with(mdw, del_2H_permil - 8*del_18O_permil)

# can also calculate LC-excess* that has better interpretation as per
# Landwehr and Coplen 2006 https://apo.ansto.gov.au/dspace/bitstream/10238/10862/1/CSP_26_web.pdf#page=153
# can then interpret in terms of SD differences
mdw$LCe_star <- with(mdw, del_2H_permil - 8.2*del_18O_permil-11.27)/1.15


#############################################################
# 3) USE TWO END-MEMBER MIXING MODEL APPROACH
# sigma_mix = Fraction_A*sigma_A + Fraction_B*sigma_B
# 1 = Fraction_A + Fraction_B
# FA = (sigma_mix - sigma_B) / (sigma_A - sigma_B)

# read in mean-annual amount-weighted precipitation from https://wateriso.utah.edu/waterisotopes/pages/data_access/ArcGrids.html 
d2H_precip <- raster('.\\USPrecip\\d2h_MA.tif')
d18O_precip <- raster('.\\USPrecip\\d18o_MA.tif')
d2H_precip_se <- raster('.\\USPrecip\\d2h_se_MA.tif')  
d18O_precip_se <- raster('.\\USPrecip\\d18o_se_MA.tif')

# read in USA river water values from https://wateriso.utah.edu/waterisotopes/pages/data_access/ArcGrids.html 
d2H_rw <- raster('.\\USSw\\d2h.tif')
d18O_rw <- raster('.\\USSw\\d18o.tif')

# read in groundwater data from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0261651
# use band 5 that is relatively complete and deep
# where band depths correspond with (in m)
# ud = c(1, 10, 25, 50, 100, 200, 500)
# ld = c(10, 25, 50, 100, 200, 500, 2000)
d2H_gw <- raster('.\\gwMapCode-v0.3\\SPATIAL-Lab-gwMapCode-0472cdd\\out\\isoscape_d2H.nc',band=5)
d18O_gw <- raster('.\\gwMapCode-v0.3\\SPATIAL-Lab-gwMapCode-0472cdd\\out\\isoscape_d18O.nc',band=5)
d2H_gw_se <- raster('.\\gwMapCode-v0.3\\SPATIAL-Lab-gwMapCode-0472cdd\\out\\isovar_d2H.nc',band=5)^2
d18O_gw_se <- raster('.\\gwMapCode-v0.3\\SPATIAL-Lab-gwMapCode-0472cdd\\out\\isovar_d18O.nc',band=5)^2

# convert points to spatial object but first put everything in the same projection
d2H_gw <- projectRaster(d2H_gw,crs=crs(d2H_precip))
d18O_gw <- projectRaster(d18O_gw,crs=crs(d2H_precip))
d2H_gw_se <- projectRaster(d2H_gw_se,crs=crs(d2H_precip))
d18O_gw_se <- projectRaster(d18O_gw_se,crs=crs(d2H_precip))
d2H_rw <- projectRaster(d2H_rw,crs=crs(d2H_precip))
d18O_rw <- projectRaster(d18O_rw,crs=crs(d2H_precip))
spdf <- SpatialPointsDataFrame(coords = mdw[,c('long','lat')], data=mdw, proj4string = CRS("+init=epsg:4326"))
spdf <- spTransform(spdf, crs(d2H_gw))

# extract isotope values for points
mdw$d2H_precip <- extract(d2H_precip,spdf)
mdw$d2H_gw <- extract(d2H_gw,spdf)
mdw$d2H_precip_se <- extract(d2H_precip_se,spdf)
mdw$d2H_gw_se <- extract(d2H_gw_se,spdf)
mdw$d2H_rw <- extract(d2H_rw,spdf)
mdw$d18O_precip <- extract(d18O_precip,spdf)
mdw$d18O_gw <- extract(d18O_gw,spdf)
mdw$d18O_precip_se <- extract(d18O_precip_se,spdf)
mdw$d18O_gw_se <- extract(d18O_gw_se,spdf)
mdw$d18O_rw <- extract(d18O_rw,spdf)

# use simple one-isotope models
mdw$FA_d2H <- with(mdw,(del_2H_permil - d2H_gw) / (d2H_precip - d2H_gw))
mdw$FA_d18O <- with(mdw,(del_18O_permil - d18O_gw) / (d18O_precip - d18O_gw))

# many points lie outside the GW/P bounds so have to fit proper model that constrains
stanFWmix <- "data{                                                             // HERE WE DEFINE THE DATA THAT WE INPUT INTO STAN
   int<lower=1> N; 
//   int<lower=1> Nsite;    
   matrix[N,2] yHO;                                                             
   matrix[2,2] dResourcemu[N];                                                  // matrix of means for each isotope * corresponding resource
   matrix[2,2] dResourcesd[N];                                                  // matrix of means for each isotope * corresponding resource
//   int<lower=1> site[N];                                                      // site ID
 }
  parameters{                                                                   // HERE WE DEFINE THE PARAMETERS THAT WE ESTIMATE IN STAN
   real<lower=0> sigma[2];                                                      // estimate 6 SD coefficients and ensure that they are positive
   real<lower=1> lambdaP;
   real<lower=0,upper=1> phiP[N]; 
   corr_matrix[2] OmegaHO; 
   real betab;   
//   real z_site[Nsite];   
  }
transformed parameters{                                                         // HERE WE DEFINE THE PARAMETERS THAT WE CALCULATE IN STAN
   simplex[2] p[N];
   vector[N] linregP;   
   vector<lower=0>[N] phiPalpha;                                                  
   vector<lower=0>[N] phiPbeta;                                                   
   matrix[N,2] muyHO;                                                       
   vector[2] sdyHO[N];                                                   
   matrix[2,2] SigmaYHO[N];                                                   // variance-covariance matrix of consumer isotope values
    for (i in 1:N)
      linregP[i] = inv_logit(betab);// + sigma[3]*z_site[site[i]]);
    phiPalpha = linregP * lambdaP;                                       
    phiPbeta  = (1-linregP) * lambdaP; 
      
    for (i in 1:N){
      p[i,1] = phiP[i];      
      p[i,2] = (1-p[i,1]);
    }
    for (i in 1:N){                                                            // calculate  mean and SD for each consumer observation with 4 isotopes and package SDs into Cholesky-decomposed VCV matrix
      muyHO[i,1] = dot_product(p[i],dResourcemu[i,1]);
      muyHO[i,2] = dot_product(p[i],dResourcemu[i,2]);
      sdyHO[i,1] = sqrt(pow(p[i,1]*dResourcesd[i,1,1],2) + pow(p[i,2]*dResourcesd[i,1,2],2) + pow(sigma[1],2));
      sdyHO[i,2] = sqrt(pow(p[i,1]*dResourcesd[i,2,1],2) + pow(p[i,2]*dResourcesd[i,2,2],2) + pow(sigma[2],2));
      SigmaYHO[i] = quad_form_diag(OmegaHO,sdyHO[i]);
    }
  }
  model{                                                                        // HERE WE DEFINE THE PRIOR LIKELIHOODS THAT WE ESTIMATE IN STAN
     matrix[2,2] LLHO[N];                                                     // Cholesky decomposition of variance-covariance matrix of consumer isotope values
     for (i in 1:N)
      LLHO[i] = cholesky_decompose(SigmaYHO[i]);
     sigma ~ normal(0,1);
     lambdaP ~ pareto(1,2);
//	 z_site ~ normal(0,1);
     betab ~ normal(0,10);
     phiP  ~ beta(phiPalpha,phiPbeta);      
     OmegaHO ~ lkj_corr(2);
     for (i in 1:N)
      yHO[i] ~ multi_normal_cholesky(muyHO[i], LLHO[i]);
  }"


foo1 <- array(c(1), c(length(which(!is.na(mdw$d18O_rw)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip))),2, 2));  
STANdata_f <- list(  yHO = as.matrix(mdw[which(!is.na(mdw$d18O_rw)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),c('d18O_rw','d2H_rw')]),
                     dResourcemu = foo1, dResourcesd = foo1)     
STANdata_f$dResourcemu <- array(data = c(unlist(mdw[which(!is.na(mdw$d18O_rw)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),c('d18O_gw','d2H_gw')]),  
										  unlist(unlist(mdw[which(!is.na(mdw$d18O_rw)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),c('d18O_precip','d2H_precip')]))),
										  dim = c(nrow(foo1),2,2))
STANdata_f$dResourcesd <- array(data = c(unlist(mdw[which(!is.na(mdw$d18O_rw)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),c('d18O_gw_se','d2H_gw_se')]),  
										  unlist(unlist(mdw[which(!is.na(mdw$d18O_rw)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),c('d18O_precip_se','d2H_precip_se')]))),
										  dim = c(nrow(foo1),2,2))
STANdata_f$N <- nrow(STANdata_f$yHO)
STANdata_f$site <- as.numeric(as.factor(as.character(mdw[which(is.na(mdw$d18O_rw)==F&is.na(mdw$d18O_gw)==F&is.na(mdw$d18O_precip)==F),'Site'])))
STANdata_f$Nsite <- max(STANdata_f$site)

# repeat only at 1 est per site as using geospatial data with no replication in estimates
STANdata_f$yHO <- STANdata_f$yHO[match(unique(STANdata_f$site),STANdata_f$site),]
STANdata_f$dResourcemu <- STANdata_f$dResourcemu[match(unique(STANdata_f$site),STANdata_f$site),,]
STANdata_f$dResourcesd <- STANdata_f$dResourcesd[match(unique(STANdata_f$site),STANdata_f$site),,]
STANdata_f$N <- nrow(STANdata_f$yHO)
STANdata_f$site <- NULL
STANdata_f$Nsite <- NULL

its = function(){ list( sigma = rep(1,2), betab = 1, #z_site = rnorm(40), 
						phiP = rep(0.5,40), lambdaP = 10, OmegaHO = diag(2) )}
stanFWmix1_p <- stan(model_code=stanFWmix, data=STANdata_f, init=its, chains=4, iter=8000, warmup=3000, thin=20, refresh=100, cores=4, open_progress=T)

# extract results
mdw$per_gw <- mdw$per_precip <- mdw$per_gw_sd <- mdw$per_precip_sd <- NA
mdw[which(!is.na(mdw$d18O_rw)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),]$per_gw <- colMeans(rstan::extract(stanFWmix1_p,'p')[[1]][,,1])[as.numeric(as.factor(as.character(mdw[which(is.na(mdw$d18O_rw)==F&is.na(mdw$d18O_gw)==F&is.na(mdw$d18O_precip)==F),'Site'])))]
mdw[which(!is.na(mdw$d18O_rw)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),]$per_precip <- colMeans(rstan::extract(stanFWmix1_p,'p')[[1]][,,2])[as.numeric(as.factor(as.character(mdw[which(is.na(mdw$d18O_rw)==F&is.na(mdw$d18O_gw)==F&is.na(mdw$d18O_precip)==F),'Site'])))]
mdw[which(!is.na(mdw$d18O_rw)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),]$per_gw_sd <- apply(rstan::extract(stanFWmix1_p,'p')[[1]][,,1],2,sd)[as.numeric(as.factor(as.character(mdw[which(is.na(mdw$d18O_rw)==F&is.na(mdw$d18O_gw)==F&is.na(mdw$d18O_precip)==F),'Site'])))]
mdw[which(!is.na(mdw$d18O_rw)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),]$per_precip_sd <- apply(rstan::extract(stanFWmix1_p,'p')[[1]][,,2],2,sd)[as.numeric(as.factor(as.character(mdw[which(is.na(mdw$d18O_rw)==F&is.na(mdw$d18O_gw)==F&is.na(mdw$d18O_precip)==F),'Site'])))]



# re-run using actual WHONDRS data
# still all 3 values per rep are identical so can refit previous model
foo1 <- array(c(1), c(length(which(!is.na(mdw$del_2H_permil)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip))),2, 2));  
STANdata_f2 <- list(  yHO = as.matrix(mdw[which(!is.na(mdw$del_2H_permil)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),c('del_18O_permil','del_2H_permil')]),
                     dResourcemu = foo1, dResourcesd = foo1)     
STANdata_f2$dResourcemu <- array(data = c(unlist(mdw[which(!is.na(mdw$del_2H_permil)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),c('d18O_gw','d2H_gw')]),  
										  unlist(unlist(mdw[which(!is.na(mdw$del_2H_permil)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),c('d18O_precip','d2H_precip')]))),
										  dim = c(nrow(foo1),2,2))
STANdata_f2$dResourcesd <- array(data = c(unlist(mdw[which(!is.na(mdw$del_2H_permil)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),c('d18O_gw_se','d2H_gw_se')]),  
										  unlist(unlist(mdw[which(!is.na(mdw$del_2H_permil)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),c('d18O_precip_se','d2H_precip_se')]))),
										  dim = c(nrow(foo1),2,2))
STANdata_f2$site <- as.numeric(as.factor(as.character(mdw[which(is.na(mdw$del_2H_permil)==F&is.na(mdw$d18O_gw)==F&is.na(mdw$d18O_precip)==F),'Site'])))
STANdata_f2$yHO <- STANdata_f2$yHO[match(unique(STANdata_f2$site),STANdata_f2$site),]
STANdata_f2$dResourcemu <- STANdata_f2$dResourcemu[match(unique(STANdata_f2$site),STANdata_f2$site),,]
STANdata_f2$dResourcesd <- STANdata_f2$dResourcesd[match(unique(STANdata_f2$site),STANdata_f2$site),,]
STANdata_f2$N <- nrow(STANdata_f2$yHO)
STANdata_f2$site <- NULL

its = function(){ list( sigma = rep(1,2), betab = 1, 
						phiP = rep(0.5,69), lambdaP = 10, OmegaHO = diag(2) )}
stanFWmix2_p <- stan(model_code=stanFWmix, data=STANdata_f2, init=its, chains=4, iter=8000, warmup=3000, thin=20, refresh=100, cores=4, open_progress=T)

# extract results
mdw$per_gw_WH <- mdw$per_precip_WH <- mdw$per_gw_WH_sd <- mdw$per_precip_WH_sd <- NA
mdw[which(!is.na(mdw$del_2H_permil)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),]$per_gw_WH <- colMeans(rstan::extract(stanFWmix2_p,'p')[[1]][,,1])[as.numeric(as.factor(as.character(mdw[which(is.na(mdw$del_2H_permil)==F&is.na(mdw$d18O_gw)==F&is.na(mdw$d18O_precip)==F),'Site'])))]
mdw[which(!is.na(mdw$del_2H_permil)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),]$per_precip_WH <- colMeans(rstan::extract(stanFWmix2_p,'p')[[1]][,,2])[as.numeric(as.factor(as.character(mdw[which(is.na(mdw$del_2H_permil)==F&is.na(mdw$d18O_gw)==F&is.na(mdw$d18O_precip)==F),'Site'])))]
mdw[which(!is.na(mdw$del_2H_permil)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),]$per_gw_WH <- apply(rstan::extract(stanFWmix2_p,'p')[[1]][,,1],2,sd)[as.numeric(as.factor(as.character(mdw[which(is.na(mdw$del_2H_permil)==F&is.na(mdw$d18O_gw)==F&is.na(mdw$d18O_precip)==F),'Site'])))]
mdw[which(!is.na(mdw$del_2H_permil)&!is.na(mdw$d18O_gw)&!is.na(mdw$d18O_precip)),]$per_precip_WH <- apply(rstan::extract(stanFWmix2_p,'p')[[1]][,,2],2,sd)[as.numeric(as.factor(as.character(mdw[which(is.na(mdw$del_2H_permil)==F&is.na(mdw$d18O_gw)==F&is.na(mdw$d18O_precip)==F),'Site'])))]

# compare estimates
with(mdw, plot(per_gw_WH,per_gw,pch=19)); abline(0,1)

# output results
write.csv(mdw[,c('rownames','Site','Substrate','Rep','dexc','LCe_star','per_gw','per_gw_sd','per_precip','per_precip_sd','per_gw_WH','per_gw_WH_sd','per_precip_WH','per_precip_WH_sd')],file='isotope_whondrs.csv',row.names=F)

