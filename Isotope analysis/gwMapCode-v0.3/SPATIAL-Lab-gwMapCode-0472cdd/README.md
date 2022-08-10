# gwMapCode
Data analysis for 3-d groundwater isoscape

Analysis workflow is organized within the following scripts:

1. 1_get_wiDB_data.R - obtain and prepare groundwater and tap water isotopic data
2. 2_prepWellDepths.R - clean and prepare well depth data
3. 3_wd_3d.R - Build 3-d models of well presence/absence and mean isotope values within voxels
4. 4_variograms.R - Build sample semivariograms and fit semivariogram models
5. 5_kriging.R - Kriging prediction
6. 6_plots.R - Generate plots and some stats/analyses for manuscript

The folder "out" contains 4 netCDF files for the 3-d groundwater isoscapes:

1. isoscape_d2H.nc - mean isoscape predictions for hydrogen
2. isoscape_d18O.nc - mean isoscape predictions for oxygen
3. isovar_d2H.nc - total prediction variance (including kriging variance and within-voxel variance) for hydrogen
4. isovar_d18O.nc - total prediction variance (including kriging variance and within-voxel variance) for hydrogen