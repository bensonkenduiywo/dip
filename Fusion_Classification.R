#======================================================================================
# Block 1: variable definitions, data import, preparation
#======================================================================================
rm(list=ls(all=TRUE))    #Clears R memory
unlink(".RData") 
if (!require("pacman")) install.packages("pacman"); library(pacman) #package manager( loads required packages/libraries in list as below if not installed they will be installed
p_load(dismo, raster,sp, rgdal, shapefiles, openxlsx, kernlab, randomForest,fields,RStoolbox,popbio)

options(warn=1)
cat("Set variables and start processing\n")
Root 		<- 'D:/JKUAT/RESEARCH_Projects/Eswatini/Data/'
Path_out    <- paste0(Root,"Output/")
setwd(Path_out)
#======================================================================================
# End of Block 1: variable definitions, data import, preparation
#======================================================================================

#======================================================================================
# Block 2: Load Sentinel 2 and UAV images
#======================================================================================
f <- list.files(paste0(Root,'S2/interim/'),pattern = (".tif$"), recursive = TRUE, full.names = TRUE)
f
s2 <- stack(f)
names(s2) <- c("b", "g","r", "nir")
extent(s2)
crs(s2)
res(s2)
nlayers(s2)
#UAV
f <- list.files(paste0(Root,'WingtraOne/'),pattern = (".tif$"), recursive = TRUE, full.names = TRUE)
f
uv <- stack(f)
names(uv) <- c("b", "g","nir", "r")
crs(uv)
extent(uv)
res(uv)
nlayers(uv)
#======================================================================================
# End of Block 2: Load Sentinel 2 and UAV images
#======================================================================================

#======================================================================================
# Block 3: Clip Sentinel 2 to UAV extents
#======================================================================================
s2_clip <- crop(s2, extent(uv),snap='near')
x11()
par(mfrow = c(1, 2), mar = c(4, 5, 1.4, 0.2)) #c(bottom, left, top, right)
plotRGB(s2_clip, r="nir", g="r", b="g", stretch="lin", axes = TRUE, main="S2", cex.axis=0.7)
plotRGB(uv, r="nir", g="r", b="g", stretch="lin", axes = TRUE, main="UAV", cex.axis=0.7)
#======================================================================================
# End of Block 3: Clip Sentinel 2 to UAV extents
#======================================================================================

#======================================================================================
# Block 4: Resample UAV to S2 resolution and fuse
#======================================================================================
uv_r <- resample(uv,s2_clip)

library(terra)
f <- list.files(paste0(Root,'S2/interim/'),pattern = (".tif$"), recursive = TRUE, full.names = TRUE)
s <- rast(f)
res(s)
ext(s)
dim(s)
names(s) <- c("b", "g","r", "nir")
#crop
s <- crop(s, ext(v), snap="near")

f <- list.files(paste0(Root,'WingtraOne/'),pattern = (".tif$"), recursive = TRUE, full.names = TRUE)
f
v <- rast(f)
res(v)
ext(v)
dim(v)
names(v) <- c("b", "g","nir", "r")
x11()
par(mfrow = c(1, 2), mar = c(4, 5, 1.4, 0.2)) #c(bottom, left, top, right)
plotRGB(s, r="nir", g="r", b="g", stretch="lin", axes = TRUE, main="S2", cex.axis=0.7)
plotRGB(v, r="nir", g="r", b="g", stretch="lin", axes = TRUE, main="UAV", cex.axis=0.7)

# Resample UAV to S2 resolution for spectral fusion
v_r <- resample(v,s,method='bilinear')
ext(v_r)==ext(s)
#Fuse by multiplication
f <- s * v_r
f
#Display fused image alongside original UAV
x11()
par(mfrow = c(1, 2), mar = c(4, 5, 1.4, 0.2)) #c(bottom, left, top, right)
plotRGB(v, r="nir", g="r", b="g", stretch="lin", axes = TRUE, main="UAV", cex.axis=0.7)
plotRGB(f, r="nir", g="r", b="g", stretch="lin", axes = TRUE, main="UAV_Fused", cex.axis=0.7)


