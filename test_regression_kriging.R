rm(list=ls(all=TRUE))    #Clears R memory
unlink(".RData") 
if (!require("pacman")) install.packages("pacman"); library(pacman) 
p_load(raster, terra, randomForest, RStoolbox, tictoc)
options(warn=1)
cat("Set variables and start processing\n")
Root 		<- 'D:/JKUAT/RESEARCH_Projects/Eswatini/Data/'
Path_out    <- paste0(Root,"Output/")
######################################################
#Sentinel 2
path <- list.files(paste0(Root,'S2/interim/'),pattern = (".tif$"), recursive = TRUE, full.names = TRUE)
s <- rast(path)
names(s) <- c("b", "g","r", "nir")
s
#UAV
mosaicname <- paste0(Path_out,'Mpolonjeni_W1_Mosaic.tif')
if(!file.exists(mosaicname)){
  folders <- list.dirs(paste0(Root,'WingtraOne/Mpolonjeni'),recursive=TRUE)[-1]
  folders 
  for (i in 1:length(folders)) {
    path <- list.files(folders[i], pattern = (".tif$"))
    # Remove all before and up to "reflectance_" in gsub
    path <- path[order(gsub(".*reflectance_","",path))][-4]
    #reorder the bands to match those in S2
    paths <- path
    paths[3] <- path[4]
    paths[4] <- path[3]
    temp <- rast(paste0(folders[i],"/",paths))
    names(temp) <- c("b", "g", "r", "nir")
    assign(paste0("v", i), temp) 
  }
}else{
  v <- rast(mosaicname)
}
######################################################
s_r <- resample(s, v, method='bilinear')

######################################################
#k-means to aid sampling points
image <- v
image[is.na(image)] <- 0
nclass <- 4
nStart <- 2
system.time(
  E <- kmeans(as.data.frame(image, na.rm=F), nclass, iter.max = 200, 
              algorithm = "Lloyd", nstart=nStart)
)
k.map <- image[[1]]
values(k.map) <- E$cluster
k.map[is.na(v[[1]])] <- NA
x11()
plot(k.map)

#####################################################
#Sample from K-means map
set.seed(530)
points <- spatSample(k.map, 900, "stratified", as.points=T, na.rm=T, values=F, xy=TRUE)
x11()
plot(k.map,  main="K-means +sampling points", axes=T, mar = c(4, 5, 1.4, 0.2), col=topo.colors(max(values(k.map),na.rm=T)))
plot(points,add=T)
###############################################################
#Use the stratified randomly sample points to train and predict simulated high resolution S2.

s2_p <- terra::extract(s_r, points)
head(s2_p)
vr_p <- terra::extract(v, points)
head(vr_p)

data <- data.frame(S2=s2_p[,-1], UAV=vr_p[,-1])
pointsDF <- as.data.frame(as(points, "Spatial"))
#Convert to spatvector
temp <- cbind(data,pointsDF)#as(points,"Spatial")
proj4Str <- "+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs"
temp <- SpatialPointsDataFrame(coords      = temp[,c("x","y")], 
                       data        = temp,
                       proj4string = CRS(proj4Str))

##RK using GLM
GLM <- glm(formula = S2.b~UAV.b , data = temp)


glm.pred.test <- predict(GLM, newdata = temp, type="response")
evalData <- sqrt(mean((glm.pred.test - temp$S2.b)^2))

# Ordinary Kriging of GLM residuals
#

statPointsTMP <- temp
statPointsTMP@data <- cbind(statPointsTMP@data, residGLM = resid(GLM))
formMod <- residGLM ~ 1
library(gstat)
mod <- vgm(model  = "Exp", psill  = 0.15, range  = 10, nugget = 0.01)
variog <- variogram(formMod, statPointsTMP)
variogFitOLS <- fit.variogram(variog, model = mod,  fit.method = 6)

# Plot the results
x11()
plot(variog, variogFitOLS, main="Semi-variogram of GLM residuals")
y <- brick(s_r)
names(y) <- c("S2.b", "S2.g","S2.r", "S2.nir")
y_PixDF <- as(y$S2.b, "SpatialPixelsDataFrame")
crs(statPointsTMP) <- crs(y_PixDF)

startTime <- Sys.time() # Start timing
cat("Start time", format(startTime),"\n")
  residKrigMap <- krige(formula = formMod ,
                      locations = statPointsTMP, 
                      model = variogFitOLS,
                      newdata = y_PixDF)
timeDiff <- Sys.time() - startTime
cat("\nCRF Processing time", format(timeDiff), "\n")
# Converting the residual krige map from RasterLayer to SpatRaster  

residKrigRstLayer <- as(residKrigMap, "RasterLayer")
residKrigRstLayer<-as(residKrigRstLayer, "SpatRaster")

# Combining residual krige map with raster object with predicted model values 
x <- v #predictor
names(x) <- c("UAV.b", "UAV.g","UAV.r", "UAV.nir")
rstPredGLM <- predict(x, GLM, type="response") #Predict raster
glmKrigMap <- rstPredGLM + residKrigRstLayer

########################################

