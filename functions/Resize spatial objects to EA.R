
library(sf)
library(terra)
library(geodata)

# Define the countries and Download the shapefiles for all countries


EA <- st_read("data/shapefiles/EA.shp")
plot(EA)

resizetoEA <- function(x,y=EA){
  crops <- terra::crop(x, y)
  resize<- terra::aggregate(crops, fact = 5, fun="mean")
  maskit <- terra::mask(resize, y)
  return(maskit)
}

resizetoEA2 <- function(x,y=EA){
  crops <- terra::crop(x, y)
  maskit <- terra::mask(crops, y)
  return(maskit)
}

#Resize spatRaster
EA2 <- vect(EA)

VectresizetoEA <- function(x){
  cro <- terra::crop(mask, EA2, mask=TRUE)
  return(cro)
}
