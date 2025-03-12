library(dplyr)
library(terra)
library(sf)
library(h3)
library(geodata)
library(httr)
library(tidyverse)
library(osmdata)

#Region of interest
EA <- st_read("data/shapefiles/EA.shp")
plot(EA)
#Popuklation raster

plot(EA_pop$pop_2020, col=c("blue","red","green"), breaks=seq(0,10000, by=1000))

EA_pop <- rast("output/rasters/covariates/EA_pop.tif")
plot(EA_pop)
EA_pop2020<-EA_pop$pop_2020

plot(EA_pop2020)

# put hexagons around them
EA_merge <- st_union(EA) # merge the EA countries into one polygon
plot(EA_merge)

EA_region_hex <- get_h3_from_sf(
  sf_object = st_as_sf(EA_merge),
  h3_res = 3)

plot(EA_region_hex) 


# filter to only hexes populated areas 
#https://h3geo.org/docs/core-library/restable/

hexes <- EA_region_hex %>%
  sf::st_transform(crs(EA_pop2020)) %>%
  mutate(
    pop = terra::extract(EA_pop2020, ., fun = "sum", na.rm = TRUE)[, 2]
  ) %>%
  filter(
    pop > 0
  ) %>%
  mutate(
    id = match(h3_index, unique(h3_index))
  ) %>%
  dplyr::select(
    -pop
  )

plot(hexes)



# mask by populated area to make the raster need for zonal stats calculations

hex_raster <- rasterize(hexes,
                        EA_pop2020,
                        field = "id",
                        fun = min)
hex_raster <- mask(hex_raster, EA_pop2020)

plot(hex_raster)

# get the version with only populated areas
hex_raster_populated <- hex_raster * EA_pop2020
hex_raster_populated[hex_raster_populated == 0] <- NA
plot(hex_raster_populated)
writeRaster(hex_raster,
            "output/rasters/derived/hex_raster_lookup.tif",
            overwrite = TRUE)
writeRaster(hex_raster_populated,
            "output/rasters/derived/hex_raster_populated_lookup.tif",
            overwrite = TRUE)
saveRDS(hexes, "output/hexes.RDS")

# compare hex raster with populated raster and region buffer
plot(EA_pop$pop_2020, col = viridis::inferno(5)[3:4])
plot(hex_raster, col = viridis::mako(123), add = TRUE)
lines(EA_merge, lwd = 3, col = "darkgoldenrod1")





