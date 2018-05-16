# Takes shapefiles extracted from the GeoGRatis website
# and converts them into rasters which are saved for use
# by another R script.

# Script created by Kevin Bell, knb226@mun.ca and modified by
# Amy Hurford ahurford@mun.ca
#******************************************************

#The raster and rgdal packages are required
library(rgdal)
library(raster)

#Set the directory
setwd("/Users/amyhurford/Dropbox/Boundary_Conditions/")


#****************************************************************************************************
#Function to take the shapefile and trun it into a raster

makeRaster = function(shp){

#remove excess features from geogratis file
shoreline = shp[shp$dsc_en == "Coast",]

#remove unnessesary attributes from the header
keep = c()
shoreline@data = shoreline@data[keep]

#create a blank raster based on the extent of the shapefile
ext = extent(bbox(shoreline))
r = raster(ext)

# set the raster's resolution
res(r) = c(1/100)

#Create a raster of the shapefile
rastermap = rasterize(shoreline, r)

#change all NA valvues to 1 and all other values to 0
rastermap[is.na(rastermap)]=1
rastermap[rastermap@data@values>1]=0

#Resize the raster to the desired size
rastermapCrop = crop(rastermap,extent(c(-56,-52,45,48)))

return(rastermapCrop)

}

#****************************************************************************************************
#Import the shapefile and covert it to a raster with the makeRaster function

shp = readOGR(dsn="Shoreline shp",layer="shoreline_1")

rastermap = makeRaster(shp)


#****************************************************************************************************
#The original raster values need to be changed to allow the distribution to be plotted over them
rastermap[rastermap@data@values==1]=2
rastermap[rastermap@data@values==0]=1
rastermap[rastermap@data@values==2]=0

plot(rastermap)
save(rastermap, file = "Burin.RData")
