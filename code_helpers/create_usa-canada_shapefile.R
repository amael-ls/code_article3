
#### Aim of prog: create a map of USA and Canada (useful to crop data)
## *Shapefile data
# The shapefiles were downloaded at:
#	1. USA: https://geowebservices.stanford.edu:443/geoserver/wfs?outputformat=SHAPE-ZIP&request=GetFeature&service=wfs&srsName=EPSG%3A4326&typeName=druid%3Avt021tk4894&version=2.0.0
#	2. CAN: https://www12.statcan.gc.ca/census-recensement/alternative_alternatif.cfm?l=eng&dispext=zip&teng=lpr_000b16a_e.zip&k=%20%20%20%2027960&loc=http://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/2016/lpr_000b16a_e.zip

#### Load libraries
library(terra)

#### Common variables
savingPath = "../raw_data/"
if (!dir.exists(savingPath))
	dir.create(savingPath)

#### Load shapefiles
usa = vect("/home/amael/shapefiles/NorthAmerica_shapefiles/usa/usa.shp")
canada = vect("/home/amael/shapefiles/NorthAmerica_shapefiles/canada/canada.shp")

#### Subset USA
## Remove the small islands (based on area)
ls_states = unique(usa$state)
sizePolygons = expanse(x = usa, unit = "km")
threshold = sizePolygons[ls_states == "Hawaii"]

ls_states = ls_states[sizePolygons > threshold]

usa = usa[usa$state %in% ls_states]

## Crop to remove islands beyond -180 (i.e. around +179) that belongs to Alaska
extent180 = ext(c(xmin = -180, xmax = 0, ymin = 24, ymax = 90))
usa = crop(x = usa, y = extent180)

#### Join USA and Canada
## Aggregate
usa = aggregate(usa)
canada = aggregate(canada)

## Project
usa = project(usa, canada)

## Join
northAmerica = rbind(usa, canada) # Yeah, sorry Mexico...

## Simplify
northAmerica = simplifyGeom(x = northAmerica, tolerance = 100)

#### Save 
writeVector(x = northAmerica, filename = paste0(savingPath, "usa-canada.shp"), filetype = "ESRI Shapefile")
