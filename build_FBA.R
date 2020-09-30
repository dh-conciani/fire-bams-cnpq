## dissolve neighbor points to polys

## packages
library (raster)
library (rgdal)
library (rgeos)
library (igraph)
library (sp)
library(RANN)
library(spatialEco)
library (sp)
library (spdep)
library (maptools)

## read lim
#lim <- readOGR ("./esec_serra_geral/shp_lim","esec23s")
## read points
omission_points <- readOGR ("./esec_serra_geral/post_proc", "2006_ommit_points")
## read poly
BA <- readOGR ("./esec_serra_geral/BA_BAMS_RAW", "2006_BA_L7L5_ESEC")
## read reference raster to create cellsize
ref_raster <- raster ("./esec_serra_geral/maxNBRidx/2018_L8L7_MAXNBR_INDEX_ESEC.tif")

## reproject BAMs output to same CRS of points
BA <- spTransform(BA, crs(omission_points))

## calc masks
mask_ref <- raster (crs = projection(ref_raster), ext = extent (ref_raster))
res(mask_ref) = res(ref_raster)
rm(ref_raster);gc()

## extract only polygons to dissolve by neighbour
dis_neigh <- subset (omission_points, POINTID== 0, drop =TRUE)
dis_neigh$POINTID <- as.numeric (1)

## rasterize only to neighbours
r_neigh <- rasterize (dis_neigh, mask_ref, dis_neigh$POINTID)
r_neigh <- clump (r_neigh, directions=8)

## raster to pol
pol_neigh <- rasterToPolygons(r_neigh, dissolve=TRUE)
cols <- colnames(as.data.frame(BA[1:4]))
pol_neigh@data[,2:4] <- c(NA,NA,NA)
colnames(pol_neigh@data) <- cols


## create search points and associate with polygon attributes
rp <- spsample(BA,n=25000,type="regular")
rp2 <- point.in.poly(rp,BA)

# search for nearest point (with radius)
nn <- nn2(coordinates(rp2),coordinates(pol_neigh),k=1,searchtype="radius",radius=5000)$nn.idx
nearestCantons2 <- rp2$date[nn]

## insert dates into created polygons
pol_neigh$date <- nearestCantons2

## create 1,11m buffer into new polygons (to solve blank spaces from datum)
pol_neigh <- gBuffer (pol_neigh, byid=TRUE, width = 0.00001, capStyle= "FLAT", joinStyle = "MITRE",
                mitreLimit=99999)

## bind pols
gBA <- bind (BA, pol_neigh)
gBA@data <- gBA@data[1]

# create empty recipe
fPol <- SpatialPolygons(list())
projection(fPol) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

## convert date to levels
gBA$date <- as.factor(gBA$date)

## start function
list_count = length (levels(gBA$date))
for (i in 1:list_count) {
  ## subset date [i]
  tPol <- subset (gBA, date==levels(gBA$date)[i], drop=FALSE)
  ## dissolve neighbor polygons
  polUnion<- aggregate(tPol, dissolve= TRUE); polUnion<- disaggregate (polUnion)
  ## insert date value into subset
  polUnion$date <- rep(levels(gBA$date)[i], length(polUnion@polygons))
  ## calc area to filter invalid geometries
  polUnion$area <- area (polUnion) ## calc area in ~meters
  ## exclude areas less than 2 pixels 
  polUnion <- subset (polUnion, area>900) ## exclude burns less than 1 pixel
  ## bind to temp* file
  fPol <- bind (fPol, polUnion)
}

## work with single polygons (without date)
## extract single points to dissolve
omission_points$POINTID <- as.numeric (omission_points$POINTID)
dis <- subset (omission_points, POINTID> 1, drop =TRUE)

## rasterize neigbors
r <- rasterize (dis, mask_ref, dis$POINTID)
r <- clump (r, directions=8)

## raster to pol
pol <- rasterToPolygons(r, dissolve=TRUE)
cols <- colnames(as.data.frame(BA[1:4]))
pol@data[,2:4] <- c(NA,NA,NA)
colnames(pol@data) <- cols

## insert date
pol@data <- pol@data[1]
pol$date <- NA
pol$area <- area (pol) ## calc area in ~meters

# bind 
BAf <- bind (fPol, pol)
BAf@data$year <- 2006

## export final BA
writeOGR(BAf, dsn=".", layer= "2006_ESEC_finalBA", driver= "ESRI Shapefile")
