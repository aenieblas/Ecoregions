## plot_bgcp.R
## ae nieblas
## 25/7/2019
## description: plot the existing biogeographical classifications of the Indian Ocean
## inputs : shapefiles of existing classification: longhurst, lme, ppow, meow, tbp, iotc 5x5, FAO 51 and 57?

rm(list=ls())
home<-'/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/'
options(stringsAsFactors=FALSE)

library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(rgdal)
library(sf)
world <- ne_countries(scale = "medium", returnclass = "sf")

## short fin mako



## rfmos
iotc_shp=paste0(home,'indian_ocean/iotc.shp')
iotc <- readOGR(dsn = iotc_shp, stringsAsFactors = F)

swiofc_shp=paste0(home,'indian_ocean/FAO_RFB_SWIOFC/RFB_SWIOFC.shp')
swiofc <- readOGR(dsn = swiofc_shp, stringsAsFactors = F)

siofa_shp='/home/ae/Documents/PERSONAL/COOOL/projects/SIOFA_toothfish_scoping/data/siofa_shapefiles/siofa_subareas_final.shp'
siofa <- readOGR(dsn = siofa_shp, stringsAsFactors = F)

## eezs
eez_shp=paste0(home,'indian_ocean/IOTC-2018-TCAC04-DATA05_-_EEZ_shapefile (1)/IOTC-2018-TCAC04-DATA05 - EEZ shapefile/IOTC_TCAC_EEZ_shapefile.shp')
eez <- readOGR(dsn = eez_shp, stringsAsFactors = F)

## bgcps
long_shp<-paste0(home,'longhurst_shapefiles/Longhurst_world_v4_2010.shp')
long_bgcp<-readOGR(dsn=long_shp, stringsAsFactors = F)

lme_shp<-paste0(home,'LME66/LMEs66.shp')
lme_bgcp<-readOGR(dsn=lme_shp, stringsAsFactors = F)

# meow_shp<-paste0(home,'/MEOW_FINAL/MEOW/meow_ecos.shp')
# meow_shp<-'/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/WCMC036_MEOW_PPOW_2007_2012_v1/DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012.shp'
meow_shp<-'/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/MEOW_FINAL/MEOW/meow_ecos.shp'
meow_bgcp<-readOGR(dsn=meow_shp,stringsAsFactors = F)

meow_ppow_nocoast_shp<-paste0(home,'WCMC036_MEOW_PPOW_2007_2012_v1/DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp')
meow_ppow_nocoast_bgcp<-readOGR(dsn=meow_ppow_nocoast_shp,stringsAsFactors = F)

## mpas
mpa_shp<-paste0(home,'/WDPA_Aug2019_marine-shapefile/')#WDPA_Aug2019_marine-shapefile-points.shp')
layers <- sf::st_layers(mpa_shp)
mpa <- sf::read_sf(mpa_shp, layer = layers$name[1])
part_take<-mpa[which(mpa$NO_TAKE=='Part'),]
no_take<-mpa[which(mpa$NO_TAKE=='All'),]


## Longhurst and LME
map <- ggplot(data = world) +
  # geom_point(data=df3, aes(lon,lat, size = CATCH))+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "grey", fill = "lightgrey")+
  # geom_polygon(data = long_bgcp, aes(x = long, y = lat, group = group), colour = "dark blue", fill = NA)+
  geom_polygon(data = lme_bgcp, aes(x = long, y = lat, group = group), colour = "dark blue", fill = NA)+
  coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)#+
  # annotate("text",label='A',x=20, y=30, size=12)
map + theme_void()
ggsave(paste0('../figures/presentation/biogeography/LME.png'), width = 8, height = 5, dpi=200)
# ggsave(paste0('../figures/presentation/biogeography/Lonhurst.png'), width = 8, height = 5, dpi=200)


# meow_bgcp<-meow_bgcp[meow_bgcp@data$TYPE=='MEOW',]
## MEOW and PPOW
map <- ggplot(data = world) +
  # geom_point(data=df3, aes(lon,lat, size = CATCH))+
  geom_polygon(data = eez, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  # geom_polygon(data = meow_bgcp, aes(x = long, y = lat, group = group), colour = "dark blue", fill = NA)+
  
  geom_polygon(data = meow_ppow_nocoast_bgcp, aes(x = long, y = lat, group = group), colour = "dark blue", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)#+
  # annotate("text",label='B',x=20, y=30, size=12)
map + theme_void()
# ggsave(paste0('../figures/presentation/biogeography/MEOW.png'), width = 8, height = 5, dpi=200)
# ggsave(paste0('../figures/presentation/biogeography/PPOW.png'), width = 8, height = 5, dpi=200)
ggsave(paste0('../figures/presentation/biogeography/PPOW_MEOW_combined.png'), width = 8, height = 5, dpi=200)


## plot no-take MPAs
map <- ggplot(data = world) +
  # geom_point(data=df3, aes(lon,lat, size = CATCH))+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(data = part_take, colour = "cyan", fill = NA)+ 
  geom_sf(data = no_take, colour = "dark blue", fill = NA)+ 
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)+
  annotate("text",label='Chagos',x=64, y=0, size=8)
  # annotate("text",label='A',x=20, y=30, size=12)
map + theme_void()
ggsave(paste0('../figures/presentation/sociopoliticaldata/MPAs.png'), width = 8, height = 5, dpi=200)


## plot rfmos
## force box around swiofc
coords_swiofc<-as.data.frame(matrix(c(30, 10,
                                 65, 10,
                                 65, 0,
                                 80, 0,
                                 80, -45,
                        30,-45,
                        30,10), 
                               ncol = 2, byrow = TRUE))
# names(coords_swiofc)<-c('lon','lat')
P1 = Polygon(coords_swiofc)
Ps1 = SpatialPolygons(list(Polygons(list(P1), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# plot(Ps1, axes = TRUE)

map <- ggplot(data = world) +
  # geom_point(data=df3, aes(lon,lat, size = CATCH))+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_polygon(data = Ps1, aes(x = long, y = lat, group = group), colour = "dark blue", fill = NA)+
  # geom_polygon(data = siofa, aes(x = long, y = lat, group = group), colour = "dark blue", fill = NA)+
  geom_polygon(data = swiofc, aes(x = long, y = lat, group = group), colour = "dark blue", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)#+
  # annotate("text",label='A',x=20, y=30, size=12)
map + theme_void()
ggsave(paste0('../figures/presentation/sociopoliticaldata/SWIOFC.png'), width = 8, height = 5, dpi=200)


## siofa with management subunits
map <- ggplot(data = world) +
  # geom_point(data=df3, aes(lon,lat, size = CATCH))+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_polygon(data = siofa, aes(x = long, y = lat, group = group), colour = "dark blue", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)#+
# annotate("text",label='B',x=20, y=30, size=12)
map + theme_void()
ggsave(paste0('../figures/presentation/sociopoliticaldata/SIOFA.png'), width = 8, height = 5, dpi=200)


## plot ebsas
## ebsas

ebsa_files<-list.files(paste0(home,'/EBSA'))
ebsa_shp<-ebsa_files[grep('.shp',ebsa_files)]

ebsa.df<-c(0,0,0)
names(ebsa.df)<-c('id','Longitude','Latitude')

for(e in 1:length(ebsa_shp)){
  eshp<-paste0(home,'/EBSA/',ebsa_shp[e])
  eval(parse(text=paste0('e',e,'<-readOGR(dsn=eshp,stringsAsFactors = F)')))
  
  eval(parse(text=paste0('ebsa.fort <- fortify(e',e,', region = "EBSA_ID")
  centroids.df<- as.data.frame(coordinates(e',e,'))
                         idList <- e',e,'@data$EBSA_ID
name<- e',e,'@data$NAME
 ebsa.df[',e,',] <- data.frame(name=name,ie = idList, lon= centroids.df$V1, lat=centroids.df$V2)
                         ')))
  
}

map <- ggplot(data = world) +
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_polygon(data = e1, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[1,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e2, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[2,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e3, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[3,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e4, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[4,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e5, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[5,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e6, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[6,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e7, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[7,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e8, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[8,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e9, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[9,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e10, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[10,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e11, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[11,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e12, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[12,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e13, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[13,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e14, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[14,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_polygon(data = e15, aes(x = long, y = lat, group = group), colour = "blue", fill = NA)+
  geom_text(data=ebsa.df[15,],aes(label =ie , x = lon, y = lat)) + #add labels at centroids
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)
  # annotate("text",label='B',x=20, y=30, size=12)
map + theme_void()
# 




worldMap <- readOGR(dsn="data", layer="TM_WORLD_BORDERS_SIMPL-0.3")
# Change "data" to your path in the above!
worldMap.fort <- fortify(world.map, region = "ISO3")
# Fortifying a map makes the data frame ggplot uses to draw the map outlines.
# "region" or "id" identifies those polygons, and links them to your data. 
# Look at head(worldMap@data) to see other choices for id.
# Your data frame needs a column with matching ids to set as the map_id aesthetic in ggplot. 
idList <- worldMap@data$ISO3
# "coordinates" extracts centroids of the polygons, in the order listed at worldMap@data
centroids.df <- as.data.frame(coordinates(worldMap))
names(centroids.df) <- c("Longitude", "Latitude")  #more sensible column names
# This shapefile contained population data, let's plot it.
popList <- worldMap@data$POP2005

pop.df <- data.frame(id = idList, population = popList, centroids.df)


  