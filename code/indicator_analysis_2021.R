## calculates the specificity and fidelity indicators that the ecoregions are based on

## extract_ppow


rm(list=ls())
options(stringsAsFactors = FALSE)

analysis<-'PPOWMEOW' #options are : PPOW,PPOWMEOW,LONG
reassign<-FALSE # for PPOW analysis - TRUE indicates that some MEOW provinces will be reassigned as PPOW
                # for PPOWMEOW analysis - TRUE indicates that PPOW provinces will be reassinged as MEOW.
layer<-'fishery' # options are : species, fishery, size
start_year<-2005
end_year<-2019

source('~/Documents/PERSONAL/COOOL/projects/Ecoregions_2021/code/iotc_ca_regrid.R')
source('~/Documents/PERSONAL/COOOL/projects/Ecoregions_2021/code/iotc_regrid.R')
source('~/Documents/PERSONAL/COOOL/projects/Ecoregions_2021/code/regrid_datalayers.R')
home<-'/home/ae/Documents/PERSONAL/COOOL/projects/Ecoregions_2021/data/'
setwd(home)
library(pacman)
p_load('rgdal','ggplot2','sf','rnaturalearth','rnaturalearthdata','reshape')
# p_load('')
library(rgdal)
library(rgeos)
library(plyr)
library(sp)
library(tidyr)
library(raster)
library(plotly)
library("forcats")
theme_set(theme_bw())

# world data for plotting
world <- ne_countries(scale = "medium", returnclass = "sf")

# iotc convention area 
iotc_shp='/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/indian_ocean/iotc.shp'
iotc <- readOGR(dsn = iotc_shp, stringsAsFactors = F)


# STEP 1 : REGRID DATA LAYERS TO IOTC 5X5 GRID and keep only those within the IOTC convention area
iotc_5deg_shp=paste0(home,'/IOTC-2018-TCAC04-DATA04_-_5x5_grids_shapefile/IOTC_TCAC_GRID_5x5.shp')
shp <- readOGR(dsn = iotc_5deg_shp, stringsAsFactors = F)

############################################## read in bcp shapefile #######################################################################

if(analysis=='PPOW'|analysis=='PPOWMEOW'){
## put ppow on 5x5 iotc grid
ppow_shp<-paste0(home,'WCMC036_MEOW_PPOW_2007_2012_v1/DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp')
# ppow_shp<-'/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/WCMC036_MEOW_PPOW_2007_2012_v1/DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012.shp'
ppow<-readOGR(dsn=ppow_shp,stringsAsFactors = F)
}

if(analysis=='PPOWMEOW'){
  meow_shp<-'/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/MEOW_FINAL/MEOW/meow_ecos.shp'
  meow<-readOGR(dsn=meow_shp,stringsAsFactors = F)
}

# if(analysis=='LONG'){
# long_shp<-paste0(home,'longhurst_shapefiles/Longhurst_world_v4_2010.shp')
# ppow<-readOGR(dsn=long_shp, stringsAsFactors = F)
# }

wgs<-proj4string(ppow) #this shows the projection, it's important that it's the same as projection of our data that we convert in shapefile later


########################### DATA FROM MJJJ METHODS ######################################
# rm(list=ls())

## regrid all 1 degrees to 5 degrees <-------------------------------------------------------- AE HERE UPDATE THE DATA FILE --------###!!!!!!!!!!!!
d1<- iotc_1_to_5_regrid(iotc_5deg_shp=iotc_5deg_shp,
                        iotc_df=paste0(home,'EAFM_data_v2.csv'),data=2017)

## data filtering
# outliers (points >95th percentile of total catch)
raised_catch<-d1
raised_catch<-raised_catch[-which(raised_catch$RAISED_CATCHES_MT>quantile(raised_catch$RAISED_CATCHES_MT,probs=0.95,na.rm=T)),]

# tropical species found <45
levels(as.factor(raised_catch$SPECIES))
raised_catch$Species2<-as.factor(raised_catch$SPECIES)
levels(raised_catch$Species2)<-c("Temperate","Tropical","Tropical","Subtropical","Tropical")

raised_catch[which(raised_catch$Species2=='Tropical' & raised_catch$lat<(-45)),'RAISED_CATCHES_MT']<-NA

raised_catch[which(raised_catch$Species2=='Tropical' & raised_catch$lat<(-45)),'RAISED_CATCHES_MT']<-NA
# which()

library(plyr)

## sum catch within each grid cell
d2<-ddply(raised_catch,.(SPECIES, YEAR, GEAR,OPERATION_TYPE,FLEET_CODE_HASHED,size,quadrant,lat,lon), function(x) sum(x$RAISED_CATCHES_MT,na.rm=T))
# m2<-ddply(d1,.(SPECIES, MONTH, GEAR,OPERATION_TYPE,FLEET_CODE_HASHED,size,quadrant,lat,lon), function(x) sum(x$RAISED_CATCHES_MT,na.rm=T))

## adapt to rtunatlas df format
d2$geographic_identifier<-paste0(d2$size,d2$quadrant,sprintf("%02d",(abs(d2$lat)-2.5)),sprintf("%03d",(d2$lon-2.5)))
d2$unit<-'MT'
d2$value<-d2$V1

## remove package plyr
detach("package:plyr", unload=TRUE) 

library(dplyr)

##### Examine the catch species compostion/location of the different fleets/type of fisheries

head(d2) ## average across all years available (2006-2016)

datos3<-d2 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(SPECIES, lon,lat, FLEET_CODE_HASHED,GEAR,OPERATION_TYPE,geographic_identifier) %>% summarize (catch=median(value,na.rm=TRUE))
obs3<-d2 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(SPECIES, lon,lat) %>% tally()
head(datos3)

names(datos3)<-c("species","lon", "lat","fleet","fishery", "operation",'geographic_identifier',"catch")
datos3$unit<-'MT'

# convert to spdf and over function
coordinates(datos3)<-~lon+lat
proj4string(iotc) <- CRS(wgs)
iotc.polygons <- spTransform(iotc, CRS("+proj=longlat"))
proj4string(datos3) <- CRS(proj4string(iotc))

# now do OVER and bind the data to the dataframe
io<-sp::over(datos3,iotc)
datos3@data = cbind(as.data.frame(datos3@data),io) 
summary(datos3)### SEE THE RESULTS FROM SUMMARY AND COLUMN TYPE
datos_m<-subset(datos3,!is.na(datos3$RFB))


datos_m$Latitude<-datos_m@coords[,2]
datos_m$Longitude<-datos_m@coords[,1]

datos<-as.data.frame(datos_m)
datos<-datos[,c('species','lon','lat','fleet','fishery','operation','catch')]
names(datos)<-c('species','Long','Lat','fleet','fishery','operation','tCatch')


############################################# bind ppow to catch data ####################################################

if(analysis=='PPOWMEOW'){
  # convert to spdf and over function
  coordinates(datos)<-~Long+Lat
  proj4string(meow) <- CRS(wgs)
  meow.polygons <- spTransform(meow, CRS("+proj=longlat"))
  proj4string(datos) <- CRS(proj4string(meow))
  
  # now do OVER and bind the data to the dataframe
  o<-over(datos,meow)
  datos@data = cbind(datos@data,o) 
  summary(datos)### SEE THE RESULTS FROM SUMMARY AND COLUMN TYPE
  datos_m<-subset(datos,!is.na(datos$PROVINCE))
  datos_m$TYPE<-'MEOW'
  datos_p<-subset(datos,is.na(datos$PROVINCE))
  # datos_p$TYPE<-'PPOW'
  
  # coordinates(datos_p)<-~Long+Lat
  proj4string(ppow) <- CRS(wgs)
  ppow.polygons <- spTransform(ppow, CRS("+proj=longlat"))
  proj4string(datos_p) <- CRS(proj4string(ppow))
  
  # now do OVER and bind the data to the dataframe
  o<-over(datos_p,ppow)
  datos_p@data = cbind(datos_p@data,o) 
  summary(datos_p)### SEE THE RESULTS FROM SUMMARY AND COLUMN TYPE
  
  library(maptools)
  dropsM <-c('ECO_CODE','PROV_CODE','RLM_CODE','ALT_CODE','ECO_CODE_X','Lat_Zone') # list of col names
  dropsP<-c('ECO_CODE','PROV_CODE','RLM_CODE','ALT_CODE','ECO_CODE_X','Lat_Zone','BIOME')
  # dropsP<-c('ECOREGION','REALM')
  datos_m <- datos_m[,!(names(datos_m) %in% dropsM)] #remove columns "AREA" and "PERIMETER"
  datos_p <- datos_p[,!(names(datos_p) %in% dropsP)] #remove columns "AREA" and "PERIMETER"
  
  dropsPP<-c('ECOREGION','PROVINCE','REALM')
  
  datos_p <- datos_p[,!(names(datos_p) %in% dropsPP)] #remove columns "AREA" and "PERIMETER"
  colnames(datos_p@data) <-colnames(datos_m@data)
  
  datos<-spRbind(datos_m,datos_p)
}else{
# convert to spdf and over function
coordinates(datos)<-~Long+Lat
proj4string(ppow) <- CRS(wgs)
ppow.polygons <- spTransform(ppow, CRS("+proj=longlat"))
proj4string(datos) <- CRS(proj4string(ppow))

# now do OVER and bind the data to the dataframe
o<-over(datos,ppow)
datos@data = cbind(datos@data,o) 
summary(datos)### SEE THE RESULTS FROM SUMMARY AND COLUMN TYPE
}
# if(analysis=='PPOW'){write.csv(datos, "CatchesbyPPOWpoly_v4.csv", row.names=FALSE)}
if(analysis=='PPOWMEOW'){write.csv(datos, "../data/CatchesbyPPOWMEOWpoly_v2021.csv", row.names=FALSE)}
# if(analysis=='LONG'){write.csv(datos, "CatchesbyLonghurstpoly_v4.csv", row.names=FALSE)}


############################################# identify provinces ##############################################
#  regroup the MEOWs and PPOWs according to these rules
# 1) isolated pixels will join the larger grouping (e.g. Amsterdam St Paul)
# 2) coastal MEOWs will be combined to the coastal PPOW
# 3) provinces where no catch exists will be removed


if(analysis=='PPOWMEOW'){
  catchesbyPPOW <-read.csv("../data/CatchesbyPPOWMEOWpoly_v2021.csv")
  
  catchesbyPPOW$lon<-catchesbyPPOW$Long
  catchesbyPPOW$lat<-catchesbyPPOW$Lat
  
  
  ## 1 - single pixel regroupings ------------------------------------------------------------
  # Sunda Shelf > Andaman
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Sunda Shelf'),'PROVINCE']<-unique(catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Andaman'),'PROVINCE'])
  catchesbyPPOW[which(catchesbyPPOW$ECOREGION=='Sunda Shelf'),'ECOREGION']<-unique(catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Andaman'),'ECOREGION'])[1]
  catchesbyPPOW[which(catchesbyPPOW$ECOREGION=='Sunda Shelf'),'REALM']<-unique(catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Andaman'),'REALM'])[2]
  catchesbyPPOW[which(catchesbyPPOW$ECOREGION=='Sunda Shelf'),'TYPE']<-unique(catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Andaman'),'TYPE'])
  
  # Amsterdam St Paul > Indian Ocean Gyre
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Amsterdam-St Paul'),'PROVINCE']<-'Indian Ocean Gyre'
  catchesbyPPOW[which(catchesbyPPOW$ECOREGION=='Sunda Shelf'),'ECOREGION']<-'Indian Ocean Gyre'
  catchesbyPPOW[which(catchesbyPPOW$ECOREGION=='Sunda Shelf'),'REALM']<-'Indian Ocean Gyre'
  catchesbyPPOW[which(catchesbyPPOW$ECOREGION=='Sunda Shelf'),'TYPE']<-'PPOW'
  
  
  # 2 - coastal MEOWs will be combined to the nearest coastal PPOW -----------------------------------
  # agulhas + agulhas current = agulhas current
  ac_LL<-data.frame(lat=c(-37.5,-37.5,-37.5,
                          -32.5,-32.5,-32.5,-32.5,
                          -27.5,-27.5),
                    lon=c(22.5,27.5,32.5,
                          22.5,27.5,32.5,37.5,
                          32.5,37.5))
  ac<-NULL
  for(i in 1:dim(ac_LL)[1]){
    ac<-c(ac,which(catchesbyPPOW$lon==ac_LL$lon[i] & catchesbyPPOW$lat==ac_LL$lat[i] ))
  }
  catchesbyPPOW[ac,'ECOREGION']<-'Agulhas Current'
  catchesbyPPOW[ac,'PROVINCE']<-'Agulhas Current'
  catchesbyPPOW[ac,'TYPE']<-'PPOW'
  
  ## somali/arabia + somali current = somali current
  som_LL<-data.frame(lat=c(2.5,
                           7.5, 7.5,
                           12.5,12.5,12.5,12.5,
                           17.5,17.5,17.5,17.5,17.5,17.5,
                           22.5,22.5,22.5,22.5,
                           27.5,27.5,27.5,27.5,27.5,27.5),
                     lon=c(47.5,
                           47.5,52.5,
                           42.5,47.5,52.5,57.5,
                           37.5,42.5,52.5,57.5,62.5,67.5,
                           37.5,52.5,57.5,62.5,
                           32.5,37.5,47.5,52.5,57.5,62.5
                     ))
  somali<-NULL
  for(i in 1:dim(som_LL)[1]){
    somali<-c(somali,which(catchesbyPPOW$lon==som_LL$lon[i] & catchesbyPPOW$lat==som_LL$lat[i] ))
  }
  catchesbyPPOW[somali,'ECOREGION']<-'Somali Current'
  catchesbyPPOW[somali,'PROVINCE']<-'Somali Current'
  catchesbyPPOW[somali,'TYPE']<-'PPOW'
  
  ## java transitional + IOMG > Java Transitional
  java<-c(which(catchesbyPPOW$lon==102.5 & catchesbyPPOW$lat==(-12.5) ),
          which(catchesbyPPOW$lon==97.5 & catchesbyPPOW$lat==(-7.5) ))
  catchesbyPPOW[java,'ECOREGION']<-'Java Transitional'
  catchesbyPPOW[java,'PROVINCE']<-'Java Transitional'
  catchesbyPPOW[java,'TYPE']<-'MEOW'
  
  
  # indonesian throughflow
  itf<-c(which(catchesbyPPOW$lon==112.5 & catchesbyPPOW$lat==(-12.5) ),
         which(catchesbyPPOW$lon==117.5 & catchesbyPPOW$lat==(-12.5) ),
         which(catchesbyPPOW$lon==122.5 & catchesbyPPOW$lat==(-12.5) ),
         which(catchesbyPPOW$lon==127.5 & catchesbyPPOW$lat==(-12.5) ))
  catchesbyPPOW[itf,'ECOREGION']<-'Indonesian Through-Flow'
  catchesbyPPOW[itf,'PROVINCE']<-'Indonesian Through-Flow'
  catchesbyPPOW[itf,'TYPE']<-'PPOW'
  
  lc<-c(which(catchesbyPPOW$lon==112.5 & catchesbyPPOW$lat==(-17.5) ),
        which(catchesbyPPOW$lon==117.5 & catchesbyPPOW$lat==(-17.5) ),
        which(catchesbyPPOW$lon==117.5 & catchesbyPPOW$lat==(-22.5) ),
        which(catchesbyPPOW$lon==112.5 & catchesbyPPOW$lat==(-22.5) ),
        which(catchesbyPPOW$lon==112.5 & catchesbyPPOW$lat==(-27.5) ),
        which(catchesbyPPOW$lon==112.5 & catchesbyPPOW$lat==(-32.5),
              which(catchesbyPPOW$lon==112.5 & catchesbyPPOW$lat==(-37.5)  )))
  catchesbyPPOW[lc,'ECOREGION']<-'Leeuwin Current'
  catchesbyPPOW[lc,'PROVINCE']<-'Leeuwin Current'
  catchesbyPPOW[lc,'TYPE']<-'PPOW'
  
  
  ## reassigning province names to realm names for N and S
  iog<-which(catchesbyPPOW$PROVINCE=='Indo-Pacific Warm Water' & catchesbyPPOW$lat<(-12.5) )
  catchesbyPPOW[iog,'ECOREGION']<-'Indian Ocean Gyre'
  catchesbyPPOW[iog,'PROVINCE']<-'Indian Ocean Gyre'
  catchesbyPPOW[iog,'TYPE']<-'PPOW'
  
  iomg<-which(catchesbyPPOW$PROVINCE=='Indo-Pacific Warm Water' & catchesbyPPOW$lat>=(-12.5) )
  catchesbyPPOW[iomg,'ECOREGION']<-'Indian Ocean Monsoon Gyre'
  catchesbyPPOW[iomg,'PROVINCE']<-'Indian Ocean Monsoon Gyre'
  catchesbyPPOW[iomg,'TYPE']<-'PPOW'
  
  
  sant<-c(which(catchesbyPPOW$lon<60 & catchesbyPPOW$lat==(-42.5) ),
          which(catchesbyPPOW$lon>97 & catchesbyPPOW$lon<140  & catchesbyPPOW$lat==(-47.5) ))
  catchesbyPPOW[sant,'ECOREGION']<-'Subantarctic'
  catchesbyPPOW[sant,'PROVINCE']<-'Subantarctic'
  catchesbyPPOW[sant,'TYPE']<-'PPOW'
  
  stc<-c(which(catchesbyPPOW$lon>58 & catchesbyPPOW$lat==(-42.5) ),
         which(catchesbyPPOW$lon>115 & catchesbyPPOW$lat==(-37.5) ),
         which(catchesbyPPOW$lon>140 & catchesbyPPOW$lat==(-47.5) ))
  catchesbyPPOW[stc,'ECOREGION']<-'Subtropical Convergence'
  catchesbyPPOW[stc,'PROVINCE']<-'Subtropical Convergence'
  catchesbyPPOW[stc,'TYPE']<-'PPOW'
  
  apf<-c(which(catchesbyPPOW$lon<93 & catchesbyPPOW$lat==(-47.5) ),
         which(catchesbyPPOW$lon>93& catchesbyPPOW$lat==(-52.5) ))
  catchesbyPPOW[apf,'ECOREGION']<-'Antarctic Polar Front'
  catchesbyPPOW[apf,'PROVINCE']<-'Antarctic Polar Front'
  catchesbyPPOW[apf,'TYPE']<-'PPOW'
  
  
  # 3) remove provinces with no catch ---------------------------------------------------------
  ant<-c(which(catchesbyPPOW$lon>80 & catchesbyPPOW$lon<97 & catchesbyPPOW$lat==(-52.5) ))
  if(sum(catchesbyPPOW[ant,'tCatch'])==0){catchesbyPPOW<-catchesbyPPOW[-ant,]}
  
  
  # 4) order the resulting provinces for better visualisation
  unique(catchesbyPPOW[,'PROVINCE'])
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Agulhas Current'),'clockwise']<-1
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Western Indian Ocean'),'clockwise']<-2
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Somali Current'),'clockwise']<-3
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='West and South Indian Shelf'),'clockwise']<-4
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Central Indian Ocean Islands'),'clockwise']<-5
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Bay of Bengal'),'clockwise']<-6
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Andaman'),'clockwise']<-7
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Java Transitional'),'clockwise']<-8
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Indonesian Through-Flow'),'clockwise']<-9
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Leeuwin Current'),'clockwise']<-10
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Southwest Australian Shelf'),'clockwise']<-11
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Subtropical Convergence'),'clockwise']<-12
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Subantarctic'),'clockwise']<-13
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Indian Ocean Monsoon Gyre'),'clockwise']<-14
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Indian Ocean Gyre'),'clockwise']<-15
  
  
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Agulhas Current'),'NSLR']<-11
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Western Indian Ocean'),'NSLR']<-1
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Somali Current'),'NSLR']<-2
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='West and South Indian Shelf'),'NSLR']<-3
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Central Indian Ocean Islands'),'NSLR']<-4
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Bay of Bengal'),'NSLR']<-5
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Andaman'),'NSLR']<-6
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Java Transitional'),'NSLR']<-7
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Indonesian Through-Flow'),'NSLR']<-8
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Leeuwin Current'),'NSLR']<-15
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Southwest Australian Shelf'),'NSLR']<-14
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Subtropical Convergence'),'NSLR']<-13
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Subantarctic'),'NSLR']<-12
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Indian Ocean Monsoon Gyre'),'NSLR']<-9
  catchesbyPPOW[which(catchesbyPPOW$PROVINCE=='Indian Ocean Gyre'),'NSLR']<-10
  

  # group by coastal/offshore and north/south
coastalN<-c('Western Indian Ocean','Somali Current','West and South Indian Shelf','Central Indian Ocean Islands','Bay of Bengal','Andaman',
            'Java Transitional','Indonesian Through-Flow')
coastalS<-c('Agulhas Current','Leeuwin Current','Southwest Australian Shelf')
offshoreN<-c('Indian Ocean Monsoon Gyre')
offshoreS<-c('Subtropical Convergence','Subantarctic','Indian Ocean Gyre')
  
catchesbyPPOW$PROVINCE2<-catchesbyPPOW$PROVINCE

catchesbyPPOW$PROVINCE2<-fct_collapse (catchesbyPPOW$PROVINCE2,
                                      coastalN=coastalN,
                                      coastalS=coastalS,
                                      offshoreN=offshoreN,
                                      offshoreS=offshoreS)
  
  
  catchesbyPPOW<-catchesbyPPOW[order(catchesbyPPOW$NSLR),]
  
  datos<-catchesbyPPOW[,c('species','fleet','fishery','operation','tCatch','ECOREGION','PROVINCE','REALM','TYPE','Long','Lat','optional','clockwise','NSLR','PROVINCE2')] 
  names(datos)<-c('Species','Fleet','Fishery','Operation','tCatch','ECOREGION','PROVINC','REALM','Type','Long','Lat','optional','clockwise','NSLR','PROVINCE2')
  
  
  
  write.csv(datos,file='../data/PPOWMEOW_data.csv',row.names = FALSE)
  
}

## plot resulting groupings

datos<-read.csv('../data/PPOWMEOW_data.csv',stringsAsFactors = FALSE)

##TOTALCATCHES
detach("package:plyr", unload=TRUE)
library(dplyr)
# datos<-datos[,1:11]
Tcatcheslonlat<-datos %>% group_by(Long,Lat,NSLR,PROVINC,PROVINCE2) %>% summarize (tCatch=sum(tCatch,na.rm=TRUE))

## N/S L/R starting from WIO
ggplot(data = world) +
  # geom_tile(data = Tcatcheslonlat, aes(x=Long,y=Lat,fill=as.factor(NSLR)))+
  geom_tile(data = Tcatcheslonlat, aes(x=Long,y=Lat,fill=as.factor(PROVINC)))+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_polygon(data = meow, aes(x = long, y = lat, group = group), colour = "black", fill = NA)+
  geom_polygon(data = ppow, aes(x = long, y = lat, group = group), colour = "black", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  geom_point(data=Tcatcheslonlat, aes(Long,Lat, size = tCatch))+xlab('Lon')+ylab('Lat')+
  coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)


## Coastal/Offshore N/S
library(RColorBrewer)
cols<-brewer.pal(8,'YlGnBu')
library(viridis) 
ggplot(data = world) +
  geom_tile(data = Tcatcheslonlat, aes(x=Long,y=Lat,fill=as.factor(PROVINCE2)))+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_polygon(data = meow, aes(x = long, y = lat, group = group), colour = "black", fill = NA)+
  geom_polygon(data = ppow, aes(x = long, y = lat, group = group), colour = "black", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  # geom_point(data=Tcatcheslonlat, aes(Long,Lat, size = tCatch))+xlab('Lon')+ylab('Lat')+
  guides(color=guide_legend("my title"))+
  coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)
ggsave(paste0('../figures/spatial/catch_dist_',layer,'_',analysis,'_regrouped_v2021.png'), width = 8, height = 5, dpi=200)




########################################### fisheries codes and groupings ########################################################
#group_by_major_fisheries
LL<-c('ELL','FLL','LL','LLCO','LLEX','SLL','LG')
PS<-c('PS','PSS','RIN','RNOF')
GILL<-c('GILL','GIOF','GL')
LINE<-c('HAND','HLOF','SPOR','TROL','TROLM','RR')
BB<-c('BB','BBOF')
OT<-c('BS','CN','FN','LIFT','TRAP','TRAW','HARP')
DSEI<-c('DSEI')

#group by large fishery and type of operation
LL_IND<-c('FLL','LL','LLEX','LG')# SLL??
LL_IND_SWO<-'ELL'
LL_IND_SHK<-'SLL'
LL_ART<-c('LLCO')

PS_IND<-c('PS','RNOF')
PS_ART<-c('PSS','RIN')

GILL_IND<-c('GIOF')
GILL_ART<-c('GILL','GL')

LINE_IND<-c('HLOF')
LINE_ART<-c('HAND','SPOR','TROL','TROLM','RR')

BB_ART<-c('BB')
BB_IND<-c('BBOF')

OT_ART<-c('BS','CN','FN','LIFT','TRAP','TRAW','HARP')

DSEI_ART<-c('DSEI')

head(datos)

library(dplyr)
library(ggplot2)
fishery<-datos %>% group_by(Long,Lat,Fishery,PROVINC) %>% summarize (catch=sum(tCatch,na.rm=TRUE)) %>% arrange(-catch)
# head(fishery)

library("forcats")

fishery$fishery2<-fishery$Fishery

fishery$fishery2<-fct_collapse (fishery$fishery2,
                                LL=LL,
                                PS=PS,
                                PS=PS,
                                GILL=GILL,
                                LINE=LINE,
                                BB=BB,
                                OT=OT,
                                DSEI=DSEI)


fishery$fishery3<-fishery$Fishery

fishery$fishery3<-fct_collapse (fishery$fishery3,
                                LL_IND =LL_IND,
                                LL_IND_SWO =LL_IND_SWO,
                                LL_ART =LL_ART,
                                PS_IND=PS_IND,
                                PS_ART=PS_ART,
                                GILL_IND=GILL_IND,
                                GILL_ART=GILL_ART,
                                LINE_IND=LINE_IND,
                                LINE_ART=LINE_ART,
                                BB_ART=BB_ART,
                                BB_IND=BB_IND,
                                OT_ART=OT_ART,
                                DSEI_ART=DSEI_ART)



fisherylevels<-fishery %>% group_by(fishery3) %>% summarize(catch=sum(catch)) %>%arrange(-catch)
fisherylevels
unique(fisherylevels$fishery3)

fishery$fishery3<-factor(fishery$fishery3,levels=unique(fisherylevels$fishery3))

g<-ggplot(fishery,aes(fishery3))
g+geom_bar(aes(weight=catch/10000, fill=fishery2))+ylab('catch (MT/10,000)')+coord_flip()


#plot for LL

fisheryLL<- filter(fishery, fishery2=="LL")


fisheryLL[which(fisheryLL$catch==0),'catch']<-NA

ggplot(data = world) + 
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ 
  geom_tile(data=fisheryLL, aes(Long,Lat, fill = catch)) +
  facet_wrap(~Fishery,ncol=4)+
  ggtitle(" Total raised catches - major fisheryLL and operation")+
  scale_fill_gradient(low = "lightyellow",
                      high = "steelblue")

ggsave("totalcatches_fisheryLL.png", width = 10, height = 5, dpi=200)

#plot for PS

fisheryPS<- filter(fishery, fishery2=="PS")


fisheryPS[which(fisheryPS$catch==0),'catch']<-NA

summary(fisheryPS)

ggplot(data = world) + 
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ 
  geom_tile(data=fisheryPS, aes(Long,Lat, fill = catch)) +
  facet_wrap(~Fishery,ncol=4)+  ggtitle(" Total raised catches - major fisheryPS and operation")+scale_fill_gradient(low = "lightyellow",
                                                                                                                     high = "steelblue")

ggsave("totalcatches_fisheriesPS.png", width = 10, height = 5, dpi=200)



#plot for GILL

fisheryGILL<- filter(fishery, fishery2=="GILL")


fisheryGILL[which(fisheryGILL$catch==0),'catch']<-NA

ggplot(data = world) + 
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ 
  geom_tile(data=fisheryGILL, aes(Long,Lat, fill = catch)) +
  facet_wrap(~Fishery,ncol=4)+
  ggtitle(" Total raised catches - major fisheryGILL and operation")+scale_fill_gradient(low = "lightyellow",high = "steelblue")

ggsave("totalcatches_fisheryGILL.png", width = 10, height = 5, dpi=200)



#plot for LINE

fisheryLINE<- filter(fishery, fishery2=="LINE")


fisheryLINE[which(fisheryLINE$catch==0),'catch']<-NA

ggplot(data = world) + 
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ 
  geom_tile(data=fisheryLINE, aes(Long,Lat, fill = catch)) +
  facet_wrap(~Fishery,ncol=3)+
  ggtitle(" Total raised catches - major fisheryLINE and operation")+
  scale_fill_gradient(low = "lightyellow",
                      high = "steelblue")

ggsave("totalcatches_fisheryLINE.png", width = 10, height = 5, dpi=200)



#plot for Baitboat

fisheryBB<- filter(fishery, fishery2=="BB")


fisheryBB[which(fisheryBB$catch==0),'catch']<-NA

ggplot(data = world) + 
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ 
  geom_tile(data=fisheryBB, aes(Long,Lat, fill = catch)) +
  facet_wrap(~Fishery,ncol=4)+
  ggtitle(" Total raised catches - major fisheryBB and operation")+
  scale_fill_gradient(low = "lightyellow",
                      high = "steelblue")

ggsave("totalcatches_fisheryBB.png", width = 10, height = 5, dpi=200)


#plot for OT_ART

fisheryOT<- filter(fishery, fishery2=="OT")


fisheryOT[which(fisheryOT$catch==0),'catch']<-NA

ggplot(data = world) + 
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ 
  geom_tile(data=fisheryOT, aes(Long,Lat, fill = catch)) +
  facet_wrap(~Fishery,ncol=4)+
  ggtitle(" Total raised catches - major fisheryOT and operation")+
  scale_fill_gradient(low = "lightyellow",
                      high = "steelblue")

ggsave("totalcatches_fisheryOT.png", width = 10, height = 5, dpi=200)


#plot for DSEI_ART

fisheryDSEI<- filter(fishery, fishery2=="DSEI")


fisheryDSEI[which(fisheryDSEI$catch==0),'catch']<-NA

ggplot(data = world) + 
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ 
  geom_tile(data=fisheryDSEI, aes(Long,Lat, fill = catch)) +
  facet_wrap(~Fishery,ncol=4)+
  ggtitle(" Total raised catches - major fisheryDSEI and operation")+
  scale_fill_gradient(low = "lightyellow",
                      high = "steelblue")

ggsave("totalcatches_fisheryDSEI.png", width = 10, height = 5, dpi=200)




######################################### species distribution overlap plot #################################################
detach("package:plyr", unload=TRUE)
library(dplyr)
library(tidyr)
datos_na<-subset(datos,!is.na(datos$PROVINC))

fishery_na<-subset(fishery,!is.na(fishery$PROVINC))
if(layer=='species'){col_plot_names<-unique(datos$Species)
  datos_total_catches_species<-datos_na %>% group_by(Long, Lat, Species)  %>% summarise (tcatch = sum(tCatch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE))
  datos_total_catches_species_wide<-datos_total_catches_species %>% select(-c(tcatch)) %>% spread(Species,freq)
  
  }
if(layer=='fishery'){
  # col_plot_names<-c('LL','FLL','ELL','PS','HAND','TROL','GILL','LLCO','GIOF','DSEI')
  col_plot_names<-as.character(unique(fishery_na$fishery3))
  # col_plot_names<-c('BB','BBOF','BS','CN','DSEI','ELL','FLL','FN','GILL','GIOF','GL','HAND','HARP','HLOF','LG','LIFT','LL','LLCO','LLEX','PS','PSS','RIN','RNOF','RR')
  datos_total_catches_species<-fishery_na %>% group_by(Long, Lat, fishery3)  %>% summarise (catch = sum(catch,na.rm=TRUE)) %>% mutate(freq = catch / sum(catch), Totalcatch=sum(catch,na.rm=TRUE))
  datos_total_catches_species_wide<-datos_total_catches_species %>% select(-c(catch)) %>% spread(fishery3,freq)
}
head(datos_total_catches_species)

n <- nrow(datos_total_catches_species_wide)
datos_total_catches_species_wide$region <- factor(1:n)
datos_total_catches_species_wide<-as.data.frame(datos_total_catches_species_wide)
head(datos_total_catches_species_wide)
## need to replace the NA by zeros if i want to plot the pie.
datos_total_catches_species_wide[is.na(datos_total_catches_species_wide)] <- 0
# 
library(scatterpie)

  # ## unity 
  ggplot(data = world) + 
    geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
    geom_sf(color = "darkgrey", fill = "lightgrey")+xlab('')+ylab('')+
    coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) +
    geom_scatterpie(aes(x=Long, y=Lat, group=region, r=Totalcatch/Totalcatch*2),
    # geom_scatterpie(aes(x=Long, y=Lat, group=region, r=Totalcatch/100),
                    data=datos_total_catches_species_wide, cols=col_plot_names,
                    color="black", alpha=.8)+geom_scatterpie_legend(2, x=130, y=20)
    # geom_polygon(data = meow, aes(x = long, y = lat, group = group), colour = "black", fill = NA,size=1.2)+
    # geom_polygon(data = ppow, aes(x = long, y = lat, group = group), colour = "black", fill = NA,size=1)+
    # geom_scatterpie_legend(datos_total_catches_species_wide$Totalcatch/100, x=130, y=20) +

  ggsave(paste0('../figures/spatial/',layer,'_dist_map_',analysis,'_v2021.png'), width = 8, height = 5, dpi=200)
  
## would be nice to fill polygons by values of a variable:
## for later : https://stackoverflow.com/questions/19791210/r-ggplot2-merge-with-shapefile-and-csv-data-to-fill-polygons


# calculate the catch per spp over all provinces
# divide

# calculate total catch per province to see where most of the catches  occur
library(plyr)
TcatchP<-ddply(datos, .(PROVINC), summarize, Tcatch = sum(tCatch,na.rm=T))
summary(TcatchP)

TcatchP[with(TcatchP, order(-Tcatch)),]

# ## remove provinces with <5th percentile of catch (MEOW/PPOW combined analysis)
# datos1<-as.data.frame(datos)
# datos1[which(datos1$tCatch==0),'tCatch']<-NA # replace 0s with NAs to calculate quantile
# quantile(datos1$tCatch,probs=0.05,na.rm=T)
# rm_prov<-TcatchP[which(TcatchP$Tcatch<quantile(datos1$tCatch,probs=0.05,na.rm=T)),'PROVINC']
# 
write.csv(datos,file=paste0(home,'datos_',layer,'_',analysis,'_regroup_v2021.csv'),row.names=FALSE)




##### LOAD DATA FROM HERE!!!! #### <-----------------------------------------------------------------------------------------
datos<-read.csv(paste0(home,'datos_',layer,'_',analysis,'_regroup_v2021.csv'))

# if(length(rm_prov)==2){datos<-datos[-which(datos$PROVINCE==rm_prov[1] |datos$PROVINCE==rm_prov[2]),]}
# if(length(rm_prov)==3){datos<-datos[-which(datos$PROVINCE==rm_prov[1] |datos$PROVINCE==rm_prov[2]|datos$PROVINCE==rm_prov[3]),]}

if(layer=='species'){
datos$Species<-as.factor(datos$Species)
datos_sp<-datos
# if(analysis=='PPOW'){
#   # remove outliers from boxplots
# datos_sp<-datos[-which(datos$PROVINC=='Leeuwin Current' & datos$tCatch>150),]
# datos_sp<-datos_sp[-which(datos_sp$PROVINC=='Subtropical Convergence' & datos_sp$tCatch>10),]
# }

datos_sp$logCatch<-log(datos_sp$tCatch)
ggplot(datos_sp,aes(Species,logCatch))+geom_boxplot()+facet_wrap(~PROVINC,ncol=3,scales="free_y")+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave(paste0('../figures/spatial/Tcatch_boxplots_byprovince_for_',layer,'_',analysis,'_unforced_v2021.png'), width = 8, height = 5, dpi=200)

ggplot(datos_sp,aes(PROVINC,logCatch))+geom_boxplot()+facet_wrap(~Species,ncol=2,scales="free_y")+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave(paste0('../figures/spatial/Tcatch_boxplots_by',layer,'_',analysis,'_unforced_v2021.png'), width = 8, height = 5, dpi=200)
# ggplot(datos,aes(Species,tCatch))+geom_boxplot()+facet_wrap(~PROVINC,ncol=4,scales="free_y")+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

}
if(layer=='fishery'){
  datos_sp<-datos
  
  # datos$Fishery<-as.factor(datos$Fishery)
  datosf<-datos_sp[grepl(paste0(col_plot_names,collapse='|'),datos_sp$Fishery),]
  datosf<-datosf[-which(datosf$Fishery=='SLL' |datosf$Fishery=='LLEX'|datosf$Fishery=='PSS'),]
  # datosf<-datosf[-which(datosf$Fishery=='PSS' |datosf$Fishery=='LLEX'),]
  
  datosf$logCatch<-log(datosf$tCatch)
  ggplot(datosf,aes(Fishery,logCatch))+geom_boxplot()+facet_wrap(~PROVINCE2,ncol=3,scales="free_y")+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(paste0('../figures/spatial/Tcatch_boxplots_byprovince_for_',layer,'_',analysis,'_unforced_v4.png'), width = 8, height = 5, dpi=200)
  
  ggplot(datosf,aes(PROVINC,logCatch))+geom_boxplot()+facet_wrap(~Fishery,ncol=3,scales="free_y")+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(paste0('../figures/spatial/Tcatch_boxplots_by',layer,'_',analysis,'_unforced_v4.png'), width = 8, height = 5, dpi=200)
}
# datos$PROVINC<-as.factor(datos$PROVINC)

# explore data
# ggsave(paste0('../figures/spatial/Tcatch_boxplots_by',layer,'_',analysis,'_unforced_v3.png'), width = 8, height = 5, dpi=200)



#################################################### specifity analysis #######################################################
col_plot_names<-c('LL','FLL','ELL','PS','HAND','TROL','GILL','LLCO','GIOF','DSEI')

library(plyr)

step1sp<-ddply(datos, .(Species,PROVINC), function(x) sum(x$tCatch)/sum(datos$tCatch[which(datos$Species==unique(x$Species))]))
step2<-ddply(step1sp,.(Species), function(x) sum(x$V1))

if(layer=='species'){
###### spp specificity
step1sp$Species2<-as.factor(step1sp$Species)
levels(step1sp$Species2)<-c("TEMPERATE","TROPICAL","TROPICAL","SUBTROPICAL","TROPICAL")
step1sp$Species <- factor(step1sp$Species, levels=c("SKJ", "YFT","BET","SWO","ALB"))
names(step1sp)<-c('Species','PROVINCE','specificity','Species2')
# step1sp$Species <- factor(step1sp$Species, levels=c("SKJ", "YFT","BET","SWO","ALB"))

ggplot(step1sp,aes(Species,specificity))+geom_bar(stat="Identity", aes(fill=Species2))+
  facet_wrap(~PROVINCE,ncol=4,scales = 'free_y')+ theme(axis.text.x = element_text(angle = 90, hjust = 1))#+
ggplot(step1sp,aes(Species,specificity))+geom_bar(stat="Identity", aes(fill=Species2))+
  facet_wrap(~PROVINCE,ncol=4)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))#+
}

if(layer=='fishery'){
###### fisheries specificity
col_plot_names<-c('LL','FLL','ELL','PS','HAND','TROL','GILL','LLCO','GIOF','DSEI')
dg<-datos[which(eval(parse(text=paste0('datos$Fishery=="',col_plot_names,'"',collapse='|')))),]
step1f<-ddply(dg, .(Fishery,PROVINC), function(x) sum(x$tCatch)/sum(dg$tCatch[which(dg$Fishery==unique(x$Fishery))]))
step2<-ddply(dg,.(Fishery), function(x) sum(x$V1))
names(step1f)<-c('Species','PROVINCE','specificity')
# 
ggplot(step1f,aes(Species,specificity))+geom_bar(stat="Identity", aes(fill=Species))+facet_wrap(~PROVINCE,ncol=4,scales="free_y")+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Specificity")
ggplot(step1f,aes(Species,specificity))+geom_bar(stat="Identity", aes(fill=Species))+facet_wrap(~PROVINCE,ncol=4)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Specificity")
}



#the specificity Ai,j is the ratio of the mean abundance of species i in the geographical cells of group j (Ni,j) to the sum of the mean abundance of species i in all the groups (Ni).

if(layer=='fishery'){
  datos_spare<-datos
colnames(datos)<-c('Species1','Fleet','Species',colnames(datos)[4:11])
datos<-datos[grepl(paste0(col_plot_names,collapse='|'),datos$Species),]
datos<-datos[-which(datos$Species=='SLL'),]
datos<-datos[-which(datos$Species=='LLEX'),]
datos<-datos[-which(datos$Species=='PSS'),]
datos<-datos[-which(datos$Species=='TROLM'),]
}
# calculate the average catch by province and species
library(plyr)
AcatchPsp<-ddply(datos, .(PROVINC,Species), summarize, meancatchSp = mean(tCatch,na.rm=T))

summary(AcatchPsp)

levels(factor(AcatchPsp$PROVINC))

#calculate the sum of the mean catches of species by province  #old WRONG code!!!!
#SumCatchPro<-ddply(AcatchPsp, .(PROVINC), transform, SumCatchPro = sum(meancatchSp))
#summary(SumCatchPro)

#calculate the sum of the mean catches of all species in province (by SPECIES)
SumCatchSpecies1<-ddply(AcatchPsp, .(Species), summarize, SumCatchSp = sum(meancatchSp,na.rm=T))
summary(SumCatchSpecies1)

#calculate the sum of the mean catches of species (by Species) and add it as a column!
SumCatchSpecies<-ddply(AcatchPsp, .(Species), transform, SumCatchSp = sum(meancatchSp,na.rm=T))

#Calculate Specificity

Specificity_Sp_per_Prov<-ddply(SumCatchSpecies, .(PROVINC,Species), transform, specificity = meancatchSp/SumCatchSp)

summary(Specificity_Sp_per_Prov)

ordered<-Specificity_Sp_per_Prov[with(Specificity_Sp_per_Prov, order(PROVINC,-specificity)),]

ordered$specificity<-round(ordered$specificity,2)
ordered[,c("PROVINC","Species","specificity")]

#visualize data --------

#add groups by climate to improve visualization
levels(as.factor(Specificity_Sp_per_Prov$Species))
Specificity_Sp_per_Prov$Species2<-as.factor(Specificity_Sp_per_Prov$Species)
levels(Specificity_Sp_per_Prov$Species2)<-c("Temperate tuna","Tropical tuna","Tropical tuna","Subtropical billfish","Tropical tuna")

summary(Specificity_Sp_per_Prov)

library(ggplot2)
if(layer=='species'){
  Specificity_Sp_per_Prov$Species <- factor(Specificity_Sp_per_Prov$Species, levels=c("SKJ", "YFT","BET","SWO","ALB"))

ggplot(Specificity_Sp_per_Prov,aes(x=Species,y=specificity))+geom_bar(stat="Identity", aes(fill=Species2))+
    facet_wrap(~PROVINC,ncol=4,scales="free_y")+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Specificity") #+ scale_fill_manual(values = c("Blue", "Red", "Green"))
  ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_specificity_',analysis,'_free_v2021.png'), width = 8, height = 5, dpi=200)

  ggplot(Specificity_Sp_per_Prov,aes(Species,specificity))+geom_bar(stat="Identity", aes(fill=Species2))+
    facet_wrap(~PROVINC,ncol=4)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Specificity") #+ scale_fill_manual(values = c("Blue", "Red", "Green"))
    ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_specificity_',analysis,'_fixed_v2021.png'), width = 8, height = 5, dpi=200)

}
if(layer=='fishery'){
  ## free y
  ggplot(Specificity_Sp_per_Prov,aes(Species,specificity))+
    geom_bar(stat="Identity", aes(fill=Species))+facet_wrap(~PROVINC,ncol=3,scales="free_y")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Specificity")
  ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_specificity_',analysis,'_free.png'), width = 8, height = 5, dpi=200)
  ## fixed y
  ggplot(Specificity_Sp_per_Prov,aes(Species,specificity))+
    geom_bar(stat="Identity", aes(fill=Species))+facet_wrap(~PROVINC,ncol=3)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Specificity")
  ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_specificity_',analysis,'_fixed.png'), width = 8, height = 5, dpi=200)
  #+ scale_fill_manual(values = c("Blue", "Red", "Green"))
}

# fixed the scale of Y -better!
# ggplot(Specificity_Sp_per_Prov,aes(species,specificity))+geom_bar(stat="Identity", aes(fill=Species2))+facet_wrap(~PROVINCE,ncol=4)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Specificity") + scale_fill_manual(values = c("Green", "Red", "Blue"))


#polar plot
ggplot(Specificity_Sp_per_Prov,aes(Species,specificity))+geom_bar(stat="Identity", aes())+facet_wrap(~PROVINC,ncol=4)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ coord_polar()
## these are swamped by the large signals, so split into signal brackets

ddply(Specificity_Sp_per_Prov,.(PROVINC), function(x) sum(x$specificity))
ordered$specificity<-round(ordered$specificity,2)
ordered[,c("PROVINC","Species","specificity")]

ggplot(Specificity_Sp_per_Prov,aes(PROVINC,specificity))+geom_bar(stat="Identity", aes(fill=Species))+facet_grid(Species~.,scales="free_y")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "Specificity")+ theme(strip.text.y = element_text(size = 10, colour = "black", angle = 0))
ggplot(Specificity_Sp_per_Prov,aes(PROVINC,specificity))+geom_bar(stat="Identity", aes(fill=Species))+facet_grid(Species~.)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "Specificity")+ theme(strip.text.y = element_text(size = 10, colour = "black", angle = 0))

# Fidelity ---------------------------------
# The fidelity Bi,j is the ratio of the number of geographical cells in group j where species i is present (Si,j) to the total number of pixels in this group (Sj)
summary(datos)

# calculate total number of pixels (so number of unique lons and lats) in each province

#prepare the data. I need to collapse the column of Species. So I can calculate number of Long+Lat per province

TcatchProvLongLat<-ddply(datos, .(PROVINC, Long,Lat), function(x) Tcatch = sum(x$tCatch,na.rm=T))
summary(TcatchProvLongLat)

#now I can count the number of long+lat per province

Number_of_cells_Prov<-ddply(TcatchProvLongLat, .(PROVINC), summarize, count = length(PROVINC))

summary(Number_of_cells_Prov)

Number_of_cells_Prov[with(Number_of_cells_Prov, order(-count)),]
par(mar=c(12,5,1,1))
mp<-barplot(Number_of_cells_Prov[with(Number_of_cells_Prov, order(-count)),'count'],ylab='Number of pixels per province')
axis(side = 1, at = mp, labels = Number_of_cells_Prov[with(Number_of_cells_Prov, order(-count)),'PROVINC'], las=2,tcl = -0.2)

# calculate total number of pixels (so number of unique lons and lats) of each species (where the species is present) in each province
summary(datos)
#Remove those cells where the total catch is ZERO. We assume the species is not present there.

dim(datos)
range(datos$tCatch,na.rm=T) # there are not cells where the catch is zero. So we keep it all.

# remove cells where catch is zero or NA - only cells with catch data remain
datos_na<-subset(datos,!is.na(datos$tCatch))
datos_na<-subset(datos_na,datos$tCatch>0)

#Calculate the number of cells per province and species where species is present
if(layer=='species'){
Number_of_cells_Prov_Sp<-ddply(datos_na, .(PROVINC,Species), function(x) count = dim(unique(x[,c('Lat','Long')]))[1])
}
if(layer=='fishery'){
  Number_of_cells_Prov_Sp<-ddply(datos_na, .(PROVINC,Species), function(x) count = dim(unique(x[,c('Lat','Long')]))[1])
}
head(Number_of_cells_Prov_Sp)	

Number_of_cells_Prov_Sp[with(Number_of_cells_Prov_Sp, order(PROVINC,-V1)),]


# I need to add both data.frames into one

head(Number_of_cells_Prov)
head(Number_of_cells_Prov_Sp)
dim(Number_of_cells_Prov)
dim(Number_of_cells_Prov_Sp)


Fidelity_Sp_per_Prov<-merge(Number_of_cells_Prov,Number_of_cells_Prov_Sp,by.x="PROVINC",by.y="PROVINC")
if(layer=='species'){names(Fidelity_Sp_per_Prov)<-c("PROVINCE","Number_of_cells_Prov","Species","Number_of_cells_Prov_Sp")}
if(layer=='fishery'){names(Fidelity_Sp_per_Prov)<-c("PROVINCE","Number_of_cells_Prov","Species","Number_of_cells_Prov_Sp")}

# plot(xx)

head(Fidelity_Sp_per_Prov)


Fidelity_Sp_per_Prov$fidelity<-Fidelity_Sp_per_Prov$Number_of_cells_Prov_Sp/Fidelity_Sp_per_Prov$Number_of_cells_Prov
# plot(Fidelity_Sp_per_Prov$Number_of_cells_Prov,Fidelity_Sp_per_Prov$fidelity,type='p',xlab='Number of pixels per province',
#      ylab='Fidelity')

levels(as.factor(Fidelity_Sp_per_Prov$Species))
Fidelity_Sp_per_Prov$Species2<-as.factor(Fidelity_Sp_per_Prov$Species)
# levels(Fidelity_Sp_per_Prov$Species2)<-c(unique(Fidelity_Sp_per_Prov$Species))
levels(Fidelity_Sp_per_Prov$Species2)<-c("Temperate","Tropical","Tropical","Subtropical","Tropical")

#add groups by climate to improve visualization
if(layer=='species'){
  Fidelity_Sp_per_Prov$Species <- factor(Fidelity_Sp_per_Prov$Species, levels=c("SKJ", "YFT","BET","SWO","ALB"))
  
  ggplot(Fidelity_Sp_per_Prov,aes(Species,fidelity))+geom_bar(stat="Identity", aes(fill=Species2))+
    facet_wrap(~PROVINCE,ncol=3,scales='free_y')+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
    labs(title = "Fidelity")#+ scale_fill_manual(values = c("Blue", "Red", "Green"))
  ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_fidelity_',analysis,'_free_y.png'), width = 8, height = 5, dpi=200)
  
ggplot(Fidelity_Sp_per_Prov,aes(Species,fidelity))+geom_bar(stat="Identity", aes(fill=Species2))+
    facet_wrap(~PROVINCE,ncol=3)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
    labs(title = "Fidelity")#+ scale_fill_manual(values = c("Blue", "Red", "Green"))
ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_fidelity_',analysis,'_fixed_y.png'), width = 8, height = 5, dpi=200)

xx<-lm(Fidelity_Sp_per_Prov$Number_of_cells_Prov_Sp ~ Fidelity_Sp_per_Prov$Number_of_cells_Prov)
summary(xx)
plot(Fidelity_Sp_per_Prov$Number_of_cells_Prov,Fidelity_Sp_per_Prov$Number_of_cells_Prov_Sp,xlab='Number of cells per Province',
     ylab='Number of cells per species per province')
abline(xx,col='red')

xx<-lm(Fidelity_Sp_per_Prov$fidelity ~ Fidelity_Sp_per_Prov$Number_of_cells_Prov)
summary(xx)

plot(Fidelity_Sp_per_Prov$Number_of_cells_Prov,Fidelity_Sp_per_Prov$fidelity,xlab='Number of cells per Province',
     ylab='Fidelity')
abline(xx,col='red')
}

if(layer=='fishery'){
  ggplot(Fidelity_Sp_per_Prov,aes(Species,fidelity))+geom_bar(stat="Identity", aes(fill=Species2))+
    facet_wrap(~PROVINCE,ncol=3)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
    labs(title = "Fidelity")#+ scale_fill_manual(values = c("Blue", "Red", "Green"))
  ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_fidelity_',analysis,'_free_y.png'), width = 8, height = 5, dpi=200)
  
  ggplot(Fidelity_Sp_per_Prov,aes(Species,fidelity))+geom_bar(stat="Identity", aes(fill=Species2))+
    facet_wrap(~PROVINCE,ncol=3,scales='free_y')+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
    labs(title = "Fidelity")#+ scale_fill_manual(values = c("Blue", "Red", "Green"))
  ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_fidelity_',analysis,'_fixed_y.png'), width = 8, height = 5, dpi=200)
  
  
  }

#Species indicator value=specificity *fidelity *100
#The indicator value (Vi,j) is calculated by multiplying the specificity and fidelity indices, because these two quantities represent independent information.

head(Specificity_Sp_per_Prov)
head(Fidelity_Sp_per_Prov)

 # as.factor(step1sp$Specie)

# risky and messy to do this..but given time...i am doing it.... we need to make it more clean and save
# FinalDF<-left_join(Specificity_Sp_per_Prov[c('PROVINC','Species','meancatchSp','SumCatchSp','specificity')],Fidelity_Sp_per_Prov,by=c('PROVINC','Species'))

if(layer=='fishery'){FinalDF<-left_join(step1f[c('PROVINCE','Species','specificity')],Fidelity_Sp_per_Prov,by=c('PROVINCE','Species'))}
if(layer=='species'){FinalDF<-left_join(step1sp[c('PROVINCE','Species','specificity')],Fidelity_Sp_per_Prov,by=c('PROVINCE','Species'))}


# if(layer=='species'){
#   names(step1sp)<-c('Species','PROVINC','specificity','Species2')
#   FinalDF<-left_join(step1sp[c('PROVINC','Species','specificity')],Fidelity_Sp_per_Prov,by=c('PROVINC','Species'))}
# 
# if(layer=='fishery'){
#   names(step1f)<-c('Fishery','PROVINC','specificity')
#   FinalDF<-left_join(step1f[c('PROVINC','Fishery','specificity')],Fidelity_Sp_per_Prov,by=c('PROVINC','Fishery'))}
# FinalDF<-cbind(Specificity_Sp_per_Prov,Fidelity_Sp_per_Prov)


FinalDF$SpIndicator<-FinalDF$specificity*FinalDF$fidelity*100
FinalDF<-FinalDF[-which(FinalDF$PROVINC=='Antarctic'|FinalDF$PROVINC=='Antarctic Polar Front'),]


ggplot(FinalDF,aes(Number_of_cells_Prov,fidelity))+
  geom_point()+
  facet_wrap(~Species)

ggplot(FinalDF,aes(Number_of_cells_Prov,specificity))+
  geom_point()+
  facet_wrap(~Species)


summary(FinalDF)

# if(layer=='species'){
  FinalDF$Species <- factor(FinalDF$Species, levels=c("SKJ", "YFT","BET","SWO","ALB"))
  ggplot(FinalDF,aes(Species,SpIndicator))+geom_bar(stat="Identity", aes(fill=Species))+
    facet_wrap(~PROVINCE,ncol=4,scales = 'free_y')+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
    labs(title = "Species Indicator Value")#+ scale_fill_manual(values = c("Blue", "Red", "Green"))
  ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_indicator_',analysis,'_facetwrap_by',layer,'_free_y_v2021.png'), width = 8, height = 5, dpi=200)
  
  ggplot(FinalDF,aes(Species,SpIndicator))+geom_bar(stat="Identity", aes(fill=Species2))+
    facet_wrap(~PROVINCE,ncol=4)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
    labs(title = "Species Indicator Value")#+ scale_fill_manual(values = c("Blue", "Red", "Green"))
  ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_indicator_',analysis,'_facetwrap_by',layer,'_fixed_y_v2021.png'), width = 8, height = 5, dpi=200)
  # ggplot(FinalDF,aes(PROVINC,SpIndicator))+geom_bar(stat="Identity", aes(fill=Species))+facet_wrap(Species~.,nrow=5,scales = 'free_y')+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "Species Indicator Value")#+ scale_fill_manual(values = c("Blue", "Red", "Green"))
# }
if(layer=='fishery'){
  ggplot(FinalDF,aes(Species,SpIndicator))+
    geom_bar(stat="Identity", aes(fill=Species))+facet_wrap(~PROVINCE,ncol=4,scales = 'free_y')+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "Species Indicator Value")#+ scale_fill_manual(values = c("Blue", "Red", "Green"))
  # ggplot(FinalDF,aes(PROVINC,SpIndicator))+geom_bar(stat="Identity", aes(fill=Species))+facet_wrap(~Species,nrow=10,scales = 'free_y')+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "Species Indicator Value")#+ scale_fill_manual(values = c("Blue", "Red", "Green"))
  # ggplot(FinalDF,aes(PROVINC,specificity))+geom_bar(stat="Identity", aes(fill=Species))+facet_grid(Species~.,scales="free_y")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "SF Indicator")+ theme(strip.text.y = element_text(size = 10, colour = "black", angle = 0))
  # ggplot(FinalDF,aes(PROVINC,specificity))+   geom_bar(stat="Identity", aes(fill=Species))+ facet_grid(Species~.,scales="free_y")+     theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "SF Indicator")+     theme(strip.text.y = element_text(size = 10, colour = "black", angle = 0))
  ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_indicator_',analysis,'_facetwrap_by',layer,'_free_y.png'), width = 8, height = 5, dpi=200)

  ggplot(FinalDF,aes(Species,SpIndicator))+
    geom_bar(stat="Identity", aes(fill=Species))+facet_wrap(~PROVINCE,ncol=4)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title = "Species Indicator Value")#+ scale_fill_manual(values = c("Blue", "Red", "Green"))
  ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_indicator_',analysis,'_facetwrap_by',layer,'_fixed.png'), width = 8, height = 5, dpi=200)

}
# # ggsave(paste0('../figures/presentation/ecologicaldata/',layer,'_indicator_',analysis,'_facetwrap_by',layer,'.png'), width = 8, height = 5, dpi=200)

## oops - overwrote v3 with v4 (v4==fewer gears)
write.csv(FinalDF,file=paste0(home,layer,'_final_indicator_df_',analysis,'__v2021.csv'),row.names=FALSE)

