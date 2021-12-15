## visualise_catch.R
## AE Nieblas
## DESCRIPTION:  based off MJJJ's code IOTC_read_catch_data_july2019.R and adapted for the updated IOTC raised catch for SKJ, ALB, YFT, BET, and SWO.

rm(list=ls())
source('~/Documents/PERSONAL/COOOL/projects/Ecoregions_2021/code/iotc_regrid.R')
home<-'/home/ae/Documents/PERSONAL/COOOL/projects/Ecoregions_2021/' ## ae remove data 
setwd(home)

type<-'species'


library(ggplot2)
theme_set(theme_bw())
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgdal)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

iotc_shp=paste0(home,'data/indian_ocean/iotc.shp')
iotc <- readOGR(dsn = iotc_shp, stringsAsFactors = F)

iotc_5deg_shp=paste0(home,'data/IOTC-2018-TCAC04-DATA04_-_5x5_grids_shapefile/IOTC_TCAC_GRID_5x5.shp')

## regrid all 1 degrees to 5 degrees
d1<- iotc_1_to_5_regrid(iotc_5deg_shp=iotc_5deg_shp,
                        iotc_df=paste0(home,'data/EAFM_data_v2.csv'),data=2017)

## regrid main spp/fisheries/fleets
# raised_catch<- iotc_1_to_5_regrid(iotc_5deg_shp=iotc_5deg_shp,iotc_df=paste0(home,'EAFM_data_v1.csv'),data=2017) ## remove???

## data filtering
# outliers (points >95th percentile of total catch)
raised_catch<-d1
raised_catch<-raised_catch[-which(raised_catch$RAISED_CATCHES_MT>quantile(raised_catch$RAISED_CATCHES_MT,probs=0.95,na.rm=T)),]

# tropical species found <45
levels(as.factor(raised_catch$SPECIES))
raised_catch$Species2<-as.factor(raised_catch$SPECIES)
levels(raised_catch$Species2)<-c("Temperate","Tropical","Tropical","Subtropical","Tropical")

raised_catch[which(raised_catch$Species2=='Tropical' & raised_catch$lat<(-45)),'RAISED_CATCHES_MT']<-NA


library(plyr)
## sum catch within each grid cell
d2<-ddply(raised_catch,.(SPECIES, YEAR, GEAR,OPERATION_TYPE,FLEET_CODE_HASHED,size,quadrant,lat,lon), function(x) sum(x$RAISED_CATCHES_MT,na.rm=T))
# m2<-ddply(d1,.(SPECIES, MONTH, GEAR,OPERATION_TYPE,FLEET_CODE_HASHED,size,quadrant,lat,lon), function(x) sum(x$RAISED_CATCHES_MT,na.rm=T))

# remove spurious catch or catch outside the IOTC CA

## adapt to rtunatlas df format
d2$geographic_identifier<-paste0(d2$size,d2$quadrant,sprintf("%02d",(abs(d2$lat)-2.5)),sprintf("%03d",(d2$lon-2.5)))
d2$unit<-'MT'
d2$value<-d2$V1

## remove package plyr
detach("package:plyr", unload=TRUE) 

library(dplyr)

##### Examine the catch species composition/location of the different fleets/type of fisheries

head(d2) 
start_year<-2005
end_year<-2019


#before we start the analysis, get an average across all the years.
# datos3 <- aggregate(d2$value, list(d2$SPECIES, d2$lon, datos2$lat,datos2$FLEET_CODE,datos2$FISHERY_CODE), median)

datos3<-d2 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(SPECIES, lon,lat, FLEET_CODE_HASHED,GEAR,OPERATION_TYPE,geographic_identifier) %>% summarize (catch=median(value,na.rm=TRUE))
obs<-d2 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(lon,lat) %>% tally()

if(type=='species'){obs3<-d2 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(SPECIES, lon,lat) %>% tally()}
if(type=='fishery'){obs3<-d2 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(GEAR, lon,lat) %>% tally()}
head(datos3)

names(datos3)<-c("species","lon", "lat","fleet","fishery", "operation",'geographic_identifier',"catch")
datos3$unit<-'MT'


##TOTALCATCHES
Tcatcheslonlat<-datos3 %>% group_by(lon,lat) %>% summarize (catch=sum(catch,na.rm=TRUE))

library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(RColorBrewer)
world <- ne_countries(scale = "medium", returnclass = "sf")
# class(world)

## PLOT BY POINTS OR BY TILES
cols<-brewer.pal(8,'YlGnBu')
# Tcatcheslonlat[which(Tcatcheslonlat$catch==0),'catch']<-NA
ggplot(data = world) + 
  # geom_point(data=Tcatcheslonlat, aes(lon,lat, size = catch)) +
  geom_tile(data=Tcatcheslonlat, aes(lon,lat, fill = catch)) +
  scale_fill_gradientn(colours=cols,limits = c(0,max(Tcatcheslonlat$catch)))+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ 
  ggtitle(" Total raised catches of open ocean species")
ggsave("figures/presentation/ecologicaldata/totalcatches_mainspp_v2021_tiles.png", width = 8, height = 5, dpi=200)

library(RColorBrewer)
cols<-brewer.pal(8,'YlOrRd')
ggplot(data = world) +
  geom_tile(data = obs, aes(x=lon,y=lat,fill=n))+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  # facet_wrap(~SPECIES,ncol=2)+
  scale_fill_gradientn(colours=cols)+
  coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
  ggtitle(paste0("Number of observations  "))
ggsave("figures/presentation/ecologicaldata/numberobservations_all.png", width = 8, height = 5, dpi=200)

##TOTALCATCHES by species - tiles v points
Tcatcheslonlatspp<-datos3 %>% group_by(lon,lat,species) %>% summarize (catch=sum(catch,na.rm=TRUE))

cols<-brewer.pal(8,'YlGnBu')
# Tcatcheslonlatspp[which(Tcatcheslonlatspp$catch==0),'catch']<-NA
ggplot(data = world) + 
  geom_tile(data = Tcatcheslonlatspp, aes(x=lon,y=lat,fill=catch))+
  scale_fill_gradientn(colours=cols)+
    geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ 
  # geom_point(data=Tcatcheslonlatspp, aes(lon,lat, size = catch)) +
  facet_wrap(~species,ncol=2)+
  ggtitle(" Total raised catches of open ocean species")
ggsave("figures/presentation/ecologicaldata/totalcatches_facetwrap_mainspp_v2021_tiles.png", width = 8, height = 5, dpi=200)


  library(RColorBrewer)
  cols<-brewer.pal(8,'YlOrRd')
  ggplot(data = world) +
    geom_tile(data = obs3, aes(x=lon,y=lat,fill=n))+
    geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
    geom_sf(color = "darkgrey", fill = "lightgrey")+
    facet_wrap(~SPECIES,ncol=2)+
    scale_fill_gradientn(colours=cols)+
    coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
    ggtitle(paste0("Number of observations  "))
  ggsave("figures/presentation/ecologicaldata/numberobservations_facetwrap_mainspp.png", width = 8, height = 5, dpi=200)
  # ggsave("figures/presentation/ecologicaldata/numberobservations_mainspp.png", width = 8, height = 5, dpi=200)
  


##TOTALCATCHES by gear
  summary(unique(datos3$fishery))
  # total catches by main fisheries
  Tcatchbyfishery<-datos3 %>% group_by(fishery) %>% summarize (catch=sum(catch,na.rm=TRUE))%>% arrange(-catch)
  Tcatchbyfishery$fishery<-as.factor(Tcatchbyfishery$fishery)

  g<-ggplot(Tcatchbyfishery[1:10,],aes(fishery))
  g+geom_bar(aes(weight=catch/10000))+ylab('catch (MT/10,000)')

  mp<-barplot(Tcatchbyfishery$catch[1:10]/10000,ylab='catch (MT/10,000)')
  axis(side = 1, at = mp, labels = Tcatchbyfishery$fishery[1:10], tcl = -0.2)

  
  
Tcatcheslonlatfish<-datos3 %>% group_by(lon,lat,fishery) %>% summarize (catch=sum(catch,na.rm=TRUE))

library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")
# class(world)

Tcatcheslonlatfish[which(Tcatcheslonlatfish$catch==0),'catch']<-NA
# ggplot(data = world) + 
#   geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
#   geom_sf(color = "darkgrey", fill = "lightgrey")+
#   coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ 
#   geom_point(data=Tcatcheslonlatfish, aes(lon,lat, size = catch)) +
#   facet_wrap(~fishery,ncol=2)+
#   ggtitle(" Total raised catches of open ocean species")
# ggsave("../figures/totalcatches_mainfishery.png", width = 8, height = 5, dpi=200)

  LL<-c('ELL','FLL','LL','LLCO')
  PS<-c('PS')
  GILL<-c('GILL','GIOF')
  OT<-c('HAND','TROL','DSEI')
  
  gear='OT'
  if(gear=='LL'){    datosg<-Tcatcheslonlatfish[grepl(paste0(LL,collapse='|'),Tcatcheslonlatfish$fishery),]
  datosg<-datosg[-which(datosg$fishery=='LLEX'),]
  datosg<-datosg[-which(datosg$fishery=='SLL'),]
  datosg<-datosg[-which(datosg$fishery=='GILL'),]
  }
  if(gear=='PS'){datosg<-Tcatcheslonlatfish[which(Tcatcheslonlatfish$fishery=='PS'),] }
  if(gear=='GILL'){datosg<-Tcatcheslonlatfish[grepl(paste0(GILL,collapse='|'),Tcatcheslonlatfish$fishery),] }
  if(gear=='OT'){
    # datosg<-Tcatcheslonlatfish[grepl(paste0(OT,collapse='|'),Tcatcheslonlatfish$fishery),] 
  datosg<- Tcatcheslonlatfish[which(Tcatcheslonlatfish$fishery=='TROL'|
                             Tcatcheslonlatfish$fishery=='DSEI'|
                               Tcatcheslonlatfish$fishery=='HAND'),]
  }
  
  library(RColorBrewer)
  cols<-brewer.pal(8,'YlGnBu')
  ggplot(data = world) +
    geom_tile(data=datosg, aes(lon,lat, fill = catch)) +
    geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
    # geom_polygon(data = ppow, aes(x = long, y = lat, group = group), colour = "black", fill = NA)+
    geom_sf(color = "darkgrey", fill = "lightgrey")+
    # geom_point(data=datosg, aes(lon,lat, size = catch)) +
    scale_fill_gradientn(colours=cols)+
    facet_wrap(~fishery,ncol=3)+
    coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+
    # facet_wrap(~fishery,ncol=3)+
    ggtitle(paste0("Sum of raised catch (MT) by ",gear))
# ggsave(paste0("../figures/fisheries/totalcatches_facetwrap_",gear,"_v2.png"), width = 8, height = 5, dpi=200)
ggsave(paste0("../figure/presentation/ecologicaldata/totalcatches_facetwrap_",gear,"_sum_v2021.png"), width = 8, height = 5, dpi=200)

# ##<----------------------------- AE check the gear groupings
# ## check ELL distribution - surface nighttime setting around Reunion
# ## check FLL distribution - targets ALB and SWO (deep?)
# unique(datos3$fishery)
# ## propose LL industrial
# LL<-c('ELL','FLL','LL')
# ## and LL artisanal
# LLCO<-c('LLCO')
# 
# # propose PS industrial
# PS<-c('PS')
# # and PS artisanal (to check distribution, otherwise ignore)
# PSCO<-c('PSS')
# 
# GILL<-c('GILL','GIOF')
# OT<-c('HAND','TROL','BB')
# 
# if(type=='fishery'){
#   gear='PS'
#   if(gear=='LL'){    datosg<-obs3[grepl(paste0(LL,collapse='|'),obs3$GEAR),]
#   datosg<-datosg[-which(datosg$GEAR=='LLEX'),]
#   datosg<-datosg[-which(datosg$GEAR=='SLL'),]
#   datosg<-datosg[-which(datosg$GEAR=='GILL'),]
#   }
#   if(gear=='PS'){datosg<-obs3[which(obs3$GEAR=='PS'),] }
#   if(gear=='GILL'){datosg<-obs3[grepl(paste0(GILL,collapse='|'),obs3$GEAR),] }
#   if(gear=='OT'){datosg<-obs3[grepl(paste0(OT,collapse='|'),obs3$GEAR),] }
#   
#   ## -------------------------------- change lines to see different gears one by one
#   library(RColorBrewer)
#   cols<-brewer.pal(8,'YlOrRd')
#   # datos1<-as.data.frame(datos)
#   ggplot(data = world) +
#     # geom_tile(data = datosg, aes(x=lon,y=lat,fill=n))+
#     geom_tile(data = obs3, aes(x=lon,y=lat,fill=n))+
#     geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
#     geom_sf(color = "darkgrey", fill = "lightgrey")+
#     # geom_point(data=datos1, aes(Long,Lat, size = tCatch)) +
#     facet_wrap(~GEAR,ncol=5)+
#     # facet_wrap(~SPECIES,ncol=2)+
#     scale_fill_gradientn(colours=cols)+
#     coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
#     ggtitle(paste0("Number of observations  "))
#   # ggsave("../figures/presentation/ecologicaldata/numberobservations_facetwrap_mainspp.png", width = 8, height = 5, dpi=200)
#   ggsave(paste0("../figures/presentation/ecologicaldata/numberobservations_total_facetwrap",gear,"_grouped_v2021.png"), width = 8, height = 5, dpi=200)
# }
# 



#### TOTAL CATCHES BY SPECIES - MJJJ
library(dplyr)
library(tidyr)
datos_total_catches_species<-datos3 %>% group_by(lon, lat, species)  %>% summarise (tcatch = sum(catch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
head(datos_total_catches_species)

datos_total_catches_species_wide<-datos_total_catches_species %>% select(-c(tcatch)) %>% spread(species,freq)
n <- nrow(datos_total_catches_species_wide)
datos_total_catches_species_wide$region <- factor(1:n)
datos_total_catches_species_wide<-as.data.frame(datos_total_catches_species_wide)
head(datos_total_catches_species_wide)
## need to replace the NA by zeros if i want to plot the pie.
datos_total_catches_species_wide[is.na(datos_total_catches_species_wide)] <- 0

library(scatterpie)
xx<-datos_total_catches_species_wide[,c('lon','lat','ALB','BET','SKJ','SWO','YFT')]
library(reshape2)
ddspp<-melt(xx,id=c('lon','lat'))

# ggplot(data = world) +
#   geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
#   geom_sf(color = "darkgrey", fill = "lightgrey")+
#   coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+
#   geom_point(data=ddspp, aes(lon,lat, size = value)) +
#   facet_wrap(~variable)+
#   ggtitle(" Total raised catches of open ocean species")
# ggsave("../figures/presentation/ecologicaldata/totalcatches_facetwrap_mainspp.png", width = 8, height = 5, dpi=200)


## proportional : /10000; unity : /Totalcatch*2
# ggplot(data = world) +
#   geom_sf(color = "black", fill = "lightgrey")+
#   geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
#   coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) +
#   geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/10000), data=datos_total_catches_species_wide,
#                   cols=c("ALB","BET","SKJ","SWO","YFT"),  color="black", alpha=.8)+
#   geom_scatterpie_legend(datos_total_catches_species_wide$Totalcatch/10000, x=130, y=20) +
#   ggtitle("Total catches by species")


## unity
ggplot(data = world) + 
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + 
  geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), 
                  data=datos_total_catches_species_wide, cols=c("BET","SKJ","SWO","YFT",'ALB'),  
                  color="black", alpha=.8)+geom_scatterpie_legend(2, x=130, y=20) +
  ggtitle("Spatial distribution of open ocean species")
ggsave("figures/presentation/ecologicaldata/totalcatches_pie_unity_mainspp_v2021.png", width = 8, height = 5, dpi=200)




#### TOTAL CATCHES BY FISHERY - MJJJ
library(dplyr)
library(tidyr)
library(scatterpie)
datos_total_catches_fishery<-datos3 %>% group_by(lon, lat, fishery)  %>% summarise (tcatch = sum(catch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
head(datos_total_catches_fishery)

## SELECT the TOP FISHERIES: 
FISH<-c('ELL','FLL','LLCO','GILL','GIOF','HAND','DSEI')
OT<-c('TROL')
datos_total_catches_fishery1<-datos_total_catches_fishery[grep(paste0(FISH,collapse='|'),datos_total_catches_fishery$fishery),]
datos_total_catches_fishery2<-datos_total_catches_fishery[which(datos_total_catches_fishery$fishery=='TROL'|
                                                                  datos_total_catches_fishery$fishery=='LL'|
                                                                  datos_total_catches_fishery$fishery=='PS'),]
datos_total_catches_fishery<-rbind(datos_total_catches_fishery1,datos_total_catches_fishery2)

datos_total_catches_fishery_wide<-datos_total_catches_fishery %>% select(-c(tcatch)) %>% spread(fishery,freq)
n <- nrow(datos_total_catches_fishery_wide)
datos_total_catches_fishery_wide$region <- factor(1:n)
datos_total_catches_fishery_wide<-as.data.frame(datos_total_catches_fishery_wide)
head(datos_total_catches_fishery_wide)
## need to replace the NA by zeros if i want to plot the pie.
datos_total_catches_fishery_wide[is.na(datos_total_catches_fishery_wide)] <- 0

## proportional : /10000; unity : /Totalcatch*2
# ggplot(data = world) + 
#   geom_sf(color = "black", fill = "lightgrey")+
#   geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
#   coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + 
#   geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/10000), data=datos_total_catches_fishery_wide, 
#                   cols=c("ALB","BET","SKJ","SWO","YFT"),  color="black", alpha=.8)+
#   geom_scatterpie_legend(datos_total_catches_fishery_wide$Totalcatch/10000, x=130, y=20) +
#   ggtitle("Total catches by fishery")

# col<-c('LL','ELL','FLL','GIOF','PS','GILL','HAND','TROL','LLCO','BB')
col<-c('DSEI', 'ELL','FLL','GILL','GIOF','HAND','LL', 'LLCO', 'PS', 'TROL')


## unity
ggplot(data = world) + 
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + 
  geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), 
                  data=datos_total_catches_fishery_wide, cols=col,  
                  color="black", alpha=.8)+geom_scatterpie_legend(2, x=130, y=20) 
ggsave("figures/presentation/ecologicaldata/piechart_unity_fisheries_v2021.png", width = 8, height = 5, dpi=200)

# +ggtitle("TROL - all coastal fleets")
## facet wrap plots


#### TOTAL CATCHES BY FISHERY AND SPECIES 
library(dplyr)
library(tidyr)
############################# LONGLINES #########################
gear='LL'
datos4<-datos3[which(datos3$fishery=='LL'|datos3$fishery=='ELL'|datos3$fishery=='FLL'|datos3$fishery=='LLCO'),]
datos_total_catches_fspp<-datos4 %>% group_by(lon, lat, species)  %>% summarise (tcatch = sum(catch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
head(datos_total_catches_fspp)

# LL<-c('ELL','FLL','LL','LG','LLCO')
# PS<-c('PS','PSS','RIN')
# GILL<-c('GILL','GIOF','GL')
# OT<-c('HAND','TROL','BB')

# datos_total_catches_fspp$fishery<-NULL
datos_total_catches_fspp_wide<-datos_total_catches_fspp %>% select(-c(tcatch)) %>% spread(species,freq)
n <- nrow(datos_total_catches_fspp_wide)
datos_total_catches_fspp_wide$region <- factor(1:n)
datos_total_catches_fspp_wide<-as.data.frame(datos_total_catches_fspp_wide)
head(datos_total_catches_fspp_wide)
## need to replace the NA by zeros if i want to plot the pie.
datos_total_catches_fspp_wide[is.na(datos_total_catches_fspp_wide)] <- 0

library(scatterpie)

# ## proportional : /10000; unity : /Totalcatch*2
# ggplot(data = world) +
#   geom_sf(color = "black", fill = "lightgrey")+
#   geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
#   coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) +
#   geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/10000), data=datos_total_catches_fspp_wide,
#                   cols=c("ALB","BET","SKJ","SWO","YFT"),  color="black", alpha=.8)+
#   geom_scatterpie_legend(datos_total_catches_fspp_wide$Totalcatch/10000, x=130, y=20) +
#   ggtitle("LL catch")

# <-c('LL','ELL','FLL','GIOF','PS','GILL','HAND','TROL','PSS','LLCO','BB','GL','LG','RIN')
## unity
ggplot(data = world) + 
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + 
  geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), 
                  data=datos_total_catches_fspp_wide, cols=c('ALB','BET','SKJ','SWO','YFT'),  
                  color="black", alpha=.8)+geom_scatterpie_legend(2, x=130, y=20) +
  ggtitle("LL - all industrial and coastal fleets")


######################## facet wrap fishery by spp -------------------- AE doesn't work for 2021....
LL<-c('ELL','FLL','LL','LLCO')
PS<-c('PS')
GILL<-c('GILL','GIOF')
OT<-c('HAND','TROL','BB')

gear='LL'
if(gear=='LL'){    datosg<-Tcatcheslonlatfish[grepl(paste0(LL,collapse='|'),Tcatcheslonlatfish$fishery),]
datosg<-datosg[-which(datosg$fishery=='LLEX'),]
datosg<-datosg[-which(datosg$fishery=='SLL'),]
datosg<-datosg[-which(datosg$fishery=='GILL'),]
}
if(gear=='PS'){datosg<-Tcatcheslonlatfish[which(Tcatcheslonlatfish$fishery=='PS'),] }
if(gear=='GILL'){datosg<-Tcatcheslonlatfish[grepl(paste0(GILL,collapse='|'),Tcatcheslonlatfish$fishery),] }
if(gear=='OT'){datosg<-Tcatcheslonlatfish[grepl(paste0(OT,collapse='|'),Tcatcheslonlatfish$fishery),] }

gear='LL'
datosg<-datos_total_catches_fspp_wide[,c('lon','lat','Totalcatch',gear)] 
library(reshape)
# for(g in 1:length(LL)){
  # xx<-melt(datos_total_catches_fishery_wide,by=c('lon','lat','region'))
  ggplot(data = world) + 
    geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
    geom_sf(color = "darkgrey", fill = "lightgrey")+
    coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + 
    facet_wrap(~fishery,ncol=2)+
    geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), 
                    data=datosg,cols=c('ALB','BET','SKJ','SWO','YFT'),  
                    color="black", alpha=.8)+geom_scatterpie_legend(2, x=130, y=20) 
  ggsave(paste0("../figures/presentation/ecologicaldata/piechart_unity_facetwrap_",gear,".png"), width = 8, height = 5, dpi=200)
# }



################### PS large scale 
datosPS_largePS_plot<-datos3 %>% filter(fishery=="PS" ) %>% group_by(lon, lat) %>% summarize (Tcatch=sum(catch,na.rm=TRUE))
summary(datosPS_largePS_plot)

## AE - checkdata on land africa
ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ geom_point(data=datosPS_largePS_plot, aes(lon,lat, size = Tcatch)) +ggtitle(" PS high seas  - all fleets")

library(dplyr)
datosPS_largePS<-datos3 %>% filter(fishery=="PS")
datosPS_largePS_species<-datosPS_largePS %>% group_by(lon, lat, species)  %>% summarise (tcatch = sum(catch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
head(datosPS_largePS_species)

datosPS_largePS_species_wide<-datosPS_largePS_species %>% select(-c(tcatch)) %>% spread(species,freq)
n <- nrow(datosPS_largePS_species_wide)
datosPS_largePS_species_wide$region <- factor(1:n)
datosPS_largePS_species_wide<-as.data.frame(datosPS_largePS_species_wide)
head(datosPS_largePS_species_wide)
## need to replace the NA by ceros if i want to plot the pie.
datosPS_largePS_species_wide[is.na(datosPS_largePS_species_wide)] <- 0

library(scatterpie)

# ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/10000), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT"),  color="black", alpha=.8)+geom_scatterpie_legend(datosPS_largePS_species_wide$Totalcatch/10000, x=130, y=20) +ggtitle("Large PS - PSLS and PSFS")

#unity
ggplot(data = world) + geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT"),  color="black", alpha=.8)+geom_scatterpie_legend(2, x=130, y=20) +ggtitle("Industrial PS ")


################### PS coastal scale 
datosPS_largePS_plot<-datos3 %>% filter(fishery=="PSS"|fishery=="RIN" ) %>% group_by(lon, lat) %>% summarize (Tcatch=sum(catch,na.rm=TRUE))
summary(datosPS_largePS_plot)


ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ geom_point(data=datosPS_largePS_plot, aes(lon,lat, size = Tcatch)) +ggtitle(" PS high seas  - all fleets")

library(dplyr)
datosPS_largePS<-datos3 %>% filter(fishery=="PSS"|fishery=='RIN')
datosPS_largePS_species<-datosPS_largePS %>% group_by(lon, lat, species)  %>% summarise (tcatch = sum(catch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
head(datosPS_largePS_species)

datosPS_largePS_species_wide<-datosPS_largePS_species %>% select(-c(tcatch)) %>% spread(species,freq)
n <- nrow(datosPS_largePS_species_wide)
datosPS_largePS_species_wide$region <- factor(1:n)
datosPS_largePS_species_wide<-as.data.frame(datosPS_largePS_species_wide)
head(datosPS_largePS_species_wide)
## need to replace the NA by ceros if i want to plot the pie.
datosPS_largePS_species_wide[is.na(datosPS_largePS_species_wide)] <- 0

library(scatterpie)

# ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/10000), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT"),  color="black", alpha=.8)+geom_scatterpie_legend(datosPS_largePS_species_wide$Totalcatch/10000, x=130, y=20) +ggtitle("Large PS - PSLS and PSFS")

#unity
ggplot(data = world) + geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT",'SWO'),  color="black", alpha=.8)+geom_scatterpie_legend(2, x=130, y=20) +ggtitle("Coastal PS (PSS, RIN) ")


################### Industrial Longline
gear='ELL'
datosPS_largePS_plot<-datos3 %>% filter(fishery==gear ) %>% group_by(lon, lat) %>% summarize (Tcatch=sum(catch,na.rm=TRUE))

# datosPS_largePS_plot<-datos3 %>% filter(fishery=="LL"|fishery=="ELL" |fishery=="FLL" ) %>% group_by(lon, lat) %>% summarize (Tcatch=sum(catch,na.rm=TRUE))
summary(datosPS_largePS_plot)


ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ geom_point(data=datosPS_largePS_plot, aes(lon,lat, size = Tcatch)) +ggtitle(" LL high seas  - all fleets")

library(dplyr)
gear='TROL'
datosPS_largePS<-datos3 %>% filter(fishery==gear)
# datosPS_largePS<-datos3 %>% filter(fishery=="LL"|fishery=="ELL" |fishery=="FLL")
datosPS_largePS_species<-datosPS_largePS %>% group_by(lon, lat, species)  %>% summarise (tcatch = sum(catch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
head(datosPS_largePS_species)

datosPS_largePS_species_wide<-datosPS_largePS_species %>% select(-c(tcatch)) %>% spread(species,freq)
n <- nrow(datosPS_largePS_species_wide)
datosPS_largePS_species_wide$region <- factor(1:n)
datosPS_largePS_species_wide<-as.data.frame(datosPS_largePS_species_wide)
head(datosPS_largePS_species_wide)
## need to replace the NA by ceros if i want to plot the pie.
datosPS_largePS_species_wide[is.na(datosPS_largePS_species_wide)] <- 0

library(scatterpie)

# ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/10000), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT"),  color="black", alpha=.8)+geom_scatterpie_legend(datosPS_largePS_species_wide$Totalcatch/10000, x=130, y=20) +ggtitle("Industrial Longline")

#unity
ggplot(data = world) + geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + 
  geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), data=datosPS_largePS_species_wide, 
                  cols=c("ALB","BET","SKJ","YFT",'SWO'),  color="black", alpha=.8)+
                  # cols=c("ALB","BET","SKJ","YFT"),  color="black", alpha=.8)+
  geom_scatterpie_legend(2, x=130, y=20) +
  ggtitle(paste("Troll (",gear,") "))
ggsave(paste0("../figures/presentation/ecologicaldata/piechart_unitty_fisheries_byspp",gear,".png"), width = 8, height = 5, dpi=200)


################### Coastal Longline
datosPS_largePS_plot<-datos3 %>% filter(fishery=="LLCO") %>% group_by(lon, lat) %>% summarize (Tcatch=sum(catch,na.rm=TRUE))
summary(datosPS_largePS_plot)


ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ geom_point(data=datosPS_largePS_plot, aes(lon,lat, size = Tcatch)) +ggtitle(" LL high seas  - all fleets")

library(dplyr)
datosPS_largePS<-datos3 %>% filter(fishery=="LLCO")
datosPS_largePS_species<-datosPS_largePS %>% group_by(lon, lat, species)  %>% summarise (tcatch = sum(catch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
head(datosPS_largePS_species)

datosPS_largePS_species_wide<-datosPS_largePS_species %>% select(-c(tcatch)) %>% spread(species,freq)
n <- nrow(datosPS_largePS_species_wide)
datosPS_largePS_species_wide$region <- factor(1:n)
datosPS_largePS_species_wide<-as.data.frame(datosPS_largePS_species_wide)
head(datosPS_largePS_species_wide)
## need to replace the NA by ceros if i want to plot the pie.
datosPS_largePS_species_wide[is.na(datosPS_largePS_species_wide)] <- 0

library(scatterpie)

# ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/10000), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT"),  color="black", alpha=.8)+geom_scatterpie_legend(datosPS_largePS_species_wide$Totalcatch/10000, x=130, y=20) +ggtitle("Coastal Longline")

#unity
ggplot(data = world) + geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT",'SWO'),  color="black", alpha=.8)+geom_scatterpie_legend(2, x=130, y=20) +
  ggtitle("Coastal Longline (LLCO) ")



################### High seas gillnets
datosPS_largePS_plot<-datos3 %>% filter(fishery=="GIOF") %>% group_by(lon, lat) %>% summarize (Tcatch=sum(catch,na.rm=TRUE))
summary(datosPS_largePS_plot)


ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ geom_point(data=datosPS_largePS_plot, aes(lon,lat, size = Tcatch)) +ggtitle(" LL high seas  - all fleets")

library(dplyr)
datosPS_largePS<-datos3 %>% filter(fishery=="GIOF")
datosPS_largePS_species<-datosPS_largePS %>% group_by(lon, lat, species)  %>% summarise (tcatch = sum(catch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
head(datosPS_largePS_species)

datosPS_largePS_species_wide<-datosPS_largePS_species %>% select(-c(tcatch)) %>% spread(species,freq)
n <- nrow(datosPS_largePS_species_wide)
datosPS_largePS_species_wide$region <- factor(1:n)
datosPS_largePS_species_wide<-as.data.frame(datosPS_largePS_species_wide)
head(datosPS_largePS_species_wide)
## need to replace the NA by ceros if i want to plot the pie.
datosPS_largePS_species_wide[is.na(datosPS_largePS_species_wide)] <- 0

library(scatterpie)

# ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/10000), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT"),  color="black", alpha=.8)+geom_scatterpie_legend(datosPS_largePS_species_wide$Totalcatch/10000, x=130, y=20) +ggtitle("Industrial Longline")

#unity
ggplot(data = world) + geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT",'SWO'),  color="black", alpha=.8)+geom_scatterpie_legend(2, x=130, y=20) +ggtitle("High seas Gillnets (GIOF)")



################### Coastal gillnets
datosPS_largePS_plot<-datos3 %>% filter(fishery=="GILL") %>% group_by(lon, lat) %>% summarize (Tcatch=sum(catch,na.rm=TRUE))
summary(datosPS_largePS_plot)


ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ geom_point(data=datosPS_largePS_plot, aes(lon,lat, size = Tcatch)) +ggtitle(" LL high seas  - all fleets")

library(dplyr)
datosPS_largePS<-datos3 %>% filter(fishery=="GILL")
datosPS_largePS_species<-datosPS_largePS %>% group_by(lon, lat, species)  %>% summarise (tcatch = sum(catch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
head(datosPS_largePS_species)

datosPS_largePS_species_wide<-datosPS_largePS_species %>% select(-c(tcatch)) %>% spread(species,freq)
n <- nrow(datosPS_largePS_species_wide)
datosPS_largePS_species_wide$region <- factor(1:n)
datosPS_largePS_species_wide<-as.data.frame(datosPS_largePS_species_wide)
head(datosPS_largePS_species_wide)
## need to replace the NA by ceros if i want to plot the pie.
datosPS_largePS_species_wide[is.na(datosPS_largePS_species_wide)] <- 0

library(scatterpie)

# ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/10000), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT"),  color="black", alpha=.8)+geom_scatterpie_legend(datosPS_largePS_species_wide$Totalcatch/10000, x=130, y=20) +ggtitle("Industrial Longline")

#unity
ggplot(data = world) + geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + 
  geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT",'SWO'),  color="black", alpha=.8)+geom_scatterpie_legend(2, x=130, y=20) +
  ggtitle("Coastal Gillnets (GILL)")

ggsave()


################### other coastal gears
datosPS_largePS_plot<-datos3 %>% filter(fishery=="HAND"|fishery=="TROL"|fishery=="BB") %>% group_by(lon, lat) %>% summarize (Tcatch=sum(catch,na.rm=TRUE))
summary(datosPS_largePS_plot)


ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+ geom_point(data=datosPS_largePS_plot, aes(lon,lat, size = Tcatch)) +ggtitle("Other major coastal fisheries (HAND, TROL, BB)")

library(dplyr)
datosPS_largePS<-datos3 %>% filter(fishery=="HAND"|fishery=="TROL"|fishery=="BB")
datosPS_largePS_species<-datosPS_largePS %>% group_by(lon, lat, species)  %>% summarise (tcatch = sum(catch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
head(datosPS_largePS_species)

datosPS_largePS_species_wide<-datosPS_largePS_species %>% select(-c(tcatch)) %>% spread(species,freq)
n <- nrow(datosPS_largePS_species_wide)
datosPS_largePS_species_wide$region <- factor(1:n)
datosPS_largePS_species_wide<-as.data.frame(datosPS_largePS_species_wide)
head(datosPS_largePS_species_wide)
## need to replace the NA by ceros if i want to plot the pie.
datosPS_largePS_species_wide[is.na(datosPS_largePS_species_wide)] <- 0

library(scatterpie)

# ggplot(data = world) + geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/10000), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT"),  color="black", alpha=.8)+geom_scatterpie_legend(datosPS_largePS_species_wide$Totalcatch/10000, x=130, y=20) +ggtitle("Industrial Longline")

#unity
ggplot(data = world) + geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "black", fill = "lightgrey")+coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) + 
  geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), data=datosPS_largePS_species_wide, cols=c("ALB","BET","SKJ","YFT",'SWO'),  color="black", alpha=.8)+
  geom_scatterpie_legend(2, x=130, y=20) +
  ggtitle("Other major coastal fisheries (HAND, TROL, BB)")















## TOTAL CATCHES BY SPECES - TUNA ATLAS
library(rtunaatlas)

source('/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/code/pie_map_AE.R')
con<-rtunaatlas::db_connection_tunaatlas_world()
rtunaatlas::pie_map(con,
                    df_input=d2 %>% filter (YEAR>=2007 & YEAR<=2017 & unit=="MT" & SPECIES %in% c("SWO","BET","YFT","SKJ","ALB")),
                    dimension_group_by="SPECIES",
                    df_spatial_code_list_name="areas_tuna_rfmos_task2",
                    function_pie_size = "unique",
                    number_of_classes=6
)


## TOTAL CATCHES BY GEAR - TUNA ATLAS
library(rtunaatlas)

source('/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/code/pie_map_AE.R')
con<-rtunaatlas::db_connection_tunaatlas_world()
rtunaatlas::pie_map(con,
                    df_input=d2 %>% filter (YEAR>=2007 & YEAR<=2017 & unit=="MT" & SPECIES %in% c("SWO","BET","YFT","SKJ","ALB")),
                    dimension_group_by="GEAR",
                    df_spatial_code_list_name="areas_tuna_rfmos_task2",
                    function_pie_size = "unique",
                    number_of_classes=6
)
