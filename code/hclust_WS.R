# hclust_on_provinces.
# aenieblas
# 9/8/19


rm(list=ls())
options(stringsAsFactors = FALSE)

analysis<-'PPOWMEOW' #options are : PPOW,PPOWMEOW,LONG
reassign<-FALSE # for PPOW analysis - TRUE indicates that some MEOW provinces will be reassigned as PPOW
# for PPOWMEOW analysis - TRUE indicates that PPOW provinces will be reassinged as MEOW.
layer<-'fishery' # options are : species, fishery, size

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
theme_set(theme_bw())

# world data for plotting
world <- ne_countries(scale = "medium", returnclass = "sf")

# iotc convention area 
iotc_shp='/home/ae/Documents/PERSONAL/COOOL/projects/Ecoregions_2021/data/indian_ocean/iotc.shp'
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
  meow_shp<-'/home/ae/Documents/PERSONAL/COOOL/projects/Ecoregions_2021/data/MEOW_FINAL/MEOW/meow_ecos.shp'
  meow<-readOGR(dsn=meow_shp,stringsAsFactors = F)
}

if(analysis=='LONG'){
  long_shp<-paste0(home,'longhurst_shapefiles/Longhurst_world_v4_2010.shp')
  ppow<-readOGR(dsn=long_shp, stringsAsFactors = F)
}

wgs<-proj4string(ppow) #this shows the projection, it's important that it's the same as projection of our data that we convert in shapefile later





#################################### cluster analysis of indicators ##################################################
## clustering by provinces, not by lat/lon to reveal if provinces already existing can be used as IOTC ecoregions

type='species' # options are 'species','fishery','combined'
doplot=TRUE
order='second'   # options are 'first', 'second'

prov='MEOWPPOW'
speciesdf<-read.csv(paste0(home,'species_final_indicator_df_',analysis,'__v2021.csv'))
speciesdf$Species2<-NULL
speciesdf$scale<-NULL
fisherydf<-read.csv(paste0(home,'fishery_final_indicator_df_',analysis,'__v2021.csv'))
fisherydf$Species2<-NULL
names(fisherydf)<-names(speciesdf)

if(type=='species'){layer='species'}
if(type=='fishery'|type=='combined'){layer='fishery'}
# options are 'species' or 'fishery' (from the indicator analysis)
## read in analysis data
datos<-read.csv(paste0(home,'datos_',layer,'_PPOWMEOW_unforced_v2021.csv'))

if(type=='combined'){sfdf<-rbind(speciesdf,fisherydf)}
if(type=='species'){sfdf<-speciesdf}
if(type=='fishery'){sfdf<-fisherydf}
# names(sfdf)<-names(speciesdf)}

PROVINC=data.frame(id=seq(1:length(unique(sfdf$PROVINC))), name=unique(sfdf$PROVINC))

## unforced assignement
PROVINC$TYPE<-c('MEOW','PPOW','MEOW','MEOW','MEOW','MEOW','PPOW','PPOW','PPOW','MEOW','PPOW','MEOW','MEOW','MEOW','PPOW','MEOW','MEOW','PPOW','PPOW','MEOW','MEOW','MEOW','MEOW','MEOW')
sfdf$TYPE<-NULL
for(i in 1:26){    sfdf[which(sfdf$PROVINC==PROVINC$name[i]),'TYPE']<-PROVINC$TYPE[i]  }

## create data frame of species/fidelity indicators
### hclust the indicator value?s
library(dplyr)
sp<-unique(sfdf$Species)
kDF<-data.frame(PROVINCE=unique(sfdf$PROVINC))
for(s in 1:length(sp)){
  xx<-sfdf[which(sfdf$Species==sp[s]),c("PROVINCE","SpIndicator")]
  kDF<-left_join(kDF,xx,by='PROVINCE')
  kDF[which(is.na(kDF[,s+1])==TRUE),s+1]<- 0
}
colnames(kDF)<-c('PROVINC',as.character(sp))
kDF_scaled<-as.data.frame(scale(kDF[,2:(length(sp)+1)]))
colnames(kDF_scaled)<-as.character(sp)
kDF_scaled[is.na(kDF_scaled)]<-0

kk<-kmeans(kDF_scaled,7,iter.max = 30)
kDF$cluster<-kk$cluster



## hclust 
dist_mat <- dist(kDF_scaled, method = 'euclidean')
if(type=='species'){c_method='average'}
if(type=='fishery'|type=='combined'){c_method='complete'}
hclust_avg <- hclust(dist_mat, method = c_method)
hclust_avg$labels<-PROVINC$name
plot(hclust_avg)

if(prov=='MEOWPPOW' & type=='species'){
  
  PROVINC1<-PROVINC
  if(order=='first'){
    PROVINC1$cluster<-NA
    PROVINC1[hclust_avg$order[1:2],'cluster']<-1
    PROVINC1[hclust_avg$order[3],'cluster']<-2
    PROVINC1[hclust_avg$order[4:24],'cluster']<-3
  }
  if(order=='second'){
    PROVINC1$cluster<-NA
    PROVINC1[hclust_avg$order[1:2],'cluster']<-1
    PROVINC1[hclust_avg$order[3],'cluster']<-2
    PROVINC1[hclust_avg$order[4:8],'cluster']<-3
    PROVINC1[hclust_avg$order[9:15],'cluster']<-4
    PROVINC1[hclust_avg$order[16:24],'cluster']<-5
    
  }
}


if(prov=='PPOWMEOW' & type=='fishery'){
  PROVINC1<-PROVINC
if(order=='first'){
  PROVINC1$cluster<-NA
  PROVINC1[hclust_avg$order[1],'cluster']<-1
  PROVINC1[hclust_avg$order[2:19],'cluster']<-2
  PROVINC1[hclust_avg$order[20:21],'cluster']<-3
  PROVINC1[hclust_avg$order[22:24],'cluster']<-4
}
if(order=='second'){
  PROVINC1$cluster<-NA
  # PROVINC1[hclust_avg$order[1],'cluster']<-1
  # PROVINC1[hclust_avg$order[2:19],'cluster']<-2
  # PROVINC1[hclust_avg$order[20:21],'cluster']<-3
  # PROVINC1[hclust_avg$order[22:24],'cluster']<-4
  
}
  
}


if(prov=='MEOWPPOW' & type=='combined'){
  
    PROVINC1<-PROVINC
    if(order=='first'){
      PROVINC1$cluster<-NA
      PROVINC1[hclust_avg$order[1],'cluster']<-1
      PROVINC1[hclust_avg$order[2:3],'cluster']<-2
      PROVINC1[hclust_avg$order[4],'cluster']<-3
      PROVINC1[hclust_avg$order[5:22],'cluster']<-4
      PROVINC1[hclust_avg$order[23:24],'cluster']<-5
    }
    if(order=='second'){
      PROVINC1$cluster<-NA
      # PROVINC1[hclust_avg$order[1],'cluster']<-1
      # PROVINC1[hclust_avg$order[2:19],'cluster']<-2
      # PROVINC1[hclust_avg$order[20:21],'cluster']<-3
      # PROVINC1[hclust_avg$order[22:24],'cluster']<-4
    }
}
    

if(doplot==TRUE){
h4_custom<-PROVINC1
colnames(h4_custom)<-c('id','PROVINC','TYPE','cluster')

datos$cluster<-NULL
for(p in 1:dim(h4_custom)[1]){
  ## find provinces back in original shapefile and add cluster
  datos[which(datos$PROVINC==h4_custom$PROVINC[p]),'cluster']<-h4_custom$cluster[p]
}

# load meow shapefile
meow_shp<-'/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/MEOW_FINAL/MEOW/meow_ecos.shp'
meow_bgcp<-readOGR(dsn=meow_shp,stringsAsFactors = F)

# find bathy
library(marmap)
bath<-getNOAA.bathy(20,150,-55,30,resolution=60,keep=T)
bathy<-fortify(bath)
bathy[which(bathy$z>0),'z']<-0


# plot clusters
library(RColorBrewer)
cols<-brewer.pal(8,'Pastel2')[1:8]
datos1<-as.data.frame(datos)

datos1<-datos1[-grep('Antarctic',datos$PROVINC),]

ggplot(data = world) +
  geom_tile(data = datos1, aes(x=Long,y=Lat,fill=as.factor(cluster)))+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_polygon(data = ppow, aes(x = long, y = lat, group = group), colour = "black", fill = NA)+
  geom_polygon(data = meow_bgcp, aes(x = long, y = lat, group = group), colour = "black", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  geom_contour(data=bathy,aes(x=x,y=y,z=z),
               breaks=c( -200, -1000),
               colour="white", size=0.7
  )+
  # geom_point(data=datos1, aes(Long,Lat, size = tCatch)) +
  scale_fill_discrete(cols)+
  coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
  ggtitle(paste(type,order,"order clusters"))
}
# ggsave(paste0('../figures/spatial/cluster_',type,'_',analysis,'_PPOWMEOW_v5_1st.png'), width = 8, height = 5, dpi=200)


