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
  meow_shp<-'/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/MEOW_FINAL/MEOW/meow_ecos.shp'
  meow<-readOGR(dsn=meow_shp,stringsAsFactors = F)
}

if(analysis=='LONG'){
  long_shp<-paste0(home,'longhurst_shapefiles/Longhurst_world_v4_2010.shp')
  ppow<-readOGR(dsn=long_shp, stringsAsFactors = F)
}

wgs<-proj4string(ppow) #this shows the projection, it's important that it's the same as projection of our data that we convert in shapefile later


#################################### cluster analysis of indicators ##################################################
## clustering by provinces, not by lat/lon to reveal if provinces already existing can be used as IOTC ecoregions
prov='MEOWPPOW'
speciesdf<-read.csv(paste0(home,'species_final_indicator_df_',analysis,'__v2021.csv'))
speciesdf$Species2<-NULL
speciesdf$scale<-NULL
fisherydf<-read.csv(paste0(home,'fishery_final_indicator_df_',analysis,'__v2021.csv'))
fisherydf$Species2<-NULL

type='species' # options are 'species','fishery','combined'
layer='species'# options are 'species' or 'fishery' (from the indicator analysis)
## read in analysis data
datos<-read.csv(paste0(home,'datos_',layer,'_PPOWMEOW_unforced_v2021.csv'))

if(type=='combined'){sfdf<-rbind(speciesdf,fisherydf)}
if(type=='species'){sfdf<-speciesdf}
if(type=='fishery'){sfdf<-fisherydf}

PROVINC=data.frame(id=seq(1:length(unique(sfdf$PROVINC))), name=unique(sfdf$PROVINC))

## manual assignment...!!
## forced assignment
# PROVINC$TYPE<-c('MEOW','PPOW','MEOW','MEOW','MEOW','MEOW','PPOW','PPOW','PPOW','MEOW','PPOW','MEOW','MEOW','MEOW','MEOW','MEOW','PPOW','MEOW',
                # 'MEOW','MEOW','MEOW')
## unforced assignement
PROVINC$TYPE<-c('MEOW','PPOW','MEOW','MEOW','MEOW','MEOW',
                'PPOW','PPOW','PPOW','MEOW','PPOW',
                'MEOW','MEOW','MEOW','PPOW','MEOW',
                'MEOW','PPOW','PPOW','MEOW','MEOW',
                'MEOW','MEOW','MEOW')

write.csv(PROVINC,file=paste0(home,'PROVINC_cluster_info_unforced_v2021.csv'),row.names = FALSE)

sfdf$TYPE<-NULL
for(i in 1:26){    sfdf[which(sfdf$PROVINC==PROVINC$name[i]),'TYPE']<-PROVINC$TYPE[i]  }

## perform clustering on separate PPOW and MEOW


if(prov=='MEOW'){sfdf<-subset(sfdf,sfdf$TYPE=='MEOW')
PROVINC<-subset(PROVINC,PROVINC$TYPE=='MEOW')
rownames(PROVINC)<-seq(1:length(unique(PROVINC$name)))}
if(prov=='PPOW'){sfdf<-subset(sfdf,sfdf$TYPE=='PPOW')
PROVINC<-subset(PROVINC,PROVINC$TYPE=='PPOW')
rownames(PROVINC)<-seq(1:length(unique(PROVINC$name)))}

## combine the DFs
# FinalDF<-left_join(Specificity_Sp_per_Prov[c('PROVINC','Species','meancatchSp','SumCatchSp','specificity')],Fidelity_Sp_per_Prov,by=c('PROVINC','Species'))

# sfdf<-left_join(fisherydf[c('PROVINC','Species')],speciesdf,by=c('PROVINC'))

## create data frame of species/fidelity indicators
### hclust the indicator value?s
library(dplyr)
sp<-unique(sfdf$Species)
# kDF<-NULL
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

kk<-kmeans(kDF_scaled,5,iter.max = 30)
kDF$cluster<-kk$cluster



## hclust 
dist_mat <- dist(kDF_scaled, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
hclust_avg$labels<-PROVINC$name
plot(hclust_avg)
# cut_avg <- cutree(hclust_avg, k = 5)
# plot(hclust_avg)
# rect.hclust(hclust_avg , k = 7, border = 2:6)
# abline(h = 5, col = 'red')

if(prov=='MEOWPPOW' & type=='species'){## complete
  
  C1<-PROVINC[hclust_avg$order[c(1:2)],'id']
  C2<-PROVINC[hclust_avg$order[3],'id']
  C3<-PROVINC[hclust_avg$order[4:10],'id']
  
  C4<-PROVINC[hclust_avg$order[11:15],'id']
  C5<-PROVINC[hclust_avg$order[16:24],'id']
  
  C11<-PROVINC[hclust_avg$order[1:2],'id']
  C22<-PROVINC[hclust_avg$order[3],'id']
  C33<-PROVINC[hclust_avg$order[4:24],'id']
  
  
}


if(prov='PPOWMEOW' & type='fishery'){## single

  C1<-PROVINC[hclust_avg$order[c(1:5)],'id']
  C2<-PROVINC[hclust_avg$order[6:11],'id']
  C3<-PROVINC[hclust_avg$order[12:15],'id']
  
  C4<-PROVINC[hclust_avg$order[16:19],'id']
  C5<-PROVINC[hclust_avg$order[20:24],'id']
  
  C11<-PROVINC[hclust_avg$order[1:5],'id']
  C22<-PROVINC[hclust_avg$order[6:15],'id']
  C33<-PROVINC[hclust_avg$order[16:24],'id']

  
}


if(prov='MEOWPPOW' & type='combined'){## ward.D2
  # abline(h=1.7,col='red')
  # text(13,3.5,'1',cex=2)
  # text(4.5,-.5,'2',cex=2)
  # legend('topright','A',bty='n',cex=3)
  
  C1<-PROVINC[hclust_avg$order[c(1:6)],'id']
  C2<-PROVINC[hclust_avg$order[7:9],'id']
  C3<-PROVINC[hclust_avg$order[10:14],'id']
  C4<-PROVINC[hclust_avg$order[15:17],'id']
  C5<-PROVINC[hclust_avg$order[18:20],'id']
  C6<-PROVINC[hclust_avg$order[21:24],'id']
  
  C11<-PROVINC[hclust_avg$order[1:6],'id']
  C22<-PROVINC[hclust_avg$order[7:11],'id']
  C33<-PROVINC[hclust_avg$order[12:24],'id']

}



PROVINC<-read.csv(paste0(home,'PROVINC_cluster_info_unforced_v2021.csv'))

PROVINC$cluster<-NA
PROVINC[!is.na(match(PROVINC$id,C1)),'cluster']<-1
PROVINC[!is.na(match(PROVINC$id,C2)),'cluster']<-2
PROVINC[!is.na(match(PROVINC$id,C3)),'cluster']<-3
PROVINC[!is.na(match(PROVINC$id,C4)),'cluster']<-4
PROVINC[!is.na(match(PROVINC$id,C5)),'cluster']<-5
PROVINC[!is.na(match(PROVINC$id,C6)),'cluster']<-6
PROVINC[!is.na(match(PROVINC$id,C7)),'cluster']<-7
PROVINC[!is.na(match(PROVINC$id,C8)),'cluster']<-8

PROVINC[!is.na(match(PROVINC$id,C11)),'cluster']<-1
PROVINC[!is.na(match(PROVINC$id,C22)),'cluster']<-2
PROVINC[!is.na(match(PROVINC$id,C33)),'cluster']<-3
PROVINC[!is.na(match(PROVINC$id,C44)),'cluster']<-4

# write.csv(PROVINC,file=paste0(home,'hclust/hclust_single_k3_custom_',type,'_PPOWMEOW_v5_2nd.csv'),row.names=FALSE)
write.csv(PROVINC,file=paste0(home,'hclust/hclust_single_k3_custom_',type,'_PPOWMEOW_v2021.csv'),row.names=FALSE)

# hclust_avg$labels<-PROVINC$name
# plot(hclust_avg)
# abline(3.5,0,col='red')

h4_custom<-read.csv(paste0(home,'hclust/hclust_single_k3_custom_',type,'_PPOWMEOW_v2021.csv'))
colnames(h4_custom)<-c('id','PROVINC','TYPE','cluster')
# h4_custom$PROVINC <- factor(h4_custom$PROVINC, levels=levels(datos$PROVINC))
# h4_custom$cluster <- factor(h4_custom$cluster, levels=levels(datos$cluster))

# datos<-datos[-grep('Antarctic',datos$PROVINC),]
datos$cluster<-NULL
for(p in 1:dim(h4_custom)[1]){
  ## find provinces back in original shapefile and add cluster
  datos[which(datos$PROVINC==h4_custom$PROVINC[p]),'cluster']<-h4_custom$cluster[p]
}


meow_shp<-'/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/MEOW_FINAL/MEOW/meow_ecos.shp'
meow_bgcp<-readOGR(dsn=meow_shp,stringsAsFactors = F)

library(marmap)
# resolution=300 (5 deg); resolution=60 (1 deg)
bath<-getNOAA.bathy(20,150,-55,30,resolution=60,keep=T)

# bath[which(bath>0)]<-0
bathy<-fortify(bath)
bathy[which(bathy$z>0),'z']<-0



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
  ggtitle(paste0(layer," second order clusters "))
# ggsave(paste0('../figures/spatial/cluster_',type,'_',analysis,'_PPOWMEOW_v5_2nd.png'), width = 8, height = 5, dpi=200)
ggsave(paste0('../figures/spatial/cluster_',type,'_',analysis,'_PPOWMEOW_v2021.png'), width = 8, height = 5, dpi=200)




ppow@data$cluster<-as.factor(as.character(ppow@data$cluster))

ppow@data$id<-seq(1:269)
ppow@data$id = rownames(ppow@data)
ppow.points = fortify(ppow, region="id")
ppow.df = join(ppow.points, ppow@data, by="id")

# ppow_io_ids<-c(seq(87,111,by=1),seq(117,120),seq(138,141),144,145,seq(204,211))
ppow_io_ids<-c(seq(87,111,by=1),seq(118,120),132,141,144,145,seq(205,210),192,193)

# datos_p<-subset(datos,is.na(datos$PROVINCE))
# dropsM <-c('ECO_CODE','PROV_CODE','RLM_CODE','ALT_CODE','ECO_CODE_X','Lat_Zone') # list of col names
# dropsP<-c('ECO_CODE','PROV_CODE','RLM_CODE','ALT_CODE','ECO_CODE_X','Lat_Zone','TYPE','BIOME')
# # dropsP<-c('ECOREGION','REALM')
# datos_m <- datos_m[,!(names(datos_m) %in% dropsM)] #remove columns "AREA" and "PERIMETER"
# datos_p <- datos_p[,!(names(datos_p) %in% dropsP)] #remove columns "AREA" and "PERIMETER"
ppow_io<-subset(ppow.df,ppow.df$long>20)
ppow_io<-subset(ppow_io,ppow_io$long<150)

ppow_io<-ppow.df[ppow.df$TYPE=='MEOW',]

library(marmap)
# download for every grid cell of IOTC CA
# resolution=300 (5 deg); resolution=60 (1 deg)
bath<-getNOAA.bathy(20,150,-55,30,resolution=60,keep=T)

# bath[which(bath>0)]<-0
bathy<-fortify(bath)
bathy[which(bathy$z>0),'z']<-0



ggplot(data=world) + 
  geom_polygon(data=ppow_io,aes(long,lat,group=group,fill=PROVINC),color='black') + 
  geom_contour(data=bathy,aes(x=x,y=y,z=z),
               breaks=c( -200, -1000),
               colour="white", size=0.7
  )+
  # geom_polygon() +
  # geom_path(color="black") +
  # coord_equal() +
  # scale_fill_brewer("Viridis")+
  # geom_sf(data=world,color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)


# ggplot(ppow_io) + 
#   aes(long,lat,group=group,fill=cluster) + 
#   geom_polygon() +
#   geom_path(color="black") +
#   coord_equal() +
#   scale_fill_brewer("Viridis")+
#   geom_sf(data=world,color = "darkgrey", fill = "lightgrey")+
#   coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)


## plot with cluster as fill
ggplot(data = world) +
  # geom_point(data=df3, aes(lon,lat, size = CATCH))+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  # geom_polygon(data = meow_bgcp, aes(x = long, y = lat, group = group), colour = "dark blue", fill = NA)+
  geom_polygon(data = ppow, aes(x = long, y = lat, group = group,fill = BIOME), colour = "cyan")+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)+
  annotate("text",label='B',x=20, y=30, size=12)
map + theme_void()



