#neritic.R
# ae nieblas
# 4/8/19
# description : visualises data layers of neritic species catches, using the publicly-available data. 

rm(list=ls())
options(stringsAsFactors=FALSE)
home<-'/home/ae/Documents/PERSONAL/COOOL/projects/Ecoregions_2021/'
setwd(home)

layer='shark' ## other option is "neritic" or "shark"
start_year=2005
end_year=2019

main_fisheries<-c('PS','FLL','LL','LLCO','GIOF','ELL','GILL','TROL','HAND','DSEI','LG','BB')

library(pacman)
p_load('rgdal','ggplot2','sf','rnaturalearth','rnaturalearthdata','reshape','rgeos','plyr','sp','tidyr','raster','plotly','sp','maps','rgeos','plotrix','dplyr','data.table')

theme_set(theme_bw())

# world data for plotting
world <- ne_countries(scale = "medium", returnclass = "sf")

source(paste0(home,'/code/pie_map_AE.R'))
source(paste0(home,'/code/iotc_ca_regrid.R'))
source(paste0(home,'/code/iotc_regrid.R'))
source(paste0(home,'/code/regrid_datalayers.R'))


# iotc convention area 
iotc_shp=paste0(home,'/data/indian_ocean/iotc.shp')
iotc <- readOGR(dsn = iotc_shp, stringsAsFactors = F)


# STEP 1 : REGRID DATA LAYERS TO IOTC 5X5 GRID and keep only those within the IOTC convention area
iotc_5deg_shp=paste0(home,'/data/IOTC-2018-TCAC04-DATA04_-_5x5_grids_shapefile/IOTC_TCAC_GRID_5x5.shp')
shp <- readOGR(dsn = iotc_5deg_shp, stringsAsFactors = F)

# 1) read the IOTC CE data for all gears/vessels
ll<-read.csv(paste0(home,'/data/IOTC-2020-DATASETS-CE_All_Gear/IOTC-2020-WPTT22-DATA04-CELongline.csv'))
sfc<-read.csv(paste0(home,'/data/IOTC-2020-DATASETS-CE_All_Gear/IOTC-2020-WPTT22-DATA05-CESurface.csv'))
othr<-read.csv(paste0(home,'/data/IOTC-2020-DATASETS-CE_All_Gear/IOTC-2020-WPTT22-DATA06-CEOther.csv'))

# neritic spp list of the IOTC
if(layer=='neritic'){ner<-c('LOT', 'KAW','FRI','BLT','COM','GUT','FRZ')}
if(layer=='shark'){ner<-c('BSH','BTH','FAL','OCS','POR','PSK','SMA','MSK','RSK','SPY','THR','SKH')}

if(!file.exists(paste0(home,'data/catch_',layer,'_2021.csv'))){
# extract spp
# ll
ll_col<-colnames(ll)[grep(paste0(ner,collapse='|'),colnames(ll))]
if(length(ll_col)>0){catch_ll<-ll[,c(colnames(ll)[1:11],ll_col[grep('MT',ll_col)])]
catch_ll<-reshape2::melt(catch_ll,id.vars=colnames(catch_ll)[1:11])
}else{catch_ll<-NULL}
# sfc
sfc_col<-colnames(sfc)[grep(paste0(ner,collapse='|'),colnames(sfc))]
if(length(sfc_col)>0){catch_sfc<-sfc[,c(colnames(sfc)[1:11],sfc_col)]
catch_sfc<-reshape2::melt(catch_sfc,id.vars=colnames(catch_sfc)[1:11])
}else{catch_sfc<-NULL}
# othr
othr_col<-colnames(othr)[grep(paste0(ner,collapse='|'),colnames(othr))]
if(length(othr_col)>0){catch_othr<-othr[,c(colnames(othr)[1:11],othr_col[grep('MT',othr_col)])] 
catch_othr<-reshape2::melt(catch_othr,id.vars=colnames(catch_othr)[1:11])
}else{catch_othr<-NULL}

catch<-rbind(catch_ll,catch_sfc,catch_othr)
catch$Species<-apply(catch,1,function(x) unlist(strsplit(as.character(x['variable']),'\\.'))[1])
write.csv(catch,file=paste0(home,'data/catch_',layer,'_2021.csv'),row.names = FALSE)
}else{catch<-read.csv(paste0(home,'data/catch_',layer,'_2021.csv'))}

## sum catch within each grid cell for every YEAR
library(plyr)
if(!file.exists(paste0(home,'data/IOTC_',layer,'_catch_annual_sum.csv'))){
d2<-ddply(catch,.(Fleet, Year, Gear,Grid,Species), function(x) catch=sum(x$value,na.rm=T))
colnames(d2)<-c("FLEET_CODE_HASHED","YEAR","GEAR","GRID_CODE","SPECIES","RAISED_CATCHES_MT")
# d1<-na.omit(d2$GRID_CODE)
write.csv(d2,paste0(home,'data/IOTC_',layer,'_catch_annual_sum.csv'),row.names = FALSE)
}else{d2<-read.csv(paste0(home,'data/IOTC_',layer,'_catch_annual_sum.csv'))}

d1<- iotc_1_to_5_regrid(iotc_5deg_shp=iotc_5deg_shp,
                        iotc_df=paste0(home,'data/IOTC_',layer,'_catch_annual_sum.csv'),data=2017)

detach("package:plyr", unload=TRUE)
library(dplyr)
library(tidyr)

# start_year<-2005

## TOTAL SUM
nerc_sum<-d1 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(lon,lat) %>% summarize (tCatch=sum(RAISED_CATCHES_MT,na.rm=TRUE))
names(nerc_sum)<-c("lon", "lat","catch")

# nerc[which(nerc$catch>60000),'catch']<-90000
library(RColorBrewer)
cols<-brewer.pal(8,'YlGnBu')
ggplot(data = world) +
  geom_tile(data=nerc_sum, aes(lon,lat, fill = catch)) +
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  # geom_point(data=nerc, aes(lon,lat, size = catch)) +
  scale_fill_gradientn(colours=cols,limits = c(0,max(nerc_sum$catch)))+
  coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
  ggtitle(paste0("Total catches of ",layer," species  ",start_year,'-',end_year))
ggsave(paste0('figures/presentation/ecologicaldata/totalcatches_',layer,'_sum',start_year,'-',end_year,'.png'), width = 8, height = 5, dpi=200)


start_year<-2015
## TOTAL MEDIAN
nerc_sum<-d1 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(lon,lat) %>% summarize (tCatch=median(RAISED_CATCHES_MT,na.rm=TRUE))
names(nerc_sum)<-c("lon", "lat","catch")


nerc_sum[which(nerc_sum$catch>quantile(nerc_sum$catch,probs=0.99)),'catch']<-quantile(nerc_sum$catch,probs=0.99)
# nerc_sum[which(nerc_sum$catch==0),'catch']<-NA

# nerc[which(nerc$catch>60000),'catch']<-90000
library(RColorBrewer)
cols<-brewer.pal(8,'YlGnBu')
ggplot(data = world) +
  geom_tile(data=nerc_sum, aes(lon,lat, fill = catch)) +
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  # geom_point(data=nerc, aes(lon,lat, size = catch)) +
  scale_fill_gradientn(colours=cols,limits = c(0,max(nerc_sum$catch,na.rm=T)))+
  coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
  ggtitle(paste0("Total catches of ",layer," species  ",start_year,'-',end_year))
ggsave(paste0('figures/presentation/ecologicaldata/totalcatches_',layer,'_sum',start_year,'-',end_year,'.png'), width = 8, height = 5, dpi=200)


## MEDIAN BY SPP
nerc_spp_median<-d1 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(SPECIES, lon,lat) %>% summarize (tCatch=median(RAISED_CATCHES_MT,na.rm=TRUE))
names(nerc_spp_median)<-c("species","lon", "lat","catch")


cols<-brewer.pal(8,'YlGnBu')
# nerc_spp_median[which(nerc_spp_median$catch>quantile(nerc_spp_median$catch,probs=0.999)),'catch']<-quantile(nerc_spp_median$catch,probs=0.999)
ggplot(data = world) +
  geom_tile(data=nerc_spp_median, aes(lon,lat, fill = catch)) +
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  # geom_point(data=nerc_spp_sum, aes(lon,lat, size = catch)) +
  facet_wrap(~species)+
  scale_fill_gradientn(colours=cols,limits = c(0,max(nerc_spp_median$catch)))+
  scale_x_continuous(breaks=seq(20,140,by=50))+
  coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
  ggtitle(paste0("Median of catches of ",layer," species  "))
ggsave(paste0('figures/presentation/ecologicaldata/totalcatches_facetwrap_spp_',layer,'_median.png'), width = 8, height = 5, dpi=200)

## SUM BY SPP
nerc_spp_sum<-d1 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(SPECIES, lon,lat) %>% summarize (tCatch=sum(RAISED_CATCHES_MT,na.rm=TRUE))
names(nerc_spp_sum)<-c("species","lon", "lat","catch")

target_spp<-c('BSH','MSK','FAL')
nerc_spp_sum<-nerc_spp_sum[grep(paste0(target_spp,collapse='|'),nerc_spp_sum$species),]
nerc_spp_sum[which(nerc_spp_sum$catch==0),'catch']<-NA

nerc_spp_sum[which(nerc_spp_sum$catch>quantile(nerc_spp_sum$catch,probs=0.999,na.rm=T)),'catch']<-quantile(nerc_spp_sum$catch,probs=0.999,na.rm=T)

cols<-brewer.pal(8,'YlGnBu')
ggplot(data = world) +
  geom_tile(data=nerc_spp_sum, aes(lon,lat, fill = catch)) +
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  # geom_point(data=nerc_spp_sum, aes(lon,lat, size = catch)) +
  facet_wrap(~species)+
  scale_fill_gradientn(colours=cols,limits = c(0,max(nerc_spp_sum$catch)))+
  scale_x_continuous(breaks=seq(20,140,by=50))+
  coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
  ggtitle(paste0("Sum of catches of ",layer," species  "))
ggsave(paste0('figures/presentation/ecologicaldata/totalcatches_facetwrap_spp_',layer,'_sum.png'), width = 8, height = 5, dpi=200)


## SUM BY SPP
nerc_spp_sum<-d1 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(SPECIES) %>% summarize (tCatch=sum(RAISED_CATCHES_MT,na.rm=TRUE))%>% arrange(-tCatch)
names(nerc_spp_sum)<-c("species","catch")

mp<-barplot(nerc_spp_sum$catch/10000,ylab='catch (MT/10,000)')
axis(side = 1, at = mp, labels = nerc_spp_sum$species, tcl = -0.2)



## MEDIAN BY GEAR
nerc_gear_median<-d1 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(GEAR, lon,lat) %>% summarize (tCatch=median(RAISED_CATCHES_MT,na.rm=TRUE))
names(nerc_gear_median)<-c("gear","lon", "lat","catch")

cols<-brewer.pal(8,'YlGnBu')
ggplot(data = world) +
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  # geom_point(data=nerc_gear_median, aes(lon,lat, size = catch)) +
  # geom_tile(data=nerc_gear_sum, aes(lon,lat, fill = lnerc)) +
  geom_tile(data=nerc_gear_median, aes(lon,lat, fill = catch)) +
  facet_wrap(~gear)+
  scale_fill_gradientn(colours=cols,limits = c(0,max(nerc_gear_median$catch)))+
  scale_x_continuous(breaks=seq(20,150,by=50))+
  coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
  ggtitle(paste0("Median catches of ",layer," gear  "))
ggsave(paste0('figures/presentation/ecologicaldata/totalcatches_facetwrap_gear_',layer,'_median.png'), width = 8, height = 5, dpi=200)

## SUM BY GEAR
nerc_gear_sum<-d1 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(GEAR, lon,lat) %>% summarize (tCatch=sum(RAISED_CATCHES_MT,na.rm=TRUE))
names(nerc_gear_sum)<-c("gear","lon", "lat","catch")

nerc_gear_sum[which(nerc_gear_sum$catch>quantile(nerc_gear_sum$catch,probs=0.999)),'catch']<-quantile(nerc_gear_sum$catch,probs=0.999)

cols<-brewer.pal(8,'YlGnBu')
ggplot(data = world) +
  geom_tile(data=nerc_gear_sum, aes(lon,lat, fill = catch)) +
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  # geom_point(data=nerc_gear_median, aes(lon,lat, size = catch)) +
  # geom_tile(data=nerc_gear_sum, aes(lon,lat, fill = lnerc)) +
  facet_wrap(~gear)+
  scale_fill_gradientn(colours=cols,limits = c(0,max(nerc_gear_sum$catch)))+
  scale_x_continuous(breaks=seq(20,150,by=50))+
  coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
  ggtitle(paste0("Sum catches of ",layer," gear  "))
ggsave(paste0('figures/presentation/ecologicaldata/totalcatches_facetwrap_gear_',layer,'_sum.png'), width = 8, height = 5, dpi=200)


## CATCH BY YEAR
nerc_year_sum<-d1 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(YEAR) %>% summarize (tCatch=sum(RAISED_CATCHES_MT,na.rm=TRUE))
names(nerc_year_sum)<-c("year","catch")

ggplot(data=nerc_year_sum,aes(x=year,y=catch))+
  geom_bar(stat='identity')+
  ylab('Total catch (MT) by year')
ggsave(paste0('figures/presentation/ecologicaldata/totalcatch_by_year_',layer,'.png'), width = 8, height = 5, dpi=200)


if(sublayer=='observations'){
  
  ner_yr_obs<-d1 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(YEAR) %>% tally()

  ggplot(data=ner_yr_obs,aes(x=YEAR,y=n))+
    geom_bar(stat='identity')+
    # facet_wrap(~SPECIES)+
    ylab('Number of observations')
    ggsave(paste0('figures/presentation/ecologicaldata/numberobservations_year_',layer,'.png'), width = 8, height = 5, dpi=200)

    
    obs<-d1 %>% filter(YEAR<=end_year & YEAR>=start_year) %>% group_by(lon,lat) %>% tally()
    
    library(RColorBrewer)
    cols<-brewer.pal(8,'YlOrRd')
    ggplot(data = world) +
      geom_tile(data = obs, aes(x=lon,y=lat,fill=n))+
      geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
      geom_sf(color = "darkgrey", fill = "lightgrey")+
      # facet_wrap(~SPECIES,ncol=2)+
      scale_fill_gradientn(colours=cols)+
      coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
      ggtitle(paste0("Number of observations  ",start_year,'-',end_year))
    ggsave(paste0("figures/presentation/ecologicaldata/numberobservations_",layer,"_all.png"), width = 8, height = 5, dpi=200)
    
}


### SCATTERPIE
library(dplyr)
library(tidyr)
# datos_na<-subset(datos,!is.na(datos$PROVINC))
datos_total_catches_species<-d1 %>% group_by(lon, lat, SPECIES)  %>% summarise (tcatch = sum(RAISED_CATCHES_MT,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE))
head(datos_total_catches_species)

datos_total_catches_species<-datos_total_catches_species[grep(paste0(target_spp,collapse='|'),datos_total_catches_species$SPECIES),]

datos_total_catches_species_wide<-datos_total_catches_species %>% select(-c(tcatch)) %>% spread(SPECIES,freq)
n <- nrow(datos_total_catches_species_wide)
datos_total_catches_species_wide$region <- factor(1:n)
datos_total_catches_species_wide<-as.data.frame(datos_total_catches_species_wide)
head(datos_total_catches_species_wide)
## need to replace the NA by zeros if i want to plot the pie.
datos_total_catches_species_wide[is.na(datos_total_catches_species_wide)] <- 0
# 
library(scatterpie)

# datos_df<-as.data.frame(datos)
col_plot_names<-unique(datos_total_catches_species$SPECIES)
# ## unity 
ggplot(data = world) + 
  # geom_tile(data = datos_df, aes(x=Long,y=Lat,fill=ProvCode))+
  # scale_fill_discrete(cols)+
  geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
  geom_sf(color = "darkgrey", fill = "lightgrey")+
  geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2),
                  data=datos_total_catches_species_wide, cols=as.character(col_plot_names),
                  color="black", alpha=.8)+geom_scatterpie_legend(2, x=130, y=20)+
  # geom_polygon(data = ppow, aes(x = long, y = lat, group = group), colour = "black", fill = NA)+
  coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE) 
ggsave(paste0('../figures/presentation/ecologicaldata/piechart_unity_',layer,'_',type,'.png'), width = 8, height = 5, dpi=200)






# 
# 
# 
# 
# 
# library(RColorBrewer)
# cols<-brewer.pal(8,'YlOrRd')
# datos1<-as.data.frame(datos)
# 
# ################################## catch distribution #############################################
# if(sublayer=='observations'){
#   names(datos1)<-c('Long','Lat','Species','Obs','lat','lon')
#   # datos1$Obs<-ddply
#   ggplot(data = world) +
#     geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
#     geom_sf(color = "darkgrey", fill = "lightgrey")+
#     geom_tile(data = datos1, aes(x=Long,y=Lat,fill=Obs))+
#         # geom_point(data=datos1, aes(Long,Lat, size = tCatch)) +
#     # facet_wrap(~Species)+
#     scale_fill_gradientn(colours=cols)+
#     coord_sf(xlim = c(15,154), ylim = c(-70,30), expand = FALSE)+
#     ggtitle(paste0("Number of observations  "))
#   ggsave(paste0('../figures/presentation/ecologicaldata/numberobservations_facetwrap_',layer,'.png'), width = 8, height = 5, dpi=200)
#   ggsave(paste0('../figures/presentation/ecologicaldata/numberobservations_',layer,'.png'), width = 8, height = 5, dpi=200)
#   
# }else{
#   ggplot(data = world) +
#   geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
#   geom_sf(color = "darkgrey", fill = "lightgrey")+
#   geom_point(data=datos1, aes(Long,Lat, size = tCatch)) +
#   coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+
#   ggtitle(paste0("Total catch (MT) ",layer))
# ggsave(paste0('../figures/presentation/ecologicaldata/totalcatch_',layer,'.png'), width = 8, height = 5, dpi=200)
# }
# 
# ggplot(data = world) +
#   geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA)+
#   geom_sf(color = "darkgrey", fill = "lightgrey")+
#   geom_point(data=datos1, aes(Long,Lat, size = tCatch)) +
#   coord_sf(xlim = c(15,154), ylim = c(-60,30), expand = FALSE)+facet_wrap(~Species,ncol=3)+
#   ggtitle(paste0("Total catch (MT) of ",layer))
# ggsave(paste0('../figures/presentation/ecologicaldata/totalcatch_facetwrap_byspp_',layer,'.png'), width = 8, height = 5, dpi=200)
# 
# 
# 
# 
# ########## catches by fishery #########
# ######## catches by fishery ###########
# #### TOTAL CATCHES BY SPECIES - MJJJ
# library(dplyr)
# library(tidyr)
# datos_total_catches_fishery<-d1 %>% group_by(lon, lat, GEAR)  %>% summarise (tcatch = sum(RAISED_CATCHES_MT,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
# head(datos_total_catches_fishery)
# 
# datos_total_catches_fishery_wide<-datos_total_catches_fishery %>% select(-c(tcatch)) %>% spread(GEAR,freq)
# n <- nrow(datos_total_catches_fishery_wide)
# datos_total_catches_fishery_wide$region <- factor(1:n)
# datos_total_catches_fishery_wide<-as.data.frame(datos_total_catches_fishery_wide)
# head(datos_total_catches_fishery_wide)
# ## need to replace the NA by zeros if i want to plot the pie.
# datos_total_catches_fishery_wide[is.na(datos_total_catches_fishery_wide)] <- 0
# 
# library(scatterpie)
# 
# ## UNITY
# # neri<-c('COM','KAW','LOT')
# 
# fish<-unique(d1$GEAR)
# ggplot(data = world) + 
#   geom_sf(color = "darkgrey", fill = "lightgrey")+
#   coord_sf(xlim = c(15,154), ylim = c(-50,30), expand = FALSE) + 
#   geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), data=datos_total_catches_fishery_wide,
#                   cols=as.character(fish),  color="black", alpha=.8)+
#   geom_scatterpie_legend(2, x=130, y=20)+
#   coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)+
#   # annotate("text",label='B',x=20, y=30, size=12)+xlab('')+ylab('')+
#   ggtitle("Total catches of sharks by fishery")
# ggsave("../figures/totalcatches_",layer,"_byfishery_unity.png", width = 8, height = 5, dpi=200)
# 
# library(plyr)
# 
# 
# 
# 
# 
# detach("package:plyr", unload=TRUE) 
# 
# library(dplyr)
# ##### Examine the catch species compostion/location of the different fleets/type of fisheries
# datos3<-d2 %>% filter(year<2017 & year>=2003) %>% group_by(species, longitude,latitude, flag,gear,gear_group) %>% summarize (catch=median(V1,na.rm=TRUE))
# names(datos3)<-c("species","lon", "lat","fleet","fishery", "operation","catch")
# 
# ## combine RIN and PSS
# datos3[which(datos3$fishery=='RIN'),'fishery']<-'PSS'
# 
# #total catches by main fisheries for neritic tunas
# neritics<-NULL
# for(n in 1:length(ner)){neritics<-rbind(neritics,datos3[which(datos3$species==ner[n]),])}
# 
# Tcatchbyfishery<-neritics %>% group_by(fishery) %>% summarize (catch=sum(catch,na.rm=TRUE))%>% arrange(-catch)
# Tcatchbyfishery$fishery<-as.factor(Tcatchbyfishery$fishery)
# 
# library(ggplot2)
# g<-ggplot(Tcatchbyfishery[1:15,],aes(fishery))
# g+geom_bar(aes(weight=catch))+ylab('catch (MT)')
# 
# # par(mar=c(1,1,1,1,1))
# mp<-barplot(Tcatchbyfishery$catch[1:13],ylab='Catch (MT/10,000)',xlab='Fishery')
# axis(side = 1, at = mp, labels = Tcatchbyfishery$fishery[1:13], tcl = -0.2)
# 
# Tcatcheslonlat<-neritics %>% group_by(lon,lat) %>% summarize (catch=sum(catch,na.rm=TRUE))
# 
# # library("ggplot2")
# theme_set(theme_bw())
# library("sf")
# library("rnaturalearth")
# library("rnaturalearthdata")
# world <- ne_countries(scale = "medium", returnclass = "sf")
# # class(world)
# 
# ############################### catch distribution #############################################
# ggplot(data = world) + 
#   geom_sf(color = "darkgrey", fill = "lightgrey")+
#   coord_sf(xlim = c(15,154), ylim = c(-50,30), expand = FALSE)+ 
#   geom_point(data=Tcatcheslonlat, aes(lon,lat, size = catch)) +  
#   coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)+
#   # annotate("text",label='A',x=20, y=30, size=12)+
#   ggtitle(" Total catches")
# ggsave("../figures/species/totalcatches_neritic.png", width = 8, height = 5, dpi=200)
# 
# ggplot(data = world) + 
#   geom_sf(color = "darkgrey", fill = "lightgrey")+
#   coord_sf(xlim = c(15,154), ylim = c(-50,30), expand = FALSE)+ 
#   geom_point(data=Tcatcheslonlat, aes(lon,lat, size = catch)) +  
#   coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)+
#   facet_wrap(~species)
#   # annotate("text",label='A',x=20, y=30, size=12)+
#   ggtitle(" Total catches")
# ggsave("../figures/species/totalcatches_neritic.png", width = 8, height = 5, dpi=200)
# 
# #### TOTAL CATCHES BY SPECIES - MJJJ
# library(dplyr)
# library(tidyr)
# datos_total_catches_fishery<-neritics %>% group_by(lon, lat, species)  %>% summarise (tcatch = sum(catch,na.rm=TRUE)) %>% mutate(freq = tcatch / sum(tcatch), Totalcatch=sum(tcatch,na.rm=TRUE)) 
# head(datos_total_catches_fishery)
# 
# datos_total_catches_fishery_wide<-datos_total_catches_fishery %>% select(-c(tcatch)) %>% spread(species,freq)
# n <- nrow(datos_total_catches_fishery_wide)
# datos_total_catches_fishery_wide$region <- factor(1:n)
# datos_total_catches_fishery_wide<-as.data.frame(datos_total_catches_fishery_wide)
# head(datos_total_catches_fishery_wide)
# ## need to replace the NA by zeros if i want to plot the pie.
# datos_total_catches_fishery_wide[is.na(datos_total_catches_fishery_wide)] <- 0
# 
# library(scatterpie)
# 
# ## UNITY
# neri<-c('COM','KAW','LOT')
# ggplot(data = world) + 
#   geom_sf(color = "darkgrey", fill = "lightgrey")+
#   coord_sf(xlim = c(15,154), ylim = c(-50,30), expand = FALSE) + 
#   geom_scatterpie(aes(x=lon, y=lat, group=region, r=Totalcatch/Totalcatch*2), data=datos_total_catches_fishery_wide,
#                   cols=as.character(neri),  color="black", alpha=.8)+
#   geom_scatterpie_legend(2, x=130, y=20)+
#   coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)+
#   annotate("text",label='B',x=20, y=30, size=12)+xlab('')+ylab('')+
#   ggtitle(paste0("Total catches of ",layer," by species"))
# ggsave("/figures/totalcatches_mainspp_byfishery_unity.png", width = 8, height = 5, dpi=200)
# 
# 
# 
# #### BATHYMETRY ####
# # install.packages('marmap')
# library(marmap)
# # download for every grid cell of IOTC CA
# # resolution=300 (5 deg); resolution=60 (1 deg)
# bath<-getNOAA.bathy(20,150,-55,30,resolution=60,keep=T)
# 
# # bath[which(bath>0)]<-0
# bathy<-fortify(bath)
# bathy[which(bathy$z>0),'z']<-0
# 
# # autoplot(bath, geom=c("r", "c"), colour="white", size=0.5) +  scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="darkgreen")
# 
# library(rgdal)
# iotc_shp='/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/indian_ocean/iotc.shp'
# iotc <- readOGR(dsn = iotc_shp, stringsAsFactors = F)
# 
# ## read in neritic spp dist shapefiles (IUCN?)
# kawa_shp='/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/habitat_shapefiles/kawakawa/Euthynnus_affinis.shp'
# kawa <- readOGR(dsn = kawa_shp, stringsAsFactors = F)
# lot_shp='/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/habitat_shapefiles/longtail/Thunnus_tonggol.shp'
# lot <- readOGR(dsn = lot_shp, stringsAsFactors = F)
# com_shp='/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/habitat_shapefiles/Commerson/Scomberomorus_commerson.shp'
# com <- readOGR(dsn = com_shp, stringsAsFactors = F)
# 
# 
# map <- ggplot(data = world) +
#   # geom_tile(data = bathy, aes(x=x,y=y,fill=z))+
#   # geom_contour(data=bathy,aes(x=x,y=y,z=z),
#   #              breaks=c( -200, -1000),
#   #              colour="white", size=0.7
#   # )+
#   geom_polygon(data = iotc, aes(x = long, y = lat, group = group), colour = "red", fill = NA,size=0.7)+
#   geom_polygon(data = kawa, aes(x = long, y = lat, group = group), colour = "orange", fill = NA,size=0.7)+
#   geom_polygon(data = lot, aes(x = long, y = lat, group = group), colour = "green", fill = NA,size=0.7)+
#   geom_polygon(data = com, aes(x = long, y = lat, group = group), colour = "pink", fill = NA,size=0.7)+
#   geom_sf(color = "darkgrey", fill = "lightgrey")+
#   coord_sf(xlim = c(15,154), ylim = c(-60,35), expand = FALSE)
#   # annotate("text",label='B',x=20, y=30, size=12)
# map + theme_void()
# ggsave(paste0('../figures/species/neritic_habitat_shp_with_bathy.png'), width = 8, height = 5, dpi=200)
