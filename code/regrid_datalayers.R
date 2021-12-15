## iotc_regrid.R
## AE Nieblas
## 13/7/2019

## DESCRIPTION : Uses the 5x5 IOTC grid shapefile to regrid the raised catch IOTC data. Note: this function assigns the corresponding 5x5 gridcode to
# each data entry, and is appropriate for use in rescaling any resolution higher than 5 degree, but not lower.
# Inputs : iotc_5deg_shp: char, the path to the 5x5 IOTC grid shapefile, available: https://www.iotc.org/documents/TCAC04/shp/5x5
#          iotc_df : char, the path to the IOTC dataframe. 
# useage : 
# home<-'/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/'
# iotc_1_to_5_regrid(iotc_5deg_shp=paste0(home,'/IOTC-2018-TCAC04-DATA04_-_5x5_grids_shapefile/IOTC_TCAC_GRID_5x5.shp'),
#                    iotc_df=paste0(home,'EAFM_data_v1.csv'))

data_regrid<-function(iotc_5deg_shp,df){
  
  require(rgdal)
  require(ggplot2)
  library(rgdal) 
  library(sp) 
  library(tidyr)
  library(maps)
  
  # df<-BETlonlat
  ## ggplot the IOTC 5x5 grid
  shp <- readOGR(dsn = iotc_5deg_shp, stringsAsFactors = F)
  
  iotc_shp='/home/ae/Documents/PERSONAL/COOOL/projects/IOTC_2019_ecoregions/data/indian_ocean/iotc.shp'
  iotc_CA <- readOGR(dsn = iotc_shp, stringsAsFactors = F)
  
  inside.grid <- !is.na(over(iotc_CA, as(shp, "SpatialPolygons")))
  shp_CA<-as.data.frame(shp[inside.grid,])
  # shp_out<-as.data.frame(shp[!inside.grid,])
  # reassigns all data points to the 5x5 grid
  # shp_CA<-separate(shp_CA, grid, c("size", "quadrant", "lat", "lon"), sep = c(1, 2, 4)) 
  
  # 
  # iotc<-separate(shp_CA,CODE,c("size", "quadrant", "lat", "lon"), sep = c(1, 2, 4))
  # iotc$lat<-as.numeric(iotc$lat)
  # iotc$lon<-as.numeric(iotc$lon)
  # 
  # iotc$lat[iotc$quadrant==1] <- iotc$lat[iotc$quadrant==1]+2.5
  # iotc$lon[iotc$quadrant==1] <- iotc$lon[iotc$quadrant==1]+2.5
  # 
  # #Q2
  # iotc$lat[iotc$quadrant==2] <- (iotc$lat[iotc$quadrant==2]+ 2.5)*-1
  # iotc$lon[iotc$quadrant==2] <- iotc$lon[iotc$quadrant==2]+2.5
  
  # library(reshape2)
  # xx<-acast(df, lon~lat, value.var="Tcatch")
  
  # iotc_df<-iotc[,c('lon','lat')]
  # iotc_df$Tcatch<-NA
  # 
  # df2<-as.data.frame(df)
  # tt<-merge(iotc_df,df2, by=c('lon','lat'))
  
  
  # map <- ggplot() + geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "grey", fill = NA)
  # map + theme_void()

  
  # convert to data.frame>spatialpointsdataframe
  coordinates(df)<- c("lon","lat")
  
  # assign same projection as IOTC shapefile
  proj4string(df) <- proj4string(shp)
  
  inside.grid <- !is.na(over(df,as(shp, "SpatialPolygons")))
  
  df2<-NULL
  df2<-df[inside.grid,]
  df2$CODE <- over(df[inside.grid,], shp)$CODE  
  
  df2<-as.data.frame(df2)
  
  shp_CA$value<-NA
  shp_CA$NAME_EN<-NULL
  
  df3<-df2[,c('CODE','value')]
  
  # combo_df<-rbind(df2[,c('CODE','Tcatch')],shp_CA[,c('CODE','Tcatch')])
  
  # length(which(duplicated(combo_df$CODE))==TRUE)
  # shp_CA$uno<-shp_CA$CODE
  
  combo<-left_join(shp_CA[1],df3,by='CODE')
  
  df4<-separate(combo, CODE, c("size", "quadrant", "lat", "lon"), sep = c(1, 2, 4)) 
  
  
  df4$lat<-as.numeric(df4$lat)
  df4$lon<-as.numeric(df4$lon)
  
  ## centers data in the grid
  #Q1
  df4$lat[df4$quadrant==1] <- df4$lat[df4$quadrant==1]+2.5
  df4$lon[df4$quadrant==1] <- df4$lon[df4$quadrant==1]+2.5
  
  #Q2
  df4$lat[df4$quadrant==2] <- (df4$lat[df4$quadrant==2]+ 2.5)*-1
  df4$lon[df4$quadrant==2] <- df4$lon[df4$quadrant==2]+2.5
  
 return(df4) 
}
  
  
  
  
  
  
  
#   
#   dfout<-df2[!inside.grid,]
#   
#   if(dim(dfout)[1]>1){
#     catch_coords<-as.data.frame(unique(dfout@coords))
#     
#     dfout<-as.data.frame(dfout)
#     
#     # sum(dfout$RAISED_CATCHES_MT)+sum(df3$RAISED_CATCHES_MT)
#     library(rgeos)
#     
#     coordinates(catch_coords)<-c('lon','lat')
#     proj4string(catch_coords)<-proj4string(shp)
#     
#     grid<-NULL
#     for (i in 1:length(catch_coords)) {
#       grid[i] <- shp$CODE[which.min(gDistance(catch_coords[i,], shp, byid=TRUE))]
#       
#       dfout[which(dfout$lat==catch_coords@coords[i,2] & dfout$lon==catch_coords@coords[i,1]),'newgrid']<-grid[i]
#     }
#     # sum(dfout$RAISED_CATCHES_MT)+sum(df3$RAISED_CATCHES_MT)
#     
#     # nearest<-as.data.frame(dfout)
#     # colnames(nearest)<-'GRID_CODE'
#     nearest_coords<-NULL
#     if(data==2017|data=='bathy'){ nearest_coords<-separate(dfout, newgrid, c("size", "quadrant", "lat", "lon"), sep = c(1, 2, 4)) }
#     # if(data==2016){df2<-separate(df, FISHING_GROUND_CODE, c("size", "quadrant", "lat", "lon"), sep = c(1, 2, 4)) }
#     # sum(nearest_coords$RAISED_CATCHES_MT)+sum(df3$RAISED_CATCHES_MT)
#     
#     nearest_coords$lat<-as.numeric(nearest_coords$lat)
#     nearest_coords$lon<-as.numeric(nearest_coords$lon)
#     
#     nearest_coords$lat[nearest_coords$quadrant==1] <- nearest_coords$lat[nearest_coords$quadrant==1]+2.5
#     nearest_coords$lon[nearest_coords$quadrant==1] <- nearest_coords$lon[nearest_coords$quadrant==1]+2.5
#     
#     #Q2
#     nearest_coords$lat[nearest_coords$quadrant==2] <- (nearest_coords$lat[nearest_coords$quadrant==2]+ 2.5)*-1
#     nearest_coords$lon[nearest_coords$quadrant==2] <- nearest_coords$lon[nearest_coords$quadrant==2]+2.5
#     
#     nearest_coords$grid<-NULL
#     
#     # sum(nearest_coords$RAISED_CATCHES_MT)+sum(df3$RAISED_CATCHES_MT)
#     
#     plot(catch_coords, type="n")
#     map("world", add=TRUE)
#     plot(shp, border="green", add=TRUE)
#     legend("topright", cex=0.85,
#            c("Point not in grid", "Reallocated to grid", "Grid boundary"),
#            pch=c(3,16, NA), lty=c(NA, NA, 1),
#            col=c("black","grey", "green"), bty="n")
#     title(expression(paste("Observations with respect to IOTC grid")))
#     
#     points(nearest_coords$lon,nearest_coords$lat, pch=16, col="gray")
#     # points(df2[inside.grid, ], pch=16, col="red")
#     
#     df3<-rbind(df3,nearest_coords)
#   }
#   # sum(df3$RAISED_CATCHES_MT)
#   
#   if(data==2017){
#     # map <- ggplot() + geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "grey", fill = NA)+
#     #   geom_point(data=df3, aes(lon,lat, size = RAISED_CATCHES_MT))
#     # map + theme_void()
#   }
#   
#   if(data==2016){
#     
#     # library("ggplot2")
#     # theme_set(theme_bw())
#     # library("sf")
#     # library("rnaturalearth")
#     # library("rnaturalearthdata")
#     # world <- ne_countries(scale = "medium", returnclass = "sf")
#     # 
#     # map <- ggplot(data = world) + 
#     #   geom_point(data=df3, aes(lon,lat, size = CATCH))+
#     #   geom_sf(color = "black", fill = "lightgrey")+
#     #   geom_polygon(data = shp@polygons[[1]], aes(x = long, y = lat, group = group), colour = "green", fill = NA)+
#     #   coord_sf(xlim = c(15,154), ylim = c(-50,30), expand = FALSE)
#     # map + theme_void()
#   }
#   
#   return(df3)
# }
# 
# 
# 
# 
# 
