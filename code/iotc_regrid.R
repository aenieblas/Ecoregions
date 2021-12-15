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

iotc_1_to_5_regrid<-function(iotc_5deg_shp,iotc_df,data=2017){
  
  require(rgdal)
  require(ggplot2)
  library(rgdal) 
  library(sp) 
  library(tidyr)
  library(maps)
  
  
  ## ggplot the IOTC 5x5 grid
  shp  <- readOGR(dsn = iotc_5deg_shp, stringsAsFactors = F)
  
  # map <- ggplot() + geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "grey", fill = NA)
  # map + theme_void()
  
  ## load iotc data
  df   <-  read.csv(iotc_df)
  # sum(df$RAISED_CATCHES_MT)
  
  # separate grid code into its component parts
  if(data==2017){df2<-separate(df, GRID_CODE, c("size", "quadrant", "lat", "lon"), sep = c(1, 2, 4)) }
  if(data==2016){df2<-separate(df, FISHING_GROUND_CODE, c("size", "quadrant", "lat", "lon"), sep = c(1, 2, 4)) }
  if(data=='size'){df2<-separate(df, Grid, c("size", "quadrant", "lat", "lon"), sep = c(1, 2, 4)) }
  

  if(grepl('bath',data)==FALSE){
  df2$lat<-as.numeric(df2$lat)
  df2$lon<-as.numeric(df2$lon)
  
  df2$lat[df2$quadrant==1] <- df2$lat[df2$quadrant==1]+2.5
  df2$lon[df2$quadrant==1] <- df2$lon[df2$quadrant==1]+2.5
  
  #Q2
  df2$lat[df2$quadrant==2] <- (df2$lat[df2$quadrant==2]+ 2.5)*-1
  df2$lon[df2$quadrant==2] <- df2$lon[df2$quadrant==2]+2.5
  }
  
  if(data=='bathy1deg'){
    # df2<-fortify(df)
  df[which(df$V3>0),'V3']<-NA
  df2<-df
  names(df2)<-c('lon','lat','depth')
  }
  if(data=='bathy'){
    df[which(df$V3>0),'V3']<-0
    df2<-df
    names(df2)<-c('lon','lat','depth')
  }
  
  # convert to data.frame>spatialpointsdataframe
  coordinates(df2)<- c("lon","lat")
  
  # assign same projection as IOTC shapefile
  proj4string(df2) <- proj4string(shp)
  
  # sum(df2$RAISED_CATCHES_MT)
  
  # removes the data points outside the IOTC grid
  # if(remove_points==TRUE){
  # find where catch points lie within the IOTC grid
  inside.grid <- !is.na(over(df2, as(shp, "SpatialPolygons")))
  
  # proportion of points within the grid
  # mean(inside.grid) 
  
  # assign the 5 degree grid code for each data point
  df2$grid <- over(df2, shp)$CODE
  
  # sum(df2$RAISED_CATCHES_MT)
  
  ## plot of the points inside and outside the IOTC grid.
  # plot(coordinates(df2), type="n")
  # map("world", add=TRUE)
  # plot(shp, border="green", add=TRUE)
  # legend("topright", cex=0.85,
  #        c("Point in grid", "Point not in grid", "Grid boundary"),
  #        pch=c(16, 16, NA), lty=c(NA, NA, 1),
  #        col=c("red", "grey", "green"), bty="n")
  # title(expression(paste("Catch observations with respect to IOTC grid")))
  # 
  # points(df2[!inside.grid, ], pch=16, col="gray")
  # points(df2[inside.grid, ], pch=16, col="red")
  
  df3<-as.data.frame(df2[inside.grid,])
  # reassigns all data points to the 5x5 grid
  df3<-separate(df3, grid, c("size", "quadrant", "lat", "lon"), sep = c(1, 2, 4)) 
  
  df3$lat<-as.numeric(df3$lat)
  df3$lon<-as.numeric(df3$lon)
  
  ## centers data in the grid
  #Q1
  df3$lat[df3$quadrant==1] <- df3$lat[df3$quadrant==1]+2.5
  df3$lon[df3$quadrant==1] <- df3$lon[df3$quadrant==1]+2.5
  
  #Q2
  df3$lat[df3$quadrant==2] <- (df3$lat[df3$quadrant==2]+ 2.5)*-1
  df3$lon[df3$quadrant==2] <- df3$lon[df3$quadrant==2]+2.5
  
  
  dfout<-df2[!inside.grid,]
  
  if(dim(dfout)[1]>1){
    catch_coords<-as.data.frame(unique(dfout@coords))
    
    dfout<-as.data.frame(dfout)
    
    # sum(dfout$RAISED_CATCHES_MT)+sum(df3$RAISED_CATCHES_MT)
    library(rgeos)
    
    coordinates(catch_coords)<-c('lon','lat')
    proj4string(catch_coords)<-proj4string(shp)
    
    grid<-NULL
    for (i in 1:length(catch_coords)) {
      grid[i] <- shp$CODE[which.min(gDistance(catch_coords[i,], shp, byid=TRUE))]
      
      dfout[which(dfout$lat==catch_coords@coords[i,2] & dfout$lon==catch_coords@coords[i,1]),'newgrid']<-grid[i]
    }
    # sum(dfout$RAISED_CATCHES_MT)+sum(df3$RAISED_CATCHES_MT)
    
    # nearest<-as.data.frame(dfout)
    # colnames(nearest)<-'GRID_CODE'
    nearest_coords<-NULL
    if(data==2017|data=='bathy'){ nearest_coords<-separate(dfout, newgrid, c("size", "quadrant", "lat", "lon"), sep = c(1, 2, 4)) }
    # if(data==2016){df2<-separate(df, FISHING_GROUND_CODE, c("size", "quadrant", "lat", "lon"), sep = c(1, 2, 4)) }
    # sum(nearest_coords$RAISED_CATCHES_MT)+sum(df3$RAISED_CATCHES_MT)
    
    nearest_coords$lat<-as.numeric(nearest_coords$lat)
    nearest_coords$lon<-as.numeric(nearest_coords$lon)
    
    nearest_coords$lat[nearest_coords$quadrant==1] <- nearest_coords$lat[nearest_coords$quadrant==1]+2.5
    nearest_coords$lon[nearest_coords$quadrant==1] <- nearest_coords$lon[nearest_coords$quadrant==1]+2.5
    
    #Q2
    nearest_coords$lat[nearest_coords$quadrant==2] <- (nearest_coords$lat[nearest_coords$quadrant==2]+ 2.5)*-1
    nearest_coords$lon[nearest_coords$quadrant==2] <- nearest_coords$lon[nearest_coords$quadrant==2]+2.5
    
    nearest_coords$grid<-NULL
    
    # sum(nearest_coords$RAISED_CATCHES_MT)+sum(df3$RAISED_CATCHES_MT)
    
    # plot(catch_coords, type="n")
    # map("world", add=TRUE)
    # plot(shp, border="green", add=TRUE)
    # legend("topright", cex=0.85,
    #        c("Point not in grid", "Reallocated to grid", "Grid boundary"),
    #        pch=c(3,16, NA), lty=c(NA, NA, 1),
    #        col=c("black","grey", "green"), bty="n")
    # title(expression(paste("Observations with respect to IOTC grid")))
    # 
    # points(nearest_coords$lon,nearest_coords$lat, pch=16, col="gray")
    # points(df2[inside.grid, ], pch=16, col="red")
    
    df3<-rbind(df3,nearest_coords)
  }
  # sum(df3$RAISED_CATCHES_MT)
  
  if(data==2017){
    # map <- ggplot() + geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "grey", fill = NA)+
    #   geom_point(data=df3, aes(lon,lat, size = RAISED_CATCHES_MT))
    # map + theme_void()
  }
  
  if(data==2016){
    
    # library("ggplot2")
    # theme_set(theme_bw())
    # library("sf")
    # library("rnaturalearth")
    # library("rnaturalearthdata")
    # world <- ne_countries(scale = "medium", returnclass = "sf")
    # 
    # map <- ggplot(data = world) + 
    #   geom_point(data=df3, aes(lon,lat, size = CATCH))+
    #   geom_sf(color = "black", fill = "lightgrey")+
    #   geom_polygon(data = shp@polygons[[1]], aes(x = long, y = lat, group = group), colour = "green", fill = NA)+
    #   coord_sf(xlim = c(15,154), ylim = c(-50,30), expand = FALSE)
    # map + theme_void()
  }
  
  return(df3)
}





