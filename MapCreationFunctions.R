
###HR####
# 1) HR_Metadata: Create a metadata with home-range of 95% and 50% KDE
# It will need "SubsetByTime" function as well to work
#Data - the movement data of your animal, with columns "X", "Y", "date", "TIME", "Flag" (if there is another ID type need to change the function)
#Subset - subset the data to every 60 minutes. if False the data will not be subset.

HR_Metadata <- function(Data,Subset = TRUE) {
  
  if (Subset == TRUE) {
    HRsData <- SubsetByTime(Data, SubsetTime = "60 mins")
  } else {
    HRsData <- Data
  }
  
  xyt<-subset(HRsData, select = c(X,Y))
  id<-subset(HRsData, select = Flag)
  locs1<-id
  coordinates(locs1)<-xyt
  
  UD3 <- kernelUD(locs1[,1], h = "href", grid = 500,same4all = FALSE, hlim = c(0.1, 2), kern = "bivnorm", extent = 1,boundary = NULL)
  homerange1 <- as.data.frame(getverticeshr(UD3, percent = 95))
  
  colnames(homerange1)[which(names(homerange1) == "area")] <- "area_95"
  
  homerange1$area_95 <- round(homerange1$area_95,2)
  
  homerange1_Core <- as.data.frame(getverticeshr(UD3, percent = 50))
  
  colnames(homerange1_Core)[which(names(homerange1_Core) == "area")] <- "area_50"
  
  homerange1_Core$area_50 <- round(homerange1_Core$area_50,2)
  
  #Add days:
  
  agg1 <- ddply(FullDataHR60_EveryHour,~Flag,summarise,Number_Of_Days=length(unique(date)))
  
  #Combine all by flag:
  
  AllHR <- merge(homerange1, homerange1_Core, by= "id")
  
  colnames(AllHR)[which(names(AllHR) == "id")] <- "Flag"
  
  AllHR <- merge(AllHR, agg1, by= "Flag")
  
  return(AllHR)
}

# SubsetByTime - this function is needed for HR_Metadata to work. It will subset the locations every X minutes
#Data - the dataframe of the ATLAS
#SubsetTime - how much time to subset (for example, "60 mins", "15 mins", ext.)

SubsetByTime <- function(Data, SubsetTime) {
  library(dplyr)
  library(data.table)
  library(lubridate)
  Flags <- unique(Data$Flag)
  
  Listone <- list()
  for (i in Flags) {
    LoopTag <- subset(Data, Flag == i)
    daysfilter <- as.factor(unique(LoopTag$date))
    listtwo <- list()
    for (ii in daysfilter) {
      LoopDay <- subset(LoopTag, date == ii)
      LoopDay <- setDT(LoopDay)[order(LoopDay)]
      
      output <- LoopDay[, .(DateTime = dateTime[1],date = date[1], Flag = Flag[1], X = X[1], Y = Y[1], TIME = TIME[1]) ,
                        by = .(Group = floor_date(dateTime, SubsetTime))]
      listtwo[[ii]] <- output
    }
    dataframetwo <- do.call(rbind.data.frame, listtwo)
    Listone[[i]] <- dataframetwo
    print(paste("Tag", i, "finish"))
  }
  NewLocations <- do.call(rbind.data.frame, Listone)
  NewLocations<-NewLocations[order(NewLocations$Flag,NewLocations$TIME),] #make sure data is sorted chronologically (per tag)
  
  return(NewLocations)
}

####Habitat Rasters with HR:####

# (1) #AddBufferToRaster: Create a buffer around a specific Raster layer
#RasterLayer - A raster layer uploaded
#BufferRadius - Radius around the layer to add to the new raster

AddBufferToRaster <- function(RasterLayer,BufferRadius) {
  
  library(terra)
  
  
  raster1.pts <- rasterToPoints(RasterLayer)
  raster1.df <- data.frame(raster1.pts)
  RasterRaRaster <- subset(raster1.df, raster1.df[,3] > 0 )
  RasterRaRaster <- rasterFromXYZ(RasterRaRaster) 
  crs(RasterRaRaster) <- itm
  z <- rast(RasterRaRaster)
  pe <- as.polygons(ext(z))
  pr <- as.polygons(z > -Inf)
  
  b <- raster::buffer(pr, width= BufferRadius) 
  r <- rast(b, ncols=ncol(RasterRaRaster), nrows=nrow(RasterRaRaster))
  xv <- rasterize(b, r)
  
  xvb <- raster(xv)
  crs(xvb) <- itm
  
  #See how much buffer was added on a map:
  
  leaflet() %>% addProviderTiles('Esri.WorldImagery')  %>%
    addRasterImage(xvb, colors = "Blue", opacity = 1, project=FALSE) %>%
    addRasterImage(RasterRaRaster, opacity = 0.9, project=FALSE) %>%
    addMeasure(
      position = "bottomleft",
      primaryLengthUnit = "meters",
      primaryAreaUnit = "sqmeters",
      activeColor = "#3D535D",
      completedColor = "#7D4479") %>%
    addScaleBar(position = "bottomleft", options = scaleBarOptions(maxWidth = 200,imperial = FALSE))
  
  return(xvb)
  
}

# (2) Adding new layer to a raster:
#NewLayer: the new tif file to add to the raster
#NumberName: Unique number to give to the new raster.
#Full_Raster: The original raster layer
#Type = the type of layer added, can be "Tif" for a tif file or "shp" for a shape file
#This function need "AddRasterTif" and "AddRasterShapeFile" to work

AddNewRaster <- function(NewLayer,NumberName,Full_Raster,Type = "Tif") {
  if (Type == "Tif") {
    AddRasterTif(NewLayer,NumberName,Full_Raster)
  } if (Type == "Shp") {
    AddRasterShapeFile(NewLayer,NumberName,Full_Raster)
  }
}

#Adding a Tif file:
#TifName: the new tif file to add to the raster
#NumberName: Unique number to give to the new raster.
#Full_Raster: The original raster layer

AddRasterTif <- function(TifName,NumberName,Full_Raster) {
  
  RasterRa=raster(TifName)
  RasterRa <- TifName
  #plot(GadashRa)
  
  raster1.pts <- rasterToPoints(RasterRa)
  raster1.df <- data.frame(raster1.pts)
  raster1.df[,3]
  
  RasterRaRaster <- subset(raster1.df, raster1.df[,3] > 0 )
  
  RasterRaRaster[,3] <- NumberName
  
  RasterRaRaster <- rasterFromXYZ(RasterRaRaster) 
  # RasterRaRaster@legend@colortable <- ColorName
  # RasterRaRaster@legend@names <- unique(TifName)
  
  crs(RasterRaRaster) <- itm
  
  #See the new raster in leaflet:
  
  Full_Raster <- merge(RasterRaRaster,Full_Raster)
  
  m <- leaflet() %>% addProviderTiles('Esri.WorldImagery')  %>%
    addRasterImage(Full_Raster, opacity = 0.9, project=FALSE)
  print(m)
  return(Full_Raster)
}

#Adding a shp file:

AddRasterShapeFile <- function(ShpName,NumberName,Full_Raster) {
  
  FirstPoly <- spTransform(ShpName,
                           crs(itm))
  r <- raster(resolution = 3)
  
  extent(r) <- extent(FirstPoly)
  res(r) <- 3
  FirstPolyB<- rasterize(FirstPoly, r,fun='first')
  values(FirstPolyB)[values(FirstPolyB) > 0] = NumberName
  
  
  #RasterRa=raster(ShpName)
  RasterRa <- FirstPolyB
  #plot(GadashRa)
  
  raster1.pts <- rasterToPoints(RasterRa)
  raster1.df <- data.frame(raster1.pts)
  raster1.df[,3]
  
  RasterRaRaster <- subset(raster1.df, raster1.df[,3] > 0 )
  
  RasterRaRaster[,3] <- NumberName
  
  RasterRaRaster <- rasterFromXYZ(RasterRaRaster) 
  #RasterRaRaster@legend@colortable <- ColorName
  #RasterRaRaster@legend@names <- unique(ShpName)
  
  crs(RasterRaRaster) <- itm
  
  #See in leaflet:
  
  origin(RasterRaRaster) <- origin(Full_Raster)
  
  Full_Raster <- merge(RasterRaRaster,Full_Raster)
  
  m <- leaflet() %>% addProviderTiles('Esri.WorldImagery')  %>%
    addRasterImage(Full_Raster, opacity = 0.9, project=FALSE)
  print(m)
  return(Full_Raster)
} 

# (3) Create a metadata combining HR and habitats from raster layer:

# Data = The movement dataframe
# RasName = The raster layer to create a metadata from
# Subset = if the data need subseting write "TRUE"

MetadataHabitatHR <- function(Data,RasName, Subset = FALSE) {
    
  if (Subset == TRUE) {
    Data <- SubsetByTime(Data, SubsetTime = "60 mins")
  }
  
  Harod_Habitats2=RasName
  
  xyt<-subset(Data, select = c(X,Y))
  id<-subset(Data, select = Flag)
  locs1<-id
  coordinates(locs1)<-xyt
  
  
  UD3 <- kernelUD(locs1[,1], h = "href", grid = 500,same4all = FALSE, hlim = c(0.1, 2), kern = "bivnorm", extent = 1,boundary = NULL)
  homerange1 <- getverticeshr(UD3, percent = 95)
  homerange1_Core <- getverticeshr(UD3, percent = 50)
  
  homerange2 <- homerange1
  homerange1_Core2 <- homerange1_Core
  
  proj4string(homerange2)<-CRS(itm)
  proj4string(homerange1_Core2)<-CRS(itm)
  
  Raster2 <- projectRaster(Harod_Habitats2, crs=itm)
  #Raster1@legend
  Combined2 <- raster::intersect(Raster2,homerange2)
  r3 <- mask(Combined2, homerange2) 
  plot(r3)
  
  FlagID <- unique(homerange2@data$id)
  BuiltAreas <- list()
  BuiltAreasFull <- list()
  print("Making Metadata")
  # i <- 4
  # Loop for every tag chosen:
  for (i in FlagID) {
    #Taking the data from the home-range:
    FullData <- subset(Data,Flag == i)
    DaysCalculated <- length(unique(FullData$date))
    
    animalnum95 <- homerange2[i,]
    animalnum50 <- homerange1_Core2[i,]
    
    #Adding area for later:
    size95 <- animalnum95$area
    size50 <- animalnum50$area
    #Combining the raster with the data:
    comination1 <- raster::intersect(Raster2,animalnum95)
    comination2 <- raster::intersect(Raster2,animalnum50)
    #Removing unwanted raster locations:
    r4 <- mask(comination1, animalnum95)
    r5 <- mask(comination2, animalnum50)
    # plot(r5)
    
    #Making a new dataframe from 95 and 50 home-range:
    tabtab <- table(round(r4@data@values,0))
    tabtab2 <- table(round(r5@data@values,0))
    
    tab2 <- data.frame(Color = names(tabtab), Count95 = as.integer(tabtab))
    
    tab2 <- tab2[1:9,]
    tab2$Color <- 3:11
    tab2$Count95 <- NA
    
    
    #3 = Main roads
    #4 = dirt road
    #5 = Fields 1
    #6 = Fields 2
    #7 = Fields 3
    #8 = mataim
    #9 = Natural
    #10 = WaterPonds
    #11 = Urban
    
    tab4 <- data.frame(Color = as.integer(names(tabtab)), Count95 = as.integer(tabtab))
    #Adding the habitats. Those lines need to be changed depending what your raster value is:
    tab2$Habitat[tab2$Color == "3"] <- "Main roads"
    tab2$Habitat[tab2$Color == "4"] <- "dirt road"
    tab2$Habitat[tab2$Color == "5"] <- "Fields 1"
    tab2$Habitat[tab2$Color == "6"] <- "Fields 2"
    tab2$Habitat[tab2$Color == "7"] <- "Fields 3"
    tab2$Habitat[tab2$Color == "8"] <- "mataim"
    tab2$Habitat[tab2$Color == "9"] <- "Natural"
    tab2$Habitat[tab2$Color == "10"] <- "WaterPonds"
    tab2$Habitat[tab2$Color == "11"] <- "Urban"
    tab2$Flag <- i
    
    
    tab2$Count95[match(tab4$Color,tab2$Color)] <- tab4$Count95
    
    tab2$Percent95 <- round(tab2$Count95/sum(tab2$Count95, na.rm=TRUE)*100,digits = 2)
    
    tab3 <- tab2[ -c(1) ]
    
    tab4 <- dcast(tab3, Flag~Habitat, value.var='Percent95')
    tab4$HR <- 95
    tab4$Size <- round(size95,digits = 2)
    
    #habitat 50:
    
    tab2b <- data.frame(Color = names(tabtab2), Count50 = as.integer(tabtab2))
    
    tab2b <- tab2b[1:9,]
    tab2b$Color <- 3:11
    tab2b$Count50 <- NA
    
    tab4b <- data.frame(Color = as.integer(names(tabtab2)), Count50 = as.integer(tabtab2))
    
    
    tab2b$Habitat[tab2b$Color == "3"] <- "Main roads"
    tab2b$Habitat[tab2b$Color == "4"] <- "dirt road"
    tab2b$Habitat[tab2b$Color == "5"] <- "Fields 1"
    tab2b$Habitat[tab2b$Color == "6"] <- "Fields 2"
    tab2b$Habitat[tab2b$Color == "7"] <- "Fields 3"
    tab2b$Habitat[tab2b$Color == "8"] <- "mataim"
    tab2b$Habitat[tab2b$Color == "9"] <- "Natural"
    tab2b$Habitat[tab2b$Color == "10"] <- "WaterPonds"
    tab2b$Habitat[tab2b$Color == "11"] <- "Urban"
    tab2b$Flag <- i
    
    tab2b$Count50[match(tab4b$Color,tab2b$Color)] <- tab4b$Count50
    sum(tab2b$Count50, na.rm=TRUE)
    
    tab2b$Percent50 <- round(tab2b$Count50/sum(tab2b$Count50, na.rm=TRUE)*100,digits = 2)
    
    tab3b <- tab2b[ -c(1) ]
    
    tab4b <- dcast(tab3b, Flag~Habitat, value.var='Percent50')
    tab4b$HR <- 50
    tab4b$Size <- round(size50, digits = 2)
    
    #Combine tab3 with tab3b, and tab4 + tab4b:
    
    tab3_new <- merge(tab3,tab3b,by=c("Habitat","Flag"))
    tab4_new <- rbind(tab4,tab4b)
    tab4_new$Days <- DaysCalculated
    
    BuiltAreas[[i]] <- tab4_new
    BuiltAreasFull[[i]] <- tab3_new
    print(paste("Finished flag number", i))
  }
  
  NumberOfArea <- do.call(rbind.data.frame, BuiltAreas)
  NumberOfAllArea <- do.call(rbind.data.frame, BuiltAreasFull)
  
  # Arranging the data:
  NumberOfArea$Flag <- as.factor(NumberOfArea$Flag)
  NumberOfAllArea$Flag <- as.factor(NumberOfAllArea$Flag)
  NumberOfAllArea$Habitat <- factor(NumberOfAllArea$Habitat, levels = c("Main roads", "dirt road", "Fields 1","Fields 2","Fields 3","mataim","Natural","WaterPonds","Urban"))
  
    return(NumberOfArea) #Single row to 95 and 50
  
}

# (4) Combining leaflet map with raster and home range of each animal individually
# This function create for each animal a HR of 95 and 50% KDE, combine it with the raster layer and present
# each animal individually in turns. Press [ENTER] to move to the next animal
#data: the full dataframe of the chosen animals
#FlagName: What animal ID to present
#RasterHabitats: The raster with the habitats
#MetadataHR: the metadata file created in "MetadataHabitatHR" function

HR_Habitat_plotting <- function(data,FlagName,RasterHabitats,MetadataHR) {
  
  for (i in FlagName) {
    Flag1 <- subset(data, Flag == i)
    xyt<-subset(Flag1, select = c(X,Y))
    id<-subset(Flag1, select = Flag)
    locs1<-id
    coordinates(locs1)<-xyt
    
    
    UD3 <- kernelUD(locs1, h = "href", grid = 500,same4all = FALSE, hlim = c(0.1, 2), kern = "bivnorm", extent = 1,boundary = NULL)
    homerange1 <- getverticeshr(UD3, percent = 95)
    homerange1_Core <- getverticeshr(UD3, percent = 50)
    
    homerange2 <- homerange1
    homerange1_Core2 <- homerange1_Core
    
    proj4string(homerange2)<-CRS(itm)
    proj4string(homerange1_Core2)<-CRS(itm)
    
    #Raster2 <- projectRaster(RasterHabitats, crs=itm)
    #Raster1@legend
    comination1 <- raster::intersect(RasterHabitats,homerange2)
    comination2 <- raster::intersect(RasterHabitats,homerange1_Core2)
    #Removing unwanted raster locations:
    r4 <- mask(comination1, homerange2)
    r5 <- mask(comination2, homerange1_Core2)
    plot(r4)
    plot(r5,add = T)
    proj4string(r4)<-CRS(itm)
    proj4string(r5)<-CRS(itm)
    pal <- colorNumeric(c("black", "brown", "green","green2","green3","white","green4","blue","grey"), values(RasterHabitats),
                        na.color = "transparent")
    
    m <- leaflet() %>% addProviderTiles('Esri.WorldImagery')  %>%
      addRasterImage(r4,colors = pal, opacity = 0.5, project=FALSE) %>%
      addRasterImage(r5,colors = pal, opacity = 1, project=FALSE)
    
    print(m)
    
    FlagMetadata <- subset(MetadataHR, Flag == i)
    
    print(FlagMetadata)
    readline(prompt="Press [enter] to continue")
  }
}

