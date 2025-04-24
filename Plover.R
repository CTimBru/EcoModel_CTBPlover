# Clears all memory in R
rm(list=ls())
require(data.table)
require(sf)
#install.packages("lwgeom") #Due to 'great Circles' in the Cumbria shapefile
require(raster)
require(dplyr)
require(terra)
require(stars)
require(randomForest)
require(plyr)
require(viridis)
require(ggplot2)

#Set working directory
wd <- "/Users/LilC/Desktop/WLAC/CTB-Plover"
setwd(wd)

#Set seed for reproducibility
set.seed(42)

#Years of Data
years_selected <- c(2019,2020,2021)

#intake species occurrence data from GBIF
#GBIF.org (22 April 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.s5cssx
input_species_data <- fread(input="Pluvialis_Apricaria.csv",sep="\t")
input_species_data <- input_species_data[input_species_data$year<2001]
input_species_data <- input_species_data[input_species_data$year>1969]

#Changes Lat|Long into a single variable, held as a column
# 4236 is WGS84: https://spatialreference.org/ref/epsg/4326/
input_species_points <- st_as_sf(input_species_data, coords=c("decimalLongitude","decimalLatitude"), crs = 4326)

#Relevant shapefile data for all gb administrative areas
gb_cumbria <- st_read("gb.shp")

#Get ID number for Cumbria in the data (GBCMA is the relevant ID.)
GBCMA_id <- which(gb_cumbria$name == 'Cumbria')

# Grab the shapefile poly for Cumbria. Use to crop & mask.
GBCMA_poly <- gb_cumbria$geometry[GBCMA_id] %>% st_transform(4326)
#xmin: -3.634185 ymin: 54.06379 xmax: -2.150328 ymax: 55.19083

#Generate a set of random background points, equal in number to actual occurrences.
num_occurrences <- nrow(input_species_data)

# From within boundary, create random not-occurrence points. 
background_points = sf::st_sample(GBCMA_poly, size=num_occurrences)
#Random forest prefers even occurances to not-occurances, other models vary.

#Convert single column coordinates to standard longitude/latitude columns
background_points <- sf::st_coordinates(background_points)

#Convert background points object to a data table
background_points <- as.data.table(background_points)

#Convert from longitude (X) and latitude (Y) columns to sf object.
background_points <- sf::st_as_sf(background_points,coords = c("X","Y"),  crs = 4326)

#Retain species occurrence data from within the SCB
GBCMA_species_points <- st_intersection(input_species_points, GBCMA_poly) 

#Historical Bioclim?
#map_layers <- list.files(path="wclim2019",pattern = "\\.tif$")
#for(map_layer in map_layers){
#  loadmap <- 
#}

#Crop & Mask cannot work with a sfc_MULTIPOLYGON
#Convert to a spatial for the crop & mask 
spd <- sf::as_Spatial(st_geometry(GBCMA_poly), IDs = as.character(1:nrow(GBCMA_poly[[1]][[1]][[1]])))

#Set to TRUE to rewrite all of the weather data
bool_PrepRasters <- FALSE

if(bool_PrepRasters == TRUE) {
  
  #One-time loop to crop&mask world climate data 2019-2021
  #CRU-TS 4.06 (Harris et al., 2020) downscaled with WorldClim 2.1 (Fick and Hijmans, 2017).
  #loops through the years
  #https://www.worldclim.org/data/monthlywth.html#
  for(selected_year in years_selected) {
    #Create a path for the relevant world climate data
    wclim_path <- paste("wclim",selected_year,sep="")
    #Create a list of files under that path 
    wclim_list <- list.files(path=wclim_path,pattern = "\\.tif$")
    #Create a path for reduced files
    GBCMAclim_path <- paste("GBCMAclim",selected_year,sep="")
    for(wclim in wclim_list){
      #Create a new filename
      GBCMA_name <- paste("GBCMA_",wclim,sep="")
      #Load raster for wclim data
      tempmap <- raster(paste(wclim_path,"/",wclim,sep=""))
      #Crop & Raster wclim data
      tempmap <- crop(tempmap,spd)
      tempmap <- mask(tempmap,spd)
      #Save cropped&masked Cumbria wclim data
      writeRaster(tempmap,paste(GBCMAclim_path,"/",GBCMA_name,sep=""),overwrite=TRUE)
    }
  }
  
  #crop wclim historical bio data
  #Fick, S.E. and R.J. Hijmans, 2017. WorldClim 2: new 1km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37 (12): 4302-4315.
  #https://www.worldclim.org/data/worldclim21.html
  bioclim_list <- list.files(path="wclim_bio1970-2000",pattern = "\\.tif$")
  for (bioclim in bioclim_list){
    #Load Raster for bioclim
    tempmap <- raster(paste("wclim_bio1970-2000/",bioclim,sep=""))
    #Crop & Raster bioclim data
    tempmap <- crop(tempmap,spd)
    tempmap <- mask(tempmap,spd)
    GBCMA_name <- paste("GBCMA_",bioclim,sep="")
    writeRaster(tempmap,paste("GBCMAclim_bio1970-2000","/",GBCMA_name,sep=""),overwrite=TRUE)
  }
  
  #crop HadGEM3 ssp885 2061-2080 data
  #
  #https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html
  hadgem3_list <- list.files(path="HadGEM3_SSP885_2061-2080",pattern = "\\.tif$")
  for (hadgem3 in hadgem3_list){
    #Load Raster Stacks for hadgem3
    tempmaps <- stack(paste("HadGEM3_SSP885_2061-2080/",hadgem3,sep=""))
    #Split the stacks, and save each layer
    for(tempmap in unstack(tempmaps)) {
      tempmap <- crop(tempmap,spd)
      tempmap <- mask(tempmap,spd)
      GBCMA_name <- paste("GBCMA_",names(tempmap),"_",hadgem3,sep="")
      writeRaster(tempmap,paste("GBCMAHadGEM3_SSP885_2061-2080","/",GBCMA_name,sep=""),overwrite=TRUE)
    }
  }
}

#https://naturalengland-defra.opendata.arcgis.com/datasets/Defra::peaty-soils-location-england/about
#BGS, Cranfield University (NSRI) and OS must be acknowledged
#Grab peat shapefile
gb_peat <- st_read("Peaty_Soils_Location_(England)___BGS_&_NSRI.shp") %>% st_transform(4326)
#CHECK THIS HACK WITH LEVI
sf_use_s2(FALSE)
#Intersect the shapefile and the Cumbria poly
GBCMA_peat <- st_intersection(gb_peat,GBCMA_poly)
#Result looks good, turning s2 back on!
sf_use_s2(TRUE)

# Dropping GBCMA_peat data except the depth|abundance categories
GBCMA_peat <- GBCMA_peat[,names(GBCMA_peat)%in% c("PCLASSDESC","geometry")]

#Empty container for the variables
PloverClim_variables <- c()

# First model Will use the wclim 1970-2000 bioclim data to predict occurrence of Golden Plover
#Set to TRUE to run the first model
bool_modelOne <- TRUE

if(bool_modelOne == TRUE) {
  #Creata a container for our rasters
  PloverClim_variables <- c()
  #Select the data source(s) for this model
  data_sources <- c("GBCMAclim_bio1970-2000")
  for(data_source in data_sources){
    #List the tif files in the chosen directory
    data_layers <- list.files(path=data_source,pattern = "\\.tif$")
    #stack together all of the tifs into the variables container
    PloverClim_variables <- stack(paste(data_source,"/",data_layers,sep=""))
    #name all of the layers in accordance to the filenames
    names(PloverClim_variables) <- data_layers
  }
  
  # We are going to convert all of the three peat layers into their own separate raster layers
  
  #Get the unique categories of peat
  unique_peat_cats <- unique(GBCMA_peat$PCLASSDESC)
  #Creat an empty container for our rasters
  peat_layers <- c()
  
  i=1
  for(peat_cat in unique_peat_cats){
    #Rasterize the peat shapefile, for one particular category of peat.
    temp <- st_rasterize(GBCMA_peat[GBCMA_peat$PCLASSDESC==peat_cat,] %>% dplyr::select(PCLASSDESC, geometry))
    
    #Rast & Raster the output
    temp <- rast(temp)
    temp <- raster(temp)
    
    #Crop & Mask to Cumbria
    temp <- crop(temp, spd)
    temp <- mask(temp, spd)
    #Resample to match spd layer, this data is categorical so is ngb better than bilinear
    temp <- resample(temp, PloverClim_variables[[1]], method="ngb")
    #Places the raster in our peat raster container
    peat_layers[[i]] <- temp
    i=i+1
  }
  
  #The data frames later prefer no spaces in their columns, so eliminate the spaces from the peat categories
  unique_peat_cats <- gsub(" ",".",unique_peat_cats)
  
  #Stack the peat layers in one raster stack
  GBCMA_peat_layers <- stack(peat_layers)
  #Name the layers in the raster stack after the peat categories
  names(GBCMA_peat_layers) <- unique_peat_cats
  #Add the peat layers onto of the existing bioClim data
  PloverClim_variables <- stack(PloverClim_variables,GBCMA_peat_layers)
  
  #extract the data from our raster layers for the observation points
  GBCMA_extracted <-raster::extract(PloverClim_variables,GBCMA_species_points)
  #copy over the column names to our extracted observations & associated variables
  colnames(GBCMA_extracted) <- names(PloverClim_variables)
  #Create a dataframe from the extracted data
  GBCMA_extracted <- as.data.frame(GBCMA_extracted)
  
  #Factorize the peat data. Current data was just IDs of shapes from the underlaying dataset, replace with 1 (peat present), & 0 (peat absent)
  for(peat_cat in unique_peat_cats){
    GBCMA_extracted[[paste(peat_cat)]][!is.na(GBCMA_extracted[[paste(peat_cat)]])] <- 1
    GBCMA_extracted[[paste(peat_cat)]][is.na(GBCMA_extracted[[paste(peat_cat)]])] <- 0
    GBCMA_extracted[[paste(peat_cat)]] <- as.factor(GBCMA_extracted[[paste(peat_cat)]])
  }
  
  #Remove any extracted point which has any remaining N/As
  GBCMA_extracted <- GBCMA_extracted[complete.cases(GBCMA_extracted ),]
  
  #Create a new column called prescence, set it to one (present)
  GBCMA_extracted$presence <- 1
  #Make the prescence column into a factor
  GBCMA_extracted$presence <- as.factor(GBCMA_extracted$presence)
  
  #COunt the number of occurances, within time-frame and having sufficient data
  num_occurrences <- nrow(GBCMA_extracted)
  
  #Create background points with variable data
  Background_extracted <-raster::extract(PloverClim_variables,background_points)
  colnames(Background_extracted) <- names(PloverClim_variables)
  Background_extracted <- as.data.frame(Background_extracted)
  
  #Again factorize peat data.
  for(peat_cat in unique_peat_cats){
    Background_extracted[[paste(peat_cat)]][!is.na(Background_extracted[[paste(peat_cat)]])] <- 1
    Background_extracted[[paste(peat_cat)]][is.na(Background_extracted[[paste(peat_cat)]])] <- 0
    Background_extracted[[paste(peat_cat)]] <- as.factor(Background_extracted[[paste(peat_cat)]])
  }
  
  #Drop any remaining N/As
  Background_extracted <- Background_extracted[complete.cases(Background_extracted ),]
  
  #Create a new column called prescence, set it to zero (not present)
  Background_extracted$presence <- 0
  #Make the prescence column into a factor
  Background_extracted$presence <- as.factor(Background_extracted$presence)
  
  #Empty containers for our model data
  raster_predict_list <- c()
  importance_list <- c()
  accuracy_list <- c()
  partial_plot_list <- c()
  
  j <- 1
  i <- 1
  i_max <- 50
  for(i in 1:i_max){
    #Get a subset of 80% of present points, and an equal number of absent points
    subset_extracted <- rbind(GBCMA_extracted[sample(nrow(GBCMA_extracted),0.8*nrow(GBCMA_extracted)),],Background_extracted[sample(nrow(Background_extracted),0.8*nrow(GBCMA_extracted)),])
    
    #Run a random forest model over this data subset.
    rf1 <- suppressWarnings(tuneRF(x=subset_extracted[,!(colnames(subset_extracted) %in% "presence")],y=subset_extracted$presence,stepFactor=1,plot=FALSE,doBest=TRUE))
    
    #Make a prediction raster from the random forest model and store it in a list.
    raster_predict_list[[i]] <- dismo::predict(PloverClim_variables,rf1,progress='text')
    
    #ERROR ON LINE 272 plot is all N/A rather than 0,1
    
    #Plot predicted raster
    plot(raster_predict_list[[i]])
    
    #Store relative importance of variable outputs as a temporary data frame.
    tmp <- as.data.frame(rf1$importance)
    #Set one column to store the variable names from the row names.
    tmp$VariableName <- rownames(tmp)
    #Store this importance data frame in the importance list.
    importance_list[[i]] <- tmp
    
    #Calculate the true skill statistic TSS to evaluate model accuracy.
    sensitivity <- rf1$confusion[[1]] / (rf1$confusion[[1]]+rf1$confusion[[2]])
    specificity <- rf1$confusion[[4]] / (rf1$confusion[[4]]+rf1$confusion[[3]])
    TSS <- sensitivity+specificity-1
    #Store TSS results
    accuracy_list[i] <- TSS
    
    #Currently Broken
    #Loop through each environmental variable and store the partial response outputs in a temporary data frame.
    #for(PloverClim_variable in PloverClim_variables){
    #  #Store partial plot chart data in a temporary data frame.
    #  tmp <- as.data.frame(partialPlot(rf1,subset_extracted[,!(colnames(subset_extracted) %in% "presence")],x.var=c(PloverClim_variable),plot=F))
    #  #Transform logistic probabilities to regular probabilities.
    #  tmp$y <- exp(tmp$y) / (1+exp(tmp$y))
    #  #Rename probability column
    #  colnames(tmp) <- c(PloverClim_variable,"Detection_Probability")
    #  #Store partial plot data in a list of data frames.
    #  partial_plot_list[[j]] <- tmp
    #  j <- j+1
    #}
    print(paste(i,Sys.time()))
  }
  
  importance_total <- rbind.fill(importance_list)
  importance_total <- aggregate(x=importance_total$MeanDecreaseGini,by = list(importance_total$VariableName),FUN = mean)
  colnames(importance_total) <- c("VariableName","Importance")
}



# Second model Will use the wclim 2019,2020, & 2021 data to predict occurrence of Golden Plover
bool_modelTwo <- TRUE
