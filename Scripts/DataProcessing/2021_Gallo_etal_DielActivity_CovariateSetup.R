##### Process Site-level Covariates at Species Specific Scale #####

# load packages and functions
# load packages and functions
source("./Scripts/2021_Gallo_etal_DielActivity_UtilityScripts.R")

# load dataset
data_full <- readRDS("./Data/2021_Gallo_etal_DielActivity_cleaned_photo_dataset_9_species.rds")

## Calculate the proportion of time available
# Day - Reference
data_full$avail_cat1 <-  (as.numeric(data_full$dusk) - as.numeric(data_full$sunriseEnd))/3600
# Dawn
data_full$avail_cat2 <-  (as.numeric(data_full$sunriseEnd) - as.numeric(data_full$dawn))/3600
# Dusk
data_full$avail_cat3 <-  (as.numeric(data_full$night) - as.numeric(data_full$dusk))/3600

#  night hour subtracted from midnight plus midnight to dawn
# then subtract the dark night chunk - 14400/60
data_full$avail_cat4 <-  (((86400 - as.numeric(data_full$night)) + as.numeric(data_full$dawn)) -
                            7200)/3600

# dark night is just 2 hours long
data_full$avail_cat5 <- 7200/3600

#check out mean avaliable time for manuscript
apply(data_full[,21:25], 2, mean)

## get daily average temperature

# need to set up your .Rprofile with a NOAA token (get your own)
  # https://www.ncdc.noaa.gov//cdo-web/token

options(noaakey = "Need your own token")
# make a table of all the noaa stations
stations <- ghcnd_stations()

# function to make a lat/lon dataframe for extracting weather data
# had to use Topeka, KS and Cedar Rapids, IA
city_latlon <-  data.frame(rbind(c("CHIL", 41.8781, -87.6298),
                                 c("ICIA",41.88, -91.72),
                                 c("MAWI", 43.0731, -89.4012),
                                 c("WIDE", 39.7447, -75.5484),
                                 c("ININ", 39.7684, -86.1581),
                                 c("AUTX", 30.2672, -97.7431),
                                 c("MAKS", 39.07, -95.63),
                                 c("FOCO", 40.5853, -105.0844),
                                 c("DECO", 39.7392, -104.9903),
                                 c("LBCA", 33.7701, -118.1937)))
colnames(city_latlon) <- c("id", "latitude", "longitude")

# find the stations avaliable - 50k was enough distance to grab a station for each city
city_stations <- meteo_nearby_stations(lat_lon_df = city_latlon,
                                       station_data = stations,
                                       radius = 50,
                                       year_min = 2017,
                                       var = "TAVG")

# just take the nearest
city_stations <- lapply(city_stations, function(x) {
  x[which(x$distance == min(x$distance)),]
})

# extract city and date
city_date <- data_full[,c("city","date")]
# only keep unique city and date combination
city_date <- unique(city_date)

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)

# loop through each row so that we can get a daily average temperature for each day
tavg <- foreach(i=1:nrow(city_date), .packages = c("rnoaa"), .verbose = TRUE) %dopar% {

  # the city for this data point
  station <- city_stations[[as.character(city_date[i,"city"])]]
  
  # get TAVG for each data point
  meteo_pull_monitors(station$id, date_min = city_date[i,"date"],
                      date_max = city_date[i,"date"], var = "TAVG")$tavg/10
                }
stopCluster(cl)

# change empty elements to NA
emptyToNA <- function(x){
  ifelse(length(x) == 0, NA, x)
}
# run function
tavg2 <- sapply(tavg, emptyToNA)
# add covariate to dataframe
city_date$tavg <- tavg2

# only March and Novemeber dates so hard coded the monthly average
# all but one were march dates
city_date[which(is.na(city_date$tavg))[-5],"tavg"] <- 4.06
# once March dates are fixed only the single November date is left
city_date[is.na(city_date$tavg),"tavg"] <- 6.17

# join with larger data and remove all the time data
# join data_full with temp data
data_full2 <- left_join(data_full[,-c(3:4,11:20)], city_date, c("city", "date")) 

## Extract spatial covariates
# extract just Lat/Lon - ALSO KEEP SITE NAME
sites <- data_full2[,c("city", "LocationName", "UTM_E", "UTM_N")]

# fix MAKS and FOCO reversed UTM's
sites_corrected <- sites
sites_corrected[which(sites_corrected[,"city"] == "MAKS"),"UTM_E"] <- 
  sites[which(sites[,"city"] == "MAKS"),"UTM_N"]
sites_corrected[which(sites_corrected[,"city"] == "MAKS"),"UTM_N"] <- 
  sites[which(sites[,"city"] == "MAKS"),"UTM_E"]
sites_corrected[which(sites_corrected[,"city"] == "FOCO"),"UTM_E"] <- 
  sites[which(sites[,"city"] == "FOCO"),"UTM_N"]
sites_corrected[which(sites_corrected[,"city"] == "FOCO"),"UTM_N"] <- 
  sites[which(sites[,"city"] == "FOCO"),"UTM_E"]

# remove duplicate site names to isolate just one record for each site
city_sites <- sites_corrected[!duplicated(sites_corrected[,"LocationName"]),]

# fix Austin sites
city_sites[which(city_sites$LocationName == "AUTX-C01-WBP"),"UTM_E"] <- -97.81902
city_sites[which(city_sites$LocationName == "AUTX-C01-WBP"),"UTM_N"] <- 30.31024
city_sites[which(city_sites$LocationName == "AUTX_CF"),"UTM_E"] <- -97.89156
city_sites[which(city_sites$LocationName == "AUTX_CF"),"UTM_N"] <- 30.336876

# load Landcover map
lc_map <- raster("~/Documents/GIS/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img")
# load impervious map
imp_map <- raster("~/Documents/GIS/nlcd_2011_impervious_2011_edition_2014_10_10/nlcd_2011_impervious_2011_edition_2014_10_10.img")

# Function that first extracts spatial data at a given scale around each site for each city
# Function then joins that information with each data point and scales values
scaleDataExtract <- function(buffer_size){
  # function to extract site level data
  sitecovsCitySites <- function(city, buff){
    
    # extract data one city at a time
    city_data <- city_sites[which(city_sites$city == city),]
    
    # fix Madison data
    if(city == "MAWI"){
      city_data$UTM_E <- city_data$UTM_E * -1
    }
    
    # run function for correct projection
    crs_ <- cityCRS(city)
    
    # create spatial points
    city_points <- st_as_sf(city_data, coords = c("UTM_E","UTM_N"),
                            crs = crs_)
    
    # run function that calls the correct population data shape file from Silvis Lab
    file <- stateFile(city)
    
    # load population data
    pop_data <- st_read("~/Documents/GIS/UWIN_GIS/state_pop10_shapefiles", 
                        layer=file, crs = 5070)
    
    # reproject city data to population data layer
    city_points_RP <- city_points %>% st_transform(st_crs(pop_data))
    
    # create buffer around sites
    city_points_buffered <- st_buffer(city_points_RP, dist=buff)
    
    # crop the data to speed up process
    pop_data_crop <- st_crop(pop_data, st_bbox(city_points_buffered))
    rm(pop_data)
    
    buffer_intersect <- st_intersection(pop_data_crop,
                                        city_points_buffered)
    # get sum
    # grab site names
    location_names <- unique(city_points_buffered$LocationName)
    
    # calculate the mean housing & population densities for each site
    hd <- rep(0, length(location_names))
    pd <- rep(0, length(location_names))
    for(i in 1:length(location_names)){
      hu <- sum(buffer_intersect[which(buffer_intersect$LocationName == 
                                         location_names[i]),"HU10"]$HU10)
      hd[i] <- hu/(st_area(city_points_buffered[i,])/1e6)
      pu <- sum(buffer_intersect[which(buffer_intersect$LocationName == 
                                         location_names[i]),"POP10"]$POP10)
      pd[i] <- pu/(st_area(city_points_buffered[i,])/1e6)
    }
    # remove shapefile to save memory
    rm(pop_data_crop)
    
    ## Extract mean impervious cover
    
    # reproject to match impervious raster
    city_points_RP <- st_transform(city_points, crs = st_crs(imp_map))
    # extract
    imp <- extract(imp_map, city_points_RP, fun=mean, buffer=buff, df=TRUE)
    
    ## Extract available habitat
    
    # reproject points to mach landcover raster
    city_points_RP <- st_transform(city_points, crs = st_crs(lc_map))
    
    # extract land cover data for each point, given buffer size
    # extract with sf objects makes a list, so use lapply to create a table of each category
    lc_extract <- extract(lc_map, city_points_RP, buffer = buff)
    
    summarizeAvaliableHabitat <- function(x){
      # create a table of proportions for each land cover class
      proportions <- prop.table(table(x))
      
      # summarize each site's data by proportion of each cover type
      # convert to data frame
      landcover <- data.frame(cover = names(proportions), percent = as.numeric(proportions))
      
      # sum across categories that we are considering habitat                       
      habitat <- sum(landcover[which(landcover$cover == 21),2], # developed open space
                     landcover[which(landcover$cover == 41),2], # forest
                     landcover[which(landcover$cover == 42),2], # forest
                     landcover[which(landcover$cover == 43),2], # forest
                     landcover[which(landcover$cover == 51),2], # shrub
                     landcover[which(landcover$cover == 52),2], # shrub
                     landcover[which(landcover$cover == 71),2], # herbaceous
                     landcover[which(landcover$cover == 72),2], # herbaceous
                     landcover[which(landcover$cover == 73),2], # herbaceous
                     landcover[which(landcover$cover == 74),2], # herbaceous
                     landcover[which(landcover$cover == 90),2], # wetland
                     landcover[which(landcover$cover == 95),2]) # wetland
      
      return(habitat)
      
    }
    
    habitat <- do.call(rbind, lapply(lc_extract, summarizeAvaliableHabitat))
    
    ## Extract mean NDVI
    
    # load LandSat rasters
    # Band 4 and Band 5 have already been isolated from the LandSat layer
    # load Band 4 raster
    b4 <- raster(paste0("./Data/Landsat/",city,"/B4.TIF"))
    # load Band 5 raster
    b5 <- raster(paste0("./Data/Landsat/",city,"/B5.TIF"))
    
    # calculate NDVI
    # (NIR - Red)/(NIR + Red)
    ndvi_raster <- (b5-b4)/(b5+b4)
    
    # reproject to match ndvi raster
    city_points_RP <- st_transform(city_points, crs = st_crs(ndvi_raster))
    # extract the mean NDVI around each point
    mean_ndvi <- extract(ndvi_raster, city_points_RP, fun=mean, na.rm = TRUE, 
                         buffer=buff, df=TRUE)
    
    # extract the proportion of the buffer that has an NDVI greater than 0 (vegation cover)
    # extract with sf objects makes a list, so use lapply to create a table of each category
    ndvi_extract <- extract(ndvi_raster, city_points_RP, buffer = buff)
    
    # calculate the proportion of a site that is covered in vegetation
    prop_ndvi_greater0 <- lapply(ndvi_extract, function (x){
      x[which(x > 0.2)] <- 1
      x[which(x <= 0.2)] <- 0
      sum(x, na.rm = TRUE)/length(x) })
    # collapse list
    ndvi_greater0 <- do.call(rbind, prop_ndvi_greater0)
    
    ## combine data
    city_site_covs <- data.frame(LocationName = city_data$LocationName, 
                                 city = city,
                                 avail_habitat = habitat,
                                 mean_imp = imp$nlcd_2011_impervious_2011_edition_2014_10_10,
                                 hd = hd,
                                 pd = pd,
                                 mean_ndvi = mean_ndvi$layer,
                                 prop_veg = ndvi_greater0)
    
    return(city_site_covs)
    
  }
  
  # run for each city
  # vector of cities
  cities <- as.character(unique(data_full2$city))
  
  # create empty list to fill
  sitecovs <- vector("list", length(cities))
  # loop through each city
  for(i in 1:length(sitecovs)){
    sitecovs[[i]] <- sitecovsCitySites(city = cities[i], buff = buffer_size)
  }
  
  # collapse the list into a data frame
  sitecovs_df <- do.call(rbind, sitecovs)
  
  # create urban index using housing density, impervious cover, and ndvi
  urb_pca <- prcomp(sitecovs_df[,c("mean_imp","hd","prop_veg")], scale = TRUE)
  
  # see how much variance is explained by each component
  expl.var <- round(urb_pca$sdev^2/sum(urb_pca$sdev^2)*100)
  # use first principle component
  sitecovs_df$urban_index <- urb_pca$x[,1]
  
  return(list(sitecovs = sitecovs_df,
              pca_results = urb_pca,
              pca_var = expl.var))
}

# extract the covariates at 3 scales
buffers <- c(1500, 1000, 500)
sitecovsBuffered <- vector("list", length(buffers))
for(i in 1:length(sitecovsBuffered)){
  sitecovsBuffered[[i]] <- scaleDataExtract(buffers[i])
}
# give each list a name
names(sitecovsBuffered) <- c("buffer1500", "buffer1000", "buffer500")

saveRDS(sitecovsBuffered, "2021_Gallo_etal_DielActivity_sitecovsBuffered.rds")

speciesCovs <- function(species){
  # combine the photo dataset for each species with the spatial covariates
  # species specific scales
  if(species == "bobcat" | species == "coyote"| species == "redfox"){
    df <- left_join(data_full2[which(data_full2$ShortName == species),], 
                    sitecovsBuffered[["buffer1500"]][["sitecovs"]][,-2], 
                    by="LocationName")
  }
  if(species == "raccoon" | species == "striped skunk"){
    df <- left_join(data_full2[which(data_full2$ShortName == species),], 
                    sitecovsBuffered[["buffer1000"]][["sitecovs"]][,-2], 
                    by="LocationName")
  }
  if(species == "e. cottontail" | species == "v. opossum"| species == "w. t. deer"){
    df <- left_join(data_full2[which(data_full2$ShortName == species),], 
                    sitecovsBuffered[["buffer500"]][["sitecovs"]][,-2], 
                    by="LocationName")
  }
  
  # scale univariate covariates for only the sites in which we detected the species
  df$habitat_scaled <- with(df, scale(avail_habitat))
  df$imp_scaled <- with(df, scale(mean_imp))
  df$hd_scaled <- with(df, scale(hd))
  df$pd_scaled <- with(df, scale(pd))
  df$ndvi_scaled <- with(df, scale(mean_ndvi))
  df$veg_scaled <- with(df, scale(prop_veg))
  df$urb_scaled <- with(df, scale(urban_index))
  df$tavg_scaled <- with(df, scale(tavg))

  return(df)
  
}

# vector of species
species_list <- c("bobcat", "coyote", "e. cottontail", "raccoon", "red fox", 
                  "striped skunk", "v. opossum", "w. t. deer")

# calculated at a species specific scale
data_list <- vector("list", length(species_list))
# run function to extract data at a species specific buffer
for(i in 1:length(data_list)){
  data_list[[i]] <- speciesCovs(species = species_list[i])
}

# name each list to stay organized
names(data_list) <- c("bobcat", "coyote", "e_cottontail", "raccoon", "red_fox", 
                      "striped_skunk", "v_opossum", "wt_deer")

lapply(data_list, function(x) { table(x$time_cat) })
  
saveRDS(data_list, "./Data/2021_Gallo_etal_DielActivity_data_list_8sp.rds")
