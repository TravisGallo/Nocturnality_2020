## Process and summarize photo data 

# load packages and functions
source("./Scripts/2021_Gallo_etal_DielActivity_UtilityScripts.R")

# list of files in the directory
city_data <- list.files("./Data/City_Files", "*.csv")

# loop through each city file and process data
# do in parallel and rbind at the end to make large data frame
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl) # register backend

data_full <- foreach(i=1:length(city_data), .packages = c("suncalc","lubridate","dplyr"), 
                     .combine = "rbind") %dopar% {
        
                      # read in city specific dataset               
                      data_read <- read.csv(paste0("./Data/City_Files/",city_data[i]), 
                                            stringsAsFactors = FALSE)
                      
                      # fix only the species we will use
                      # fix v. opossum
                      data_read[grep('posss?um', data_read$ShortName),"ShortName"] <- "v. opossum"
                      # fix w. t. deer
                      data_read[grep('deer?', data_read$ShortName),"ShortName"] <- "w. t. deer"
                      # fix striped skunk
                      data_read[grep('skunk?', data_read$ShortName),"ShortName"] <- "striped skunk"
                      # fix e. cottontail
                      data_read[which(data_read$ShortName == "rabbit"),"ShortName"] <- "e. cottontail"
        
                      # grab only species (including misspelled) we will work with
                      data <- data_read[which(data_read$ShortName == "bobcat" |
                                                data_read$ShortName == "coyote" |
                                                data_read$ShortName == "w. t. deer" |
                                                data_read$ShortName == "v. opossum" |
                                                data_read$ShortName == "striped skunk" |
                                                data_read$ShortName == "e. cottontail" |
                                                data_read$ShortName == "raccoon" |
                                                data_read$ShortName == "red fox" |
                                                data_read$ShortName == "gray fox" ),]
                                                
                      # convert time to a time object
                      # Central time cities
                      if(substr(city_data[i],1,4) == "ICIA" | substr(city_data[i],1,4) == "MAWI" |
                         substr(city_data[i],1,4) == "AUTX" | substr(city_data[i],1,4) == "MAKS" |
                         substr(city_data[i],1,4) == "CHIL"){
                        data$ImageDate <- strptime(data$ImageDate,format="%m/%d/%Y %H:%M", 
                                                   tz = "US/Central")
                      }
                      # Eastern time cities
                      if(substr(city_data[i],1,4) == "WIDE" | substr(city_data[i],1,4) == "ININ"){
                        data$ImageDate <- strptime(data$ImageDate,format="%m/%d/%Y %H:%M", 
                                                   tz = "US/Eastern")
                      }
                      # Mountain Standard time cities
                      if(substr(city_data[i],1,4) == "DECO" | substr(city_data[i],1,4) == "FOCO"){
                        data$ImageDate <- strptime(data$ImageDate,format="%m/%d/%Y %H:%M", 
                                                   tz = "US/Mountain")
                      }
                      # Pacific time cities
                      if(substr(city_data[i],1,4) == "LBCA"){
                        data$ImageDate <- strptime(data$ImageDate,format="%m/%d/%Y %H:%M", 
                                                   tz = "US/Pacific")
                      }
        
                      # remove any data points that have NA for ImageDate
                      # need this info for analysis so no good to us otherwise
                      data <- data[!is.na(data[,"ImageDate"]),]
        
                      # remove incorrect time stamps
                      data <- data[which(data$ImageDate > as.Date('2017-01-01') & 
                                           data$ImageDate < as.Date('2019-01-01')),]
        
                      # remove messed up WIDE sites
                      if(substr(city_data[i],1,4) == "WIDE"){
                      # remove the April-May season
                      data <- data[-which(data$ImageDate > as.Date("2018-04-01") & 
                                         data$ImageDate < as.Date('2018-05-31')),]
                      }
                      
                      # remove pictures that fall within 15 minutes of each other
                      # seperate by species and then by site
                      list_by_species <- split(data, data$ShortName)
                      species_list_by_site <- lapply(list_by_species, function(x) { 
                        split(x, x$LocationName) })
                      # keep photos that are greater than 15 minutes apart
                      photo_list <- species_list_by_site
                      for(spec in 1:length(list_by_species)){
                        for(site in 1:length(species_list_by_site[[spec]])){
                          to_keep <- c(0,which(diff(species_list_by_site[[spec]][[site]]$ImageDate) > 900)) + 1
                          photo_list[[spec]][[site]] <- species_list_by_site[[spec]][[site]][to_keep,]
                        }
                      }
                      
                      # bring data back together
                      collapse1 <- lapply(photo_list, function(x) do.call(rbind, unname(x)))
                      data <- do.call("rbind", unname(collapse1))
                
                      # extract just the date
                      data$date <- as.Date(data$ImageDate, "%m/%d/%Y")
                      # extract just the time
                      data$time <- hms(format(data$ImageDate, "%H:%M:%S"))
        
                      # get the different sun times based on dates
                      sun_times <- as.data.frame(matrix(NA, nrow=nrow(data), ncol = 10))
                      colnames(sun_times) <- c("dawn","sunrise","sunriseEnd","dusk",
                                               "solarNoon","sunsetStart","sunset",
                                               "night","nadir","nightEnd")
        
                      # function to decide the city and therefore the Lat Lon
                      cityLatLon <-  function(city_code){switch(city_code,
                                                                "CHIL" = c(41.8781, -87.6298),
                                                                "ICIA" = c(41.6611, -91.5302),
                                                                "MAWI" = c(43.0731, -89.4012),
                                                                "WIDE" = c(39.7447, -75.5484),
                                                                "ININ" = c(39.7684, -86.1581),
                                                                "AUTX" = c(30.2672, -97.7431),
                                                                "MAKS" = c(39.1836, -96.5717),
                                                                "FOCO" = c(40.5853, -105.0844),
                                                                "DECO" = c(39.7392, -104.9903),
                                                                "LBCA" = c(33.7701, -118.1937))
                      }
                      
                      # function to get time zone
                      cityTZ <-  function(city_code){switch(city_code,
                                                                "CHIL" = "US/Central",
                                                                "ICIA" ="US/Central",
                                                                "MAWI" = "US/Central",
                                                                "WIDE" = "US/Eastern",
                                                                "ININ" = "US/Eastern",
                                                                "AUTX" = "US/Central",
                                                                "MAKS" = "US/Central",
                                                                "FOCO" = "US/Mountain",
                                                                "DECO" = "US/Mountain",
                                                                "LBCA" = "US/Pacific")
                      }
                      
                      # run functions to get lat/lon and time zones
                      latlon <- cityLatLon(substr(city_data[i],1,4))
                      tz <- cityTZ(substr(city_data[i],1,4))
        
                      for(j in 1:nrow(data)){
              
                        # calculate the sun times for each data point
                        sun_times[j,] <- format(getSunlightTimes(date=data$date[j],
                                                                 lat=latlon[1],
                                                                 lon=latlon[2],
                                                                 keep=c("dawn","sunrise",
                                                                        "sunriseEnd","dusk",
                                                                        "solarNoon",
                                                                        "sunsetStart"
                                                                        ,"sunset","night",
                                                                        "nadir","nightEnd")
                                                                 , tz=tz)[4:13], "%H:%M:%S")
                      }
        
                      # convert all to hours minutes seconds (period object)
                      sun_times2 <- sun_times %>% mutate_all(hms)
                      
                      # add only the data needed to the list element
                      full_df <- data.frame(city=substr(city_data[i],1,4),
                                            data[,c("ShortName","ImageDate","time","date",
                                                    "LocationName","UTM_E","UTM_N","UTMZone")], 
                                            sun_times2, 
                                            time_cat = rep(0, nrow(data)))
        
                      # label time categories for each photo
                      # Day - Reference
                      full_df$time_cat[which(as.numeric(full_df$time) >= as.numeric(full_df$sunriseEnd) & 
                                                    as.numeric(full_df$time) <= as.numeric(full_df$dusk))]=1
                      # Dawn
                      full_df$time_cat[which(as.numeric(full_df$time) >= as.numeric(full_df$dawn) & 
                                                  as.numeric(full_df$time) <= as.numeric(full_df$sunriseEnd))]=2
                      # Dusk
                      full_df$time_cat[which(as.numeric(full_df$time) >= as.numeric(full_df$dusk) & 
                                                  as.numeric(full_df$time) <= as.numeric(full_df$night))]=3
                      # Night - can change all after night and before dawn to night time_utc
                      # next steps will overwrite these night categories to dark night
                      full_df$time_cat[which(as.numeric(full_df$time) >= as.numeric(full_df$night) |
                                                  as.numeric(full_df$time) < as.numeric(full_df$dawn))]=4
                      
                      time_check <- function(x){
                        x <- as.numeric(x)
                        ifelse(x> 86400, x - 86400, x)
                      }
                      
                      # Dark Night - when darkest hour is before midnight
                      for(n_row in 1:nrow(full_df)){
                        if(as.numeric(full_df[n_row,"nadir"]) > (3600*12)){
                          if(as.numeric(full_df[n_row,"time"]) >= (as.numeric(full_df[n_row,"nadir"]) - 3600)){
                            full_df[n_row,"time_cat"] <- 5
                          }
                          if(as.numeric(full_df[n_row,"time"]) <= time_check(full_df[n_row,"nadir"] + 3600)){
                            full_df[n_row,"time_cat"] <- 5
                          }
                        }
                      }
                      
                      # Dark Night - when darkest hour is after midnight
                      for(n_row in 1:nrow(full_df)){
                        if(as.numeric(full_df[n_row,"nadir"]) < (3600*12)){
                          if(as.numeric(full_df[n_row,"time"]) <= (as.numeric(full_df[n_row,"nadir"]) + 3600)){
                            full_df[n_row,"time_cat"] = 5
                          }
                          if(as.numeric(full_df[n_row,"time"]) >= (86400 - (3600 - as.numeric(full_df[n_row,"nadir"])))){
                            full_df[n_row,"time_cat"] <- 5
                          }
                        }
                      }
  
  return(full_df)

                     }

stopCluster(cl)


# save the full processed dataset
saveRDS(data_full, "./Data/2019-03-26_cleaned_photo_dataset_9_species.rds")

