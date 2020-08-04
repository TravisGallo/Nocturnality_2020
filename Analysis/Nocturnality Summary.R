data <- readRDS("~/Downloads/2019-03-26_data_list_8sp.rds")

# matrix of total species detections
species_detections <- do.call(rbind, lapply(data, nrow))

# median and mean detections per species
med_detections <- median(species_detections[,1])
mean_detections <- mean(species_detections[,1])

# top species and where they were found

species_in_cities <- lapply(data, function(x) { unique(x$city)})



names(data$bobcat)
