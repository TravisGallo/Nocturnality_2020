#load packages and functions
source("nocturnality_utility.R")

# load cleaned data set
data_list <- readRDS("./Data/2019-03-26_data_list_8sp.rds")

# divide data up by city
# raccoons are in all cities
cities <- unique(data_list[[4]]$city)

# creat a list within a list (species then city)
species_city_list <- lapply(data_list, function(x){
  city_data <- vector("list", length(cities))
  for(city in 1:length(cities)){
    city_data[[city]] <- x[which(x$city == cities[city]),]
  }
  names(city_data) <- cities
  return(city_data)
})

# caclulate a risk ratio for each species and each city
rr_list <- vector("list", length(species_city_list))
for(i in 1:length(rr_list)){
  rr_list[[i]] <- lapply(species_city_list[[i]], function(x){
    # seperate out photos above and below mean urbanization
    urb_high <- x[which(x$urban_index > mean(x$urban_index)),]
    urb_low <- x[which(x$urban_index < mean(x$urban_index)),]
    
    # calculate proportion of nightime photos for each urbanization level
    x_high <- (nrow(urb_high[which(urb_high$time_cat == 4 | 
                                          urb_high$time_cat == 5),]) + 0.0001)/
      (nrow(urb_high) + 0.0001)
    
    x_low <- (nrow(urb_low[which(urb_low$time_cat == 4 | 
                                    urb_low$time_cat == 5),]) + 0.0001)/
      (nrow(urb_low) + 0.0001)
    
    log(x_high/x_low) })
}

# Calculate risk ratio variance from Gaynor et al. 2018
variance_list <- vector("list", length(species_city_list))
for(i in 1:length(rr_list)){
  variance_list[[i]] <- lapply(species_city_list[[i]], function(x){
    # seperate out photos above and below mean urbanization
    urb_high <- x[which(x$urban_index > mean(x$urban_index)),]
    urb_low <- x[which(x$urban_index < mean(x$urban_index)),]
    # variance calculation
    ((1/nrow(urb_high[which(urb_high$time_cat == 4 | urb_high$time_cat == 5),]) + 0.0001) -
      (1/nrow(urb_high) + 0.0001)) + 
      ((1/nrow(urb_low[which(urb_low$time_cat == 4 |  urb_low$time_cat == 5),]) + 0.0001) -
      (1/nrow(urb_low) + 0.0001))
     })
}


# collapse the lists down to one list
rr_species <- lapply(rr_list, function(x){ do.call("c", x) })
variance_species <- lapply(variance_list, function(x){ do.call("c", x) })
# collapse list down to a matrix
rr_mat <- do.call("rbind", rr_species)
var_mat <- do.call("rbind", variance_species)
# give rownames
rownames(rr_mat) <- names(species_city_list)
rownames(var_mat) <- names(species_city_list)

# find where each species had at least 1 detection
det_list <- vector("list", length(species_city_list))
for(i in 1:length(det_list)){
  det_list[[i]] <- lapply(species_city_list[[i]], function(x){nrow(x)})
}

det_species <- lapply(det_list, function(x){ do.call("c", x) })
det_mat <- do.call("rbind", det_species)
rownames(det_mat) <- names(species_city_list)

# make risk ratios NA if we never detected the species in the city
rr_mat[det_mat == 0] <- NA
var_mat[det_mat == 0] <- NA

# if we could not calculate variance we will not report the RR
rr_mat[var_mat == Inf | var_mat == "NaN"] <- NA

# find average housing density for each city
sitecovs <- readRDS("2019-03-18_sitecovsBuffered.rds")
mean_urb <- sitecovs$buffer1000$sitecovs %>% 
  group_by(city) %>% 
  summarise(mean_urb = mean(urban_index)) %>% 
  as.data.frame()
# order least urban to urban
mean_urb <- mean_urb[order(mean_urb$mean_urb),]

# order the rr_mat in the same way and make it only two decimal points
ord <- as.character(mean_urb$city)
rr_mat_ordered <- round(rr_mat[,ord],2)

## Plot a heat map
library(gplots)
library(RColorBrewer)

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("blue", "white", "black"))(n = 100)

tiff("test_rr.tif", width = 10, height = 10, units = "in", res = 300, compression = "lzw")

heatmap.2(rr_mat_ordered,
          
          # dendrogram control
          Rowv = FALSE,
          Colv = FALSE,
          dendrogram="none",
          
          # data scaling
          scale = "none",
          
          # colors
          col = my_palette,
          symbreaks = TRUE,
          
          # seperations
          rowsep = 0:8, colsep = 0:11, sepwidth = c(0.02,0.02),
          sepcolor = "white",
          
          # cell labeling
          cellnote = rr_mat_ordered,
          notecex =1 ,
          notecol = "hotpink",
          na.color = "antiquewhite2",
          
          # level trace
          trace="none",
          
          # row and column labels
          margins = c(10,10),
          labRow = c("Bobcat", "Coyote", "E. cottontail", 
                     "Raccoon", "Red fox", "Striped skunk",
                     "V. opossum", "White-tailed deer"),
          srtRow = 0,
          adjRow = c(0.5,0.6),
          offsetRow = -56,
          labCol = c("Manhattan, KS", "Iowa City, IA", "Indianapolis, IN",
                     "Ft. Collins, CO", "Austin, TX", "Wilmington, DE", "Madison, WI",
                     "Denver, CO", "Chicago, IL", "Long Beach, CA"),
          srtCol = 45,
          
          # color key + density info
          key = FALSE,
          
          # plot labels
          main = "",
          xlab = "",
          ylab = "")

dev.off()
 