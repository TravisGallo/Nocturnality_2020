### update and load required packages ###
package_load<-function(packages = NULL, quiet=TRUE, verbose=FALSE, warn.conflicts=FALSE){
  
  # download required packages if they're not already
  
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload,, quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

# the packages needed for this analysis
packs <- c('dplyr','reshape2','stringr','lubridate','suncalc','nnet','rmutil','scales',
           'mcmcplots','runjags','rjags','coda',
           'raster','sp','rgdal','rgeos','maptools','sf',
           'parallel','foreach','doParallel',
           'devtools', 'rnoaa', 'gplots', 'RColorBrewer')
# load packages from utility script
package_load(packs)

## Spatial Functions
# switch function to decide the city specific CRS
cityCRS <-  function(city_code){switch(city_code,
                                       "CHIL" = "+init=epsg:26916",
                                       "ICIA" = "+init=epsg:26915",
                                       "MAWI" = "+init=epsg:4326",
                                       "WIDE" = "+init=epsg:26918",
                                       "ININ" = "+init=epsg:26916",
                                       "AUTX" = "+init=epsg:4326",
                                       "MAKS" = "+init=epsg:26914",
                                       "FOCO" = "+init=epsg:4326",
                                       "DECO" = "+init=epsg:26913",
                                       "LBCA" = "+init=epsg:26911")
}

stateFile <- function(city_code){switch(city_code,
                                        "CHIL" = "il_blk10_PLA",
                                        "ICIA" = "ia_blk10_PLA",
                                        "MAWI" = "wi_blk10_PLA",
                                        "WIDE" = "de_blk10_Census_change_1990_2010_PLA2",
                                        "ININ" = "in_blk10_PLA",
                                        "AUTX" = "tx_blk10_PLA",
                                        "MAKS" = "ks_blk10_PLA",
                                        "FOCO" = "co_blk10_PLA",
                                        "DECO" = "co_blk10_PLA",
                                        "LBCA" = "ca_blk10_Census_change_1990_2010_PLA2")
}

# Function to plot CI with polygons
polygonsFnc <- function(new_cov, lo, hi, color){
  x1 <- new_cov
  x2 <- rev(new_cov)
  y1 <- lo
  y2 <- rev(hi)
  polygon(c(x1, x2), c(y1, y2), col = alpha(color, 0.4), border = NA)
}