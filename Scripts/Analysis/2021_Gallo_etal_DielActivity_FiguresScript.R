## Script for Figures

# load packages and functions
source("./Scripts/2021_Gallo_etal_DielActivity_UtilityScripts.R")

# load cleaned data set
data_list <- readRDS("./Data/2021_Gallo_etal_DielActivity_data_list_8sp.rds")

# read in model results
# models take a long time to run so we saved for later use
#model_list <- readRDS("model_list.rds")

# same list names for book keeping
names(model_list) <- c("bobcat", "coyote", "e_cottontail", "raccoon", "red_fox",
                   "striped_skunk", "v_opposum", "w_t_deer")
names(data_list) <- names(model_list)

# read in the new group mean centered variables
gmc_variables <- readRDS("./Data/group_mean_centered_variables.rds")

# group mean centered variable created NA for the one bobcat observation in AUTX
# so we removed that observation
data_list[[1]] <- data_list[[1]][-1,]

# bring in original site covs for scaling new covariate values
site_covs <- readRDS("./Data/2021_Gallo_etal_DielActivity_sitecovsBuffered.rds")

##################################
##### Figure 4: Nocturnality #####
##################################

# setting up data to plot

# create new data range for predicting across: 

# store them in a list to loop through
new_cov <- list(
#avail habitat
new_habitat = seq(0, 1, 0.01),
# impervious cover
new_imp = seq(0, 100, 1),
# vegetation cover
new_veg = seq(0,1,0.01),
# population
new_pop = seq(log(0.01),log(14000),0.1),
# daily average temp
new_temp = seq(-24, 33, 0.5))

# function to calculate linear predictions
linearPrediction <- function(species){
  
  # call correct model
  data <- as.matrix(model_list[[species]])
  
  # call the correct covariate data for the respective species
  cov_data <- data_list[[species]]
  # log the population data
  cov_data$log_pd <- log(cov_data$pd + 0.01)
  
  # grab only covariates that we will use to make a loop easier
  covs <- cov_data[,c("avail_habitat", "mean_imp", "prop_veg", "log_pd", "tavg")]

  # [samples x category] matrix for each parameter
  b_mu_matrix <- data[,grep("b_mu", colnames(data))]
  # create a list of coeff matrices to loop through
  coeff_list <- list(b1_matrix = data[,grep("b1", colnames(data))],
                     b2_matrix = data[,grep("b2", colnames(data))],
                     b3_matrix = data[,grep("b3", colnames(data))],
                     b4_matrix = data[,grep("b4", colnames(data))],
                     b5_matrix = data[,grep("b5", colnames(data))])
  
  # calculate mean offset
  mean_offset <- apply(data_list[[species]][,9:13], 2, mean)
  
  # function to calculate the probability of each category
  # also calculates the total probability of night activity
  # this will be done one species at a time
  
  predictProbs <- function(new_cov, original_cov, beta_matrix){
    
    # scale new covariate values with original mean
    new_cov_scaled <- (new_cov - mean(new_cov)) / 
      sd(original_cov)
    
    # looping through with the linear predictor to predict probability
    lin_pred <- probs <- array(NA, dim = c(length(new_cov_scaled), 
                                           ncol(beta_matrix), 
                                           nrow(beta_matrix)))
    # linear predictor
    for(samp in 1:dim(lin_pred)[3]){
      for(cat in 1:dim(lin_pred)[2]){
        for(predictor in 1:dim(lin_pred)[1]){
          # create a linear predictor at each step of the mcmc
          lin_pred[predictor,cat,samp] <- log(mean_offset[cat]) + 
            b_mu_matrix[samp,cat] + beta_matrix[samp,cat] * 
            new_cov_scaled[predictor]
        }
      }
      
      # calculate the probabilities using softmax function
      probs[,,samp] <- exp(lin_pred[,,samp]) / 
        rowSums(exp(lin_pred[,,samp]))
    }
    
    # add a new column and calculate night and dark night combined (summed)
    new_array <- array(NA, dim = c(dim(probs)[1],
                                   dim(probs)[2] + 1,
                                   dim(probs)[3]))
    for(i in 1:dim(probs)[3]){
      new_array[,,i] <- cbind(probs[,,i],
                              probs[,4,i] + probs[,5,i])
    }
    
    # calculate the quantiles
    quants <- apply(new_array, c(1,2), quantile, 
                    probs = c(0.025, 0.5, 0.975))
    
    return(quants)
  }
  
  # loop through each covariate
  quants <- vector("list", length(new_cov))
  for(i in 1:length(new_cov)){
    quants[[i]] <- predictProbs(new_cov = new_cov[[i]], 
                         original_cov = covs[,i], 
                         beta_matrix = coeff_list[[i]])
  }
  
  # give list names for better book keeping
  names(quants) <- colnames(covs)
 
  # calculate average probability
  intercept_pred <- probs_intercept <- matrix(NA, nrow = nrow(b_mu_matrix), 
                                              ncol = ncol(b_mu_matrix))
  
  for(cat in 1:ncol(intercept_pred)){
    # create a linear predictor at each step of the mcmc
    intercept_pred[,cat] <- log(mean_offset[cat]) + b_mu_matrix[,cat]
  }
  for(cat in 1:ncol(probs_intercept)){
    probs_intercept[,cat] <- exp(intercept_pred[,cat]) / 
      rowSums(exp(intercept_pred))
  }
  
  # add a column that combines night and dark night
  new_matrix <- cbind(probs_intercept, probs_intercept[,4] + probs_intercept[,5])
  
  # calculate the median and CI's
  quants_intercept <- apply(new_matrix, 2, quantile, 
                            probs = c(0.025, 0.5, 0.975))
  
  return(list(avg = quants_intercept,
              quants = quants))
  
}

# store these results in a list, one element per species
prediction_list <- vector("list", length(model_list))

for(i in 1:length(prediction_list)){
  prediction_list[[i]] <- linearPrediction(names(model_list)[i])
}


## Plotting Figure 4

plotProbability <- function(habitat, veg, imp, pop, temp){
  
  # Plotting available habitat
  if(i == 5 | i == 7 | i == 8){
    plot(habitat[2,,6] ~ new_cov[["new_habitat"]], type = "l", col = "black", 
         xaxt="n", yaxt="n", ylim = c(0,1), ylab = "", xlab = "", lwd=3)
    axis(2, seq(0,1,0.5), labels = FALSE)
    axis(1,seq(0,1,0.25), labels = FALSE)
    polygonsFnc(new_cov[["new_habitat"]], lo = habitat[1,,6], 
                hi = habitat[3,,6], col = "black")
  } else {
    plot(habitat[2,,6] ~ new_cov[["new_habitat"]], type = "l", col = "gray60", 
         xaxt="n", yaxt="n", ylim = c(0,1), ylab = "", xlab = "", lwd=3)
    axis(2, seq(0,1,0.5), labels = FALSE)
    axis(1,seq(0,1,0.25), labels = FALSE)
    polygonsFnc(new_cov[["new_habitat"]], lo = habitat[1,,6], 
                hi = habitat[3,,6], col = "gray")
  }
  
  # Plotting veg cover
  if(i == 3){
    plot(veg[2,,6] ~ new_cov[["new_veg"]], type = "l", col = "black", xaxt="n", yaxt="n",
         ylim = c(0,1), ylab = "", xlab = "", lwd=3)
    axis(2, seq(0,1,0.5), labels = FALSE)
    axis(1,seq(0,1,0.25), labels = FALSE)
    polygonsFnc(new_cov[["new_veg"]], lo = veg[1,,6], 
                hi = veg[3,,6], col = "black")
  } else {
    plot(veg[2,,6] ~ new_cov[["new_veg"]], type = "l", col = "gray60", xaxt="n", yaxt="n",
         ylim = c(0,1), ylab = "", xlab = "", lwd=3)
    axis(2, seq(0,1,0.5), labels = FALSE)
    axis(1,seq(0,1,0.25), labels = FALSE)
    polygonsFnc(new_cov[["new_veg"]], lo = veg[1,,6], 
                hi = veg[3,,6], col = "gray")
  }
  
  # Plotting impervious cover
  if(i == 8 | i == 4){
    plot(imp[2,,6] ~ new_cov[["new_imp"]], type = "l", col = "black", xaxt="n", yaxt="n",
         ylim = c(0,1), ylab = "", xlab = "", lwd=3)
    axis(2, seq(0,1,0.5), labels = FALSE)
    axis(1,seq(0,100,25), labels = FALSE)
    polygonsFnc(new_cov[["new_imp"]], lo = imp[1,,6], 
                hi = imp[3,,6], col = "black")
  } else {
    plot(imp[2,,6] ~ new_cov[["new_imp"]], type = "l", col = "gray60", xaxt="n", yaxt="n",
         ylim = c(0,1), ylab = "", xlab = "", lwd=3)
    axis(2, seq(0,1,0.5), labels = FALSE)
    axis(1,seq(0,100,25), labels = FALSE)
    polygonsFnc(new_cov[["new_imp"]], lo = imp[1,,6], 
                hi = imp[3,,6], col = "gray")
  }
  
  # Plotting population
  if(i == 2 | i == 8 | i == 3){
    plot(pop[2,,6] ~ new_cov[["new_pop"]], type = "l", col = "black", xaxt="n", yaxt="n",
         ylim = c(0,1), ylab = "", xlab = "", las = 1, lwd=3)
    axis(2, seq(0,1,0.5), labels = FALSE)
    axis(1,seq(-5,10,2.5), labels = FALSE)
    polygonsFnc(new_cov[["new_pop"]], lo = pop[1,,6], 
                hi = pop[3,,6], col = "black")
  } else {
    plot(pop[2,,6] ~ new_cov[["new_pop"]], type = "l", col = "gray60", xaxt="n", yaxt="n",
         ylim = c(0,1), ylab = "", xlab = "", las = 1, lwd=3)
    axis(2, seq(0,1,0.5), labels = FALSE)
    axis(1,seq(-5,10,2.5), labels = FALSE)
    polygonsFnc(new_cov[["new_pop"]], lo = pop[1,,6], 
                hi = pop[3,,6], col = "gray")
  }
  
  # Plotting avg daily temp
  if(i == 3 | i == 4 | i == 6 | i == 7 | i == 8){
    plot(temp[2,,6] ~ new_cov[["new_temp"]], type = "l", col = "black", xaxt="n", 
         yaxt="n", ylim = c(0,1), ylab = "", xlab = "", las = 1, lwd=3)
    axis(2, seq(0,1,0.5), labels = FALSE)
    axis(1,seq(-25,35,10), labels = FALSE)
    polygonsFnc(new_cov[["new_temp"]], lo = temp[1,,6], 
                hi = temp[3,,6], col = "black")
  } else {
    plot(temp[2,,6] ~ new_cov[["new_temp"]], type = "l", col = "gray60", xaxt="n", 
         yaxt="n", ylim = c(0,1), ylab = "", xlab = "", las = 1, lwd=3)
    axis(2, seq(0,1,0.5), labels = FALSE)
    axis(1,seq(-25,35,10), labels = FALSE)
    polygonsFnc(new_cov[["new_temp"]], lo = temp[1,,6], 
                hi = temp[3,,6], col = "gray")
  }
  
}

# species order to match the average plot
new_order <- c(1,2,5,4,6,7,3,8)

layout_mat <- matrix(1:40, nrow = 8, ncol = 5, byrow = TRUE)
{tiff("2022-03-11_ProbabilityPlot_NightCombined.tiff", width = 17, height = 20, units = "in", 
     compression = "lzw", res = 300)
layout(layout_mat)
par(mar=c(1,1,1,1))
for(i in new_order){
  plotProbability(prediction_list[[i]][[2]][["avail_habitat"]],
                  prediction_list[[i]][[2]][["prop_veg"]],
                  prediction_list[[i]][[2]][["mean_imp"]],
                  prediction_list[[i]][[2]][["log_pd"]],
                  prediction_list[[i]][[2]][["tavg"]])
}

dev.off()}

#####################################
##### Figure 3: selection plots #####
####### Written by M. Fidino ########
#####################################

# turn each into a matrix obect

model_list <- lapply(
  model_list,
  function(x) as.matrix(as.mcmc.list(x), chains = TRUE)
)

# calculate quantiles
model_sum <- lapply(
  model_list,
  function(x) t(
    apply(
      x,
      2,
      function(y) quantile(exp(y), probs = c(0.025,0.5,0.975))
    )
  )
)

# see which ones are significant
msign <- lapply(
  model_sum,
  function(x) t(apply(x, 1, function(x) all(x>1)| all(x<1)))
)


sigs <- model_sum
for(i in 1:length(model_sum)){
  sigs[[i]] <- cbind(model_sum[[i]], as.numeric(msign[[i]]))
  colnames(sigs[[i]])[4] <- "significant"
}

# order them by time category. Specifically the slope terms

slopes <- lapply(
  sigs,
  function(x) x[grep("b1|b2|b3|b4|b5", row.names(x)),]
)

# order by time category
slopes <- lapply(
  slopes,
  function(x) x[order(rep(1:5, 5)),]
)

# add the names of the covariates
terms <- c(
  "Avail. greenspace",
  "Impervious cover",
  "Vegetation cover",
  "Population",
  "Daily avg. temp"
)
cats <- c("Day", "Dawn", "Dusk", "Night", "Deep night")

for(i in 1:length(slopes)){
  slopes[[i]] <- data.frame(slopes[[i]])
  slopes[[i]]$covar <- rep(terms, 5)
  slopes[[i]]$cat <- rep(cats, each = 5)
}

# reorder slopes to match other plots
slopes <- slopes[new_order]

# set up plots for each species

xl <- 1

for(i in 1:length(slopes)){
  
  tiff(
    paste0(names(slopes)[i], "_selection_estimates.tiff"),
    height = 4,
    width = 3,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  
  
  {tiff(
    paste0("selection_plot2.tiff"),
    height = 8,
    width = 8,
    units = "in",
    res = 600,
    compression = "lzw"
  )
    
    #windows(8,8)
    mr <- 3
    no <- 2
    m <- matrix(
      c(rep(0,no), rep(1:4,each = mr),rep(0,no), rep(5:8, each = mr)),
      ncol = mr*4 + no,
      nrow = 2,
      byrow = TRUE
    )
    
    pnames <- c(
      "Bobcat",
      "Coyote",
      "Red fox",
      "Raccoon",
      "Striped skunk",
      "Virginia opossum",
      "Eastern cottontail",
      "White-tailed deer"
    )
    
    layout(m)
    #par(mar = c(4,1,1,1))
    for(i in 1:length(slopes)){
      df <- slopes[[i]]
      
      df <- df[-which(df$cat == "Day"),]
      df
      if(i %in% c(1:4)){
        par(mar = c(1,1,4,1))
      } else {
        par(mar = c(4,1,1,1))
      }
      
      plot(1~1, type = "n", xlim = c(0,4), ylim = c(1,24), xlab = "", ylab = "",
           xaxt = "n", yaxt = "n", bty = "l")
      ur <- par("usr")
      my_y <- 1:24
      my_y_names <- c(terms, "", terms, "", terms, "", terms, "")
      
      if(i %in% c(1,5)){
        mtext(text = my_y_names, 2, line = 0.5, at = my_y, las = 1 , cex = 0.85)
      }
      lines(x = c(1,1), y = c(0, 24), lty = 2)
      par(xpd = NA)
      text(x = ur[1]-0.1, y = ur[4] + 0.5,
           labels = 
             paste0(
               LETTERS[i],
               ") ",
               pnames[i]
             ), pos = 4, cex = 1.4
      )
      
      
      for(j in 1:4){
        rect(
          xleft = ur[1],
          xright = ur[2],
          ybottom = which(my_y_names == "")[j] - 0.5,
          ytop = which(my_y_names == "")[j] + 0.5,
          col = "gray80",
          border = "black"
        )
      }
      axis(2, at = c(1:24)[my_y_names != ""], labels = FALSE, tck = -0.0123)
      axis(2, at =  1:25, labels = FALSE, tck = 0)
      text(x = 2, y = which(my_y_names == ""),
           labels = c("Dawn", "Dusk", "Night","Deep night"))
      
      ys <- my_y[-which(my_y_names == "")]
      
      for(k in 1:20){
        lines(x = df[k, c(1,3)], y = rep(ys[k], 2), lwd = 2)
      }
      
      points(x = df[,2], y = ys, pch = 21, cex = 1.2,
             bg = ifelse(df$significant == 1, "black", "white"))
      axis(1, at = seq(0, 4, 0.5), labels = FALSE, tck = -0.0125)
      axis(1, at = seq(0, 0, 0.25), labels = FALSE, tck = -0.0125/2)
      if(i %in% c(5:8)){
        mtext(
          sprintf("%.0f",seq(0,4,1)),1, line = 0.4, at = seq(0,4,1), cex= 0.85
        )
      }
      
      if(i == 6){
        mtext("Relative selection", 1, at = 5, line = 2.5, cex = 1.2)
      }
      
      #if(i == 1){
      #legend(x = -3.825, y = ur[4] + 3.6, c("no", "yes"), pch = 21, pt.bg = c("black", "white"),
      #        title = "95% CI overlaps 1?", horiz = TRUE, cex = 1.2, bty = "n")
      #}
      if(i %in% c(4,8)){
        rect(
          xleft = 2.25,
          xright = ur[2],
          ybottom = ur[3],
          ytop = 5,
          col = "white",
          border = "black"
        )
        legend("bottomright", c("yes", "no"), pch = 21, pt.bg = c("white", "black"),
               title = "95% CI\noverlaps 1?", cex = 1, bty = "n")
        if(i == 8){
          dev.off()
        }
      }
      #dev.off()
      
    }
  }
  dev.off()
}

################################
##### Figure 2: Plasticity #####
################################

# function to calculate linear predictions
linearPrediction <- function(species){
  
  # call correct model
  data <- as.matrix(model_list[[species]])
  
  # call the correct covariate data for the respective species
  cov_data <- data_list[[species]]
  
  ## Process coefficents
  
  # get city specific intercepts
  b0_matrix <- data[,grep("b0", colnames(data))]
  
  # get number of cities for the respective cities
  ncity <- n_distinct(cov_data$city)
  
  # convert to a list - each city is an element
  b0_city <- matrix(NA, nrow = ncity, ncol = 5)
  for(i in 1:ncity){
    b0_city[i,] <- apply(b0_matrix[,grep(paste0("^b0\\[",i,","), colnames(b0_matrix))],
                         2, median)
  }
  
  # create a list of coeff matrices to loop through
  coeff_mat <- matrix(NA, nrow = 5, ncol = 5)
  coeff_mat[1,] <- apply(data[,grep("b1", colnames(data))], 2, median)
  coeff_mat[2,] <- apply(data[,grep("b2", colnames(data))], 2, median)
  coeff_mat[3,] <- apply(data[,grep("b3", colnames(data))], 2, median)
  coeff_mat[4,] <- apply(data[,grep("b4", colnames(data))], 2, median)
  coeff_mat[5,] <- apply(data[,grep("b5", colnames(data))], 2, median)
  
  # calculate mean offset
  mean_offset <- apply(data_list[[species]][,9:13], 2, mean)
  
  ## Process covariate data
  
  # grab only covariates that we will use
  covs_tmp <- cov_data[,c("city", "LocationName", "habitat_scaled", "imp_scaled", 
                          "veg_scaled", "pd_scaled")]
  # remove duplicates to just get the data for each site the species was detected
  site_covs <- covs_tmp[!duplicated(covs_tmp),]
  
  # calculate avg tmp for city
  site_temp <- cov_data %>% group_by(LocationName) %>%
    summarize(avg_temp = mean(tavg_scaled))
  # join to covariate data
  site_covs_final <- left_join(site_covs, site_temp, by = "LocationName")
  
  # create a list to separate out cities
  # we will loop and predict at each site using the city-specific intercept
  cities <- unique(site_covs_final$city)
  city_sites <- list()
  for(i in 1:length(cities)){
    city_sites[[i]] <- site_covs_final[which(site_covs_final$city == cities[i]),]
  }
  
  ## Calculate the probability of each category at each site
  city_probs <- list()
  for(i in 1:length(city_sites)){
    
    # looping through with the linear predictor to predict probability
    lin_pred <- probs <- matrix(NA, nrow = nrow(city_sites[[i]]),
                                ncol = ncol(coeff_mat))
    
    for(s in 1:nrow(lin_pred)){
      for(c in 1:ncol(lin_pred)){
        # create a linear predictor at each step of the mcmc
        lin_pred[s,c] <- log(mean_offset[c]) + 
          b0_city[i,c] + 
          coeff_mat[1,c] * city_sites[[i]][s,"habitat_scaled"] +
          coeff_mat[2,c] * city_sites[[i]][s,"imp_scaled"] +
          coeff_mat[3,c] * city_sites[[i]][s,"veg_scaled"] +
          coeff_mat[4,c] * city_sites[[i]][s,"pd_scaled"] +
          coeff_mat[5,c] * city_sites[[i]][s,"avg_temp"]
      }
    }
    
    # softmax function to convert to probability
    for(c in 1:ncol(probs)){
      # calculate the probabilities with softmax function
      probs[,c] <- exp(lin_pred[,c]) / 
        rowSums(exp(lin_pred))
    }
    
    city_probs[[i]] <- probs
  }
  
  # we did this one city at a time, now combine them for a full dataset
  species_probs <- do.call("rbind", city_probs)
  
  # order by early night
  species_probs_ordered <- species_probs[order(species_probs[,4]),]
  
  return(species_probs_ordered)
  
}

# calculate species probabilities at each site
each_species <- list()
for(i in 1:length(model_list)){
  each_species[[i]] <- linearPrediction(names(model_list)[i])
}


species_name <- c("Bobcat", "Coyote", "Eastern cottontail", "Raccoon", "Red fox",
                  "Striped skunk", "Virginia opposum", "White-tailed deer")

# look at range of probabilities for each species
species_prob_range <- lapply(each_species, function(x){
  t(apply(x, 2, range)) })
range_table <- do.call("rbind", species_prob_range)

## Plotting Figure 3

new_order <- c(1,2,5,4,6,7,3,8)

layout_mat <- matrix(1:8, nrow = 4, ncol = 2, byrow=TRUE) 

{tiff(
  "2022-03-11_plasticity_fig.tiff",
  height = 175,
  width = 114,
  units = "mm",
  res = 300,
  compression = "lzw"
)
  
  layout(layout_mat)
  par(oma = c(1.5,4,0,0))
  
  for(i in new_order){
    par(mar=c(0.5,0,2,0) + 0.2)
    barplot(t(each_species[[i]]), col = c("#E5DE44","lightpink1","#B9A0D2","#053752","black"),
            space = 0, border = NA, axes = FALSE)
    title(main = species_name[i], cex = 5)
    mtext("Individual sampling sites", side = 1, line = 0, 
          outer = TRUE, font = 2)
    if(i < 5){
      axis(2, at = seq(0,1,0.2), las = 1, outer = TRUE)
      mtext("Probability of use in each time category", side = 2, line = 2.55, 
            outer = TRUE, font = 2)
    }
    
  }
  
  dev.off()
  
}

#######################################
##### Figure 1: Random intercept ######
#######################################

#95% predictive interval

getResults <- function(species){
  
  off_set <- log(apply(data_list[[species]][,9:13], 2, mean))
  
  data <- as.matrix(model_list[[species]])
  
  # [samples x category] matrix for b_mu and sigma parameters
  b_mu_matrix <- data[,grep("b_mu", colnames(data))]
  sigma_matrix <- data[,grep("sigma", colnames(data))]
  sigma_matrix <- cbind(0, sigma_matrix)
  
  # add the offset to the beta mu parameter
  mu_offset <- b_mu_matrix + rep(off_set, each = nrow(b_mu_matrix))
  
  # get probability and credible intervals for intercept with offset
  b_mu_probs <- matrix(NA, nrow = nrow(mu_offset), ncol = ncol(mu_offset))
  for(i in 1:ncol(b_mu_probs)){
    b_mu_probs[,i] <- exp(mu_offset[,i])/
      rowSums(exp(mu_offset))
  }
  
  # get quantiles
  b_mu_quants <- apply(b_mu_probs, 2, quantile, probs = c(0.025, 0.5, 0.975))
  colnames(b_mu_quants) <- c("cat1", "cat2", "cat3", "cat4", "cat5")
  
  # posterior prediction
  # we dont end up using this
  pp_y <- matrix(NA, nrow = nrow(b_mu_matrix), ncol = ncol(b_mu_matrix))
  for(i in 1:ncol(b_mu_matrix)){
    pp_y[,i] <- rnorm(nrow(b_mu_matrix), mean = mu_offset[,i], sd = sigma_matrix[,i])
  }
  # softmax to predict probability
  pp_probs <- matrix(NA, nrow = nrow(pp_y), ncol = ncol(pp_y))
  for(i in 1:ncol(pp_probs)){
    pp_probs[,i] <- exp(pp_y[,i])/
      rowSums(exp(pp_y))
  }
  # get quantiles
  pp_quants <- apply(pp_probs, 2, quantile, probs = c(0.025, 0.975))
  colnames(pp_quants) <- c("cat1", "cat2", "cat3", "cat4", "cat5")
  
  # get city specific intercepts
  b0_matrix <- data[,grep("b0", colnames(data))]
  
  # get number of cities for the respective cities
  ncity <- n_distinct(data_list[[species]]$city)
  
  # convert to a list - each city is an element
  b0_city <- list()
  for(i in 1:ncity){
    b0_city[[i]] <- b0_matrix[,grep(paste0("^b0\\[",i), colnames(b0_matrix))]
  }
  
  # add offset
  mu_b0_city <- lapply(b0_city, function(x){
    x + rep(off_set, each = nrow(x))
  })
  
  # softmax for probs in each city
  city_probs <- lapply(mu_b0_city, function(x){
    tmp <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for(i in 1:ncol(tmp)){
      tmp[,i] <- exp(x[,i])/
        rowSums(exp(x))
    }
    return(tmp)
  })
  
  city_quants <- lapply(city_probs, function(x){
    apply(x, 2, quantile, probs = c(0.025, 0.5, 0.975))
  })
  
  names(city_quants) <- unique(data_list[[species]]$city)
  
  
  return(list(b_mu_quants = b_mu_quants,
              pp_quants = pp_quants,
              city_quants = city_quants
  ))
  
}

results_list <- vector("list", length(model_list))
for(s in 1:length(model_list)){
  results_list[[s]] <- getResults(names(model_list)[s])
}

names(results_list) <- names(model_list)

# fix and create vector of species names for y-lab
tp <- gsub("_"," ", names(results_list))

tp <- gsub("e ", "", tp)

tp <- gsub("v", "Virginia", tp, ignore.case = FALSE)

tp <- gsub("w t", "White-tailed", tp, ignore.case = FALSE)

# change names capitalization
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
tp <- firstup(tp)

# we want it to plot in the same order as our other probability figure
# so reorder the list
names(results_list) <- tp
reorder_names <- rev(c("Bobcat", "Coyote", "Red fox", "Raccoon", "Striped skunk", 
                       "Virginia opposum", "Cottontail", "White-tailed deer"))
results_list <- results_list[reorder_names]

# set up layout for multi-panel plot
l_mat <- matrix(1:5, 1, 5, byrow = TRUE)

# vector of main titles to loop through for each panel
main_title <- c("Day", "Dawn", "Dusk", "Night", "Deep Night")

# And the figure!
{tiff(
  "random_interecept_figure.tiff",
  height = 4,
  width = 12,
  units = "in",
  res = 300,
  compression = "lzw"
)
  
  #par(xpd = NA)
  par(oma = c(2,19,2,0))
  layout(l_mat)
  
  for(c in c(2,1,3:5)){
    
    par(mar=c(2,1,2,1) + 0.1)
    
    # y_min modifies the vertical placement of the legend in this
    #  plot. I used this to play around with the appropriate placement 
    
    y_min <- 16.5
    plot(1~1, type = "n", xlim = c(0,0.8), 
         ylim = c(0.5,length(results_list)+0.5),
         bty = "l", xlab = "",
         ylab = "", yaxt = "n",
         xaxt = "n",
         yaxs = "i",xaxs = "i")
    
    set.seed(90210)
    
    axis(1, at= seq(0,0.8, 0.2), labels = F, tck = -0.025, lwd = 1)
    axis(1, at= seq(0,0.8, 0.2/2), labels = F, tck = -0.0125, lwd = 1)
    
    axis(2, at = seq(0.5, length(results_list)+0.5),
         labels = F, tck = 0, lwd = 1)
    axis(2, at = seq(1, length(results_list)),
         labels = F, tck = -0.0125, lwd = 1)
    
    mtext(
      sprintf("%.1f", seq(0,0.8, 0.2)),
      1,
      at = seq(0,0.8,0.2),
      line = 0.5,
      cex = 0.9
    )
    
    if(c == 2){
      mtext(names(results_list), 2, line = 0.5, at = seq(1, length(results_list)),
            cex = 2, las = 1)
    }
    
    if(c == 3){
      mtext("Probability of activity", 1, 2.5, font = 2, cex = 2, outer=FALSE)
    }
    
    mtext(main_title[c], 3, 1, cex = 2, font = 2)
    
    for(i in 1:length(results_list)){
      
      #  The city-specific estimates for each species
      #  This is the reason we set the seed at the beginning.
      #  We are jittering the y-axis a TINY bit.
      
      city_medians <- lapply(results_list[[i]][["city_quants"]], function(x){
        x[2,c]
      })
      
      points(
        x = do.call("c", city_medians),
        y = jitter(rep(i, length(city_medians)), 0.2, amount = 0.2),
        pch = 21, bg = scales::alpha("grey", 0.5),
        col = "darkgrey",
        cex = 2
      )
      
      
      # plot credible interval
      lines(
        x = results_list[[i]][["b_mu_quants"]][c(1,3),c],
        y = rep(i,2), 
        col = "black",
        lwd = 4 
      )
      
      
      # Add among-city average occupancy for each species
      points(
        x = results_list[[i]][["b_mu_quants"]][2,c],
        y = i, 
        cex = 1.75,
        pch=21,
        bg = "black"
      )
    }
  }
  
  dev.off()}


#################################
##
## Calculating Odds Ratios
##
#################################

# function to calculate linear predictions
calculateOddsRatio <- function(species){
  
  # call correct model
  data <- as.matrix(model_list[[species]])
  
  # call the correct covariate data for the respective species
  cov_data <- data_list[[species]]
  # log the population data
  cov_data$log_pd <- log(cov_data$pd + 0.01)
  
  # grab only covariates that we will use to make a loop easier
  covs <- cov_data[,c("avail_habitat", "mean_imp", "prop_veg", "log_pd", "tavg")]
  
  # [samples x category] matrix for each parameter
  b_mu_matrix <- data[,grep("b_mu", colnames(data))]
  # create a list of coeff matrices to loop through
  coeff_list <- list(b1_matrix = data[,grep("b1", colnames(data))],
                     b2_matrix = data[,grep("b2", colnames(data))],
                     b3_matrix = data[,grep("b3", colnames(data))],
                     b4_matrix = data[,grep("b4", colnames(data))],
                     b5_matrix = data[,grep("b5", colnames(data))])
  
  # calculate mean offset
  mean_offset <- apply(data_list[[species]][,9:13], 2, mean)
  
  # function to calculate the probability of each category
  # also calculates the total probability of night activity
  # this will be done one species at a time
  
  covariateOddsRatios <- function(beta_matrix){
    
    # looping through with the linear predictor to predict probability
    pred <- probs <- array(NA, dim = c(nrow(beta_matrix), 
                                       ncol(beta_matrix),
                                       3))
    
    # create a linear predictor for each step of the mcmc at each sd
    for(st_dev in 0:2){
      for(cat in 1:ncol(pred)){
        pred[,cat,(st_dev+1)] <- log(mean_offset[cat]) + 
          b_mu_matrix[,cat] + beta_matrix[,cat] * st_dev
      }
    }
    
    # use linear predictor and softmax to calculate probability at each sd
    for(st_dev in 1:3){
      for(cat in 1:ncol(pred)){
        # calculate the probabilities
        probs[,cat,st_dev] <- exp(pred[,cat,st_dev]) / 
          rowSums(exp(pred[,,st_dev]))
      }
    }
    
    # add a new column and calculate night and dark night combined (summed)
    probs_array <- array(NA, dim = c(dim(probs)[1],
                                     dim(probs)[2] + 1,
                                     dim(probs)[3]))
    
    for(i in 1:dim(probs_array)[3]){
      probs_array[,,i] <- cbind(probs[,,i],
                                probs[,4,i] + probs[,5,i])
    }
    
    
    ## Calculate odds ratio for each category
    
    odds_ratio <- matrix(NA, ncol = 2, nrow = dim(probs_array)[1])
    
    for(i in 2:3){
      odds_ratio[,i-1] <- (probs_array[,6,i]/(1-probs_array[,6,i])) / 
        (probs_array[,6,1]/(1-probs_array[,6,1]))
    }
    
    # calculate the quantiles
    
    quants <- apply(odds_ratio, 2, quantile, 
                    probs = c(0.025, 0.5, 0.975))
    
    # make it a dataframe
    odds_ratio_quants <- data.frame(or_sd1 = quants[,1], or_sd2 = quants[,2])
    
    return(list(or_quants = odds_ratio_quants,
                odds_ratios = odds_ratio,
                probability_array = probs_array,
                linear_predictor_array = pred))
    
  }
  
  # run this to get information for each covariate
  species_odds_ratios <- vector("list", length(coeff_list))
  
  for(i in 1:length(coeff_list)){
    species_odds_ratios[[i]] <- covariateOddsRatios(coeff_list[[i]])
  }
  
  # give each list element a name for book keeping
  names(species_odds_ratios) <- c("avail_habitat", "imp_cover", "veg_cover", 
                                  "log_pop", "tavg")
  
  return(species_odds_ratios)
  
}

# run through and calculate the results for each species

species_odds_ratios <- vector("list", length(model_list))

for(i in 1:length(model_list)){
  species_odds_ratios[[i]] <- calculateOddsRatio(i)
}

names(species_odds_ratios) <- names(model_list)

saveRDS(species_odds_ratios, "./Analysis/2021-03-18_odds_ratios.rds")

species_odds_ratios <- readRDS("./Analysis/2021-03-18_odds_ratios.rds")

# get range of values for a 1 and 2 unit increase- toggle parameter
lapply(data_list, function(x) { mean(x$pd) })
lapply(data_list, function(x) { sd(x$pd) + mean(x$pd) })
lapply(data_list, function(x) { sd(x$pd)*2 + mean(x$pd) })

# create OR table for veg cover
veg_or <- lapply(species_odds_ratios, function(x){
  x$veg_cover$or_quants
})

veg_or_text <- do.call("rbind",
                       lapply(veg_or, function(x){
                         apply(x, 2, function(y){
                           paste0(round(y[2], 2), 
                                  " (", round(y[1], 2), 
                                  "-", round(y[3], 2), ")")
                         })
                       })
)

# create OR table for Available Habitat
habitat_or <- lapply(species_odds_ratios, function(x){
  x$avail_habitat$or_quants
})

habitat_or_text <- do.call("rbind",
                           lapply(habitat_or, function(x){
                             apply(x, 2, function(y){
                               paste0(round(y[2], 2), 
                                      " (", round(y[1], 2), 
                                      "-", round(y[3], 2), ")")
                             })
                           })
)

# create OR table for Impervious cover
imp_or <- lapply(species_odds_ratios, function(x){
  x$imp_cover$or_quants
})

imp_or_text <- do.call("rbind",
                       lapply(imp_or, function(x){
                         apply(x, 2, function(y){
                           paste0(round(y[2], 2), 
                                  " (", round(y[1], 2), 
                                  "-", round(y[3], 2), ")")
                         })
                       })
)

# create OR table for pop dense
pop_or <- lapply(species_odds_ratios, function(x){
  x$log_pop$or_quants
})

pop_or_text <- do.call("rbind",
                       lapply(pop_or, function(x){
                         apply(x, 2, function(y){
                           paste0(round(y[2], 2), 
                                  " (", round(y[1], 2), 
                                  "-", round(y[3], 2), ")")
                         })
                       })
)

# create OR table for temp
temp_or <- lapply(species_odds_ratios, function(x){
  x$tavg$or_quants
})

temp_or_text <- do.call("rbind",
                        lapply(temp_or, function(x){
                          apply(x, 2, function(y){
                            paste0(round(y[2], 2), 
                                   " (", round(y[1], 2), 
                                   "-", round(y[3], 2), ")")
                          })
                        })
)

full_or_table <- cbind(habitat_or_text, imp_or_text, veg_or_text, pop_or_text,
                       temp_or_text)

# reorder to match figures
or_table <- full_or_table[new_order,] 

write.csv(or_table, "OR_table.csv")
