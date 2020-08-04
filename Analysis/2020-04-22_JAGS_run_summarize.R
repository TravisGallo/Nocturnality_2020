
# load packages and functions
source("nocturnality_utility.R")

# load cleaned data set
data_list <- readRDS("./Data/2019-03-26_data_list_8sp.rds")


# function to set up data and run JAGS model
runJAGS <- function(species){
  
  # species specific data
  df <- data_list[[species]]
  
  # cottontail workaround to remove LBCA data (desert cottontail)
  if(species == "e_cottontail"){
    df <- df[-which(df$city == "LBCA"),]
  }
  
  # creating some dimensions
  ndata <- nrow(df) # number of data points
  ncat <- n_distinct(df$time_cat) # number of categories
  nbeta <- 5 # number of variables
  
  # making a named list to give to JAGS
  jags_data <- list(y = as.numeric(factor(df$time_cat)),
                    veg = df$veg_scaled[,1],
                    imp = df$imp_scaled[,1],
                    log_pd = scale(log(df$pd + 0.01))[,1],
                    habitat = df$habitat_scaled[,1],
                    temp = df$tavg_scaled[,1],
                    avail_time = as.matrix(df[,9:13]),
                    ndata = ndata,
                    ncat = ncat,
                    city = as.numeric(factor(df$city)),
                    ncity = n_distinct(df$city),
                    nbeta = nbeta)
  
  # JAGS settings
  
  # Create initial values
  # function that generates starting values for parameters
  inits <- function(chain){
    gen_list <- function(chain = chain){
      list(b_mu = c(NA, rnorm(ncat-1)),
           b1 = c(NA, rnorm(ncat-1)),
           b2 = c(NA, rnorm(ncat-1)),
           b3 = c(NA, rnorm(ncat-1)),
           b4 = c(NA, rnorm(ncat-1)),
           b5 = c(NA, rnorm(ncat-1)),
           b.lambda = runif(nbeta,1,10),
           tau = c(NA, rgamma(ncat-1,1,1)),
           .RNG.name = switch(chain,
                              "1" = "base::Wichmann-Hill",
                              "2" = "base::Marsaglia-Multicarry",
                              "3" = "base::Super-Duper",
                              "4" = "base::Mersenne-Twister",
                              "5" = "base::Wichmann-Hill",
                              "6" = "base::Marsaglia-Multicarry",
                              "7" = "base::Super-Duper",
                              "8" = "base::Mersenne-Twister",
                              "9" = "base::Wichmann-Hill",
                              "10" = "base::Marsaglia-Multicarry",
                              "11" = "base::Super-Duper",
                              "12" = "base::Mersenne-Twister",
                              "13" = "base::Wichmann-Hill",
                              "14" = "base::Marsaglia-Multicarry"),
           .RNG.seed = sample(1:1e+06, 1)
      )
    }
    return(switch(chain,           
                  "1" = gen_list(chain),
                  "2" = gen_list(chain),
                  "3" = gen_list(chain),
                  "4" = gen_list(chain),
                  "5" = gen_list(chain),
                  "6" = gen_list(chain),
                  "7" = gen_list(chain),
                  "8" = gen_list(chain),
                  "9" = gen_list(chain),
                  "10" = gen_list(chain),
                  "11" = gen_list(chain),
                  "12" = gen_list(chain),
                  "13" = gen_list(chain),
                  "14" = gen_list(chain)
    )
    )
  }
  
  # parameters to monitor
  params <- c("b0", "b_mu", "b1", "b2", "b3", "b4", "b5", "sigma", "b.lambda")
  
  # settings
  n_chains <- 14
  adapt_steps <- 8000
  burn_in <- 8000
  sample_steps <- 60000
  thin_steps <- 2
  
  # run the jags model
  mod_mcmc <- as.mcmc.list(run.jags(model = "2020-04-22_softmax_offset_ind_vars_CATS.R",
                                    modules = "glm on",
                                    monitor = params, 
                                    data = jags_data,  
                                    inits = inits, 
                                    n.chains = n_chains,
                                    adapt = adapt_steps,
                                    burnin = burn_in,
                                    sample = ceiling(sample_steps / n_chains),
                                    thin = thin_steps,
                                    summarise = FALSE,
                                    plots = FALSE,
                                    method = "parallel"))
  
  return(mod_mcmc)
}

# names for saving
species_names <- c("bobcat", "coyote", "e_cottontail", "raccoon", "red_fox", "striped_skunk",
                   "v_opossum", "w_t_deer")

jags_list <- vector("list", length(data_list))
for(i in 5:length(jags_list)){  
  jags_list[[i]] <- runJAGS(names(data_list)[i])
  names(jags_list)[i] <- names(data_list)[i]
  saveRDS(jags_list[[i]], paste0(species_names[i],"_model_run_ind_vars",Sys.Date(),".rds"))
}

# save JAGS list for later use
saveRDS(jags_list, paste0(Sys.Date(),"_jags_model_ind_vars_list.rds"))


## Plotting

# load in data from saved models
model_list <- list(bobcat = readRDS("bobcat_model_run_ind_vars2020-04-22.rds"),
                   coyote = readRDS("coyote_model_run_ind_vars2020-04-23.rds"),
                   e_cottontail = readRDS("e_cottontail_model_run_ind_vars2020-04-24.rds"),
                   racooon = readRDS("raccoon_model_run_ind_vars2020-05-02.rds"),
                   red_fox = readRDS("red_fox_model_run_ind_vars2020-04-29.rds"),
                   striped_skunk = readRDS("striped_skunk_model_run_ind_vars2020-04-29.rds"),
                   v_opposum = readRDS("v_opossum_model_run_ind_vars2020-04-30.rds"),
                   w_t_deer = readRDS("w_t_deer_model_run_ind_vars2020-05-01.rds"))

# manually check convergence
source("DBDA2E-utilities.R")
# Diagnostics
for(j in 1:8){
  for(i in 1:5){
    #diagMCMC(jags_list[[j]], paste0("b1[",i,"]"))
    diagMCMC(model_list[[j]], paste0("b1[",i,"]"))
  }
}

# plot diagnostics manually
mcmcplot(model_list[[8]])

#load data list
#data_list <- readRDS("2019-03-26_data_list_8sp.rds")

# create new data range for predicting across: 
#avail habitat
new_habitat <- seq(0, 1, 0.01)
# impervious cover
new_imp <- seq(0, 100, 1)
# vegetation cover
new_veg <- seq(0,1,0.01)
# population
new_pop <- seq(log(0.01),log(14000),0.1)
# daily average temp
new_temp <- seq(-24, 33, 0.5)

# function to calculate linear predictions
linearPrediction <- function(species){
  
  # call correct model
  data <- as.matrix(model_list[[species]])
  # for organizational purposes match names for both model and data list
  names(data_list) <- names(model_list)
  
  # [samples x category] matrix for each parameter
  b_mu_matrix <- data[,grep("b_mu", colnames(data))]
  b1_matrix <- data[,grep("b1", colnames(data))]
  b2_matrix <- data[,grep("b2", colnames(data))]
  b3_matrix <- data[,grep("b3", colnames(data))]
  b4_matrix <- data[,grep("b4", colnames(data))]
  b5_matrix <- data[,grep("b5", colnames(data))]
  
  # calculate mean offset
  mean_offset <- apply(data_list[[species]][,9:13], 2, mean)
  
  ## Predict avaliable habitat
  # scale new sequence with original mean
  new_habitat_scaled <- (new_habitat - mean(data_list[[species]]$avail_habitat)) / 
    sd(data_list[[species]]$avail_habitat)
  
  # predicting
  lin_pred_habitat <- probs_habitat <- array(NA, dim = c(length(new_habitat_scaled), 
                                                         ncol(b1_matrix), 
                                                         nrow(b1_matrix)))
  
  for(samp in 1:dim(lin_pred_habitat)[3]){
    for(cat in 1:dim(lin_pred_habitat)[2]){
      for(predictor in 1:dim(lin_pred_habitat)[1]){
        # create a linear predictor at each step of the mcmc
        lin_pred_habitat[predictor,cat,samp] <- log(mean_offset[cat]) + b_mu_matrix[samp,cat] + 
          b1_matrix[samp,cat] * new_habitat_scaled[predictor]
      }
    }
    # calculate the probabilities
    probs_habitat[,,samp] <- exp(lin_pred_habitat[,,samp]) / 
      rowSums(exp(lin_pred_habitat[,,samp]))
  }
  
  # calculate the quantiles
  quants_habitat <- apply(probs_habitat, c(1,2), quantile, probs = c(0.025, 0.5, 0.975))
  
  
  ## Predicting for Impervious covariate
  # scale new sequence with original mean
  new_imp_scaled <- (new_imp - mean(data_list[[species]]$mean_imp)) / 
    sd(data_list[[species]]$mean_imp)
  
  # predicting
  lin_pred_imp <- probs_imp <- array(NA, dim = c(length(new_imp_scaled), 
                                                 ncol(b2_matrix), 
                                                 nrow(b2_matrix)))
  
  for(samp in 1:dim(lin_pred_imp)[3]){
    for(cat in 1:dim(lin_pred_imp)[2]){
      for(predictor in 1:dim(lin_pred_imp)[1]){
        # create a linear predictor at each step of the mcmc
        lin_pred_imp[predictor,cat,samp] <- log(mean_offset[cat]) + b_mu_matrix[samp,cat] + 
          b2_matrix[samp,cat] * new_imp_scaled[predictor]
      }
    }
    # calculate the probabilities
    probs_imp[,,samp] <- exp(lin_pred_imp[,,samp]) / 
      rowSums(exp(lin_pred_imp[,,samp]))
  }
  
  # calculate the quantiles
  quants_imp <- apply(probs_imp, c(1,2), quantile, probs = c(0.025, 0.5, 0.975))
  
  
  ## Predicting for Veg covariate
  # scale new sequence with original mean
  new_veg_scaled <- (new_veg - mean(data_list[[species]]$prop_veg)) / 
    sd(data_list[[species]]$prop_veg)
  
  # predicting
  lin_pred_veg <- probs_veg <- array(NA, dim = c(length(new_veg_scaled), 
                                                 ncol(b3_matrix), 
                                                 nrow(b3_matrix)))
  
  for(samp in 1:dim(lin_pred_veg)[3]){
    for(cat in 1:dim(lin_pred_veg)[2]){
      for(predictor in 1:dim(lin_pred_veg)[1]){
        # create a linear predictor at each step of the mcmc
        lin_pred_veg[predictor,cat,samp] <- log(mean_offset[cat]) + b_mu_matrix[samp,cat] + 
          b3_matrix[samp,cat] * new_veg_scaled[predictor]
      }
    }
    # calculate the probabilities
    probs_veg[,,samp] <- exp(lin_pred_veg[,,samp]) / 
      rowSums(exp(lin_pred_veg[,,samp]))
  }
  
  # calculate the quantiles
  quants_veg <- apply(probs_veg, c(1,2), quantile, probs = c(0.025, 0.5, 0.975))
  
  
  ## Predicting for Population covariate
  # scale new sequence with original mean
  new_pop_scaled <- (new_pop - mean(log(data_list[[species]]$pd + 0.01))) / 
    sd(log(data_list[[species]]$pd + 0.01))
  
  # predicting
  lin_pred_pop <- probs_pop <- array(NA, dim = c(length(new_pop_scaled), 
                                                 ncol(b4_matrix), 
                                                 nrow(b4_matrix)))
  
  for(samp in 1:dim(lin_pred_pop)[3]){
    for(cat in 1:dim(lin_pred_pop)[2]){
      for(predictor in 1:dim(lin_pred_pop)[1]){
        # create a linear predictor at each step of the mcmc
        lin_pred_pop[predictor,cat,samp] <- log(mean_offset[cat]) + b_mu_matrix[samp,cat] + 
          b4_matrix[samp,cat] * new_pop_scaled[predictor]
      }
    }
    # calculate the probabilities
    probs_pop[,,samp] <- exp(lin_pred_pop[,,samp]) / 
      rowSums(exp(lin_pred_pop[,,samp]))
  }
  
  # calculate the quantiles
  quants_pop <- apply(probs_pop, c(1,2), quantile, probs = c(0.025, 0.5, 0.975))
  
  
  ## Predicting for Temp covariate
  # scale new sequence with original mean
  new_temp_scaled <- (new_temp - mean(data_list[[species]]$tavg)) / 
    sd(data_list[[species]]$tavg)
  
  # predicting
  lin_pred_temp <- probs_temp <- array(NA, dim = c(length(new_temp_scaled), 
                                                   ncol(b5_matrix), 
                                                   nrow(b5_matrix)))
  
  for(samp in 1:dim(lin_pred_temp)[3]){
    for(cat in 1:dim(lin_pred_temp)[2]){
      for(predictor in 1:dim(lin_pred_temp)[1]){
        # create a linear predictor at each step of the mcmc
        lin_pred_temp[predictor,cat,samp] <- log(mean_offset[cat]) + b_mu_matrix[samp,cat] + 
          b5_matrix[samp,cat] * new_temp_scaled[predictor]
      }
    }
    # calculate the probabilities
    probs_temp[,,samp] <- exp(lin_pred_temp[,,samp]) / 
      rowSums(exp(lin_pred_temp[,,samp]))
  }
  
  # calculate the quantiles
  quants_temp <- apply(probs_temp, c(1,2), quantile, probs = c(0.025, 0.5, 0.975))
  
  ## Calculate average selection
  intercept_pred <- probs_intercept <- matrix(NA, nrow = nrow(b_mu_matrix), ncol = ncol(b_mu_matrix))
  
  for(cat in 1:ncol(intercept_pred)){
    # create a linear predictor at each step of the mcmc
    intercept_pred[,cat] <- log(mean_offset[cat]) + b_mu_matrix[,cat]
  }
  for(cat in 1:ncol(probs_intercept)){
    probs_intercept[,cat] <- exp(intercept_pred[,cat]) / 
      rowSums(exp(intercept_pred))
  }
  
  quants_intercept <- apply(probs_intercept, 2, quantile, probs = c(0.025, 0.5, 0.975))
  
  return(list(avg = quants_intercept,
              habitat = quants_habitat,
              imp = quants_imp,
              veg = quants_veg,
              pop = quants_pop,
              temp = quants_temp))
  
}

prediction_list <- vector("list", length(model_list))

for(i in 1:length(prediction_list)){
  prediction_list[[i]] <- linearPrediction(names(model_list)[i])
}



## Plotting

plotPredictions <- function(intercept, habitat, imp, veg, pop, temp){
  
  # plotting average selection
  # create color list
  cat_cols <- c("blue","yellow","orange","grey","black")
  # Create empty plot
  plot(1:10, xlim=c(0.75,5.25), ylim=c(0,1), col="white", xaxt="n", yaxt="n",
       xlab="", ylab="")
  axis(2, seq(0,1,0.5), labels=FALSE)
  axis(1,seq(1,5,1),labels=FALSE)
  for(j in 1:5){
    arrows(j,intercept[1,j],j,intercept[3,j], length = 0.05, angle = 90, code = 3, 
           col = cat_cols[j], lwd=3.5)
    points(j,intercept[2,j], pch=16, col = cat_cols[j], cex = 2.5)
  }
  
  # Plotting avaliable habitat
  # plot day category listed
  plot(habitat[2,,1] ~ new_habitat, type = "l", col = "blue", 
       xaxt="n", yaxt="n", ylim = c(0,1), ylab = "", xlab = "", lwd=3)
  axis(2, seq(0,1,0.5), labels = FALSE)
  axis(1,seq(0,1,0.25), labels = FALSE)
  polygonsFnc(new_habitat, lo = habitat[1,,1], 
              hi = habitat[3,,1], col = "blue")
  # plot dawn category listed
  lines(habitat[2,,2] ~ new_habitat, col = "yellow", lwd=3)
  polygonsFnc(new_habitat, lo = habitat[1,,2], 
              hi = habitat[3,,2], col = "yellow")
  # plot dusk category listed
  lines(habitat[2,,3] ~ new_habitat, col = "orange", lwd=3)
  polygonsFnc(new_habitat, lo = habitat[1,,3], 
              hi = habitat[3,,3], col = "orange")
  # plot night category listed
  lines(habitat[2,,4] ~ new_habitat, col = "grey", lwd=3)
  polygonsFnc(new_habitat, lo = habitat[1,,4], 
              hi = habitat[3,,4], col = "grey")
  
  # plot dark night category listed
  lines(habitat[2,,5] ~ new_habitat, col = "black", lwd=3)
  polygonsFnc(new_habitat, lo = habitat[1,,5], 
              hi = habitat[3,,5], col = "black")
  
  # Plotting impervious cover
  # plot day category listed
  plot(imp[2,,1] ~ new_imp, type = "l", col = "blue", xaxt="n", yaxt="n",
       ylim = c(0,1), ylab = "", xlab = "", lwd=3)
  axis(2, seq(0,1,0.5), labels = FALSE)
  axis(1,seq(0,100,25), labels = FALSE)
  polygonsFnc(new_imp, lo = imp[1,,1], 
              hi = imp[3,,1], col = "blue")
  # plot dawn category listed
  lines(imp[2,,2] ~ new_imp, col = "yellow", lwd=3)
  polygonsFnc(new_imp, lo = imp[1,,2], 
              hi = imp[3,,2], col = "yellow")
  # plot dusk category listed
  lines(imp[2,,3] ~ new_imp, col = "orange", lwd=3)
  polygonsFnc(new_imp, lo = imp[1,,3], 
              hi = imp[3,,3], col = "orange")
  # plot night category listed
  lines(imp[2,,4] ~ new_imp, col = "grey", lwd=3)
  polygonsFnc(new_imp, lo = imp[1,,4], 
              hi = imp[3,,4], col = "grey")
  # plot dark night category listed
  lines(imp[2,,5] ~ new_imp, col = "black", lwd=3)
  polygonsFnc(new_imp, lo = imp[1,,5], 
              hi = imp[3,,5], col = "black") 
  
  
  # Plotting veg cover
  # plot day category listed
  plot(veg[2,,1] ~ new_veg, type = "l", col = "blue", xaxt="n", yaxt="n",
       ylim = c(0,1), ylab = "", xlab = "", lwd=3)
  axis(2, seq(0,1,0.5), labels = FALSE)
  axis(1,seq(0,1,0.25), labels = FALSE)
  polygonsFnc(new_veg, lo = veg[1,,1], 
              hi = veg[3,,1], col = "blue")
  # plot dawn category listed
  lines(veg[2,,2] ~ new_veg, col = "yellow", lwd=3)
  polygonsFnc(new_veg, lo = veg[1,,2], 
              hi = veg[3,,2], col = "yellow")
  # plot dusk category listed
  lines(veg[2,,3] ~ new_veg, col = "orange", lwd=3)
  polygonsFnc(new_veg, lo = veg[1,,3], 
              hi = veg[3,,3], col = "orange")
  # plot night category listed
  lines(veg[2,,4] ~ new_veg, col = "grey", lwd=3)
  polygonsFnc(new_veg, lo = veg[1,,4], 
              hi = veg[3,,4], col = "grey")
  # plot dark night category listed
  lines(veg[2,,5] ~ new_veg, col = "black", lwd=3)
  polygonsFnc(new_veg, lo = veg[1,,5], 
              hi = veg[3,,5], col = "black") 
  
  
  # Plotting pop
  # plot day category listed
  plot(pop[2,,1] ~ new_pop, type = "l", col = "blue", xaxt="n", yaxt="n",
       ylim = c(0,1), ylab = "", xlab = "", las = 1, lwd=3)
  axis(2, seq(0,1,0.5), labels = FALSE)
  axis(1,seq(-5,10,2.5), labels = FALSE)
  polygonsFnc(new_pop, lo = pop[1,,1], 
              hi = pop[3,,1], col = "blue")
  # plot dawn category listed
  lines(pop[2,,2] ~ new_pop, col = "yellow", lwd=3)
  polygonsFnc(new_pop, lo = pop[1,,2], 
              hi = pop[3,,2], col = "yellow")
  # plot dusk category listed
  lines(pop[2,,3] ~ new_pop, col = "orange", lwd=3)
  polygonsFnc(new_pop, lo = pop[1,,3], 
              hi = pop[3,,3], col = "orange")
  # plot night category listed
  lines(pop[2,,4] ~ new_pop, col = "grey", lwd=3)
  polygonsFnc(new_pop, lo = pop[1,,4], 
              hi = pop[3,,4], col = "grey")
  # plot dark night category listed
  lines(pop[2,,5] ~ new_pop, col = "black", lwd=3)
  polygonsFnc(new_pop, lo = pop[1,,5], 
              hi = pop[3,,5], col = "black") 
  
  
  # Plotting avaliable temp
  # plot day category listed
  plot(temp[2,,1] ~ new_temp, type = "l", col = "blue", xaxt="n", yaxt="n",
       ylim = c(0,1), ylab = "", xlab = "", las = 1, lwd=3)
  axis(2, seq(0,1,0.5), labels = FALSE)
  axis(1,seq(-25,35,10), labels = FALSE)
  polygonsFnc(new_temp, lo = temp[1,,1], 
              hi = temp[3,,1], col = "blue")
  # plot dawn category listed
  lines(temp[2,,2] ~ new_temp, col = "yellow", lwd=3)
  polygonsFnc(new_temp, lo = temp[1,,2], 
              hi = temp[3,,2], col = "yellow")
  # plot dusk category listed
  lines(temp[2,,3] ~ new_temp, col = "orange", lwd=3)
  polygonsFnc(new_temp, lo = temp[1,,3], 
              hi = temp[3,,3], col = "orange")
  # plot night category listed
  lines(temp[2,,4] ~ new_temp, col = "grey", lwd=3)
  polygonsFnc(new_temp, lo = temp[1,,4], 
              hi = temp[3,,4], col = "grey")
  # plot dark night category listed
  lines(temp[2,,5] ~ new_temp, col = "black", lwd=3)
  polygonsFnc(new_temp, lo = temp[1,,5], 
              hi = temp[3,,5], col = "black")
  
}


layout_mat <- matrix(1:48, nrow = 8, ncol = 6, byrow = TRUE)
tiff("test_plots3.tiff", width = 17, height = 20, units = "in", compression = "lzw", res = 300)
layout(layout_mat)
par(mar=c(4,3.5,3,1) + 0.1)
for(i in 1:length(prediction_list)){
  plotPredictions(prediction_list[[i]][["avg"]],
                  prediction_list[[i]][["habitat"]], 
                  prediction_list[[i]][["imp"]],
                  prediction_list[[i]][["veg"]],
                  prediction_list[[i]][["pop"]],
                  prediction_list[[i]][["temp"]])
}

dev.off()

# save as an SVG
layout_mat <- matrix(1:48, nrow = 8, ncol = 6, byrow = TRUE)
svg("mn_plot_base.svg", width = 17, height = 20, bg = FALSE)
layout(layout_mat)
par(mar=c(4,3.5,3,1) + 0.1)
for(i in 1:length(prediction_list)){
  plotPredictions(prediction_list[[i]][["avg"]],
                  prediction_list[[i]][["habitat"]], 
                  prediction_list[[i]][["imp"]],
                  prediction_list[[i]][["veg"]],
                  prediction_list[[i]][["pop"]],
                  prediction_list[[i]][["temp"]])
}

dev.off()



