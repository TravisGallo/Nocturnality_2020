
# load packages and functions
source("./Scripts/2021_Gallo_etal_DielActivity_UtilityScripts.R")

# load cleaned data set
data_list <- readRDS("./Data/2021_Gallo_etal_DielActivity_data_list_8sp.rds")

# Mean center the covariates by city to reduce the sensitivity to sample size
# this is in response to a reviewer comment and therefore not in the original
# data processing script

# removed previously scaled variables to eliminate confusion
data_list <- lapply(data_list, function(x){
  df <- x %>%
    dplyr::select(!(habitat_scaled:tavg_scaled))
  return(df)
})

# need to log population data first
data_list <- lapply(data_list, function(x){
  x$log_pd <- log(x$pd + 0.01)
  return(x)
})

# calculate group mean center with city as group
gmc_variables <- lapply(data_list, function(x){
  vars <- x %>% 
    dplyr::select(c(city,tavg:log_pd)) 
  gmc <- vars %>%
    group_by(city) %>%
    summarize_all(scale, scale = F)
  gmc_scaled <- gmc %>%
    mutate(tavg_scaled = tavg / sd(vars$tavg),
           habitat_scaled = avail_habitat / sd(vars$avail_habitat),
           imp_scaled = mean_imp / sd(vars$mean_imp),
           pd_scaled = log_pd / sd(vars$log_pd),
           veg_scaled = prop_veg / sd(vars$prop_veg)) %>%
    as.data.frame(.)
  return(gmc_scaled)
})

# remove rows with NaN since a single city observation will not work anymore
gmc_variables <- lapply(gmc_variables, na.omit)

saveRDS(gmc_variables, "./Data/group_mean_centered_variables.rds")

# need to make this better
data_list[[1]] <- data_list[[1]][-1,]
# function to set up data and run JAGS model
runJAGS <- function(species){
  
  # species specific data
  df <- data_list[[species]]
  
  # species specific scaled variables
  scaled_vars <- gmc_variables[[species]]
  
  # cottontail workaround to remove LBCA data (desert cottontail)
  if(species == "e_cottontail"){
    df <- df[-which(df$city == "LBCA"),]
  }
  
  # creating some dimensions for JAGS
  ndata <- nrow(df) # number of data points
  ncat <- n_distinct(df$time_cat) # number of categories
  nbeta <- 5 # number of variables
  
  # making a named data list to give to JAGS
  jags_data <- list(y = as.numeric(factor(df$time_cat)),
                    veg = scaled_vars$veg_scaled[,1],
                    imp = scaled_vars$imp_scaled[,1],
                    log_pd = scaled_vars$pd_scaled[,1],
                    habitat = scaled_vars$habitat_scaled[,1],
                    temp = scaled_vars$tavg_scaled[,1],
                    avail_time = as.matrix(df[,9:13]),
                    ndata = ndata,
                    ncat = ncat,
                    city = as.numeric(factor(df$city)),
                    ncity = n_distinct(df$city),
                    nbeta = nbeta)
  
  # JAGS settings
  
  # create initial values
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
  mod_mcmc <- as.mcmc.list(run.jags(model = "./Scripts/Analysis/2021_Gallo_etal_DielActivity_SoftmaxModel.R",
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

# loop through and run model for each species
model_list <- vector("list", length(data_list))
for(i in 1:length(jags_list)){  
  jags_list[[i]] <- runJAGS(names(data_list)[i])
  saveRDS(jags_list[[i]], paste0("species_",i,".rds"))
  names(jags_list)[i] <- names(data_list)[i]
}

# models take a long time to run so we saved for later use
#saveRDS(model_list, "model_list.rds")

