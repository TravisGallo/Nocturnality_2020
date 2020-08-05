

#######################
#
# summarising and plotting the slope terms
#
# Written by M. Fidino 
#
#######################
# load packages and functions
source("./Analysis/nocturnality_utility.R")

# load in data from saved models
model_list <- list(
  bobcat = readRDS("./Analysis/bobcat_model_run_ind_vars2020-04-22.rds"),
  coyote = readRDS("./Analysis/coyote_model_run_ind_vars2020-04-23.rds"),
  e_cottontail = readRDS("./Analysis/e_cottontail_model_run_ind_vars2020-04-24.rds"),
  racooon = readRDS("./Analysis/raccoon_model_run_ind_vars2020-05-02.rds"),
  red_fox = readRDS("./Analysis/red_fox_model_run_ind_vars2020-04-29.rds"),
  striped_skunk = readRDS("./Analysis/striped_skunk_model_run_ind_vars2020-04-29.rds"),
  v_opposum = readRDS("./Analysis/v_opossum_model_run_ind_vars2020-04-30.rds"),
  w_t_deer = readRDS("./Analysis/w_t_deer_model_run_ind_vars2020-05-01.rds")
)



# turn each into a matrix obect

model_list <- lapply(
  model_list,
  function(x) as.matrix(as.mcmc.list(x), chains = TRUE)
)

# calculate quantiles
model_sum <- lapply(
  model_list,
  function(x) t(apply(x, 2, quantile, probs = c(0.025,0.5,0.975)))
)

# see which ones are significant
msign <- lapply(
  model_sum,
  function(x) t(apply(x, 1, sign))
)

msign <- lapply(
  msign,
  rowSums
)

sigs <- model_sum
for(i in 1:length(model_sum)){
  sigs[[i]] <- cbind(model_sum[[i]], as.numeric(abs(msign[[i]]) == 3))
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
  "available habitat",
  "impervious cover",
  "vegetation",
  "population density",
  "temperature"
)
cats <- c("day", "dawn", "dusk", "night", "dark night")

for(i in 1:length(slopes)){
  slopes[[i]] <- data.frame(slopes[[i]])
  slopes[[i]]$covar <- rep(terms, 5)
  slopes[[i]]$cat <- rep(cats, each = 5)
}

# set up plots for each species

xl <- 1
for(i in 1:length(slopes)){
  
  tiff(
    paste0("./Figures/", names(slopes)[i], "_slope_estimates.tiff"),
    height = 6,
    width = 6,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  df <- slopes[[i]]
  
  df <- df[-which(df$cat == "day"),]
  df
  
  #windows(6,6)
  par(mar = c(4,8,1,1))
  plot(1~1, type = "n", xlim = c(-2,2), ylim = c(1,24), xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "l")
  my_y <- 1:24
  my_y_names <- c(terms, "", terms, "", terms, "", terms, "")
  
  
  mtext(text = my_y_names, 2, line = 0.5, at = my_y, las = 1 )
  lines(x = c(0,0), y = c(0, 24), lty = 2)
  for(j in 1:4){
    rect(
      xleft = -2.5,
      xright = 2,
      ybottom = which(my_y_names == "")[j] - 0.5,
      ytop = which(my_y_names == "")[j] + 0.5,
      col = "gray80",
      border = "black"
    )
  }
  axis(2, at = c(1:24)[my_y_names != ""], labels = FALSE, tck = -0.0123)
  axis(2, at =  1:25, labels = FALSE, tck = 0)
  text(x = 0, y = which(my_y_names == ""),
       labels = c("dawn", "dusk", "night","dark night"))
  
  ys <- my_y[-which(my_y_names == "")]

  for(k in 1:20){
    lines(x = df[k, c(1,3)], y = rep(ys[k], 2), lwd = 2)
  }
  
  points(x = df[,2], y = ys, pch = 21, cex = 1.2,
         bg = ifelse(df$significant == 1, "black", "white"))
  axis(1, at = seq(-2, 2, 0.5), labels = FALSE, tck = -0.0125)
  axis(1, at = seq(-2, 2, 0.25), labels = FALSE, tck = -0.0125/2)
  mtext(
    sprintf("%.2f",seq(-2,2,0.5)),1, line = 0.4, at = seq(-2,2,0.5)
  )
  mtext("Slope term estimates", 1, at = 0, line = 1.75)
  
  #legend("bottomright", c("yes", "no"), pch = 21, pt.bg = c("black", "white"),bty = "n",
   #      title = "95% CI overlaps 0?")
 dev.off()
}
