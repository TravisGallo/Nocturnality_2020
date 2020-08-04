model {
  
  # priors
  
  # priors for random effect on intercept
  for(r in 1:ncity){
    b0[r,1] <- 0
    for(k in 2:ncat){
      b0[r,k] ~ dnorm(b_mu[k],tau[k])
    }
  }
  
  # first category set to reference (0)
  b_mu[1] <- 0
  b1[1] <- 0
  b2[1] <- 0
  b3[1] <- 0
  b4[1] <- 0
  b5[1] <- 0
  
  # Categorical Structured Lasso (CATS)
  # http://gph.is/10o9uZf
  # lamba for each variable
  for(b in 1:nbeta){
    b.lambda[b] ~ dgamma(0.001, 0.001)
  }
  
  # priors on beta parameters
  for(k in 2:ncat){
    # random effect hyper priors
    b_mu[k] ~ dnorm(0,2)
    tau[k] ~ dgamma(1,1)
    sigma[k] <- 1/sqrt(tau[k])
    
    # LASSO on variable
    # double exponential (Laplace) distribution
    # each beta gets its own lambda, but its the same for each linear predictor
    b1[k] ~ ddexp(0, b.lambda[1])
    b2[k] ~ ddexp(0, b.lambda[2])
    b3[k] ~ ddexp(0, b.lambda[3])
    b4[k] ~ ddexp(0, b.lambda[4])
    b5[k] ~ ddexp(0, b.lambda[5])
  }
  
  # softmax function
  for(i in 1:ndata){
    y[i] ~ dcat(explambda[1:ncat,i])
    for(k in 1:ncat){
      explambda[k,i] <- exp(log(avail_time[i,k]) + b0[city[i],k] + b1[k]*habitat[i] + 
                              b2[k]*imp[i] +  b3[k]*veg[i] + b4[k]*log_pd[i] + 
                              b5[k]*temp[i])
    }
  }
  
}