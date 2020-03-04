

### Enter bootstrap:

bootfun <- function(param, Sigma_est, frac){
	###  Update Sigma:
		b_yest <- param[1:2]
		b_zest <- param[3:4]
		sigma2_est <- param[5]
		rho_est <- param[6]  


 
 #generate data: 
  randomeffect <- mvrnorm(n =D, mu = c(0,0), Sigma = Sigma_est)
  b = rep(randomeffect[,1],times = rep((1:ni)*D,each = D/ni))
  u = rep(randomeffect[,2],times = rep((1:ni)*D,each = D/ni))
 
  g1 = g_est(x, b, u, b_yest, b_zest, sigma2_est, rho_est)
  pi_y = exp(g1)/(1+exp(g1))    #pi(y = 1) 
  ####################    pi_y <- exp(cbind(1,x)%*%beta_y +b)/(1+ exp(cbind(1,x)%*%beta_y + b))
  areafac <- rep(1:length(randomeffect[,1]), times = rep((1:ni)*D,each = D/ni))
  delta <- sapply(pi_y, function(p){ rbinom(1, 1, prob = p)})

  y <- sapply(pi_y, function(x){ rbinom(1, prob = x, size= 1)})

  #sample z data
  z = rnorm(N, mean = Astar_est(x, y, b, u, b_yest, b_zest, sigma2_est, rho_est), sd = sqrt(sigma2_est))
  
  #combine all information 
  data <- data.frame(x,y,z,group,b,u,pi = pi_y,g1)
  data <- cbind(id = 1:dim(data)[1], data)
  

#time.start <- Sys.time()
 sampledata <- data %>% group_by(group) %>%  sample_frac(frac)

   paramInit <- initialvalueconsistent(sampledata)
  param <- paramInit[[1]]
  sigmaub1 <- paramInit[[2]]	  
  cormod1 <- paramInit[[3]]
 


 glm2a <- paramInit[[4]]
  glm1 <- paramInit[[5]]

  Sigma_est <- diag(c(sqrt(param[7]), param[8]))%*%matrix(c(1,cormod1 , cormod1 , 1), 2, 2)%*%diag(c(sqrt(param[7]), param[8]))

   paramEMEst1 <- MCEMEst(sampledata,2, param, Sigma_est, glm2a, glm1)

   predictBIV <- predictMC(sampledata, paramEMEst1[[2]],  paramEMEst1[[1]], fitteddata, data)
  
   zpop <-   tapply(data$z, data$group, mean)
   ypop  <-    tapply(data$y, data$group, mean)
   ratpop  <-  tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean)
   covpop  <-   areacov(c(data$y, data$z), data$group)

  
	
   list(predictBIV, zpop, ypop, ratpop, covpop)


}

msehatfun <- function(estmat, popmat){
	apply( (estmat - popmat)^2, 2, mean)
}




### Enter bootstrap:

bootfun2 <- function(param, Sigma_est, frac){
	###  Update Sigma:
		b_yest <- param[1:2]
		b_zest <- param[3:4]
		sigma2_est <- param[5]
		rho_est <- param[6]  


 
 #generate data: 
  randomeffect <- mvrnorm(n =D, mu = c(0,0), Sigma = Sigma_est)
  b = rep(randomeffect[,1],times = rep((1:ni)*D,each = D/ni))
  u = rep(randomeffect[,2],times = rep((1:ni)*D,each = D/ni))
 
  g1 = g_est(x, b, u, b_yest, b_zest, sigma2_est, rho_est)
  pi_y = exp(g1)/(1+exp(g1))    #pi(y = 1) 
  ####################    pi_y <- exp(cbind(1,x)%*%beta_y +b)/(1+ exp(cbind(1,x)%*%beta_y + b))
  areafac <- rep(1:length(randomeffect[,1]), times = rep((1:ni)*D,each = D/ni))
  delta <- sapply(pi_y, function(p){ rbinom(1, 1, prob = p)})

  y <- sapply(pi_y, function(x){ rbinom(1, prob = x, size= 1)})

  #sample z data
  z = rnorm(N, mean = Astar_est(x, y, b, u, b_yest, b_zest, sigma2_est, rho_est), sd = sqrt(sigma2_est))
  
  #combine all information 
  data <- data.frame(x,y,z,group,b,u,pi = pi_y,g1)
  data <- cbind(id = 1:dim(data)[1], data)
  

#time.start <- Sys.time()
 sampledata <- data %>% group_by(group) %>%  sample_frac(frac)

   paramInit <- initialvalueconsistent(sampledata)
  param <- paramInit[[1]]
  sigmaub1 <- paramInit[[2]]	  
  cormod1 <- paramInit[[3]]
 


 glm2a <- paramInit[[4]]
  glm1 <- paramInit[[5]]

  Sigma_est <- diag(c(sqrt(param[7]), param[8]))%*%matrix(c(1,cormod1 , cormod1 , 1), 2, 2)%*%diag(c(sqrt(param[7]), param[8]))

   paramEMEst1 <- MCEMEst(sampledata,2, param, Sigma_est, glm2a, glm1)

   predictBIV <- predictMC(sampledata, paramEMEst1[[2]],  paramEMEst1[[1]], fitteddata, data)
  
   zpop <-   tapply(data$z, data$group, mean)
   ypop  <-    tapply(data$y, data$group, mean)
   ratpop  <-  tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean)
   covpop  <-   areacov(c(data$y, data$z), data$group)

  
	
   list(predictBIV, zpop, ypop, ratpop, covpop)


}

  
