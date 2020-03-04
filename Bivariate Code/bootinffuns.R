genboot2Inf <-  function(param, lwmodel, Sigma_est, sampledata,  data, steps, frac){
  ###  Update Sigma:
  b_yest <- param[1:2]
  b_zest <- param[3:4]
  sigma2_est <- param[5]
  rho_est <- param[6]  
  #generate data: 
  randomeffect <- mvrnorm(n =D, mu = c(0,0), Sigma = Sigma_est)
  b =  randomeffect[,1] 
  names(b) <- 1:D
  u = randomeffect[,2] 
  names(u) <- 1:D

  upop <-  u[as.character(data$group)]
  bpop <- b[as.character(data$group)]
 
  g1pop = g_est(data$x, bpop, upop, b_yest, b_zest, sigma2_est, rho_est)
  pi_ypop = exp(g1pop)/(1+exp(g1pop))    #pi(y = 1) 
  ####################    pi_y <- exp(cbind(1,x)%*%beta_y +b)/(1+ exp(cbind(1,x)%*%beta_y + b))
  # areafac <- rep(1:length(randomeffect[,1]), times = rep((1:ni)*D,each = D/ni))
  deltapop <- sapply(pi_ypop, function(p){ rbinom(1, 1, prob = p)})

  ypop <- sapply(pi_ypop, function(x){ rbinom(1, prob = x, size= 1)})
  #sample z data
  zpop = rnorm(nrow(data), mean = Astar_est(data$x, ypop, bpop, upop, b_yest, b_zest, sigma2_est, rho_est), sd = sqrt(sigma2_est))
   
  popdata <- data.frame(x = data$x, group = data$group, y = ypop, z = zpop)

  sampledataboot <- popdata[sampledata$id,]

  sampledataboot <- cbind(id = 1:dim(sampledataboot)[1], sampledataboot)
  fitteddataboot <- popdata[-sampledata$id,]  
  
  mmw <- cbind(1, sampledataboot$z, sampledataboot$y, model.matrix(~as.factor(sampledataboot$group)  )[,-1])
  lw <- as.vector(mmw%*%coef(lwmodel)) + rnorm(nrow(mmw), mean = 0, sd = summary(lwmodel)$sigma^2)

#  lw <- simulate(lwmodel)
  
 list(sampledataboot, randomeffect, popdata, fitteddataboot, lw)

}



estbootSimpInitValsInf <- function(param, Sigma_est,   data, steps, sampledataboot, randomeffect,  lw, Msamps){

###  Update Sigma:
   b_yest <- param[1:2]
   b_zest <- param[3:4]
   sigma2_est <- param[5]
   rho_est <- param[6]  

   Sigma_estold <- Sigma_est
   paramold <- param

   #time.start <- Sys.time()
   # sampledata <- data %>% group_by(group) %>%  sample_frac(frac)

   paramInit <- initialvalueconsistent(sampledataboot)
   param <- paramInit[[1]]
   sigmaub1 <- paramInit[[2]]      
   cormod1 <- paramInit[[3]]
 
   if(param[7] > 0 & param[8] > 0){

        glm2a <- paramInit[[4]]
        glm1 <- paramInit[[5]]

        Sigma_est <- diag(c(sqrt(param[7]),sqrt( param[8])))%*%matrix(c(1,cormod1 , cormod1 , 1), 2, 2)%*%diag(c(sqrt(param[7]), sqrt(param[8])))
        print(Sigma_est)
    #  paramEMEst1 <- MCEMEst(sampledataboot,steps, param, Sigma_est, glm2a, glm1, Msamps)

  #   predictBIV <- predictMCM1(sampledataboot, param, Sigma_est, fitteddata, data, Msamps)

    }else{
                Sigma_est <- diag(c(sqrt(param[7]), sqrt(param[8]) ))%*%matrix(c(1,0 , 0, 1), 2, 2)%*%diag(c(sqrt(param[7]), sqrt(param[8]) ))
        print(Sigma_est)
     
        

   #      predictBIV <- predictMCM1(sampledataboot,  param[1:6], matrix(c(param[7],0,0,param[8]), nrow = 2,byrow = TRUE), fitteddata, data, Msamps) 
    }

 #       zpop <-   tapply( data$z , data$group, mean)
  #      ypop  <-    tapply(data$y, data$group, mean)
   #     ratpop  <-  tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean)
   #    covpop  <-   areacovData(c(data$y, exp(data$z)), data$group, nrow(data))

     #   predictBIVM1 <- predictMCM1(sampledataboot,   paramold, Sigma_estold, fitteddataboot,  data, Msamps)

           paramEMEst1 <- list(Sigma_est, param)
      #  Sigma_est <- paramEMEst1[[1]]; param <- paramEMEst1[[2]][1:6]

     gammamod <- lm(lw ~sampledataboot$z + sampledataboot$y+as.factor(sampledataboot$group))

     as.vector(c(as.vector(Sigma_est), param[1:6], gammamod$coef[c(2,3)]))  


}





