initialvalueconsistent <- function(sddata){
  glm1 <- lmer(z ~ x + y + (1|group), data = sddata)
  ui_glm <- ranef(glm1)$group %>% unlist()
 # u_true <- randomeffect[,2]
 # b_true <- randomeffect[,1]
  s = as.data.frame(VarCorr(glm1))[2,4]  
  
###  glm2a <- glmer(y ~ x + (1|group), offset = z*fixef(glm1)[3]/s   , data= sddata , family = binomial(link = "logit"), nAGQ = 15)

  glm2a <- glmer(y ~ x + z + (1|group),   data= sddata , family = binomial(link = "logit"), nAGQ = 15)

  bi_glm <- ranef(glm2a)$group %>% unlist()

  glm2Ind <- glmer(y~x + (1|group), data = sddata, family = binomial(link = "logit"), nAGQ = 15)


  vhat1 <- (1/summary(glm1)$sigma^4)*vcov(glm1)[3,3]  + (fixef(glm1)[3]^2)/(summary(glm1)$sigma^8)*2*summary(glm1)$sigma^4/nrow(sddata)
  vhat2 <- vcov(glm2a)[3,3]
  rhohatinit <- (fixef(glm1)[3]/s/vhat1 + fixef(glm2a)[3]/vhat2)/(1/vhat1 + 1/vhat2)

  sigmaub1 <-  (VarCorr(glm2Ind)$group[[1]] - rhohatinit^2*as.data.frame(VarCorr(glm1))[1,4] - as.data.frame(VarCorr(glm2a))[4])/(2*rhohatinit)

  cormod1 <- as.vector(sigmaub1/sqrt( as.data.frame(VarCorr(glm2a))[4])/sqrt(as.data.frame(VarCorr(glm1))[1,4]))[[1]]
  cormod2 <- cor(cbind(bi_glm,ui_glm))[1,2]

  if(cormod1 < -1 | cormod1 > 1 | is.na(cormod1)){cormod1 <- cormod2}
        

  #initial value from glmm model:   new initial values!
  param = c(fixef(glm2a)[1:2],fixef(glm1)[1:2],s,rhohatinit,
            as.data.frame(VarCorr(glm2a))[4] %>% unlist, as.data.frame(VarCorr(glm1))[1,4],
            (cor(cbind(bi_glm,ui_glm)))%>%as.vector %>%"["(2)) %>% unlist 

  list(param, sigmaub1, cormod1, glm2a, glm1)

}



MCEMEst <- function(sampledata, steps, param, Sigma_est, glm2a, glm1, Msamps = 2000){

        iterMC <- 0
        storepars <- c()

        repeat{

                iterMC <- iterMC + 1
        #######################  
                D <- length(unique(group))
                randMC <- replicate(Msamps, mvrnorm(D, mu = c(0,0), Sigma = Sigma_est) )  #n * 2 matrix
  
  
                ###  Update Sigma:
                b_yest <- param[1:2]
                b_zest <- param[3:4]
                sigma2_est <- param[5]
                rho_est <- param[6]  

                #outMC1 <- sapply(1:100, UpdateUBMean, randMC,sampledata, b_yest, b_zest, sigma2_est, rho_est, D)
                #UpdateLLMean(otherparm, randMC = randMC, data1 = sampledata, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
                outMC1  <- parSapply(cl = cl, X = 1:Msamps,  UpdateUBMean, randMC, sampledata, b_yest, b_zest, sigma2_est, rho_est, D)

                f <- outMC1[1:D,]
                fb <- outMC1[(D+1):(2*D),]
                fu <- outMC1[(2*D+1):(3*D),]
                fb2 <- outMC1[(3*D+1):(4*D),]
                fu2 <- outMC1[(4*D+1):(5*D),]
                fub <- outMC1[(5*D+1):(6*D),]

                
                den <- apply(f,1,mean)
                bmean <- apply(fb, 1, mean)/den
                umean <- apply(fu, 1, mean)/den
                b2mean <-  apply(fb2, 1, mean)/den
                u2mean <- apply(fu2, 1, mean)/den
                bumean <- apply(fub, 1, mean)/den
                llhood <- sum(log(den))
                SigUpdate <- matrix(c(mean(b2mean), mean(bumean), mean(bumean), mean(u2mean)), nrow = 2, byrow= TRUE)
                otherparm <- c(b_yest, b_zest, sigma2_est, rho_est)

                lbs <- c( summary( glm2a )$coefficients[1:2,1] - 5*summary(glm2a)$coefficients[1:2,2],  summary(glm1)$coefficients[1:2,1] -  5*summary(glm1)$coefficients[1:2,2], 0.0001, ( summary(glm1)$coefficients[3,1] -  5*summary(glm1)$coefficients[3,2])/sigma2_est )
                ubs <- c( summary( glm2a )$coefficients[1:2,1] + 5*summary(glm2a)$coefficients[1:2,2],  summary(glm1)$coefficients[1:2,1] + 5*summary(glm1)$coefficients[1:2,2], 5*sigma2_est, (summary(glm1)$coefficients[3,1] + 5*summary(glm1)$coefficients[3,2] )/sigma2_est)
                #if(iterMC == 1){
                timeStart <- Sys.time()
                parUpdate <- optimParallel(otherparm , fn = UpdateLLMean, parallel = list(cl = cl), randMC = randMC, data1 = sampledata, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D , method = "L-BFGS-B", control = list(maxit = 10, factr = 10^-3) , lower = lbs,upper = ubs) 

                ##### parUpdate <- optim(otherparm , fn = UpdateLLMean, randMC = randMC, data1 = data1, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D, control = list(maxit = 5))# method = "L-BFGS-B", control = list(trace = 1,maxit = 10, factr = 10^-3), lower = lbs,upper = ubs) 
                timeEnd <- Sys.time()

            reldiff <- abs( c(c(parUpdate$par, as.vector(SigUpdate)) - c(param[1:6], as.vector(Sigma_est)))) / abs( c(param[1:6], as.vector(Sigma_est)))

                corsigu <- cov2cor(SigUpdate)
                corsige <- cov2cor(Sigma_est)   
            reldiff2 <- abs( (c(parUpdate$par, corsigu[1,2], SigUpdate[1,1], SigUpdate[2,2]  ) - c(param[1:6], corsige[1,2], Sigma_est[1,1], Sigma_est[2,2] ) ) / abs(  c(param[1:6], corsige[1,2], Sigma_est[1,1], Sigma_est[2,2] ) ) + 0.01)

                param <- parUpdate$par

                #}

                Sigma_est <- SigUpdate

                #U <- jacobian(func = UpdateLLMean, x = otherparm,   randMC = randMC, data1 = data1, b_yest=b_yest, b_zest=b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
                #H <- hessian(func = UpdateLLMean, x = otherparm,   randMC = randMC, data1 = data1, b_yest=b_yest, b_zest=b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
           
                storepars <- rbind(storepars, c(param, as.vector(SigUpdate)))

                print(paste("MC", reldiff2, iterMC))
                print(paste("Log Lik", parUpdate$value, iterMC))

                if(max(reldiff2) < 0.01 | iterMC == steps){break}

        }

        list(SigUpdate, param, storepars)
}





MCEMEstDF <- function(sampledata, steps, param, Sigma_est, glm2a, glm1, Msamps = 2000){

        iterMC <- 0
        storepars <- c()

        repeat{

                iterMC <- iterMC + 1
        #######################  
                D <- length(unique(group))
                randMC <- replicate(Msamps, mvrnorm(D, mu = c(0,0), Sigma = Sigma_est) )  #n * 2 matrix
  
  
                ###  Update Sigma:
                b_yest <- param[1:2]
                b_zest <- param[3:4]
                sigma2_est <- param[5]
                rho_est <- param[6]  

                #outMC1 <- sapply(1:100, UpdateUBMean, randMC,sampledata, b_yest, b_zest, sigma2_est, rho_est, D)
                #UpdateLLMean(otherparm, randMC = randMC, data1 = sampledata, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
                outMC1  <- parSapply(cl = cl, X = 1:Msamps,  UpdateUBMean, randMC, sampledata, b_yest, b_zest, sigma2_est, rho_est, D)

                f <- outMC1[1:D,]
                fb <- outMC1[(D+1):(2*D),]
                fu <- outMC1[(2*D+1):(3*D),]
                fb2 <- outMC1[(3*D+1):(4*D),]
                fu2 <- outMC1[(4*D+1):(5*D),]
                fub <- outMC1[(5*D+1):(6*D),]

                
                den <- apply(f,1,mean)
                bmean <- apply(fb, 1, mean)/den
                umean <- apply(fu, 1, mean)/den
                b2mean <-  apply(fb2, 1, mean)/den
                u2mean <- apply(fu2, 1, mean)/den
                bumean <- apply(fub, 1, mean)/den
                llhood <- sum(log(den))
                SigUpdate <- matrix(c(mean(b2mean), mean(bumean), mean(bumean), mean(u2mean)), nrow = 2, byrow= TRUE)*D/(D-5)
                otherparm <- c(b_yest, b_zest, sigma2_est, rho_est)

                lbs <- c( summary( glm2a )$coefficients[1:2,1] - 5*summary(glm2a)$coefficients[1:2,2],  summary(glm1)$coefficients[1:2,1] -  5*summary(glm1)$coefficients[1:2,2], 0.0001, ( summary(glm1)$coefficients[3,1] -  5*summary(glm1)$coefficients[3,2])/sigma2_est )
                ubs <- c( summary( glm2a )$coefficients[1:2,1] + 5*summary(glm2a)$coefficients[1:2,2],  summary(glm1)$coefficients[1:2,1] + 5*summary(glm1)$coefficients[1:2,2], 5*sigma2_est, (summary(glm1)$coefficients[3,1] + 5*summary(glm1)$coefficients[3,2] )/sigma2_est)
                #if(iterMC == 1){
                timeStart <- Sys.time()
                parUpdate <- optimParallel(otherparm , fn = UpdateLLMean, parallel = list(cl = cl), randMC = randMC, data1 = sampledata, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D , method = "L-BFGS-B", control = list(maxit = 10, factr = 10^-3) , lower = lbs,upper = ubs) 

                ##### parUpdate <- optim(otherparm , fn = UpdateLLMean, randMC = randMC, data1 = data1, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D, control = list(maxit = 5))# method = "L-BFGS-B", control = list(trace = 1,maxit = 10, factr = 10^-3), lower = lbs,upper = ubs) 
                timeEnd <- Sys.time()

            reldiff <- abs( c(c(parUpdate$par, as.vector(SigUpdate)) - c(param[1:6], as.vector(Sigma_est)))) / abs( c(param[1:6], as.vector(Sigma_est)))

                corsigu <- cov2cor(SigUpdate)
                corsige <- cov2cor(Sigma_est)   
            reldiff2 <- abs( (c(parUpdate$par, corsigu[1,2], SigUpdate[1,1], SigUpdate[2,2]  ) - c(param[1:6], corsige[1,2], Sigma_est[1,1], Sigma_est[2,2] ) ) / abs(  c(param[1:6], corsige[1,2], Sigma_est[1,1], Sigma_est[2,2] ) ) + 0.01)

                param <- parUpdate$par

                #}

                Sigma_est <- SigUpdate

                #U <- jacobian(func = UpdateLLMean, x = otherparm,   randMC = randMC, data1 = data1, b_yest=b_yest, b_zest=b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
                #H <- hessian(func = UpdateLLMean, x = otherparm,   randMC = randMC, data1 = data1, b_yest=b_yest, b_zest=b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
           
                storepars <- rbind(storepars, c(param, as.vector(SigUpdate)))

                print(paste("MC", reldiff2, iterMC))
                print(paste("Log Lik", parUpdate$value, iterMC))

                if(max(reldiff2) < 0.01 | iterMC == steps){break}

        }

        list(SigUpdate, param, storepars)
}




predictMCM1 <- function(data1, param, Sigma_est, fitteddata, data, Msamps){
        
        randMC <- replicate(Msamps, mvrnorm(D, mu = c(0,0), Sigma = Sigma_est) )  #n * 2 matrix
  
#       outmean <- parSapply(cl, X = 1:Msamps, ComputeCondMeanZYFull,   randMC, data1, fitteddata,  param[c(1,2)], param[c(3,4)], param[5], param[6], D)
#       outmeanbar <- apply(outmean, 1, mean)
#       zmeanNS <- outmeanbar[(D+1):(2*D)]/outmeanbar[1:D]
#       zmeanS <- tapply(data1$z, data1$group, mean)
#       ymeanNS <- outmeanbar[(2*D+1):(3*D)]/outmeanbar[1:D]
#       ymeanS <- tapply(data1$y, data1$group, mean)
#       nis <- table(data1$group)
#       Nis <- table(data$group)
#       zmu <- nis/Nis*zmeanS + (Nis - nis)/Nis*zmeanNS/(Nis-nis)
#       ymu <- nis/Nis*ymeanS + (Nis - nis)/Nis*ymeanNS/(Nis-nis)
#       xrat <- cbind(1, data$x, 1)
#     xratfit <- as.vector(xrat%*%fixef(glm1))
#       ratBV <- tapply(xratfit, data$group, mean) + ranef(glm1)[[1]][,1]

#       msebivZ <- mean((zmu - tapply(data$z, data$group, mean))^2)
#       mseobsZ <- mean((zmeanS - tapply(data$z, data$group, mean))^2)

#       msebivY <- mean((ymu - tapply(data$y, data$group, mean))^2)
#       mseobsY <- mean((ymeanS - tapply(data$y, data$group, mean))^2)

#       ratpop <-  tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean)

        fnum <- parSapply(cl, X = 1:dim(randMC)[3], NumForImpSampFun, randMC, data1,  param[c(1,2)], param[c(3,4)], param[5], param[6], D)
        probImp <- fnum/apply(fnum, 1, sum)

        YZgen <- replicate(Msamps, SimulateYZ(data, randMC, param[c(1,2)], param[c(3,4)], param[5], param[6], probImp, D))

        Ygen <- YZgen[1:N,]; Ygen[data1$id,] <- data1$y
        Zgen <- YZgen[(N+1):(2*N),]; Zgen[data1$id,] <- data1$z

        #ymuEBPBIV <- tapply(apply(Ygen, 1, mean), data$group, mean)
        #zmuEBPBIV <- tapply(apply(Zgen, 1, mean), data$group, mean)
        #ratEBP <-  tapply(apply(Ygen*Zgen, 1, mean), data$group, mean)/ymuEBPBIV
        #covYZEBPs <- apply(rbind(Ygen, Zgen), 2, areacov, data$group)
        #covEBP <- apply(covYZEBPs, 1, mean)
        
#      zmuEBPBIVM1 <- tapply(apply(Zgen, 1, mean), data$group, var)
#      ymuEBPBIVM1 <- tapply(apply(Ygen, 1, mean), data$group, var)
#       Rdel <- Ygen*Zgen - 
#       ratEBPM1 <-  tapply(apply(Ygen*Zgen, 1, mean), data$group, mean)/ymuEBPBIV

        predebp <- computeEBPgroupmean(Ygen, Zgen, data$group)


        predebp

}



predictMCM2 <- function(data1, param, Sigma_est, fitteddata, data, Msamps){
        
        randMC <- replicate(Msamps, mvrnorm(D, mu = c(0,0), Sigma = Sigma_est) )  #n * 2 matrix
  
#       outmean <- parSapply(cl, X = 1:Msamps, ComputeCondMeanZYFull,   randMC, data1, fitteddata,  param[c(1,2)], param[c(3,4)], param[5], param[6], D)
#       outmeanbar <- apply(outmean, 1, mean)
#       zmeanNS <- outmeanbar[(D+1):(2*D)]/outmeanbar[1:D]
#       zmeanS <- tapply(data1$z, data1$group, mean)
#       ymeanNS <- outmeanbar[(2*D+1):(3*D)]/outmeanbar[1:D]
#       ymeanS <- tapply(data1$y, data1$group, mean)
#       nis <- table(data1$group)
#       Nis <- table(data$group)
#       zmu <- nis/Nis*zmeanS + (Nis - nis)/Nis*zmeanNS/(Nis-nis)
#       ymu <- nis/Nis*ymeanS + (Nis - nis)/Nis*ymeanNS/(Nis-nis)
#       xrat <- cbind(1, data$x, 1)
#     xratfit <- as.vector(xrat%*%fixef(glm1))
#       ratBV <- tapply(xratfit, data$group, mean) + ranef(glm1)[[1]][,1]

#       msebivZ <- mean((zmu - tapply(data$z, data$group, mean))^2)
#       mseobsZ <- mean((zmeanS - tapply(data$z, data$group, mean))^2)

#       msebivY <- mean((ymu - tapply(data$y, data$group, mean))^2)
#       mseobsY <- mean((ymeanS - tapply(data$y, data$group, mean))^2)

#       ratpop <-  tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean)

        fnum <- parSapply(cl, X = 1:dim(randMC)[3], NumForImpSampFun, randMC, data1,  param[c(1,2)], param[c(3,4)], param[5], param[6], D)
        probImp <- fnum/apply(fnum, 1, sum)

        YZgen <- replicate(Msamps, SimulateYZ(data, randMC, param[c(1,2)], param[c(3,4)], param[5], param[6], probImp, D))

        Ygen <- YZgen[1:N,]; Ygen[data1$id,] <- data1$y
        Zgen <- YZgen[(N+1):(2*N),]; Zgen[data1$id,] <- data1$z

        #ymuEBPBIV <- tapply(apply(Ygen, 1, mean), data$group, mean)
        #zmuEBPBIV <- tapply(apply(Zgen, 1, mean), data$group, mean)
        #ratEBP <-  tapply(apply(Ygen*Zgen, 1, mean), data$group, mean)/ymuEBPBIV
        #covYZEBPs <- apply(rbind(Ygen, Zgen), 2, areacov, data$group)
        #covEBP <- apply(covYZEBPs, 1, mean)
        
#      zmuEBPBIVM1 <- tapply(apply(Zgen, 1, mean), data$group, var)
#      ymuEBPBIVM1 <- tapply(apply(Ygen, 1, mean), data$group, var)
#       Rdel <- Ygen*Zgen - 
#       ratEBPM1 <-  tapply(apply(Ygen*Zgen, 1, mean), data$group, mean)/ymuEBPBIV

        predebp1 <- computeEBPgroupmean(Ygen, Zgen, data$group )
        predebp2 <- computeEBPgroupmeanWeights(Ygen, Zgen, data$group, rep(1,nrow(data)))


        list(predebp1, predebp2)

}



computeEBPgroupmean <- function(Ygen, Zgen, group){
        YmeanL <- apply(Ygen, 2, function(x){ tapply(x, group, mean)} )
        ZmeanL <- apply(Zgen, 2, function(x){ tapply(x, group, mean)} )
        YZmeanL <- apply(Ygen*Zgen, 2, function(x){ tapply(x, group, mean)} )

        ymuebp <- apply(YmeanL,  1, mean)
        zmuebp <- apply(ZmeanL,  1, mean)
        YZmuebp <- apply(YZmeanL,  1, mean)

        ymuebpM1 <- apply(YmeanL,  1, var)
        zmuebpM1 <- apply(ZmeanL,  1, var)
        yzmuebpM1 <- apply(YZmeanL,  1, var)
        ratmuebp <- YZmuebp/ymuebp

      Rdel <- YZmeanL - ratmuebp*YmeanL 
        ratmuM1 <- apply(Rdel,1,var)/ymuebp^2

 list(zmuebp, ymuebp, ratmuebp, zmuebpM1, ymuebpM1, ratmuM1)

}

 genboot2 <- function(param, Sigma_est, sampledata,  data, steps, frac){
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
 list(sampledataboot, randomeffect, popdata, fitteddataboot)

}



computeEBPgroupmeanrat2 <- function(Ygen, Zgen, group){
        YmeanL <- apply(Ygen, 2, function(x){ tapply(x, group, mean)} )
        ZmeanL <- apply(Zgen, 2, function(x){ tapply(x, group, mean)} )
        YZmeanL <- apply(Ygen*Zgen, 2, function(x){ tapply(x, group, mean)} )

        ymuebp <- apply(YmeanL,  1, mean)
        zmuebp <- apply(ZmeanL,  1, mean)
        Ratmuebp <- apply(YZmeanL/YmeanL,  1, mean)

        ymuebpM1 <- apply(YmeanL,  1, var)
        zmuebpM1 <- apply(ZmeanL,  1, var)
	  ratmuebpM1 <- apply(YZmeanL/YmeanL,  1, var)


 list(zmuebp, ymuebp, Ratmuebp, zmuebpM1, ymuebpM1, ratmuebpM1)

}
 



estboot2 <- function(param, Sigma_est,   data, steps, sampledataboot, randomeffect, fitteddata, Msamps){

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

        Sigma_est <- diag(c(sqrt(param[7]), param[8]))%*%matrix(c(1,cormod1 , cormod1 , 1), 2, 2)%*%diag(c(sqrt(param[7]), param[8]))
        print(Sigma_est)
         paramEMEst1 <- MCEMEst(sampledataboot,steps, param, Sigma_est, glm2a, glm1, Msamps)

        predictBIV <- predictMCM1(sampledataboot, paramEMEst1[[2]],  paramEMEst1[[1]], fitteddata, data, Msamps)

   }else{

        predictBIV <- predictMCM1(sampledataboot,  param[1:6], matrix(c(param[7],0,0,param[8]), nrow = 2,byrow = TRUE), fitteddata, data, Msamps) 
   }

        zpop <-   tapply( data$z , data$group, mean)
        ypop  <-    tapply(data$y, data$group, mean)
        ratpop  <-  tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean)
   #    covpop  <-   areacovData(c(data$y, exp(data$z)), data$group, nrow(data))

        predictBIVM1 <- predictMCM1(sampledataboot,   paramold, Sigma_estold, fitteddata,  data, Msamps)

        
        
   list(predictBIV, zpop, ypop, ratpop,   predictBIVM1, as.vector(c(as.vector(paramEMEst1[[1]]), paramEMEst1[[2]])) )


}

estboot2InitVals <- function(param, Sigma_est,   data, steps, sampledataboot, randomeffect, fitteddata, Msamps){

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

        Sigma_est <- diag(c(sqrt(param[7]), param[8]))%*%matrix(c(1,cormod1 , cormod1 , 1), 2, 2)%*%diag(c(sqrt(param[7]), param[8]))
        print(Sigma_est)
      ##paramEMEst1 <- MCEMEst(sampledataboot,steps, param, Sigma_est, glm2a, glm1, Msamps)

        predictBIV <- predictMCM1(sampledataboot, param, Sigma_est, fitteddata, data, Msamps)

   }else{

        predictBIV <- predictMCM1(sampledataboot,  param[1:6], matrix(c(param[7],0,0,param[8]), nrow = 2,byrow = TRUE), fitteddata, data, Msamps) 
   }

        zpop <-   tapply( data$z , data$group, mean)
        ypop  <-    tapply(data$y, data$group, mean)
        ratpop  <-  tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean)
   #    covpop  <-   areacovData(c(data$y, exp(data$z)), data$group, nrow(data))

        predictBIVM1 <- predictMCM1(sampledataboot,   paramold, Sigma_estold, fitteddataboot,  data, Msamps)

        
        
   list(predictBIV, zpop, ypop, ratpop,   predictBIVM1, as.vector(c(as.vector(Sigma_est), param[1:6])) )


}


estbootSimp <- function(param, Sigma_est,   data, steps, sampledataboot, randomeffect, fitteddata, Msamps){

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
      paramEMEst1 <- MCEMEst(sampledataboot,steps, param, Sigma_est, glm2a, glm1, Msamps)

  #   predictBIV <- predictMCM1(sampledataboot, param, Sigma_est, fitteddata, data, Msamps)

    }else{
    		Sigma_est <- diag(c(sqrt(param[7]), sqrt(param[8]) ))%*%matrix(c(1,0 , 0, 1), 2, 2)%*%diag(c(sqrt(param[7]), sqrt(param[8]) ))
        print(Sigma_est)
     
	   paramEMEst1 <- list(Sigma_est, param)

   #      predictBIV <- predictMCM1(sampledataboot,  param[1:6], matrix(c(param[7],0,0,param[8]), nrow = 2,byrow = TRUE), fitteddata, data, Msamps) 
    }

 #       zpop <-   tapply( data$z , data$group, mean)
  #      ypop  <-    tapply(data$y, data$group, mean)
   #     ratpop  <-  tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean)
   #    covpop  <-   areacovData(c(data$y, exp(data$z)), data$group, nrow(data))

     #   predictBIVM1 <- predictMCM1(sampledataboot,   paramold, Sigma_estold, fitteddataboot,  data, Msamps)

        
        Sigma_est <- paramEMEst1[[1]]; param <- paramEMEst1[[2]][1:6]

     as.vector(c(as.vector(Sigma_est), param[1:6]))  


}



estbootSimpInitVals <- function(param, Sigma_est,   data, steps, sampledataboot, randomeffect, fitteddata, Msamps){

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

     as.vector(c(as.vector(Sigma_est), param[1:6]))  


}

 computeEBPgroupmeanWeights <- function(Ygen, Zgen, group, weights){
	denw <- tapply( weights, group, sum)
   	 YmeanL <- apply(Ygen, 2, function(x){ tapply(x*weights, group, sum)/denw} )
        ZmeanL <- apply(Zgen, 2, function(x){ tapply(x*weights, group, sum)/denw} )
        YZmeanL <- apply(Ygen*Zgen, 2, function(x){ tapply(x*weights, group, sum)/denw} )

    


        ymuebp <- apply(YmeanL,  1, mean)
        zmuebp <- apply(ZmeanL,  1, mean)
        YZmuebp <- apply(YZmeanL,  1, mean)
         ratmuebp <- YZmuebp/ymuebp

	   ymuebpM1 <-  sapply(1:D, function(i){ cwfun( Ygen[group == i,], weights[group == i])})
	   zmuebpM1 <-  sapply(1:D, function(i){ cwfun( Zgen[group == i,], weights[group == i])})

        names(ratmuebp) <- 1:length(table(group))
         

        Rdel <- Ygen*Zgen - ratmuebp[as.character(group)]*Ygen
        ratmuM1 <-  sapply(1:D, function(i){ cwfun( Rdel[group == i,], weights[group == i])})/ymuebp^2

 list(zmuebp, ymuebp, ratmuebp, zmuebpM1, ymuebpM1, ratmuM1)

}



cwfun <- function(Ygen, weights){
	  Cyy <- cov(t(Ygen))
	  sum(diag(weights)%*%Cyy%*%diag(weights))/sum(weights)^2
}




