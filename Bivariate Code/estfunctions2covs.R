

g2covs <- function(xy, xz, betay, betaz, bi, ui, sigma2, rho, group){
  xy = cbind(1,xy)
  xz <- cbind(1,xz)
  names(bi) <- 1:length(table(group)); names(ui) <- 1:length(table(group))
  bilong <- bi[as.character(group)]
  uilong <- ui[as.character(group)]
 # print(length(uilong))
  mu = xy%*%betay+bilong
  delta = xz%*%betaz+uilong
  return(rho*delta+mu+sigma2*rho^2/2)
}


Astar2covs <-  function(xy, xz, betay, betaz,   ui, sigma2, rho, y, group){
  xy <- cbind(1, xy)
  xz <- cbind(1, xz)
  names(ui) <- 1:length(table(group))
  uilong <- ui[as.character(group)]
  
  delta = xz%*%betaz+uilong
  return(delta + sigma2*rho*y)
}
	
initialvalueconsistent2covs <-  function(sddata, xy, xz ){

  glm1 <- lmer(z ~ xz + y + (1|group), data = sddata)
  ui_glm <- ranef(glm1)$group %>% unlist()
 # u_true <- randomeffect[,2]
 # b_true <- randomeffect[,1]
  s = as.data.frame(VarCorr(glm1))[2,4]  
  
###  glm2a <- glmer(y ~ x + (1|group), offset = z*fixef(glm1)[3]/s   , data= sddata , family = binomial(link = "logit"), nAGQ = 15)

  glm2a <- glmer(y ~ xy + z + (1|group),   data= sddata , family = binomial(link = "logit"), nAGQ = 15)

  bi_glm <- ranef(glm2a)$group %>% unlist()

  glm2Ind <- glmer(y~ cbind(xy, xz) + (1|group), data = sddata, family = binomial(link = "logit"), nAGQ = 15)

  pz <- length(fixef(glm1))
  py <- length(fixef(glm2a))
  vhat1 <- (1/summary(glm1)$sigma^4)*vcov(glm1)[pz,pz]  + (fixef(glm1)[pz]^2)/(summary(glm1)$sigma^8)*2*summary(glm1)$sigma^4/nrow(sddata)
  vhat2 <- vcov(glm2a)[py,py]
  rhohatinit <- (fixef(glm1)[pz]/s/vhat1 + fixef(glm2a)[py]/vhat2)/(1/vhat1 + 1/vhat2)

  sigmaub1 <-  (VarCorr(glm2Ind)$group[[1]] - rhohatinit^2*as.data.frame(VarCorr(glm1))[1,4] - as.data.frame(VarCorr(glm2a))[4])/(2*rhohatinit)

  cormod1 <- as.vector(sigmaub1/sqrt( as.data.frame(VarCorr(glm2a))[4])/sqrt(as.data.frame(VarCorr(glm1))[1,4]))[[1]]
  cormod2 <- cor(cbind(bi_glm,ui_glm))[1,2]

  if(cormod1 < -1 | cormod1 > 1 | is.na(cormod1)){cormod1 <- cormod2}
        

  #initial value from glmm model:   new initial values!
  param = c(fixef(glm2a)[-py] ,fixef(glm1)[-pz] ,s,rhohatinit,
            as.data.frame(VarCorr(glm2a))[4] %>% unlist, as.data.frame(VarCorr(glm1))[1,4],
            (cor(cbind(bi_glm,ui_glm)))%>%as.vector %>%"["(2)) %>% unlist 

  list(param, sigmaub1, cormod1, glm2a, glm1)

}

UpdateUBMean2covs <-  function(index,randMC ,data1, b_yest, b_zest, sigma2_est, rho_est, xy, xz  ) {
  bi = randMC[,1,index]
  ui = randMC[,2,index]
  astarvecest <- as.vector(Astar2covs( xy,  xz,  b_yest , b_zest,   ui, sigma2_est , rho_est, data1$y, data1$group ))
  normglog  <- -0.5/sigma2_est*(data1$z - astarvecest)^2 
  g1log <- g2covs( xy, xz, b_yest, b_zest, bi, ui, sigma2_est, rho_est, data1$group)
  gxlog <- g1log*data1$y - log(1 + exp(g1log))
  f <- exp(tapply(normglog + gxlog, data1$group, sum))#*sapply(1:D, function(i){ dmvnorm( c(bi[i],ui[i]), mean = c(0,0), sigma = SIGMACUR)})
  as.vector(cbind(1, bi, ui, bi^2, ui^2, bi*ui)*as.vector(f))
}


 MCEMEst2cov <-  function(sampledata, steps, param, Sigma_est, glm2a, glm1, py, pz, xy, xz, Msamps = 2000){

        iterMC <- 0
        storepars <- c()

        repeat{

                iterMC <- iterMC + 1
        #######################  
                D <- length(unique(group))
                randMC <- replicate(Msamps, mvrnorm(D, mu = c(0,0), Sigma = Sigma_est) )  #n * 2 matrix
  
             
                ###  Update Sigma:
                b_yest <- param[1:(py-1)]
                b_zest <- param[(py ):(py + pz-2)]
                sigma2_est <- param[(py + pz-1)]
                rho_est <- param[py + pz]  

                #outMC1 <- sapply(1:100, UpdateUBMean, randMC,sampledata, b_yest, b_zest, sigma2_est, rho_est, D)
                #UpdateLLMean(otherparm, randMC = randMC, data1 = sampledata, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
                outMC1  <- parSapply(cl = cl, X = 1:Msamps,  UpdateUBMean2covs, randMC, sampledata, b_yest, b_zest, sigma2_est, rho_est, xy, xz)

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

              #  lbs <- c( summary( glm2a )$coefficients[1:2,1] - 5*summary(glm2a)$coefficients[1:2,2],  summary(glm1)$coefficients[1:2,1] -  5*summary(glm1)$coefficients[1:2,2], 0.0001, ( summary(glm1)$coefficients[3,1] -  5*summary(glm1)$coefficients[3,2])/sigma2_est )
              #  ubs <- c( summary( glm2a )$coefficients[1:2,1] + 5*summary(glm2a)$coefficients[1:2,2],  summary(glm1)$coefficients[1:2,1] + 5*summary(glm1)$coefficients[1:2,2], 5*sigma2_est, (summary(glm1)$coefficients[3,1] + 5*summary(glm1)$coefficients[3,2] )/sigma2_est)
                #if(iterMC == 1){
                timeStart <- Sys.time()
                parUpdate <- optimParallel(otherparm , fn = UpdateLLMean2covs, parallel = list(cl = cl), randMC = randMC, data1 = sampledata, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D , xy = xy, xz = xz, method = "L-BFGS-B", control = list(maxit = 10, factr = 10^-3) ) 

                ##### parUpdate <- optim(otherparm , fn = UpdateLLMean, randMC = randMC, data1 = data1, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D, control = list(maxit = 5))# method = "L-BFGS-B", control = list(trace = 1,maxit = 10, factr = 10^-3), lower = lbs,upper = ubs) 
                timeEnd <- Sys.time()

            reldiff <- abs( c(c(parUpdate$par, as.vector(SigUpdate)) - c(param[1:(py + pz)], as.vector(Sigma_est)))) / abs( c(param[1:(py + pz)], as.vector(Sigma_est)))

                corsigu <- cov2cor(SigUpdate)
                corsige <- cov2cor(Sigma_est)   
            reldiff2 <- abs( (c(parUpdate$par, corsigu[1,2], SigUpdate[1,1], SigUpdate[2,2]  ) - c(param[1:6], corsige[1,2], Sigma_est[1,1], Sigma_est[2,2] ) ) / abs(  c(param[1:6], corsige[1,2], Sigma_est[1,1], Sigma_est[2,2] ) ) + 0.01)

                param <- parUpdate$par

                #}

                Sigma_est <- SigUpdate

                #U <- jacobian(func = UpdateLLMean, x = otherparm,   randMC = randMC, data1 = data1, b_yest=b_yest, b_zest=b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
                #H <- hessian(func = UpdateLLMean, x = otherparm,   randMC = randMC, data1 = data1, b_yest=b_yest, b_zest=b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D, py = py, pz = pz, xy = xy, xz = xz)
           
                storepars <- rbind(storepars, c(param, as.vector(SigUpdate)))

                print(paste("MC", reldiff2, iterMC))
                print(paste("Log Lik", parUpdate$value, iterMC))

                if(max(reldiff2) < 0.01 | iterMC == steps){break}

        }

        list(SigUpdate, param, storepars)
}



UpdateLLMean2covs <- function(otherparm, randMC, data1, b_yest, b_zest, sigma2_est, rho_est, D, xy, xz) {
  outMCLL <- sapply(1:dim(randMC)[3],LLFullCond1b2covs ,otherparm, randMC, data1, b_yest, b_zest, sigma2_est, rho_est, D, ncol(xy) + 2, ncol(xz) + 2, xy, xz)
  llmeanf  <- apply(outMCLL, 1, mean)
  llmean <- -sum(llmeanf[(D+1):(2*D)]/llmeanf[ 1:D])
  llmean
}

LLFullCond1b2covs <- function(index, otherparm, randMC, data1,   b_yest, b_zest, sigma2_est, rho_est, D, py, pz, xy, xz){
   bi = randMC[,1,index]
   ui = randMC[,2,index]
   zbetsig <- otherparm[(py):(py + pz-1)]

   ybet <- otherparm[1:(py-1)]
   rho <- otherparm[py + pz] 
   names(bi) <- names(ui) <- 1:D
   bilong <- bi[as.character(data1$group)] 
   uilong <- ui[as.character(data1$group)]
  
   astarvecestPAR <- as.vector(Astar2covs( xy,  xz , ybet, zbetsig[-length(zbetsig)],  ui, zbetsig[length(zbetsig)] , rho , data1$y, data1$group ))

   normglog  <- -0.5/zbetsig[3]*(data1$z - astarvecestPAR)^2 - log(zbetsig[3])/2
   g1log <-   g2covs( xy,  xz,ybet, zbetsig[-length(zbetsig)], bi, ui, zbetsig[length(zbetsig)], rho,  data1$group)
   gxlog <- as.vector(g1log*data1$y - log(1 + exp(g1log)))
   llarea <- as.vector(tapply(normglog   + gxlog, data1$group, sum)   )

   astarvecest <- as.vector(Astar2covs( xy, xz,  b_yest , b_zest, ui, sigma2_est , rho_est,data1$y, data1$group ))
   normglog  <- -0.5/sigma2_est*(data1$z - astarvecest)^2 
   g1log <-   g2covs( xy, xz,b_yest, b_zest, bi, ui, sigma2_est, rho_est, data1$group)

   gxlog <- as.vector(g1log*data1$y - log(1 + exp(g1log)))
   f <- exp(tapply(normglog + gxlog, data1$group, sum))
  # logdensarea <- as.vector(tapply(logzdens , data1$group, sum))
   as.vector(cbind(f, f*llarea))
}


NumForImpSampFun2covs <- function(index,   randMC, data1,   b_yest, b_zest, sigma2_est, rho_est, D, xy, xz){
   bi = randMC[,1,index]
   ui = randMC[,2,index]
   zbetsig <- c(b_zest, sigma2_est)
   ybet <- b_yest
   rho <- rho_est
  
  astarvecest <- as.vector(Astar2covs(xy , xz ,  b_yest , b_zest, ui, sigma2_est , rho_est, data1$y,  data1$group ))
   normglog  <- -0.5/sigma2_est*(data1$z - astarvecest)^2 
   g1log <-   g2covs(xy, xz,  b_yest, b_zest, bi, ui, sigma2_est, rho_est, data1$group )
   gxlog <- as.vector(g1log*data1$y - log(1 + exp(g1log)))
   f <- exp(tapply(normglog + gxlog, data1$group, sum))
  # logdensarea <- as.vector(tapply(logzdens , data1$group, sum))
  f
}


SimulateYZ2covs <- function(data, randMC, b_yest, b_zest, sigma2_est, rho_est, probImp,xy, xz, D){
        ubindexes <- apply(probImp, 1,  function(vec){ sample(1:length(vec), size = 1, prob = vec)})
        ub <- sapply(1:length(ubindexes), function(i){ randMC[i,,ubindexes[i]]})
        bi <- ub[1,]; names(bi) <- 1:D
        ui <- ub[2,]; names(ui) <- 1:D
        gfunAll <- g2covs(xy, xz,  b_yest, b_zest, bi, ui, sigma2_est, rho_est, data$group )
      piyall <- exp(gfunAll)/(1+exp(gfunAll))
        Ygen <- sapply(piyall, function(x){ rbinom(1, prob = x, size= 1)})       
      muzcond <- cbind(1,xz)%*%b_zest + sigma2_est*rho_est*Ygen +  ui[as.character(data$group)]
        Zgen <- rnorm(dim(data)[1], mean = muzcond, sd = sqrt(sigma2_est))
        c(Ygen, Zgen)
}

 
              
predictMCM12covs <- function(data1, param, Sigma_est, fitteddata, data,py, pz, xy, xz, xypop,xzpop, Msamps){
        
        randMC <- replicate(Msamps, mvrnorm(D, mu = c(0,0), Sigma = Sigma_est) )  #n * 2 matrix
 
        fnum <- parSapply(cl, X = 1:dim(randMC)[3], NumForImpSampFun2covs, randMC, data1,  param[1:(py-1)],param[(py):(py + pz-2)] , param[py + pz -1], param[py + pz],   D, xy, xz)
        probImp <- fnum/apply(fnum, 1, sum)

        YZgen <- replicate(Msamps, SimulateYZ2covs(data, randMC,  param[1:(py-1)],param[(py):(py + pz-2)] , param[py + pz -1], param[py + pz], probImp,xypop, xzpop, D))

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

