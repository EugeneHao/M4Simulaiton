

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


MCEMEstFI  <- function(sampledata, steps, param, Sigma_est, glm2a, glm1){

	iterMC <- 0
	storepars <- c()
      randMCFix <- replicate(100, mvrnorm(D, mu = c(0,0), Sigma = diag(2)) ) 
	#randMCFix1 <- randMCFix2 <- qnorm(seq(0.001, 0.999, by = 0.001))
	#randMCFix <- rbind(randMCFix1, randMCFix2)
	
	repeat{

		iterMC <- iterMC + 1
	#######################  
 		D <- length(unique(group))
            SVDSig <- eigen(Sigma_est)
 	      SqrtSig <- SVDSig$vectors%*%diag(sqrt(SVDSig$values))%*%t(SVDSig$vector)
	      randMC <- sapply(as.list(1:D), function(i){t( SqrtSig%*%randMCFix[i,,])} , simplify = "array")
  	      randMC <- aperm(randMC, perm = c(3,2,1))
  		###randMC <- randMCFix

		###  Update Sigma:
		b_yest <- param[1:2]
		b_zest <- param[3:4]
		sigma2_est <- param[5]
		rho_est <- param[6]  

		#outMC1 <- sapply(1:100, UpdateUBMean, randMC,sampledata, b_yest, b_zest, sigma2_est, rho_est, D, Sigma_est)
		#UpdateLLMean(otherparm, randMC = randMC, data1 = sampledata, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
		outMC1  <- parSapply(cl = cl, X = 1:100,  UpdateUBMean, randMC, sampledata, b_yest, b_zest, sigma2_est, rho_est, D)

		f <- outMC1[1:100,]
		fb <- outMC1[101:200,]
		fu <- outMC1[201:300,]
		fb2 <- outMC1[301:400,]
		fu2 <- outMC1[401:500,]
		fub <- outMC1[501:600,]

		bmean <- apply(fb, 1, mean)/apply(f, 1, mean)
		umean <- apply(fu, 1, mean)/apply(f, 1, mean)
		b2mean <-  apply(fb2, 1, mean)/apply(f, 1, mean)
		u2mean <- apply(fu2, 1, mean)/apply(f, 1, mean)
		bumean <- apply(fub, 1, mean)/apply(f, 1, mean)

		llhood <- sum(log(den))

		SigUpdate <- matrix(c(mean(b2mean), mean(bumean), mean(bumean), mean(u2mean)), nrow = 2, byrow= TRUE)
		otherparm <- c(b_yest, b_zest, sigma2_est, rho_est)

		lbs <- c( summary( glm2a )$coefficients[,1] - 5*summary(glm2a)$coefficients[,2],  summary(glm1)$coefficients[1:2,1] -  5*summary(glm1)$coefficients[1:2,2], 0.0001, max(c(summary(glm1)$coefficients[3,1] -  5*summary(glm1)$coefficients[3,2],-1)))
		ubs <- c( summary( glm2a )$coefficients[,1] + 5*summary(glm2a)$coefficients[,2],  summary(glm1)$coefficients[1:2,1] + 5*summary(glm1)$coefficients[1:2,2], 5*sigma2_est, min(c(summary(glm1)$coefficients[3,1] + 5*summary(glm1)$coefficients[3,2],1)))
		#if(iterMC == 1){

		parUpdate <- optimParallel(otherparm , fn = UpdateLLMean, parallel = list(cl = cl), randMC = randMC, data1 =sampledata, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D, method = "L-BFGS-B", control = list(maxit = 10), lower = lbs,upper = ubs) 
		##### parUpdate <- optim(otherparm , fn = UpdateLLMean, randMC = randMC, data1 = data1, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D, control = list(maxit = 5))# method = "L-BFGS-B", control = list(trace = 1,maxit = 10, factr = 10^-3), lower = lbs,upper = ubs) 


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

		print(paste("FI", reldiff2, iterMC))

		if(max(reldiff2) < 0.01 | iterMC == steps){break}

	}

	list(SigUpdate, param, storepars)
}





