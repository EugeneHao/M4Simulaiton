
#####  Prediction functions:

predictMC <- function(data1, param, Sigma_est, fitteddata, data){
	
	randMC <- replicate(1000, mvrnorm(D, mu = c(0,0), Sigma = Sigma_est) )  #n * 2 matrix
  
	outmean <- parSapply(cl, X = 1:1000, ComputeCondMeanZYFull,   randMC, data1, fitteddata,  param[c(1,2)], param[c(3,4)], param[5], param[6], D)
	outmeanbar <- apply(outmean, 1, mean)
	zmeanNS <- outmeanbar[(D+1):(2*D)]/outmeanbar[1:D]
	zmeanS <- tapply(data1$z, data1$group, mean)
	ymeanNS <- outmeanbar[(2*D+1):(3*D)]/outmeanbar[1:D]
	ymeanS <- tapply(data1$y, data1$group, mean)
	nis <- table(data1$group)
	Nis <- table(data$group)
	zmu <- nis/Nis*zmeanS + (Nis - nis)/Nis*zmeanNS/(Nis-nis)
	ymu <- nis/Nis*ymeanS + (Nis - nis)/Nis*ymeanNS/(Nis-nis)
	ratBV <- zmu/ymu

#	msebivZ <- mean((zmu - tapply(data$z, data$group, mean))^2)
#	mseobsZ <- mean((zmeanS - tapply(data$z, data$group, mean))^2)

#	msebivY <- mean((ymu - tapply(data$y, data$group, mean))^2)
#	mseobsY <- mean((ymeanS - tapply(data$y, data$group, mean))^2)

#	ratpop <-  tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean)

	fnum <- parSapply(cl, X = 1:dim(randMC)[3], NumForImpSampFun, randMC, data1,  param[c(1,2)], param[c(3,4)], param[5], param[6], D)
	probImp <- fnum/apply(fnum, 1, sum)

	YZgen <- replicate(1000, SimulateYZ(data, randMC, param[c(1,2)], param[c(3,4)], param[5], param[6], probImp, D))

	Ygen <- YZgen[1:N,]; Ygen[data1$id,] <- data1$y
	Zgen <- YZgen[(N+1):(2*N),]; Zgen[data1$id,] <- data1$z

	ymuEBPBIV <- tapply(apply(Ygen, 1, mean), data$group, mean)
	zmuEBPBIV <- tapply(apply(Zgen, 1, mean), data$group, mean)
	ratEBP <-  tapply(apply(Ygen*Zgen, 1, mean), data$group, mean)/ymuEBPBIV
	covYZEBPs <- apply(rbind(Ygen, Zgen), 2, areacov, data$group)
	covEBP <- apply(covYZEBPs, 1, mean)

	list(zmu, ymu,  zmuEBPBIV,  ymuEBPBIV, ratEBP, covEBP)

}


MCEMEstLoop <- function(sampledata, steps, param, Sigma_est, glm2a, glm1){

	iterMC <- 0
	storepars <- c()

      loop <- loop1 <- 0

	repeat{

		iterMC <- iterMC + 1
	#######################  
 		D <- length(unique(group))
	  	randMC <- replicate(2000, mvrnorm(D, mu = c(0,0), Sigma = Sigma_est) )  #n * 2 matrix
  
  
		###  Update Sigma:
		b_yest <- param[1:2]
		b_zest <- param[3:4]
		sigma2_est <- param[5]
		rho_est <- param[6]  

		#outMC1 <- sapply(1:100, UpdateUBMean, randMC,sampledata, b_yest, b_zest, sigma2_est, rho_est, D)
		#UpdateLLMean(otherparm, randMC = randMC, data1 = sampledata, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
		outMC1  <- parSapply(cl = cl, X = 1:2000,  UpdateUBMean, randMC, sampledata, b_yest, b_zest, sigma2_est, rho_est, D)

		f <- outMC1[1:100,]
		fb <- outMC1[101:200,]
		fu <- outMC1[201:300,]
		fb2 <- outMC1[301:400,]
		fu2 <- outMC1[401:500,]
		fub <- outMC1[501:600,]
		den <- apply(f,1,mean)
		bmean <- apply(fb, 1, mean)/den
		umean <- apply(fu, 1, mean)/den
		b2mean <-  apply(fb2, 1, mean)/den
		u2mean <- apply(fu2, 1, mean)/den
		bumean <- apply(fub, 1, mean)/den

		llhood <- sum(log(den))

		SigUpdate <- matrix(c(mean(b2mean), mean(bumean), mean(bumean), mean(u2mean)), nrow = 2, byrow= TRUE) 
		otherparm <- c(b_yest, b_zest, sigma2_est, rho_est)

		lbs <- c( summary( glm2a )$coefficients[,1] - 5*summary(glm2a)$coefficients[,2],  summary(glm1)$coefficients[1:2,1] -  5*summary(glm1)$coefficients[1:2,2], 0.0001, max(c(summary(glm1)$coefficients[3,1] -  5*summary(glm1)$coefficients[3,2],-1)))
		ubs <- c( summary( glm2a )$coefficients[,1] + 5*summary(glm2a)$coefficients[,2],  summary(glm1)$coefficients[1:2,1] + 5*summary(glm1)$coefficients[1:2,2], 5*sigma2_est, min(c(summary(glm1)$coefficients[3,1] + 5*summary(glm1)$coefficients[3,2],1)))
		#if(iterMC == 1){

		if(loop == 0){

			parUpdate <- optimParallel(otherparm , fn = UpdateLLMean, parallel = list(cl = cl), randMC = randMC, data1 = sampledata, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D, method = "L-BFGS-B", control = list(maxit = 10), lower = lbs,upper = ubs) 
			##### parUpdate <- optim(otherparm , fn = UpdateLLMean, randMC = randMC, data1 = data1, b_yest = b_yest, b_zest = b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D, control = list(maxit = 5))# method = "L-BFGS-B", control = list(trace = 1,maxit = 10, factr = 10^-3), lower = lbs,upper = ubs) 
		}

            reldiff <- abs( c(c(parUpdate$par, as.vector(SigUpdate)) - c(param[1:6], as.vector(Sigma_est)))) / (abs( c(param[1:6], as.vector(Sigma_est))) + 0.01)

 		corubhatupdate <- cov2cor(SigUpdate)[1,2]
		reldiff2 <-  abs(  c(parUpdate$par,  corubhatupdate, SigUpdate[1,1], SigUpdate[2,2])  - c( param[1:6], cov2cor(Sigma_est)[1,2], Sigma_est[1,1], Sigma_est[2,2]  )   ) / (abs(  c( param[1:6], cov2cor(Sigma_est)[1,2], Sigma_est[1,1], Sigma_est[2,2]  ) ) + 0.01)

		if( max(reldiff2[-c(length(reldiff2), length(reldiff2)-1, length(reldiff2)-2)]) < 0.015 & loop1 == 0){
			loop <- 1; loop1 <- 1
		}

		if(loop == 1 & loop1 == 1){
			check <- max(reldiff2)
			if(check < 0.01){ loop <- 0 }
		}
		
		param <- parUpdate$par

		#}

		Sigma_est <- SigUpdate

		#U <- jacobian(func = UpdateLLMean, x = otherparm,   randMC = randMC, data1 = data1, b_yest=b_yest, b_zest=b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
		#H <- hessian(func = UpdateLLMean, x = otherparm,   randMC = randMC, data1 = data1, b_yest=b_yest, b_zest=b_zest, sigma2_est=sigma2_est, rho_est = rho_est, D = D)
	   
		storepars <- rbind(storepars, c(param, as.vector(SigUpdate)))

		print(paste("MC", reldiff2, iterMC))

		if(iterMC == 2){  storepars2 <- c(param, as.vector(SigUpdate))} 

		if( (max(reldiff2) < 0.015 & loop == 0) | iterMC == steps){break}

	}

	list(SigUpdate, param, storepars, storepars2)
}


