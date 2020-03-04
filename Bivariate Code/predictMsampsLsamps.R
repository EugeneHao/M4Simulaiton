
 
SimulateYZL  <- function(l, data,randMC, b_yest, b_zest, sigma2_est, rho_est, probImp, D, Msamps, Lsamps){
        ubindexes <- apply(probImp, 1,  function(vec){sample( l + (0:(Msamps-1))*Lsamps, size = 1, prob = vec[ l + (0:(Msamps-1))*Lsamps])})
        ub <- sapply(1:length(ubindexes), function(i){ randMC[i,,ubindexes[i]]})
        bi <- ub[1,]; names(bi) <- 1:D
        ui <- ub[2,]; names(ui) <- 1:D
        gfunAll <- g_est2(data$x, bi , ui , b_yest, b_zest, sigma2_est, rho_est, data$group, D)
        piyall <- exp(gfunAll)/(1+exp(gfunAll))
        Ygen <- sapply(piyall, function(x){ rbinom(1, prob = x, size= 1)})       
        muzcond <- cbind(1, data$x)%*%b_zest + sigma2_est*rho_est*Ygen +  ui[as.character(data$group)]
        Zgen <- rnorm(dim(data)[1], mean = muzcond, sd = sqrt(sigma2_est))
        c(Ygen, Zgen)
}

predictMCM1LS <- function(data1, param, Sigma_est, fitteddata, data, Msamps, Lsamps){
        
        randMC <- replicate(Msamps*Lsamps, mvrnorm(D, mu = c(0,0), Sigma = Sigma_est) )  #n * 2 matrix
  
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

###################   YZgen <- replicate(Msamps, SimulateYZ(data, randMC, param[c(1,2)], param[c(3,4)], param[5], param[6], probImp, D))
	  YZgen <- sapply(1:Lsamps, SimulateYZL, data, randMC, param[c(1,2)], param[c(3,4)], param[5], param[6], probImp, D, Lsamps, Msamps)

        Ygen <- YZgen[1:N,]; Ygen[data1$id,] <- data1$y
        Zgen <- YZgen[(N+1):(2*N),]; Zgen[data1$id,] <- data1$z

        #ymuEBPBIV <- tapply(apply(Ygen, 1, mean), data$group, mean)
        #zmuEBPBIV <- tapply(apply(Zgen, 1, mean), data$group, mean)
        #ratEBP <-  tapply(apply(Ygen*Zgen, 1, mean), data$group, mean)/ymuEBPBIV
        #covYZEBPs <- apply(rbind(Ygen, Zgen), 2, areacov, data$group)
        #covEBP <- apply(covYZEBPs, 1, mean)
        
#       zmuEBPBIVM1 <- tapply(apply(Zgen, 1, mean), data$group, var)
#       ymuEBPBIVM1 <- tapply(apply(Ygen, 1, mean), data$group, var)
#       Rdel <- Ygen*Zgen - 
#       ratEBPM1 <-  tapply(apply(Ygen*Zgen, 1, mean), data$group, mean)/ymuEBPBIV

        predebp <- computeEBPgroupmean(Ygen, Zgen, data$group)


        predebp

}

SimulateYZL2covs <- function(l, data, randMC, b_yest, b_zest, sigma2_est, rho_est, probImp,xy, xz, D, Msamps, Lsamps ){
        ubindexes <- apply(probImp, 1,  function(vec){sample( l + (0:(Msamps-1))*Lsamps, size = 1, prob = vec[ l + (0:(Msamps-1))*Lsamps])})
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



predictMCM12covsLsampsMsamps <-  function(data1, param, Sigma_est, fitteddata, data,py, pz, xy, xz, xypop,xzpop, Msamps, Lsamps){
        
        randMC <- replicate(Msamps*Lsamps, mvrnorm(D, mu = c(0,0), Sigma = Sigma_est) )  #n * 2 matrix
 
        fnum <- parSapply(cl, X = 1:dim(randMC)[3], NumForImpSampFun2covs, randMC, data1,  param[1:(py-1)],param[(py):(py + pz-2)] , param[py + pz -1], param[py + pz],   D, xy, xz)
        probImp <- fnum/apply(fnum, 1, sum)

        YZgen <- sapply(1:Lsamps,  SimulateYZL2covs, data, randMC,  param[1:(py-1)],param[(py):(py + pz-2)] , param[py + pz -1], param[py + pz], probImp,xypop, xzpop, D, Msamps, Lsamps)

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





