rm(list = ls(all = TRUE))

library(MASS)
library(dplyr)
library(tidyr)
library(parallel)
library(lme4)
library(ggplot2)
library(gridExtra)
library(optimParallel)
library(sampling)

source("bivariateSAEFunctions.R")
source("paramestfuns.R")
source("predictMCFun.R")
source("predictMsampsLsamps.R")
source("bootstrapfunction.R")
source("estpredbootfunsRevInf.R")
source("estpredbootfunsRev.R")
source("bootinffuns.R")


computeEBPgroupmean <- function(Ygen, Zgen, group ){
        YmeanL <- apply(Ygen, 2, function(x){ tapply(x, group, mean)} )
        ZmeanL <- apply(Zgen, 2, function(x){ tapply(x, group, mean)} )
        YZmeanL <- apply(Ygen*Zgen, 2, function(x){ tapply(x, group, mean)} )
        Y0ZmeanL <- apply((1 - Ygen)*Zgen, 2, function(x){ tapply(x, group, mean)} )

        ymuebp <- apply(YmeanL,  1, mean)
        zmuebp <- apply(ZmeanL,  1, mean)
        YZmuebp <- apply(YZmeanL,  1, mean)
        Y0Zmuebp <- apply(Y0ZmeanL,  1, mean)

        ymuebpM1 <- apply(YmeanL,  1, var)
        zmuebpM1 <- apply(ZmeanL,  1, var)
        yzmuebpM1 <- apply(YZmeanL,  1, var)
        ratmuebp <- YZmuebp/ymuebp
       ratmuebp0 <- Y0Zmuebp/(1-ymuebp)

      Rdel <- YZmeanL - ratmuebp*YmeanL 
        ratmuM1 <- apply(Rdel,1,var)/ymuebp^2
      Rdel0 <- Y0ZmeanL - ratmuebp0*(1-YmeanL )
 rat0muM1 <- apply(Rdel0,1,var)/((1-ymuebp)^2)

 list(zmuebp, ymuebp, ratmuebp, zmuebpM1, ymuebpM1, ratmuM1, ratmuebp0 , rat0muM1)

}



rho = 0.5
beta_y = c(-2, 0.5)
beta_z = c(0, 1)
sigma2 = 1
###########SIGMA = diag(c( sqrt(0.5) , 1))%*%matrix(c(1,.4,.4,1), nrow = 2)%*%diag(c(sqrt(0.5), 1))   #variance of random effect
SIGMA <- matrix(c(0.25, 0.05, 0.05, 0.0625), nrow = 2, byrow = TRUE)

paramInfoInit <- c()
paramUpdates <- c()
SigUpdates <- c()

#####  Store All:
zmus <- c()
ymus <- c()

zmuInds <- c()
ymuInds <- c()

zpops <- c()
ypops <- c()
ratpops <- c()
rat0pops <- c()
covpops <- c()
ymuEBPBIVs <- c()
zmuEBPBIVs <-   c()
ymuBIVs <- c()
zmuBIVs <-c()
ratBIVs <- c()
rat0BIVs <- c()
   
ratEBPs <-    c()
rat0EBPs <-    c()

covEBPs <-   c()

ymuEBPINDs  <-   c()
zmuEBPINDs <-   c()
ratINDs <-   c()
rat0INDs <- c()
covEBPINDs  <-   c()

paramEM1s <- c()
paramEM1FIs <- c()


ymuEBPBIVM1s <- c()
zmuEBPBIVM1s <-   c()
ratEBPM1s <- c()
rat0EBPM1s <- c()

mhzebps  <- c(); mhyebps <- c(); mhratebps <- c(); mhrat0ebps <- c()
mhzebpM1s  <- c(); mhyebpM1s <- c(); mhratebpM1s <-c(); mhrat0ebpM1s <- c()
mhzebpM1bs  <- c(); mhyebpM1bs <- c(); mhratebpM1bs <- c(); mhrat0ebpM1bs <- c()
mhzebpM2s  <- c();  mhyebpM2s <-c();  mhratebpM2s <-c();  mhrat0ebpM2s <- c()
mhzebpM2bs  <-c();  mhyebpM2bs <-c();  mhratebpM2bs <- c(); mhrat0ebpM2bs <- c()

vhatparboots <-	c()

ydirs <- c()
zdirs <- c()
ratdirs <- c()

iter <- 0
D <- 20
Nis <- rep(c(100, 200, 300, 400 ), each = 5)
Nisall <- rep(Nis, Nis)
group <- rep(1:D, Nis)
N <- length(group)

   ymuEBPINITs <- c()
   zmuEBPINITs <- c()
   ratEBPINITs <-  c()
   rat0EBPINITs <-  c()

   ymuEBPINITInfs <- c()
   zmuEBPINITInfs <- c()
   ratEBPINITInfs <- c() 
   rat0EBPINITInfs <- c() 
  
library("parallel")
cl <- makeCluster(getOption ("cl.cores", 2)) 
clusterExport(cl, varlist = list("Astar_est2", "g_est2", "LLFullCond1b"))

mhzmus <- c(); mhymus <- c(); mhzebps  <- c(); mhyebps <- c()
mhrats <-  c(); mhcovs <- c() 
B <- 2

repeat{


iter <- iter + 1

###  Miniature simulaton study for initial values:
  x =  runif(N, 1.5, 3.5)

  #generate data: 
  randomeffect <- mvrnorm(n =D, mu = c(0,0), Sigma = SIGMA)
  b = rep(randomeffect[,1],times = Nis)
  u = rep(randomeffect[,2],times = Nis)
 
  g1 = g(x, b, u)
  pi_y = exp(g1)/(1+exp(g1))    #pi(y = 1) 
  ####################    pi_y <- exp(cbind(1,x)%*%beta_y +b)/(1+ exp(cbind(1,x)%*%beta_y + b))
  y <- sapply(pi_y, function(x){ rbinom(1, prob = x, size= 1)})

  #sample z data
  z = rnorm(N, mean = Astar(x, y, b, u), sd = sqrt(sigma2))
  
  #combine all information 
  data <- data.frame(x,y,z,group,b,u,pi = pi_y,g1)
  data <- cbind(id = 1:dim(data)[1], data)

##### Define selection probabilities:
delta <-  rnorm(N,0,1)
delta[delta > 2] <- 2
delta[delta < -2] <- -2
piij1 <- exp( -z/3 + y/2+ delta/5 )
piij2 <- tapply(piij1, group, sum)[as.character(group)]
names(Nis) <- 1:D
piij <- (0.05*Nis[as.character(group)])*piij1/piij2

sampind  <- unlist(sapply(1:D, function(i){ UPsystematic(piij[group ==i])})) 
sampindex <- which(sampind == 1)

sampledata <- data[sampindex,]

fitteddata <- data[-sampindex,]

  paramInit <- initialvalueconsistent(sampledata)
  param <- paramInit[[1]]
  sigmaub1 <- paramInit[[2]]	  
  cormod1 <- paramInit[[3]]
 
  paramInfoInit <- rbind(paramInfoInit, c(param, param[9]*sqrt(param[7])*sqrt(param[8]), sigmaub1))

  glm2a <- paramInit[[4]]
  glm1 <- paramInit[[5]]

if(param[7] > 0 & param[8] > 0){
  Sigma_est <- diag(c(sqrt(param[7]), sqrt(param[8])))%*%matrix(c(1,cormod1 , cormod1 , 1), 2, 2)%*%diag(c(sqrt(param[7]), sqrt(param[8])))
}else{
  Sigma_est <- diag(c(sqrt(param[7]), sqrt(param[8])))%*%matrix(c(1, 0 , 0 , 1), 2, 2)%*%diag(c(sqrt(param[7]), sqrt(param[8])))
}

w <- 1/piij 
coefwmod <- lm(log(w[sampindex]) ~z + y  + as.factor(group), data = sampledata)
coefw <- coefwmod$coef

 # predictINIT <- predictMCM1(sampledata,param[1:6],  Sigma_est, fitteddata, data, 200)
  
paramInf <- param 
paramInf[1] <- param[1] + (summary(paramInit[[5]])$sigma^2)*coefw[3] 
paramInf[3] <- param[3] + coefw[2] 

  predictINITINF <- predictMCM1InfLS(sampledata, param[1:6], paramInf[1:6],  Sigma_est, fitteddata, data, 200, 100)


 #  ymuEBPINITs <- rbind(ymuEBPINITs, predictINIT[[2]]) 
 #  zmuEBPINITs <-   rbind( zmuEBPINITs , predictINIT[[1]])
 #  ratEBPINITs <-    rbind(ratEBPINITs, predictINIT[[3]]) 
 
   ymuEBPINITInfs <- rbind(ymuEBPINITInfs, predictINITINF[[2]]) 
   zmuEBPINITInfs <-   rbind( zmuEBPINITInfs , predictINITINF[[1]])
   ratEBPINITInfs <-    rbind(ratEBPINITInfs, predictINITINF[[3]]) 
   rat0EBPINITInfs <-    rbind(rat0EBPINITInfs, predictINITINF[[7]]) 
  
 # paramEMEst1 <- MCEMEst(sampledata,1, param, Sigma_est, glm2a, glm1, Msamps = 500)

 #  paramEMEstFI <- MCEMEstFI(sampledata, 5, param, Sigma_est, glm2a, glm1)

 #paramEM1s <- rbind(paramEM1s, c(paramEMEst1[[2]], as.vector(paramEMEst1[[1]]), dim(paramEMEst1[[3]])[1]))
# paramEM1FIs <- rbind(paramEM1FIs, c(paramEMEstFI[[2]], as.vector(paramEMEstFI[[1]]), dim(paramEMEstFI[[3]])[1]))
 
  paramEMEst1 <- list(Sigma_est, param[1:6])

 # predictBIV <- predictMCM1(sampledata, paramEMEst1[[2]],  paramEMEst1[[1]], fitteddata, data, 2000)
  
 # ymuEBPBIVs <- rbind(ymuEBPBIVs, predictBIV[[2]]) 
 # zmuEBPBIVs <-   rbind( zmuEBPBIVs , predictBIV[[1]])
 # ratEBPs <-    rbind(ratEBPs, predictBIV[[3]]) 

  ymuEBPBIVM1s <- rbind(ymuEBPBIVM1s, predictINITINF[[5]]) 
  zmuEBPBIVM1s <-   rbind( zmuEBPBIVM1s , predictINITINF[[4]])
  ratEBPM1s <-    rbind(ratEBPM1s, predictINITINF[[6]]) 
  rat0EBPM1s <-    rbind(rat0EBPM1s, predictINITINF[[8]]) 


  zpops <- rbind(zpops, tapply(data$z, data$group, mean))
  ypops <- rbind(ypops, tapply(data$y, data$group, mean))
  ratpops <- rbind(ratpops, tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean))
  rat0pops <- rbind(rat0pops, tapply(data$z*(1-data$y), data$group, mean)/tapply(1-data$y, data$group, mean))

  #covpops <- rbind(covpops, areacov(c(data$y, data$z), data$group))

   
  #####  Enter bootstrap:

    b <- 0
   zmuebpboot <- c(); ymuebpboot <- c(); ratboot <- c(); covboot <- c(); rat0boot <- c()
   zpopboot <- c(); ypopboot <- c(); ratpopboot <- c(); covpopboot <- c(); rat0popboot <- c()
   zmuebpbootM1 <- c(); ymuebpbootM1 <- c(); ratbootM1 <- c(); rat0bootM1 <- c()
   zbayesboot <- c(); ybayesboot <- c(); ratbayesboot <- c(); rat0bayesboot <- c()
   zbayesbootM1 <- c(); ybayesbootM1 <- c(); ratbayesbootM1 <- c(); rat0bayesbootM1 <- c()
   bootpars <- c()
   zmuebpbootM2 <- c()
   ymuebpbootM2<- c()
   ratebpbootM2 <- c()
   rat0ebpbootM2 <- c()
   zorigmat <- c()
   yorigmat <- c()
   ratorigmat <- c()
   rat0origmat <- c()

   zmuebpbootM2c <- c()
   ymuebpbootM2c<- c()
   ratebpbootM2c <- c()
   rat0ebpbootM2c <- c()
				
   zmuebpbootM2d <- c()
   ymuebpbootM2d<- c()
   ratebpbootM2d <- c()
   rat0ebpbootM2d <- c()

  repeat{
	b <- b + 1
	bootsamp <-genboot2Inf(param, coefwmod ,  Sigma_est, sampledata, data, 1, 0.05)
	bootout <- estbootSimpInitValsInf (param, Sigma_est, bootsamp[[3]],1,  bootsamp[[1]],bootsamp[[2]], bootsamp[[5]],500)
      paramb <- bootout[5:10]
      paramInfb <- paramb
      paramInfb[c(1,3)] <- paramb[c(1,3)] + bootout[c(12,11)]
      Sigma_estb <- matrix(bootout[c(1,2,3,4)],nrow=2,byrow = TRUE)
      predictBIVb <- predictMCM1InfLS(sampledata, paramb, paramInfb[1:6],  Sigma_estb, fitteddata, data, 200, 100)
      
		zmuebpbootM1 <- rbind(zmuebpbootM1, predictBIVb[[4]])
		ymuebpbootM1 <- rbind(ymuebpbootM1, predictBIVb[[5]])
		ratbootM1 <- rbind(ratbootM1, predictBIVb[[6]])
		rat0bootM1 <- rbind(rat0bootM1, predictBIVb[[8]])
			
		zmuebpbootM2 <- rbind(zmuebpbootM2, predictBIVb[[1]])
		ymuebpbootM2 <- rbind(ymuebpbootM2,  predictBIVb[[2]])
		ratebpbootM2 <- rbind(ratebpbootM2,  predictBIVb[[3]])
		rat0ebpbootM2 <- rbind(rat0ebpbootM2,  predictBIVb[[7]])

	#	zmuebpbootM2c <- rbind(zmuebpbootM2c, predictBIVc[[1]])
	#	ymuebpbootM2c <- rbind(ymuebpbootM2c,  predictBIVc[[2]])
#		ratebpbootM2c <- rbind(ratebpbootM2c,  predictBIVc[[3]])

#		zmuebpbootM2d <- rbind(zmuebpbootM2d, predictBIVd[[1]])
#		ymuebpbootM2d <- rbind(ymuebpbootM2d,  predictBIVd[[2]])
#		ratebpbootM2d <- rbind(ratebpbootM2d,  predictBIVd[[3]])

		zorigmat <- rbind(zorigmat, predictINITINF[[1]])
		yorigmat <- rbind(yorigmat, predictINITINF[[2]])
		ratorigmat <- rbind(ratorigmat,predictINITINF[[3]])
		rat0origmat <- rbind(rat0origmat,predictINITINF[[7]])

	if(b == B){break}

 }

	#mhzebp  <- msehatfun(zmuebpboot, zpopboot); mhyebp <- msehatfun(ymuebpboot, ypopboot); mhrat <- msehatfun(ratboot, ratpopboot);  
	#mhzebpM1  <- msehatfun(zbayesboot, zpopboot); mhyebpM1 <- msehatfun(ybayesboot, ypopboot); mhratM1 <- msehatfun(ratbayesboot, ratpopboot);  
	mhzebpM1b  <- apply(zmuebpbootM1,2,mean); mhyebpM1b <-apply(ymuebpbootM1,2,mean); mhratM1b <- apply(ratbootM1,2,mean);  mhrat0M1b <- apply(rat0bootM1,2,mean);  
#	mhzebpM2  <- msehatfun(zmuebpbootM2c, zmuebpbootM2d); mhyebpM2 <- msehatfun(ymuebpbootM2c,ymuebpbootM2d); mhratM2 <- msehatfun(ratebpbootM2c, ratebpbootM2d);  
	mhzebpM2b  <- msehatfun(zmuebpbootM2, zorigmat); mhyebpM2b <- msehatfun(ymuebpbootM2,yorigmat); mhratM2b <- msehatfun(ratebpbootM2, ratorigmat); mhrat0M2b <- msehatfun(rat0ebpbootM2, rat0origmat);  

      #mhzebps  <- rbind(mhzebps, mhzebp); mhyebps <- rbind(mhyebps, mhyebp); mhratebps <- rbind(mhratebps, mhrat)
      #mhzebpM1s  <- rbind(mhzebpM1s, mhzebpM1); mhyebpM1s <- rbind(mhyebpM1s, mhyebpM1); mhratebpM1s <- rbind(mhratebpM1s, mhratM1)
	mhzebpM1bs  <- rbind(mhzebpM1bs, mhzebpM1b); mhyebpM1bs <- rbind(mhyebpM1bs, mhyebpM1b); mhratebpM1bs <- rbind(mhratebpM1bs, mhratM1b);  mhrat0ebpM1bs <- rbind(mhrat0ebpM1bs, mhrat0M1b)

#	mhzebpM2s  <- rbind(mhzebpM2s, mhzebpM2); mhyebpM2s <- rbind(mhyebpM2s, mhyebpM2); mhratebpM2s <- rbind(mhratebpM2s, mhratM2)
	mhzebpM2bs  <- rbind(mhzebpM2bs, mhzebpM2b); mhyebpM2bs <- rbind(mhyebpM2bs, mhyebpM2b); mhratebpM2bs <- rbind(mhratebpM2bs, mhratM2b); mhrat0ebpM2bs <- rbind(mhrat0ebpM2bs, mhrat0M2b)
#
#	vhatparboots <- rbind(vhatparboots, apply(bootpars, 2, var))

#if(iter%%10 == 0){save.image("BootInfInitPart1SmallSig.Rdata")}
		
print(iter)

  if(iter%%25 == 0){save.image("InfBootRevRPart1.Rdata")}

if(iter == 2){break}

}


maxcnt <- iter

matbiasbc2 <- cbind(tapply(apply( ymuEBPBIVM1s[1:maxcnt,]^2/mhyebpM1bs[1:maxcnt,] +mhyebpM2bs[1:maxcnt,]  , 2, mean), Nis, mean)/tapply(apply(( ymuEBPINITInfs[1:maxcnt,] - ypops[1:maxcnt,] )^2, 2, mean), Nis, mean) ,
tapply(apply( zmuEBPBIVM1s[1:maxcnt,]^2/mhzebpM1bs[1:maxcnt,]+mhzebpM2bs[1:maxcnt,], 2, mean), Nis, mean)/tapply(apply(( zmuEBPINITInfs[1:maxcnt,] - zpops[1:maxcnt,])^2, 2, mean), Nis, mean) ,
tapply(apply( ratEBPM1s[1:maxcnt,]^2/mhratebpM1bs[1:maxcnt,] +mhratebpM2bs[1:maxcnt,], 2, mean), Nis, mean)/tapply(apply(( ratEBPINITInfs[1:maxcnt,] - ratpops[1:maxcnt,])^2, 2, mean), Nis, mean),
 tapply(apply( rat0EBPM1s[1:maxcnt,]^2/mhrat0ebpM1bs[1:maxcnt,] +mhrat0ebpM2bs[1:maxcnt,], 2, mean), Nis, mean)/tapply(apply(( rat0EBPINITInfs[1:maxcnt,] - rat0pops[1:maxcnt,])^2, 2, mean), Nis, mean)
)

matbiasbca <- cbind(tapply(apply( 2*ymuEBPBIVM1s[1:maxcnt,] - mhyebpM1bs[1:maxcnt,] +mhyebpM2bs[1:maxcnt,]  , 2, mean), Nis, mean)/tapply(apply(( ymuEBPINITInfs[1:maxcnt,] - ypops[1:maxcnt,] )^2, 2, mean), Nis, mean) ,
tapply(apply( 2*zmuEBPBIVM1s[1:maxcnt,] - mhzebpM1bs[1:maxcnt,]+mhzebpM2bs[1:maxcnt,], 2, mean), Nis, mean)/tapply(apply(( zmuEBPINITInfs[1:maxcnt,] - zpops[1:maxcnt,])^2, 2, mean), Nis, mean) ,
tapply(apply( 2*ratEBPM1s[1:maxcnt,] - mhratebpM1bs[1:maxcnt,] +mhratebpM2bs[1:maxcnt,], 2, mean), Nis, mean)/tapply(apply(( ratEBPINITInfs[1:maxcnt,] - ratpops[1:maxcnt,])^2, 2, mean), Nis, mean) ,
tapply(apply( 2*rat0EBPM1s[1:maxcnt,] - mhrat0ebpM1bs[1:maxcnt,] +mhrat0ebpM2bs[1:maxcnt,], 2, mean), Nis, mean)/tapply(apply(( rat0EBPINITInfs[1:maxcnt,] - rat0pops[1:maxcnt,])^2, 2, mean), Nis, mean) 
)

matbiasnobc <- cbind(tapply(apply( ymuEBPBIVM1s[1:maxcnt,] +mhyebpM2bs[1:maxcnt,]  , 2, mean), Nis, mean)/tapply(apply((  ymuEBPINITInfs[1:maxcnt,] - ypops[1:maxcnt,] )^2, 2, mean), Nis, mean) ,
tapply(apply( zmuEBPBIVM1s[1:maxcnt,]+mhzebpM2bs[1:maxcnt,], 2, mean), Nis, mean)/tapply(apply(( zmuEBPINITInfs[1:maxcnt,] - zpops[1:maxcnt,])^2, 2, mean), Nis, mean) ,
tapply(apply( ratEBPM1s[1:maxcnt,] +mhratebpM2bs[1:maxcnt,], 2, mean), Nis, mean)/tapply(apply(( ratEBPINITInfs[1:maxcnt,] - ratpops[1:maxcnt,])^2, 2, mean), Nis, mean) ,
tapply(apply( ratEBPM1s[1:maxcnt,] +mhrat0ebpM2bs[1:maxcnt,], 2, mean), Nis, mean)/tapply(apply(( rat0EBPINITInfs[1:maxcnt,] - rat0pops[1:maxcnt,])^2, 2, mean), Nis, mean) 
)

100*(cbind(matbiasnobc, matbiasbca, matbiasbc2) - 1)
  





