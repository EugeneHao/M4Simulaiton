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
source("bootstrapfunction.R")
source("estpredbootfunsRevInf.R")
source("predictMsampsLsamps.R")


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
SIGMA = diag(c( sqrt(0.5) , 1))%*%matrix(c(1,.4,.4,1), nrow = 2)%*%diag(c(sqrt(0.5), 1))   #variance of random effect


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
  
ratEBPs <-    c()
covEBPs <-   c()

ymuEBPINDs  <-   c()
zmuEBPINDs <-   c()
ratINDs <-   c()
covEBPINDs  <-   c()

paramEM1s <- c()
paramEM1FIs <- c()


ymuEBPBIVM1s <- c()
zmuEBPBIVM1s <-   c()
ratEBPM1s <- c()

mhzebps  <- c(); mhyebps <- c(); mhratebps <- c()
mhzebpM1s  <- c(); mhyebpM1s <- c(); mhratebpM1s <-c()
mhzebpM1bs  <- c(); mhyebpM1bs <- c(); mhratebpM1bs <- c()
	
vhatparboots <-	c()

ydirs <- c()
zdirs <- c()
ratdirs <- c()

iter <- 0
D <- 100
Nis <- rep(c(100, 200, 300, 400 ), each = 25)
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
B <- 20

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
coefw <- lm(log(w[sampindex]) ~ sampledata$z + sampledata$y  + as.factor(sampledata$group), data = sampledata)$coef

  predictINIT <- predictMCM1LS(sampledata,param[1:6],  Sigma_est, fitteddata, data, 200, 200)
  
paramInf <- param 
paramInf[1] <- param[1] + coefw[3]*summary(paramInit[[5]])$sigma^2
paramInf[3] <- param[3] + coefw[2] 

  predictINITINF <- predictMCM1InfLS(sampledata, param[1:6], paramInf[1:6],  Sigma_est, fitteddata, data, 200, 200)


   ymuEBPINITs <- rbind(ymuEBPINITs, predictINIT[[2]]) 
   zmuEBPINITs <-   rbind( zmuEBPINITs , predictINIT[[1]])
   ratEBPINITs <-    rbind(ratEBPINITs, predictINIT[[3]]) 
   rat0EBPINITs <-    rbind(rat0EBPINITs, predictINIT[[7]]) 
 
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

  zpops <- rbind(zpops, tapply(data$z, data$group, mean))
  ypops <- rbind(ypops, tapply(data$y, data$group, mean))
  ratpops <- rbind(ratpops, tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean))
  rat0pops <- rbind(rat0pops, tapply(data$z*(1-data$y), data$group, mean)/tapply(1-data$y, data$group, mean))

  #covpops <- rbind(covpops, areacov(c(data$y, data$z), data$group))

  #####################################    Independent with X (EBP for independent case...)
#  glm1Ind <- lmer(z~x + (1|group), data = sampledata)
#  glm2Ind <- glmer(y~x + (1|group), data = sampledata, famil = binomial(link = "logit"))

#  paramInd <- c(fixef(glm2Ind), fixef(glm1Ind), summary(glm1Ind)$sigma^2, 0)
#  SigmaInd <- matrix(c(VarCorr( glm2Ind)[[1]][1,1], 0, 0, VarCorr( glm1Ind)[[1]][1,1]), 2,2)

#  predInd <- predictMCM1(sampledata, paramInd,  SigmaInd, fitteddata, data, 2000)

 # ymuEBPINDs <- rbind(ymuEBPINDs, predInd[[2]])
  #zmuEBPINDs<- rbind(zmuEBPINDs, predInd[[1]])
  #ratINDs <-  rbind(ratINDs, predInd[[3]])

  #####################################    
#  zdirs <- rbind(zdirs, tapply(sampledata$z, sampledata$group, mean))
 # ydirs <- rbind(ydirs, tapply(sampledata$y, sampledata$group, mean))
  #yzdir <-  tapply(sampledata$y*sampledata$z, sampledata$group, mean) 

  #ratdirs <- rbind(ratdirs, yzdir/ tapply(sampledata$y, sampledata$group, mean))


print(iter)

#########################################    if(iter%%20 == 0){save.image("SimRho.4Update4-22-2019.Rdata")}

if(iter == 100){break}

}



maxcnt <- iter

biasfun <- function(a, b){
apply(a - b,2,mean)
}
sdfun <- function(a, b){
apply(a - b,2,sd)
}

 
 
popestlist <- vector("list",4)
popestlist[[1]] <- zpops[1:maxcnt,]
popestlist[[2]] <- ypops[1:maxcnt,]
popestlist[[3]] <- ratpops[1:maxcnt,]
popestlist[[4]] <- rat0pops[1:maxcnt,]

 

initestlist <- vector("list",4)
initestlist[[1]] <- zmuEBPINITs[1:maxcnt,]
initestlist[[2]] <- ymuEBPINITs[1:maxcnt,]
initestlist[[3]] <- ratEBPINITs[1:maxcnt,]
initestlist[[4]] <- rat0EBPINITs[1:maxcnt,]


initinfestlist <- vector("list",4)
initinfestlist[[1]] <- zmuEBPINITInfs[1:maxcnt,]
initinfestlist[[2]] <- ymuEBPINITInfs[1:maxcnt,]
initinfestlist[[3]] <- ratEBPINITInfs[1:maxcnt,]
initinfestlist[[4]] <- rat0EBPINITInfs[1:maxcnt,]

biascompare <- c()
for(i in 1:4){
 	biascompare <- cbind(biascompare, 
	tapply(biasfun(initestlist[[i]], popestlist[[i]]), Nis, mean),
	tapply(biasfun(initinfestlist[[i]], popestlist[[i]]), Nis, mean)
	)
}

sdcompare <- c()
for(i in 1:4){
 	sdcompare <- cbind(sdcompare, 
	tapply(sdfun(initestlist[[i]], popestlist[[i]]), Nis, mean),
	tapply(sdfun(initinfestlist[[i]], popestlist[[i]]), Nis, mean)
	)
}


biascompare <- c()
for(i in 1:4){
 	biascompare <- cbind(biascompare, 
	 biasfun(initestlist[[i]], popestlist[[i]]),  
	 biasfun(initinfestlist[[i]], popestlist[[i]]) 
	)
}

sdcompare <- c()
for(i in 1:4){
 	sdcompare <- cbind(sdcompare, 
	 sdfun(initestlist[[i]], popestlist[[i]]) ,
	 sdfun(initinfestlist[[i]], popestlist[[i]]) 
	)
}

biascompare2 <- rbind(biascompare[,c(1,3,5,7)], biascompare[,c(2,4,6,8)])
sdcompare2 <- rbind(sdcompare[,c(1,3,5,7)], sdcompare[,c(2,4,6,8)])

tall <- as.vector(biascompare2/sdcompare2*sqrt(maxcnt))

dftstat <- data.frame( Param = c(rep("Z",200), rep("Y",200), rep("R1",200), rep("R0",200)) , Nis = rep(rep(Nis,times = 2), times = 4) ,  Method = rep(rep( c("SRS", "Inf"), each = 100), times = 4),  tstat =tall )

ggplot(dftstat, aes(x = Nis, y = tstat, fill = Method)) +facet_grid( ~ Param )+ geom_boxplot( ) + theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"), legend.text = element_text(size = 16),legend.title = element_text(size = 16))

