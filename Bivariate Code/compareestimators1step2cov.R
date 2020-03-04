rm(list = ls(all = TRUE))

library(MASS)
library(dplyr)
library(tidyr)
library(parallel)
library(lme4)
library(ggplot2)
library(gridExtra)
library(optimParallel)

######  Set working directory to folder containing supporting files below:

source("bivariateSAEFunctions.R")
source("paramestfuns.R")
source("predictMCFun.R")
source("bootstrapfunction.R")
source("estpredbootfunsRev.R")
source("estfunctions2covs.R")
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


rho = 2
beta_y = c(-4, 0.5)
beta_z = c(-1, 1, 1.5)
sigma2 = 0.75
SIGMA = diag(c( sqrt(0.5) ,1 ))%*%matrix(c(1,.4,.4,1), nrow = 2)%*%diag(c(sqrt(0.5), 1))   #variance of random effect


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
rat0EBPs <- c()
covEBPs <-   c()

ymuEBPINDs  <-   c()
zmuEBPINDs <-   c()
ratINDs <-   c()
rat0INDs <-   c()

covEBPINDs  <-   c()

paramEM1s <- c()
paramEM1FIs <- c()


ymuEBPBIVM1s <- c()
zmuEBPBIVM1s <-   c()
ratEBPM1s <- c()

mhzebps  <- c(); mhyebps <- c(); mhratebps <- c()
mhzebpM1s  <- c(); mhyebpM1s <- c(); mhratebpM1s <-c()
mhzebpM1bs  <- c(); mhyebpM1bs <- c(); mhratebpM1bs <- c()


ymuEBPINITs <- c()   
zmuEBPINITs <- c()   
ratEBPINITs <-  c()  
rat0EBPINITs <- c()  
 

	
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

library("parallel")
cl <- makeCluster(getOption ("cl.cores", 2)) 
clusterExport(cl, varlist = list("Astar_est2", "g_est2", "LLFullCond1b", "Astar2covs", "g2covs", "LLFullCond1b2covs"))

mhzmus <- c(); mhymus <- c(); mhzebps  <- c(); mhyebps <- c()
mhrats <-  c(); mhcovs <- c() 
B <- 20

repeat{


iter <- iter + 1

###  Miniature simulaton study for initial values:
  xy1 <-   runif(N, 1.5, 3.5)
  xz1 <-   runif(N, 1.5, 3.5)
  xz2 <-   (xz1 - mean(xz1))^2 - mean((xz1 - mean(xz1))^2)
  
  xypop <- as.matrix(xy1)
  xzpop <- cbind(xz1, xz2)
  

  #generate data: 
  randomeffect <- mvrnorm(n =D, mu = c(0,0), Sigma = SIGMA)
  b = rep(randomeffect[,1],times = Nis)
  u = rep(randomeffect[,2],times = Nis)
 
  g1 = g2covs(xy1,xzpop,beta_y, beta_z, randomeffect[,1], randomeffect[,2], sigma2, rho, group)
  pi_y = exp(g1)/(1+exp(g1))    #pi(y = 1) 
  ####################    pi_y <- exp(cbind(1,x)%*%beta_y +b)/(1+ exp(cbind(1,x)%*%beta_y + b))
  y <- sapply(pi_y, function(x){ rbinom(1, prob = x, size= 1)})

  #sample z data
  z = rnorm(N, mean = Astar2covs(xy1, xzpop, beta_y, beta_z, randomeffect[,2],   sigma2, rho, y, group), sd = sqrt(sigma2))
  
  #combine all information 
  data <- data.frame( y,z,group,b,u,pi = pi_y,g1, xy1, xz1, xz2)
  data <- cbind(id = 1:dim(data)[1], data)
  

#time.start <- Sys.time()
 sampledata <- data %>% group_by(group) %>%  sample_frac(0.05)
#time.end <- Sys.time()
  #unobversed data
  
 xysamp <- as.matrix(xypop[sampledata$id,])
 xzsamp <- xzpop[sampledata$id,]
 
 fitteddata <- data[-which(data$id %in% sampledata$id),]
  
 xyfit <- as.matrix(xypop[-sampledata$id,])
 xzfit <- xzpop[-sampledata$id,]
 
  paramInit <- initialvalueconsistent2covs(sampledata, xysamp, xzsamp)
  param <- paramInit[[1]]
  sigmaub1 <- paramInit[[2]]	  
  cormod1 <- paramInit[[3]]
 
  paramInfoInit <- rbind(paramInfoInit, c(param, param[9]*sqrt(param[7])*sqrt(param[8]), sigmaub1))

  glm2a <- paramInit[[4]]
  glm1 <- paramInit[[5]]

###  number of regression parameters, including the rho term
py <- ncol(xysamp) + 2
pz <- ncol(xzsamp) + 2

if(param[py + pz + 1] > 0 & param[py + pz + 2] > 0){
  Sigma_est <- diag(c(sqrt(param[py + pz + 1]), sqrt(param[py + pz + 2])))%*%matrix(c(1,cormod1 , cormod1 , 1), 2, 2)%*%diag(c(sqrt(param[py + pz + 1]), sqrt(param[py + pz + 2])))
}else{
  Sigma_est <- diag(c(sqrt(param[py + pz + 1]), sqrt(param[py + pz + 2])))%*%matrix(c(1, 0 , 0 , 1), 2, 2)%*%diag(c(sqrt(param[py + pz + 1]), sqrt(param[py + pz + 2])))
}

 # paramEMEst1 <- MCEMEst2cov(sampledata,1, param, Sigma_est, glm2a, glm1, 3,4, xysamp, xzsamp, Msamps = 200)
 # paramEMEst1 <-  list(Sigma_est, param)

 #  paramEMEstFI <- MCEMEstFI(sampledata, 5, param, Sigma_est, glm2a, glm1)

 # paramEM1s <- rbind(paramEM1s, c(paramEMEst1[[2]], as.vector(paramEMEst1[[1]]) ))
# paramEM1FIs <- rbind(paramEM1FIs, c(paramEMEstFI[[2]], as.vector(paramEMEstFI[[1]]), dim(paramEMEstFI[[3]])[1]))
 
  paramEMEst1 <- list(Sigma_est, param[1: (py + pz)])

  predictBIV <- predictMCM12covsLsampsMsamps(sampledata, paramEMEst1[[2]],  paramEMEst1[[1]], fitteddata, data, py, pz, xysamp, xzsamp, xypop,xzpop, 200, 200)
  
  ymuEBPBIVs <- rbind(ymuEBPBIVs, predictBIV[[2]]) 
  zmuEBPBIVs <-   rbind( zmuEBPBIVs , predictBIV[[1]])
  ratEBPs <-    rbind(ratEBPs, predictBIV[[3]]) 
  rat0EBPs <-    rbind(rat0EBPs, predictBIV[[7]]) 


  ymuEBPINITs <-  rbind(ymuEBPINITs, predictBIV[[2]]) 
  zmuEBPINITs <- rbind( zmuEBPINITs , predictBIV[[1]])
  ratEBPINITs <-  rbind(ratEBPINITs, predictBIV[[3]]) 
  rat0EBPINITs <- rbind(rat0EBPINITs, predictBIV[[7]]) 
 

  ymuEBPBIVM1s <- rbind(ymuEBPBIVM1s, predictBIV[[5]]) 
  zmuEBPBIVM1s <-   rbind( zmuEBPBIVM1s , predictBIV[[4]])
  ratEBPM1s <-    rbind(ratEBPM1s, predictBIV[[6]]) 


  zpops <- rbind(zpops, tapply(data$z, data$group, mean))
  ypops <- rbind(ypops, tapply(data$y, data$group, mean))
  ratpops <- rbind(ratpops, tapply(data$z*data$y, data$group, mean)/tapply(data$y, data$group, mean))
  rat0pops <- rbind(rat0pops, tapply(data$z*(1-data$y), data$group, mean)/tapply(1 - data$y, data$group, mean))

  #covpops <- rbind(covpops, areacov(c(data$y, data$z), data$group))

  #####################################    Independent with X (EBP for independent case...)
  glm1Ind <- lmer(z~xzsamp + (1|group), data = sampledata)
  glm2Ind <- glmer(y~xysamp + (1|group), data = sampledata, family = binomial(link = "logit"))

  paramInd <- c(fixef(glm2Ind), fixef(glm1Ind), summary(glm1Ind)$sigma^2, 0)
  SigmaInd <- matrix(c(VarCorr( glm2Ind)[[1]][1,1], 0, 0, VarCorr( glm1Ind)[[1]][1,1]), 2,2)

  predInd <- predictMCM12covsLsampsMsamps(sampledata, paramInd,  SigmaInd, fitteddata, data,py, pz, xysamp, xzsamp, xypop,xzpop, 200, 200)

  ymuEBPINDs <- rbind(ymuEBPINDs, predInd[[2]])
  zmuEBPINDs<- rbind(zmuEBPINDs, predInd[[1]])
  ratINDs <-  rbind(ratINDs, predInd[[3]])
  rat0INDs <-  rbind(rat0INDs, predInd[[7]])

  #####################################    
  zdirs <- rbind(zdirs, tapply(sampledata$z, sampledata$group, mean))
  ydirs <- rbind(ydirs, tapply(sampledata$y, sampledata$group, mean))
  yzdir <-  tapply(sampledata$y*sampledata$z, sampledata$group, mean) 

  ratdirs <- rbind(ratdirs, yzdir/ tapply(sampledata$y, sampledata$group, mean))


print(iter)


 if(iter%%20 == 0){save.image("StrongCovs.Rdata")}

if(iter == 100){break}

}



save.image("TwoCovStrongInit.Rdata")

 






