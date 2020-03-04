
##############################
### base functions: 
##############################
#g function  pi(y = 1) = exp(g)/(1+exp(g))
g <- function(x, bi, ui) {
  x = cbind(1,x)
  mu = x%*%beta_y+bi
  delta = x%*%beta_z+ui
  return(rho*delta+mu+sigma2*rho^2/2)
}


Astar <-function(x, y, bi, ui) {
  x = cbind(1,x)
  mu = x%*%beta_y+bi
  delta = x%*%beta_z+ui
  return(delta + sigma2*rho*y)
}


g_est <- function(x, bi, ui, by, bz, sigma12, rho1) {
  x = cbind(1,x)
  mu = x%*%by+bi
  delta = x%*%bz+ui
  return(rho1*delta+mu+sigma12*rho1^2/2)
}


Astar_est <-function(x, y, bi, ui, by, bz, sigma12, rho1) {
  x = cbind(1,x)
  mu = x%*%by+bi
  delta = x%*%bz+ui
  return(delta + sigma12*rho1*y)
}

#g_x : return sum of log(marginal probability of y_ij =1 )
g_x <- function(x, y, bi, ui,by, bz, sigma12, rho1) {
  g1 = g_est(x, bi, ui, by, bz, sigma12, rho1)
  return (sum(y*g1-log(1+exp(g1))))
}


#record the sum of the kernel of conditional density distribution of z_ij
normkernel <- function(x, y, z, bi, ui,by, bz, sigma12, rho1) {
  temp = Astar_est(x, y, bi, ui, by, bz, sigma12, rho1)
  return (-sum((z-temp)^2)/2/sigma12)
}

##################################
#### funtions for optimization: 
##################################
##estimate (bi, ui)   
rand_est_MC <- function(data1, gp_eft_area, by, bz, sigma12, rho1) {
  bi = gp_eft_area[1]
  ui = gp_eft_area[2]
  temp = data1 %>% group_by(group) %>% summarise(
    b = exp(normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ g_x(x,y,bi,ui,by, bz, sigma12, rho1)) * bi,
    u = exp(normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ g_x(x,y,bi,ui,by, bz, sigma12, rho1)) * ui,
    f = exp(normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ g_x(x,y,bi,ui,by, bz, sigma12, rho1)))
  return(temp[,2:4])   
}

##estimate (bi^2, ui^2, bi*ui)
rand2_est_MC <- function(data1, gp_eft_area, by, bz, sigma12, rho1) {
  bi = gp_eft_area[1]
  ui = gp_eft_area[2]
   temp = data1 %>% group_by(group) %>% summarise(
    b2 = exp(normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ g_x(x,y,bi,ui,by, bz, sigma12, rho1)) * bi^2,
    u2 = exp(normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ g_x(x,y,bi,ui,by, bz, sigma12, rho1)) * ui^2,
    bu = exp(normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ g_x(x,y,bi,ui,by, bz, sigma12, rho1)) * bi * ui,
    f  = exp(normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ g_x(x,y,bi,ui,by, bz, sigma12, rho1)))
  return(temp[,2:5])  
}


#used for optim in param_est_f_EM 
#param_p is the previous parameter: 
loglike_est_MC <- function(data1, gp_eft_area, by, bz, sigma12, rho1, param_p) {
  bi = gp_eft_area[1:100]
  ui = gp_eft_area[101:200]
  data1 <- data.frame(data1, bi = rep(bi, rep((1:10)*5, each = 10)),
                 ui = rep(ui, rep((1:10)*5, each = 10)))
  temp = data1 %>% group_by(group) %>% 
    summarise(
      l = normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ 
        g_x(x,y,bi,ui,by, bz, sigma12, rho1),
      f= exp(normkernel(x,y,z,bi,ui,param_p[1:2], param_p[3:4], param_p[5], param_p[6])+ 
               g_x(x,y,bi,ui,param_p[1:2], param_p[3:4], param_p[5], param_p[6])))
  return(temp[,2:3])  
}
##for optim  Note: the order is : by, bz, sigma1, rho1
## this function is used for optim by approximate bi, ui with expection
##data contians: x, y, z, E(b), E(u)
 
##############################################################################
#Next parts are new
##############################################################################
UpdateUBMean <- function(index,randMC ,data1, b_yest, b_zest, sigma2_est, rho_est, D  ) {
  bi = randMC[,1,index]
  ui = randMC[,2,index]
  astarvecest <- as.vector(Astar_est2(data1$x, data1$y,  bi, ui, b_yest , b_zest,sigma2_est , rho_est, data1$group, D))
  normglog  <- -0.5/sigma2_est*(data1$z - astarvecest)^2 
  g1log <-   g_est2(data1$x,   bi, ui, b_yest, b_zest, sigma2_est, rho_est, data1$group, D)
  gxlog <- g1log*data1$y - log(1 + exp(g1log))
  f <- exp(tapply(normglog + gxlog, data1$group, sum))#*sapply(1:D, function(i){ dmvnorm( c(bi[i],ui[i]), mean = c(0,0), sigma = SIGMACUR)})
  as.vector(cbind(1, bi, ui, bi^2, ui^2, bi*ui)*as.vector(f))
}


##################   Fix LLFullCond1b Function ..... 	
LLFullCond1b <- function(index, otherparm, randMC, data1,   b_yest, b_zest, sigma2_est, rho_est, D){
   bi = randMC[,1,index]
   ui = randMC[,2,index]
   zbetsig <- otherparm[c(3,4,5)]
   ybet <- otherparm[c(1,2)]
   rho <- otherparm[6]
   names(bi) <- names(ui) <- 1:D
   bilong <- bi[as.character(data1$group)] 
   uilong <- ui[as.character(data1$group)]
  
   astarvecestPAR <- as.vector(Astar_est2(data1$x, data1$y,  bi, ui, ybet , zbetsig[1:2], zbetsig[3] , rho , data1$group, D))
   normglog  <- -0.5/zbetsig[3]*(data1$z - astarvecestPAR)^2 - log(zbetsig[3])/2
   g1log <-   g_est2(data1$x,   bi, ui, ybet, zbetsig[1:2], zbetsig[3], rho , data1$group, D)
   gxlog <- as.vector(g1log*data1$y - log(1 + exp(g1log)))
   llarea <- as.vector(tapply(normglog   + gxlog, data1$group, sum)   )

   astarvecest <- as.vector(Astar_est2(data1$x, data1$y,  bi, ui, b_yest , b_zest,sigma2_est , rho_est, data1$group, D))
   normglog  <- -0.5/sigma2_est*(data1$z - astarvecest)^2 
   g1log <-   g_est2(data1$x,   bi, ui, b_yest, b_zest, sigma2_est, rho_est, data1$group, D)
   gxlog <- as.vector(g1log*data1$y - log(1 + exp(g1log)))
   f <- exp(tapply(normglog + gxlog, data1$group, sum))
  # logdensarea <- as.vector(tapply(logzdens , data1$group, sum))
   as.vector(cbind(f, f*llarea))
}


UpdateLLMean <- function(otherparm, randMC, data1, b_yest, b_zest, sigma2_est, rho_est, D) {
  outMCLL <- sapply(1:dim(randMC)[3],LLFullCond1b ,otherparm, randMC, data1, b_yest, b_zest, sigma2_est, rho_est, D)
  llmeanf  <- apply(outMCLL, 1, mean)
  llmean <- -sum(llmeanf[(D+1):(2*D)]/llmeanf[ 1:D])
  llmean
}
  


LLFullCond1bZ <- function(index, otherparm, randMC, data1,   b_yest, b_zest, sigma2_est, rho_est, D){
   bi = randMC[,1,index]
   ui = randMC[,2,index]
   zbetsig <- otherparm#[c(3,4,5)]
   #ybet <- otherparm[c(1,2)]
   names(bi) <- names(ui) <- 1:D
   bilong <- bi[as.character(data1$group)] 
   uilong <- ui[as.character(data1$group)]
   logzdens <- log(dnorm(data1$z, mean = cbind(1,data1$x)%*%zbetsig[1:2]+uilong, sd = sqrt(zbetsig[3])))
   logydens <- log(dbinom(data1$y, prob = exp( cbind(1,data1$x)%*%ybet+bilong)/(1 + exp(cbind(1,data1$x)%*%ybet+bilong)), size = 1))  
   astarvecest <- as.vector(Astar_est2(data1$x, data1$y,  bi, ui, b_yest , b_zest,sigma2_est , rho_est, data1$group, D))
   normglog  <- -0.5/sigma2_est*(data1$z - astarvecest)^2 
   g1log <-   g_est2(data1$x,   bi, ui, b_yest, b_zest, sigma2_est, rho_est, data1$group, D)
   gxlog <- as.vector(g1log*data1$y - log(1 + exp(g1log)))
   f <- exp(tapply(normglog + gxlog, data1$group, sum))
   logdensarea <- as.vector(tapply(logzdens , data1$group, sum))
   as.vector(cbind(f, f*logdensarea))
}

Astar_est2 <- function(x, y, bi, ui, by, bz, sigma12, rho1, group, D) {
  x = cbind(1,x)
  names(ui) <- 1:D
  uilong <- ui[as.character(group)]
  delta = x%*%bz+uilong
  return(delta + sigma12*rho1*y)
}

normkernel2 <- function(x, y, z, bi, ui,by, bz, sigma12, rho1) {
  temp = 
  return (-sum((z-temp)^2)/2/sigma12)
}

g_est2 <- function(x, bi, ui, by, bz, sigma12, rho1,group, D) {
  x = cbind(1,x)
  names(bi) <- 1:D
  bilong <- bi[as.character(group)]
  names(ui) <- 1:D
  uilong <- ui[as.character(group)]
  mu = x%*%by+bilong
  delta = x%*%bz+uilong
  return(rho1*delta+mu+sigma12*rho1^2/2)
}



#####################################################################################################

rand2_est_MC <- function(gp_eft_area, data1,  by, bz, sigma12, rho1) {
  bi = gp_eft_area[,1]
  ui = gp_eft_area[,2]
   temp = data1 %>% group_by(group) %>% summarise(
    b2 = exp(normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ g_x(x,y,bi,ui,by, bz, sigma12, rho1)) * bi^2,
    u2 = exp(normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ g_x(x,y,bi,ui,by, bz, sigma12, rho1)) * ui^2,
    bu = exp(normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ g_x(x,y,bi,ui,by, bz, sigma12, rho1)) * bi * ui,
    f  = exp(normkernel(x,y,z,bi,ui,by, bz, sigma12, rho1)+ g_x(x,y,bi,ui,by, bz, sigma12, rho1)))
  return(temp[,2:5])  
}

dbvnfun <- function(vhatdirvec, bumatvec, Sig, bmean){
	VhatCov <- Sig + diag(vhatdirvec)
      -0.5*log(det(VhatCov)) - 0.5*t(bumatvec - bmean)%*%solve(VhatCov)%*%(bumatvec - bmean)
}


bvnormopt <- function(corpar, vhatdir, varpar, bumat, D){
      bmean <- apply(bumat, 1, mean)
	Sig <- diag(c(sqrt(varpar[1]), varpar[2]))%*%matrix(c(1,corpar,corpar,1), nrow = 2)%*%diag(c(sqrt(varpar[1]), varpar[2]))   #variance of random effect
      llterms <- sapply(1:D, function(i){ dbvnfun( vhatdir[,i], bumat[,i], Sig, bmean)})
	-sum(llterms)
}




##################  Obtain Conditional Expectation...
ComputeCondMeanZYFull <- function(index,   randMC, data1, nonsampdata,  b_yest, b_zest, sigma2_est, rho_est, D){
   bi = randMC[,1,index]
   ui = randMC[,2,index]
   zbetsig <- c(b_zest, sigma2_est)
   ybet <- b_yest
   rho <- rho_est
   names(bi) <- names(ui) <- 1:D
   bilong <- bi[as.character(nonsampdata$group)] 
   uilong <- ui[as.character(nonsampdata$group)]


   g1log <- g_est2(nonsampdata$x,   bi, ui, ybet, zbetsig[1:2], zbetsig[3], rho , nonsampdata$group, D)
   egxY1 <- exp(g1log)/(1 + exp(g1log))
   eZ <- cbind(1, nonsampdata$x)%*%zbetsig[1:2] + uilong + egxY1*(rho*zbetsig[3])
   sumnonsampZ <- tapply(eZ, nonsampdata$group, sum)
   sumnonsampY <- tapply(egxY1, nonsampdata$group, sum)

#   astarvecestPAR <- as.vector(Astar_est2(data1$x, data1$y,  bi, ui, ybet , zbetsig[1:2], zbetsig[3] , rho , data1$group, D))
#   normglog  <- -0.5/zbetsig[3]*(data1$z - astarvecestPAR)^2 - log(zbetsig[3])/2
#   gxlog <- as.vector(g1log*data1$y - log(1 + exp(g1log)))
#   llarea <- as.vector(tapply(normglog   + gxlog, data1$group, sum)   )

   astarvecest <- as.vector(Astar_est2(data1$x, data1$y,  bi, ui, b_yest , b_zest,sigma2_est , rho_est, data1$group, D))
   normglog  <- -0.5/sigma2_est*(data1$z - astarvecest)^2 
   g1log <-   g_est2(data1$x,   bi, ui, b_yest, b_zest, sigma2_est, rho_est, data1$group, D)
   gxlog <- as.vector(g1log*data1$y - log(1 + exp(g1log)))
   f <- exp(tapply(normglog + gxlog, data1$group, sum))
  # logdensarea <- as.vector(tapply(logzdens , data1$group, sum))
   as.vector(cbind(f, f*sumnonsampZ, f*sumnonsampY))
}
 
################    Implement EBP procedure for the bivariate model.....

NumForImpSampFun <- function(index,   randMC, data1,   b_yest, b_zest, sigma2_est, rho_est, D){
   bi = randMC[,1,index]
   ui = randMC[,2,index]
   zbetsig <- c(b_zest, sigma2_est)
   ybet <- b_yest
   rho <- rho_est
   names(bi) <- names(ui) <- 1:D
   bilong <- bi[as.character(data1$group)] 
   uilong <- ui[as.character(data1$group)]

  astarvecest <- as.vector(Astar_est2(data1$x, data1$y,  bi, ui, b_yest , b_zest,sigma2_est , rho_est, data1$group, D))
   normglog  <- -0.5/sigma2_est*(data1$z - astarvecest)^2 
   g1log <-   g_est2(data1$x,   bi, ui, b_yest, b_zest, sigma2_est, rho_est, data1$group, D)
   gxlog <- as.vector(g1log*data1$y - log(1 + exp(g1log)))
   f <- exp(tapply(normglog + gxlog, data1$group, sum))
  # logdensarea <- as.vector(tapply(logzdens , data1$group, sum))
  f
}
 
SimulateYZ <- function(data,randMC, b_yest, b_zest, sigma2_est, rho_est, probImp, D){
	ubindexes <- apply(probImp, 1,  function(vec){ sample(1:length(vec), size = 1, prob = vec)})
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


areacov <- function(yzvec,group){
 	yvec <- yzvec[1:N]
	zvec <- yzvec[(N+1):(2*N)]
	sapply(unique(group), function(g){ cov(yvec[group == g], zvec[group == g])})
}


GenZCondIndFun <- function( data, glm1Ind, nis){
	ranefuind <- ranef(glm1Ind)[[1]][,1]
	names(ranefuind) <- 1:D
	ranefuindlongPOP <- ranefuind[as.character(data$group)]
	mucondpop <- as.vector(cbind(1,  data$x)%*%fixef(glm1Ind) + ranefuindlongPOP )
	vcondu  <- data.frame(VarCorr(glm1Ind))[1,4]*data.frame(VarCorr(glm1Ind))[2,4]/nis/(data.frame(VarCorr(glm1Ind))[1,4]+data.frame(VarCorr(glm1Ind))[2,4]/nis)
	ugen <- rnorm(D, mean = 0, sd = sqrt(vcondu))  	
	names(ugen) <- 1:D
	egen <- rnorm(N, mean = 0, sd = sqrt(data.frame(VarCorr(glm1Ind))[2,4]))		
	mucondpop + ugen[as.character(data$group)] + egen
}
	




