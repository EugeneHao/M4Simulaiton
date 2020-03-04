# functions 
#("D:/PhD/Project2019/GauMultiModel")

## this function provides prob of each category for all data, also used for generating data of truedata and bootstrap
# delta0 : n * 1 matrix
# P0 : 2 * 1  matrix
# Ychoice : 3 * 2 matrix
# SigmaZ0 : a number \
#mu0 : n * 2 matrix
multiYprob <- function(P0, SigmaZ0, delta0, mu0) 
{
  temp = (delta0 %*% t(P0) + mu0) %*% t(Ychoice)    #n * 3 
  PSP =  SigmaZ0 * P0 %*% t(P0)/2
  rm(SigmaZ0)
  rm(delta0)
  rm(mu0)
  fixpart = lapply(1:3, FUN = function(x) t(Ychoice[x,]) %*% PSP %*% Ychoice[x,]) %>% unlist()
  prob = sweep(temp, STATS = fixpart, FUN =  "+", MARGIN = 2) %>% exp()
  rm(temp)
  rm(fixpart)
  sumprob = rowSums(prob)
  prob = sweep(prob, STATS = sumprob, FUN = "/", MARGIN = 1)
  return(prob)
}


rcppEBPonedraw <-function(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, deltanew_fix, munew_fix, weight, rand_matrix, m)
{
  onedraw = apply(weight, MARGIN = 1, FUN = function(x) sample(1:m, size = 1, prob = x))
  rm(weight)
  rand_one = sapply(1:D,  FUN = function(x) rand_matrix[x, onedraw[x], ]) %>% t()
  rm(rand_matrix)
  rm(onedraw)
  delta_new = deltanew_fix + rep(rand_one[,1], predictsize)
  mu_new = munew_fix + rep(rand_one[,2:3], rep(predictsize,2)) %>% matrix(ncol = 2)
  n = length(deltanew_fix)
  rand_vec = runif(n)
  rm(deltanew_fix)
  rm(munew_fix)
  Yprob_new = multiYprob(P_hat, SigmaZ_hat, delta_new, mu_new)
  
  rcppYZdraw(delta_new, mu_new, Yprob_new, P_hat, SigmaZ_hat, rand_vec) -> temp
  rm(delta_new)
  rm(mu_new)
  rm(Yprob_new)
  
  Z_draw = rnorm(n, mean = temp[,4], sd = sqrt(SigmaZ_hat))
  data.frame(fitting, Y1 = temp[,1], Y2 = temp[,2], Z = Z_draw, Yid = temp[,3] %>% as.factor) %>% 
    rbind(., data.frame(data)) %>% 
    group_by(group) %>% summarise(z1mean = mean(Z),
                                  y10mean = mean(1-Y1-Y2), y11mean = mean(Y1), y12mean = mean(Y2), 
                                  z1y10 = sum((1-Y1-Y2)*Z)/sum(1-Y1-Y2), 
                                  z1y11 = sum(Y1*Z)/sum(Y1), z1y12 = sum(Y2*Z)/sum(Y2)) %>% as.matrix() -> mean_draw
  return(mean_draw[,-1])
}


rcppEBPfun <- function(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat,  SigmaUV_hat, m, tt)
{
  
  PSP =  SigmaZ_hat *P_hat  %*% t(P_hat)/2
  probfix = sapply(1:3, FUN = function(x) t(Ychoice[x,]) %*% PSP %*% Ychoice[x,]) 
  
  aux = cbind(1, data$X2)
  delta_fix = aux %*% Bz_hat
  mu_fix = aux %*% By_hat
  Y = cbind(data$Y1, data$Y2)
  Z = data$Z
  A_fix = Y %*% P_hat * SigmaZ_hat - Z
  ylabel = data$Yid
  
  deltanew_fix = cbind(1,fitting$X2) %*% Bz_hat
  munew_fix = cbind(1,fitting$X2) %*% By_hat
  
  rand_matrix= rmvnorm(m*D, sigma = SigmaUV_hat) %>% array(., dim = c(D, m, 3))
  
  id = rep(1:D, samplesize)
  
  logdens = rcppgroupdens(P_hat, SigmaZ_hat, rand_matrix, as.vector(delta_fix), mu_fix, id, A_fix,
                          unlist(ylabel), probfix, m, Ychoice, D)
  dens = (logdens - median(logdens)) %>% round(3) %>% exp()
  weight = sweep(dens, STATS= rowSums(dens), FUN = "/" , MARGIN = 1)
  
  temp = sapply(1:tt, FUN = function(x) rcppEBPonedraw(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, 
                                                       deltanew_fix, munew_fix, weight, rand_matrix, m)) %>% 
    array(., dim = c(D, 7, tt))
  
  return(temp)
}

rcppSigmaUVfun <- function(data, SigmaZ_hat, P_hat, Bz_hat, By_hat, SigmaUV_hat, m)
{
  rand_origin = rmvnorm(m*D, sigma = SigmaUV_hat)
  rand_matrix  = rand_origin %>% array(., dim = c(D, m, 3))
  PSP =  SigmaZ_hat *P_hat  %*% t(P_hat)/2
  probfix = sapply(1:3, FUN = function(x) t(Ychoice[x,]) %*% PSP %*% Ychoice[x,]) 
  
  aux = cbind(1, data$X2)
  delta_fix = aux %*% Bz_hat
  mu_fix = aux %*% By_hat
  Y = cbind(data$Y1, data$Y2)
  ylabel = data$Yid
  Z = data$Z
  A_fix = Y %*% P_hat * SigmaZ_hat - Z
  id = rep(1:D, samplesize)
  rcppgroupdens(P_hat, SigmaZ_hat, rand_matrix, as.vector(delta_fix), mu_fix, id, A_fix,
                unlist(ylabel), probfix, m, Ychoice, D) -> logdens
  rm(rand_matrix)
  quant_logdens = apply(logdens, MARGIN = 1, FUN = function(x) quantile(x, 0.6)) 
  dens = sweep(logdens, MARGIN = 1, FUN = "-", STATS = quant_logdens) %>% exp()
  weight = sweep(dens, STATS= rowSums(dens), FUN = "/" , MARGIN = 1)/D
  
  return(rcppSigmaUVweightavg(rand_origin, weight))
}


###### Bootstrap part
#################################
## Simulate sample data for both true data and bootstrap
simuldata <- function(auxillary, Bz_hat, By_hat,  P_hat, SigmaZ_hat, SigmaUV_hat, groupsize)
{
  randef = rmvnorm(D, sigma = SigmaUV_hat) %>%
    rep(., rep(groupsize,3)) %>% matrix(ncol = 3) %>% as.data.frame() %>% 
    "names<-"(c("u1", "v11", "v12"))
  X = cbind(1, auxillary)
  
  delta_draw = (X %*% Bz_hat + randef[, 1]) %>% as.matrix()
  
  mu_draw = (X%*% By_hat + randef[, 2:3]) %>% as.matrix()
  
  
  group = rep(c(1:D), groupsize)   
  simudata <- data.frame(group, X, randef, delta_draw, mu_draw)
  Yprob = multiYprob(P_hat, SigmaZ_hat, delta_draw, mu_draw)
  
  Y = lapply(1:dim(Yprob)[1], FUN = function(x) sample(1:3, size = 1, prob = Yprob[x,])) %>% unlist()
  
  
  Yid = Y  %>% as.factor()
  Y = Ychoice[Y,]  %>% 'colnames<-'(paste("Y",  1:2, sep = ""))
  
  condZmean = delta_draw + SigmaZ_hat * Y %*% P_hat 
  
  Z = rnorm(sum(groupsize), mean = condZmean, sd = sqrt(SigmaZ_hat)) 
  return(data.frame(group, X2 = auxillary, Y, Z, Yid))
}

getgy_one <- function(P0, Z, mu0, ylabel, probfix) 
{
  temp = (Z %*% t(P0) + mu0) %*% t(Ychoice)  #n * 3
  prob = temp + probfix
  index = cbind(1:dim(prob)[1], ylabel)
  sumprob = exp(prob) %>% rowSums()  %>% log()
  return(prob[index] - sumprob)
}

Ydensity_onerand <- function(ylabel, Z, By_hat, P_hat, probfix,mu_fix, samplesize, V)
{
  n = sum(samplesize)
  mu = mu_fix + rep(V, each = n)
  gy = getgy_one(P_hat, Z, mu, ylabel, probfix) + 0.5
  data.frame(group = rep(1:D, samplesize), value= gy) -> temp
  temp %>% group_by(group) %>% summarise(m = sum(value)) %>% "["(,2) %>% unlist() %>% exp() %>% return()
}

##### for multinomial univariate model: 
multiYprob_uni <- function(mu0) 
{
  prob = (mu0 %*% t(Ychoice)) %>% exp    #n * 3 
  rm(mu0)
  sumprob = rowSums(prob)
  prob = sweep(prob, STATS = sumprob, FUN = "/", MARGIN = 1)
  return(prob)
}

## return the log density of Y in univariate case, return a n*m matrix
# mu0:  n * 2 * m array
getgy_uni <- function(mu0, ylabel, m) 
{
  prob = sapply(1:m, FUN = function(x) (mu0[,x,] %*% t(Ychoice))) %>% 
    array(., dim = c(dim(mu0)[1], 3, m)) #n * 3 * m
  index = cbind(rep(1:dim(prob)[1], m), rep(ylabel, m), rep(1:m, each = dim(prob)[1]))
  sumprob = apply(exp(prob), MARGIN = c(1,3), FUN = sum) %>% as.vector() %>% log()
  return(matrix(prob[index] - sumprob, ncol = m))
}

# return group log density of Y|X 
grouplogdensity_uni <- function(ylabel, By_hat, rand_sample, mu_fix, m)
{
  n = sum(samplesize)
  rand_array = rand_sample %>% rep(., rep(samplesize, m*2)) %>% array(., dim = c(n, m, 2))
  rm(rand_sample)
  mu = sweep(rand_array, MARGIN = c(1,3), FUN = "+", STATS = mu_fix)
  rm(mu_fix)
  rm(rand_array)
  gy = getgy_uni(mu, ylabel,  m)  
  rm(mu)
  data.frame(group = rep(1:D, samplesize), gy) -> temp
  result = sapply(1:D, FUN = function(x) temp[temp$group==x, -1] %>% colSums() %>% as.vector) %>% t() 
  return(result)
}

#for multivariate model:
EBPonedraw_uni <-function(data, fitting, By_hat, SigmaV_hat, munew_fix, weight, rand_matrix, m, Zdraw)
{
  onedraw = apply(weight, MARGIN = 1, FUN = function(x) sample(1:m, size = 1, prob = x))
  rm(weight)
  rand_one = sapply(1:D,  FUN = function(x) rand_matrix[x, onedraw[x], ]) %>% t()
  rm(rand_matrix)
  rm(onedraw)
  mu_new = munew_fix + rep(rand_one, rep(predictsize,2)) %>% matrix(ncol = 2)
  rm(munew_fix)
  Yprob_new = multiYprob_uni(mu_new)
  Y_draw = sapply(1:dim(Yprob_new)[1], FUN = function(x) sample(1:3, size = 1, prob = Yprob_new[x,]))
  n = dim(Yprob_new)[1]
  Yid = Y_draw %>% as.factor()
  
  Y_draw = Ychoice[Y_draw, ]
  rm(mu_new)
  rm(Yprob_new)
  data.frame(fitting, Y1 = Y_draw[,1], Y2 = Y_draw[,2], Z = Zdraw, Yid = Yid) %>% 
    rbind(., data.frame(data)) %>% 
    group_by(group) %>% summarise(z1mean = mean(Z),
                                  y10mean = mean(1-Y1-Y2), y11mean = mean(Y1), y12mean = mean(Y2), 
                                  z1y10 = sum((1-Y1-Y2)*Z)/sum(1-Y1-Y2), 
                                  z1y11 = sum(Y1*Z)/sum(Y1), z1y12 = sum(Y2*Z)/sum(Y2)) %>% as.matrix() -> mean_draw
  return(mean_draw[,-1])
}


## get the EBP entire results, a D * P * t array
#logdens is a m * D  matrix
#parameter of interests: Z, Y1, Y2, Y3, Z|Y1, Z|Y2, Z|Y3
EBPfun_uni <- function(data, fitting, By_hat, SigmaV_hat, m, tt, Zdraw)
{
  aux = cbind(1, data$X2)
  mu_fix = aux %*% By_hat
  Y = cbind(data$Y1, data$Y2)
  ylabel = data$Yid
  munew_fix = cbind(1,fitting$X2) %*% By_hat
  rand_matrix= rmvnorm(m*D, sigma = SigmaV_hat) %>% array(., dim = c(D, m, 2))
  logdens = grouplogdensity_uni(ylabel, By_hat, rand_matrix,  mu_fix, m)
  dens = (logdens - median(logdens)) %>% round(3) %>% exp()
  weight = sweep(dens, STATS= rowSums(dens), FUN = "/" , MARGIN = 1)
  
  temp = sapply(1:tt, FUN = function(x) EBPonedraw_uni(data, fitting, By_hat, SigmaV_hat, 
                                                       munew_fix, weight, rand_matrix, m, Zdraw)) %>% 
    array(., dim = c(D, 7, tt))
  
  return(temp)
}



oneboot <- function(data, fitting, Bz0, By0, SigmaZ0, P0, SigmaUV0, m, tt, numornot, EMsize)
{
  popaux <- rbind(data[,1:2], fitting) %>% arrange(., group)
  popdata0 <- simuldata(popaux$X2 %>% as.vector(), Bz0, By0, P0, SigmaZ0, SigmaUV0, popsize)
  sampledata0 <- left_join(data[,1:2], popdata0)
  lm_z = lmer(Z ~ X2 + Yid + (1|group), data = sampledata0)
  #fixef(lm_zi): Bz_i0, Bz_i1, (Sigma_Z*P)_i
  Bz_hat = fixef(lm_z)[1:2]
  U_hat = ranef(lm_z)$group
  SigmaZ_hat = VarCorr(lm_z) %>% as.data.frame() %>% '['(2,4)
  P_hat1 = as.vector(fixef(lm_z)[3:4])/SigmaZ_hat
  prec_P1 = solve(vcov(lm_z)[3:4, 3:4]/SigmaZ_hat^2)
  #########################################
  #next we use Y_i|Z, Y_{-i}  
  
  glm_noA2 = glmer(Yid ~ X2 + Z  + (1|group), data = sampledata0[sampledata0$Yid!= '2', ], 
                   family = binomial, nAGQ = 20)
  glm_noA3 = glmer(Yid ~ X2 + Z  + (1|group), data = sampledata0[sampledata0$Yid!= '3', ], 
                   family = binomial, nAGQ = 20)
  
  #glm_noA3 %>% fixef : By, Pm', W 
  By_hat = cbind(glm_noA3 %>% fixef %>% "["(1:2),
                 glm_noA2 %>% fixef %>% "["(1:2))
  
  V1 <- rep(0, 40)
  V2 <- rep(0, 40)
  
  V1[rownames(ranef(glm_noA3)$group ) %>% as.numeric()] <- ranef(glm_noA3)$group %>% unlist
  V2[rownames(ranef(glm_noA2)$group ) %>% as.numeric()] <- ranef(glm_noA2)$group %>% unlist
  
  V_hat = cbind(V1, V2)
  
  rand_hat = cbind(U_hat, V_hat)
  smallsd = sqrt(c(0.0001, 0.0008, 0.0008))/2
  if(isSingular(lm_z) == TRUE)
    rand_hat[,1] = rnorm(D, sd = smallsd[1])
  if(isSingular(glm_noA3) == TRUE)
    rand_hat[,2] = rnorm(D, sd = smallsd[2])
  if(isSingular(glm_noA2) == TRUE)
    rand_hat[,3] = rnorm(D, sd = smallsd[3])
  SigmaUV_hat = rand_hat %>% cov
  
  P_hat2 = c(as.vector(fixef(glm_noA3)[3]), 
             as.vector(fixef(glm_noA2)[3]))
  
  prec_P2 = solve(diag(c(diag(vcov(glm_noA3))[3], diag(vcov(glm_noA2))[3])))
  
  P_hat = solve(prec_P1 + prec_P2) %*% (prec_P1 %*% P_hat1 + prec_P2 %*% P_hat2) %>% as.vector()
  
  if(numornot == TRUE)
  {
    expand.grid(x = (1:10)/20, y = (1:10)/20)  %>% as.matrix-> paramlist
    Corr = cor(SigmaUV_hat)[6]
    pll <- NULL
    PSP =  SigmaZ_hat *P_hat  %*% t(P_hat)/2
    probfix = sapply(1:3, FUN = function(x) t(Ychoice[x,]) %*% PSP %*% Ychoice[x,]) 
    
    aux = cbind(1, sampledata0$X2)
    mu_fix = aux %*% By_hat
    Y = cbind(sampledata0$Y1, sampledata0$Y2)
    Z = sampledata0$Z
    ylabel = sampledata0$Yid
    
    for(pos in 1:100)
    {
      SigmaV_hat = diag(paramlist[pos,]) %*% matrix(c(1,Corr,Corr,1), 2) %*% diag(paramlist[pos,])
      L = t(chol(SigmaV_hat)) * sqrt(2)
      V = L %*% t(expand.grid(X_k, X_k))
      sapply(1:dim(V)[2], FUN = function(x) 
        Ydensity_onerand(ylabel, Z, By_hat, P_hat, probfix, mu_fix, samplesize, V[,x])) -> densY
      weight = expand.grid(W_k, W_k) %>% apply(., MARGIN = 1, FUN = prod)
      (t(densY) * weight/sqrt(pi)) %>% colMeans(na.rm = T) %>% "*"(100) %>% log %>% sum(na.rm = T) %>% "*"(-1) -> Yintegal
      pll <- c(pll, Yintegal)
    }
    paramlist[which.min(pll), ] -> newVsd
    SigmaV_hat = diag(newVsd) %*% matrix(c(1, Corr, Corr, 1), nrow = 2) %*% diag(newVsd)
    SigmaUV_hat[2:3, 2:3] <- SigmaV_hat
  }
  
  for(j in 1:EMsize)
  {
    for(k in 1:10)
      SigmaUV_hat = rcppSigmaUVfun(sampledata0, SigmaZ_hat, P_hat,  Bz_hat, By_hat, SigmaUV_hat, m)
    if(update_all == TRUE)
    {
      rcppupdatefix(sampledata0, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, 20) -> temp
      Bz_hat = temp[1:2]
      By_hat = temp[3:6] %>% matrix(., 2)
      SigmaZ_hat = temp[7]
      P_hat = temp[8:9]
      
    }
  }
  #EBP for NEW method: 
  rcppEBPfun(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, m, tt) -> temp
  bootmean = apply(temp, MARGIN = c(1,2), FUN = mean) %>% as.vector()
  bootvar = apply(temp, MARGIN = c(1,2), FUN = var) %>% as.vector()
  #EBP for traditional method: 
  rcppEBPfun(sampledata0, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, m, tt) -> temp2
  bootmean2 = apply(temp2, MARGIN = c(1,2), FUN = mean) %>% as.vector()
  popdata0 %>% group_by(group) %>% 
    summarise(z1mean = mean(Z), y10mean = mean(1-Y1-Y2), y11mean = mean(Y1), y12mean = mean(Y2), 
              z1y10 = sum((1-Y1-Y2)*Z)/sum(1-Y1-Y2), z1y11 = sum(Y1*Z)/sum(Y1), z1y12 = sum(Y2*Z)/sum(Y2)) %>% 
    "["(,-1) %>% unlist() %>% as.vector() -> bootpopmean
  return(cbind(bootmean, bootvar, bootpopmean, bootmean2))
}



############################################################
# This function is one replicate of Monte Carlo simulation
onestep <- function(m, t, b, x, EMsize, numornot)
{
  set.seed(x)
  X = rnorm(D*Pop, sd= 1)  
  truedata = simuldata(X, Bz, By, P, SigmaZ, SigmaUV, popsize)
  
  mapply(sample_n, split(truedata, truedata$group), samplesize, SIMPLIFY = FALSE) %>%
    bind_rows() -> sampledata
  
  ###############################################
  #first we use Z|Y
  lm_z = lmer(Z ~ X2 + Yid + (1|group), data = sampledata)
  #fixef(lm_zi): Bz_i0, Bz_i1, (Sigma_Z*P)_i
  Bz_hat = fixef(lm_z)[1:2]
  
  U_hat = ranef(lm_z)$group
  
  
  SigmaZ_hat = VarCorr(lm_z) %>% as.data.frame() %>% '['(2,4)
  
  
  P_hat1 = as.vector(fixef(lm_z)[3:4])/SigmaZ_hat
  
  prec_P1 = solve(vcov(lm_z)[3:4, 3:4]/SigmaZ_hat^2)
  #########################################
  #next we use Y_i|Z, Y_{-i}  
  
  glm_noA2 = glmer(Yid ~ X2 + Z  + (1|group), data = sampledata[sampledata$Yid!= '2', ], 
                   family = binomial, nAGQ = 20)
  glm_noA3 = glmer(Yid ~ X2 + Z  + (1|group), data = sampledata[sampledata$Yid!= '3', ], 
                   family = binomial, nAGQ = 20)
  
  #glm_noA3 %>% fixef : By, Pm', W 
  By_hat = cbind(glm_noA3 %>% fixef %>% "["(1:2),
                 glm_noA2 %>% fixef %>% "["(1:2))
  
  V1 <- rep(0, 40)
  V2 <- rep(0, 40)
  
  V1[rownames(ranef(glm_noA3)$group ) %>% as.numeric()] <- ranef(glm_noA3)$group %>% unlist
  V2[rownames(ranef(glm_noA2)$group ) %>% as.numeric()] <- ranef(glm_noA2)$group %>% unlist
  
  V_hat = cbind(V1, V2)
  rand_hat = cbind(U_hat, V_hat)
  smallsd = sqrt(c(0.0001, 0.0008, 0.0008))/2
  if(isSingular(lm_z) == TRUE)
    rand_hat[,1] = rnorm(D, sd = smallsd[1])
  if(isSingular(glm_noA3) == TRUE)
    rand_hat[,2] = rnorm(D, sd = smallsd[2])
  if(isSingular(glm_noA2) == TRUE)
    rand_hat[,3] = rnorm(D, sd = smallsd[3])
  SigmaUV_hat = rand_hat %>% cov
  
  P_hat2 = c(as.vector(fixef(glm_noA3)[3]), 
             as.vector(fixef(glm_noA2)[3]))
  
  prec_P2 = solve(diag(c(diag(vcov(glm_noA3))[3], diag(vcov(glm_noA2))[3])))
  
  P_hat = solve(prec_P1 + prec_P2) %*% (prec_P1 %*% P_hat1 + prec_P2 %*% P_hat2) %>% as.vector()
  #P_hat = (P_hat1/var_P1 + P_hat2/var_P2)/(1/var_P1 + 1/var_P2)
  
  if(numornot == TRUE)
  {
    expand.grid(x = (1:10)/20, y = (1:10)/20)  %>% as.matrix-> paramlist
    Corr = cor(SigmaUV_hat)[6]
    pll <- NULL
    PSP =  SigmaZ_hat *P_hat  %*% t(P_hat)/2
    probfix = sapply(1:3, FUN = function(x) t(Ychoice[x,]) %*% PSP %*% Ychoice[x,]) 
    
    aux = cbind(1, sampledata$X2)
    mu_fix = aux %*% By_hat
    Y = cbind(sampledata$Y1, sampledata$Y2)
    Z = sampledata$Z
    ylabel = sampledata$Yid
    
    for(pos in 1:100)
    {
      SigmaV_hat = diag(paramlist[pos,]) %*% matrix(c(1,Corr,Corr,1), 2) %*% diag(paramlist[pos,])
      L = t(chol(SigmaV_hat)) * sqrt(2)
      V = L %*% t(expand.grid(X_k, X_k))
      sapply(1:dim(V)[2], FUN = function(x) 
        Ydensity_onerand(ylabel, Z, By_hat, P_hat, probfix, mu_fix, samplesize, V[,x])) -> densY
      weight = expand.grid(W_k, W_k) %>% apply(., MARGIN = 1, FUN = prod)
      (t(densY) * weight/sqrt(pi)) %>% colMeans(na.rm = T) %>% "*"(100) %>% log %>% sum(na.rm = T) %>% "*"(-1) -> Yintegal
      pll <- c(pll, Yintegal)
    }
    paramlist[which.min(pll), ] -> newVsd
    SigmaV_hat = diag(newVsd) %*% matrix(c(1, Corr, Corr, 1), nrow = 2) %*% diag(newVsd)
    SigmaUV_hat[2:3, 2:3] <- SigmaV_hat
  }
  
  for(j in 1:EMsize)
  {
    for(k in 1:10)
      SigmaUV_hat = rcppSigmaUVfun(sampledata, SigmaZ_hat, P_hat, Bz_hat, By_hat, SigmaUV_hat, m)
    if(update_all == TRUE)
    {
      rcppupdatefix(sampledata, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, 20) -> temp
      Bz_hat = temp[1:2]
      By_hat = temp[3:6] %>% matrix(., 2)
      SigmaZ_hat = temp[7]
      P_hat = temp[8:9]
    }
  }
  
  ############################
  unsampledata = anti_join(truedata, sampledata) %>% "["(,1:2)
  
  #true mean 
  truedata %>% group_by(group) %>% summarise(z1mean = mean(Z), y10mean = mean(1-Y1-Y2), 
                                             y11mean = mean(Y1), y12mean = mean(Y2), 
                                             z1y10 = sum((1-Y1-Y2)*Z)/sum(1-Y1-Y2), 
                                             z1y11 = sum(Y1*Z)/sum(Y1), z1y12 = sum(Y2*Z)/sum(Y2)) -> truemean
  
  
  #only use observed data
  sampledata %>% group_by(group) %>% summarise(z1mean = mean(Z), y10mean = mean(1-Y1-Y2), 
                                               y11mean = mean(Y1), y12mean = mean(Y2), 
                                               z1y10 = sum((1-Y1-Y2)*Z)/sum(1-Y1-Y2), 
                                               z1y11 = sum(Y1*Z)/sum(Y1), z1y12 = sum(Y2*Z)/sum(Y2)) -> samplemean
  
  
  
  #univariate model 
  um_z = lmer(Z ~ X2 + (1|group), data = sampledata)
  
  um_y11 = glmer(Y1 ~ X2 + (1|group), data = sampledata[sampledata$Y2!=1,], family = binomial, nAGQ = 20)
  um_y12 = glmer(Y2 ~ X2 + (1|group), data = sampledata[sampledata$Y1!=1,], family = binomial, 
                 nAGQ = 20)
  
  zpredict = predict(um_z, unsampledata)
  
  By1_uni <- fixef(um_y11)%>%unlist()
  By2_uni <- fixef(um_y12)%>%unlist()
  
  V1_uni <- rep(0, 40)
  V2_uni <- rep(0, 40)

  V1_uni[rownames(ranef(um_y11)$group ) %>% as.numeric()] <- ranef(um_y11)$group %>% unlist
  V2_uni[rownames(ranef(um_y12)$group ) %>% as.numeric()] <- ranef(um_y12)$group %>% unlist
  
  SigmaV_hat <- cov(cbind(V1_uni, V2_uni))
  
  By_hat_uni <- cbind(By1_uni, By2_uni)
  EBPfun_uni(sampledata, unsampledata, By_hat_uni, SigmaV_hat, m, 50, zpredict) %>% 
    apply(., MARGIN = c(1,2), FUN = mean) -> unimodelmean
  
  
  #Hierchical Model: 
  #lmhier_z = lmer(Z ~ X2 + Yid + (1|group), data = sampledata)
  
  #lmhier_y11 = glmer(Y1 ~ X2 + (1|group), data = sampledata[sampledata$Y2!=1,], family = binomial, nAGQ = 20)
  #lmhier_y12 = glmer(Y2 ~ X2 + (1|group), data = sampledata[sampledata$Y1!=1,], family = binomial, 
   #                  nAGQ = 20)
  
  #Bz_hier = fixef(lmhier_z)[1:2]
  #U_hier = ranef(lmhier_z)$group
  #SigmaZ_hier = VarCorr(lmhier_z) %>% as.data.frame() %>% '['(2,4)
  #P_hier= as.vector(fixef(lmhier_z)[3:4])
  
  
  #By1_hier <- fixef(lmhier_y11)%>%unlist()
  #By2_hier <- fixef(lmhier_y12)%>%unlist()
  
  #V1_hier <- rep(0, 40)
  #V2_hier <- rep(0, 40)
  
  #V1_hier[rownames(ranef(lmhier_y11)$group ) %>% as.numeric()] <- ranef(lmhier_y11)$group %>% unlist
  #V2_hier[rownames(ranef(lmhier_y12)$group ) %>% as.numeric()] <- ranef(lmhier_y12)$group %>% unlist
  
  #SigmaUV_hier <- cov(cbind(U_hier, V1_hier, V2_hier))
  
  #By_hier <- cbind(By1_hier, By2_hier)
  
  #for(j in 1:EMsize)
  #  SigmaUV_hier = rcppSigmaUVfunhier(sampledata, SigmaZ_hier, P_hier, Bz_hier, By_hier, SigmaUV_hier, m)
  
  
  #rcppEBPfunhier(sampledata, unsampledata, Bz_hier, By_hier, SigmaZ_hier, P_hier,  SigmaUV_hier, m, t) %>% 
  #  apply(., MARGIN = c(1,2), FUN = mean) -> hiermodelmean
  
  
  # Joint Model: 
  #onegroupEBP(onegroupdata, otherX, Bz_hat, By_hat, SigmaZ_hat, P_hat, W_hat, SigmaUV_hat, m, t)
  rcppEBPfun(sampledata, unsampledata , Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, m, t)  -> temp
  
  jointmean = apply(temp, MARGIN = c(1,2), FUN = mean)
  
  #######Next calculate MSE
  MSE_sample = apply((truemean[,-1] - samplemean[,-1])^2, MARGIN = 2, FUN = mean, na.rm =T)
  MSE_unimodel = apply((truemean[,-1] - unimodelmean)^2, MARGIN = 2, FUN = mean, na.rm =T)
  MSE_jointmodel = apply((truemean[,-1] - jointmean)^2, MARGIN = 2, FUN = mean, na.rm =T)
  #MSE_hiermodel = apply((truemean[,-1] - hiermodelmean)^2, MARGIN = 2, FUN = mean, na.rm =T)
  #cbind(MSE_sample, MSE_unimodel, MSE_hiermodel, MSE_jointmodel) %>% round(7) %>% "*"(100) -> MSEresult
  
  mclapply(1:b, FUN = function(x) oneboot(sampledata, unsampledata, Bz_hat, By_hat, 
                                          SigmaZ_hat, P_hat, SigmaUV_hat, m, t, numornot, EMsize), mc.cores = cores) %>% do.call(cbind,.)-> boot
  
  result = cbind(truemean[,-1] %>% unlist, samplemean[,-1] %>% unlist, 
                 unimodelmean %>% as.vector, jointmean %>% as.vector, 
                 temp %>% matrix(.,ncol = t), boot)
  #return(MSEresult)
}

