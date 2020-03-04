#this code consider a two dimension continuous variables and two binary variables
#
# I will use C, L and K to denote the general case of dimensions of variables 


#  P: K * (C-L)
#  SigmaZ0 : K * K matrix (I think we should estimate the precision matrix)
#  GammaZ0 = solve(SigmaZ0)
#  delta0: n*K matrix
#  mu0: n * (C-L)  
#  Ychoice C * (C-L)
#  PSP : K * K
#  W0: (C-L) * (C-L) , with diagonal 0
#  return a n * C matrix
multiYprob <- function(P0, GammaZ0, delta0, mu0, W0) 
{
  temp = (delta0 %*% P0 + mu0) %*% t(Ychoice)    #n * C
  PSP =  t(P0) %*% solve(GammaZ0) * P0/2
  fixpart = diag(Ychoice %*% (PSP+W0/2) %*% t(Ychoice))
  prob = sweep(temp, STATS = fixpart, FUN =  "+", MARGIN = 2) %>% exp()
  rm(temp)
  rm(fixpart)
  sumprob = rowSums(prob)
  prob = sweep(prob, STATS = sumprob, FUN = "/", MARGIN = 1)
  return(prob)
}


# Bz_hat: p * K is the transform of the Bz_hat in the model (p = 2 here) 
# Bz_hat = [Bz1, Bz2, ...]
# By_hat = [By11, By12, ... , By_L,(C_L-1)]
simuldata <- function(auxillary, Bz_hat, By_hat,  P_hat, GammaZ_hat, SigmaUV_hat, groupsize, W_hat)
{
  randef = rmvnorm(D, sigma = SigmaUV_hat) %>%
    rep(., rep(groupsize,C)) %>% matrix(ncol = C) %>% as.data.frame() 
  X = cbind(1, auxillary)
  delta_draw = (X %*% Bz_hat + randef[, 1:K]) %>% as.matrix()
  SigmaZ0 = solve(GammaZ_hat)
  mu_draw = (X%*% By_hat + randef[, (K+1):(K+C-L)]) %>% as.matrix()
  group = rep(c(1:D), groupsize)   
  Yprob = multiYprob(P_hat, GammaZ_hat, delta_draw, mu_draw, W_hat)
  
  Y = lapply(1:dim(Yprob)[1], FUN = function(x) sample(1:categtotal, size = 1, prob = Yprob[x,])) %>% unlist()
  Yid = Y  %>% as.factor()
  name <- NULL
  for(i in 1:L)
    for(j in 1:Cl[i])
      name <- c(name, paste("Y", i, j, sep= ""))
  Y = Ychoice[Y,]  %>% 'colnames<-'(name)
  
  condZmean = delta_draw + Y %*% t(P_hat) %*% SigmaZ0 
  Z1 = rnorm(sum(groupsize), mean = condZmean[,1], sd = sqrt(SigmaZ0[1]))
  cond_Z2 = condZmean[,2] + SigmaZ0[1,2]/SigmaZ0[1,1]*(Z1 - condZmean[,1])
  var_Z2 = SigmaZ0[2,2] - SigmaZ0[1,2]^2/SigmaZ0[1,1]
  Z2 = rnorm(D*Pop,mean = cond_Z2, sd = sqrt(var_Z2))
  return(data.frame(group, X2 = auxillary, Y, Z1, Z2, Yid))
}

#getAIC("Z2", "X2", c("Z1", "Y11", "Y21"), sampledata)

getAIC <- function(response, fixedvar, targetvar, data1, glm = FALSE)
{
  AICbest <- 100000000000
  modelid <- NULL
  fmlaselect <- NULL
  
  fmla0 <- paste(response, "~", fixedvar, sep = "")
  k <- length(targetvar)
  for(i in 1:2^k -1)
  {
    fmla = fmla0   
    tobit <- intToBits(i) 
    id <- sapply(1:k, FUN= function(x) tobit[x] %>% strtoi)
    idselect <- which(id>0)
    for(j in idselect)
      fmla <- paste(fmla , "+", targetvar[j], sep = "")
    if(glm == FALSE)
    {
      fmla <- paste(fmla, "+(1|group)", sep = " ")
      aic = lmer(as.formula(fmla), data = data1, REML = FALSE) %>% AIC
    }
    if(glm == TRUE)
    {
      aic = glm(as.formula(fmla), data = data1, family = binomial) %>% AIC
      fmla <- paste(fmla, "+(1|group)", sep = " ")
    }
    if(aic < AICbest)
    {
      AICbest <- aic
      modelid <- idselect
      fmlaselect <- fmla
    }
  }
  return(list(AICbest, modelid, fmlaselect))
}

#initialest("X2", c("Z1", "Z2", "Y11", "Y21"), 2, 2, 40, sampledata)
# This functon return a list: 
#[[1]]: graphic structure  (K+C-L) * (K+C-L)
#[[2]]: fixed parameter corresponding with correlation (K+C-L):(K+C-L)
#[[3]]: fixed parameter corresponding with auxillary information (K+C-L):(1+p)
#[[4]]: initial SigmaUV_hat (K+C-L)*(K+C-L)
initialestfun <- function(aux, resp, zsize, ysize, data1, fulledge = FALSE)
{
  formulas <- NULL
  modelid <- list()
  edgematrix <- matrix(0, nrow = zsize+ysize, ncol = zsize + ysize)
  corrmatrix <- matrix(0, nrow = zsize+ysize, ncol = zsize + ysize)
  fixefmatrix <- NULL
  randefmatrix  <- NULL
  for(i in 1:length(resp))
  {
    if(i <= zsize)
      temp = getAIC(resp[i], aux, resp[-i], data1)
    if(i > zsize)
      temp = getAIC(resp[i], aux, resp[-i], data1, glm = TRUE)
    id = temp[[2]]
    modelid[[i]] <- c(id[id<i], id[id>=i] + 1)
    if (fulledge == TRUE)
      modelid[[i]] = (1:4)[-i]
    edgematrix[i,modelid[[i]]] <- 1
  }
  
  # let edge matrix symmetric
  for(i in 1:(zsize+ysize))
    for(j in 1:(zsize+ysize))
      edgematrix[i,j] = max(edgematrix[i,j], edgematrix[j,i])
  
  # fit each condtional model 
  for(i in 1:(zsize + ysize))
  {
    fmla <- paste(resp[i], "~", aux, sep = "")
    for(j in which(edgematrix[i,] > 0))
      fmla <- paste(fmla, "+", resp[j], sep = "")
    fmla <- paste(fmla, "+(1|group)", sep = "")
    if(i <= zsize)
    {
      temp = lmer(as.formula(fmla), data = data1)
      VarCorr(temp) %>% as.data.frame()%>% "["(2,4) -> corrmatrix[i,i]
    }
    
    if(i > zsize)
      temp = glmer(as.formula(fmla), data = data1, family = "binomial", nAGQ = 20)
    paramest = fixef(temp)
    fixefmatrix <- rbind(fixefmatrix, paramest[1:2])
    corrmatrix[i, which(edgematrix[i,] > 0)] <- paramest[-(1:2)]
    
    V <- rep(0, D)
    V[rownames(ranef(temp)$group ) %>% as.numeric()] <- ranef(temp)$group %>% unlist
    smallsd = sqrt(0.0008)/2
    if(isSingular(temp) == TRUE)
      V= rnorm(D, sd = smallsd)
    randefmatrix <- cbind(randefmatrix, V)
  }
  
  #let corrmatrix symmetric
  temp = matrix(0, nrow = zsize+ysize, ncol = zsize + ysize)
  for(i in 1: (zsize + ysize))
    for(j in 1:(zsize + ysize))
    {
      if(abs(sign(corrmatrix[i,j])) + abs(sign(corrmatrix[j,i])) == 2)
        temp[i,j] = temp[j,i] = (corrmatrix[i,j] + corrmatrix[j,i])/2
      else
        temp[i,j] = temp[j,i] = corrmatrix[i,j] + corrmatrix[j,i]
    }
  
  corrmatrix <- temp
  return(list(edgematrix = edgematrix, 
              corrmatrix = corrmatrix,
              auxmatrix = fixefmatrix, 
              sigmaUV = cov(randefmatrix)))
}

#this function returns the initial estimates of 
#(Bz, By, P, W, SigmaZ, SigmaUV, and edge)
firststep <- function(zsize, ysize, data1, fulledge)
{
  aux <- "X2"
  resp <- c("Z1", "Z2", "Y11", "Y21")
  iniest <- initialestfun(aux, resp, zsize, ysize, data1, fulledge)
  Lambda <- diag(diag(iniest[[2]][1:zsize, 1:zsize]))
  M <- iniest$corrmatrix[1:zsize, 1:zsize] - Lambda
  SigmaZ_hat <- solve(-solve(Lambda)%*%M + Lambda)
  
  P_hat <- iniest$corrmatrix[1:zsize, (1+zsize):(zsize+ysize)]
  
  W_hat <- iniest$corrmatrix[(1+zsize):(zsize+ysize), (1+zsize):(zsize+ysize)]
  
  Bz_hat <- t(SigmaZ_hat %*% solve(Lambda) %*% iniest$auxmatrix[1:zsize, ])
  
  By_hat <- iniest$auxmatrix[(1+zsize):(zsize+ysize),] %>% t()
  
  SigmaUV_hat <- iniest$sigmaUV
  
  return(list(edge = iniest$edgematrix, 
              Bz_hat = Bz_hat, By_hat = By_hat, 
              P_hat = P_hat, W_hat  = W_hat, SigmaZ_hat = SigmaZ_hat,
              SigmaUV_hat = SigmaUV_hat))
}

#rcppSigmaUVfun(data, SigmaZ_hat, P_hat, Bz_hat, By_hat, SigmaUV_hat, W_hat, m)
rcppSigmaUVfun <- function(data, SigmaZ_hat, P_hat, Bz_hat, By_hat, SigmaUV_hat, W_hat, m)
{  
  rand_origin = rmvnorm(m*D, sigma = SigmaUV_hat)
  rand_matrix  = rand_origin %>% array(., dim = c(D, m, zsize + ysize))
  PSP = t(P_hat) %*% SigmaZ_hat %*% P_hat/2 + W_hat/2
  probfix = sapply(1:(zsize+ysize), FUN = function(x) t(Ychoice[x,]) %*% PSP %*% Ychoice[x,]) 
  
  aux = cbind(1, data$X2)
  delta_fix = aux %*% Bz_hat
  mu_fix = aux %*% By_hat
  Y = data[,3:(2+ysize)] %>%as.matrix()
  ylabel = data$Yid
  Z = data[,(2+ysize+1):(2+ysize+zsize)]
  A_fix = Y %*% t(P_hat) %*% SigmaZ_hat - Z
  id = rep(1:D, samplesize)
  
  rcppmulgroupdens(P0 = P_hat, G_hat = solve(SigmaZ_hat), rand_matrix = rand_matrix, delta_fix = delta_fix, mu_fix = mu_fix, 
                   id = id, A_fix =A_fix %>% as.matrix(), ylabel = ylabel, probfix = probfix, m = m, Ychoice = Ychoice, 
                   zsize = zsize, ysize =ysize, categtotal = categtotal, D = D) -> logdens
  
  rm(rand_matrix)
  quant_logdens = apply(logdens, MARGIN = 1, FUN = function(x) quantile(x, 0.6)) 
  dens = sweep(logdens, MARGIN = 1, FUN = "-", STATS = quant_logdens) %>% exp()
  weight = sweep(dens, STATS= rowSums(dens), FUN = "/" , MARGIN = 1)/D
  
  return(rcppmulSigmaUVweightavg(rand_origin, weight, zsize, ysize))
}


##### for multivarnomial univariate model
multiYprob_uni <- function(mu0) 
{
  prob = (mu0 %*% t(Ychoice)) %>% exp    #n * 4 
  rm(mu0)
  sumprob = (1+prob[,2]) * (1 + prob[,3])
  prob = sweep(prob, STATS = sumprob, FUN = "/", MARGIN = 1)
  return(prob)
}

## return the log density of Y in univariate case, return a n*m matrix
# mu0:  n * 2 * m array
# sumprob: (1+exp(m1))(1+exp(m2))
getgy_uni <- function(mu0, ylabel, m) 
{
  prob = sapply(1:m, FUN = function(x) (mu0[,x,] %*% t(Ychoice))) %>% 
    array(., dim = c(dim(mu0)[1], 4, m)) #n * 4 * m
  index = cbind(rep(1:dim(prob)[1], m), rep(ylabel, m), rep(1:m, each = dim(prob)[1]))
  sumprob = (1 + exp(prob[,2,])) *(1 +exp(prob[,3,]))
  return(matrix(prob[index] - log(sumprob), ncol = m))
}

# return group log density of Y|X (this case it is two binomial models)
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
  Y_draw = sapply(1:dim(Yprob_new)[1], FUN = function(x) sample(1:categtotal, size = 1, prob = Yprob_new[x,]))
  n = dim(Yprob_new)[1]
  Yid = Y_draw %>% as.factor()
  
  Y_draw = Ychoice[Y_draw, ]
  rm(mu_new)
  rm(Yprob_new)
  data.frame(fitting, Y11 = Y_draw[,1], Y21 = Y_draw[,2], Z1 = Zdraw[,1], Z2 = Zdraw[,2], Yid = Yid) %>% 
    rbind(., data.frame(data)) %>% 
    group_by(group) %>% summarise(z1mean = mean(Z1), z2mean = mean(Z2), 
                                  y10mean = mean(1-Y11), y11mean = mean(Y11), 
                                  y20mean = mean(1-Y21), y21mean = mean(Y21), 
                                  z1y10 = sum((1-Y11)*Z1)/sum(1-Y11), 
                                  z1y11 = sum(Y11*Z1)/sum(Y11),
                                  z1y20 = sum((1-Y21)*Z1)/sum(1-Y21), 
                                  z1y21 = sum(Y21*Z1)/sum(Y21),
                                  z2y10 = sum((1-Y11)*Z2)/sum(1-Y11), 
                                  z2y11 = sum(Y11*Z2)/sum(Y11),
                                  z2y20 = sum((1-Y21)*Z2)/sum(1-Y21), 
                                  z2y21 = sum(Y21*Z2)/sum(Y21))  %>% as.matrix() -> mean_draw
  return(mean_draw[,-1])
}


## get the EBP entire results, a D * P * t array
#logdens is a m * D  matrix
#parameter of interests: Z, Y1, Y2, Y3, Z|Y1, Z|Y2, Z|Y3
EBPfun_uni <- function(data, fitting, By_hat, SigmaV_hat, m, tt, Zdraw)
{
  aux = cbind(1, data$X2)
  mu_fix = aux %*% By_hat
  Y = cbind(data$Y11, data$Y21)
  ylabel = data$Yid
  munew_fix = cbind(1,fitting$X2) %*% By_hat
  rand_matrix= rmvnorm(m*D, sigma = SigmaV_hat) %>% array(., dim = c(D, m, 2))
  logdens = grouplogdensity_uni(ylabel, By_hat, rand_matrix,  mu_fix, m)
  dens = (logdens - median(logdens)) %>% round(3) %>% exp()
  weight = sweep(dens, STATS= rowSums(dens), FUN = "/" , MARGIN = 1)
  
  temp = sapply(1:tt, FUN = function(x) EBPonedraw_uni(data, fitting, By_hat, SigmaV_hat, 
                                                       munew_fix, weight, rand_matrix, m, Zdraw)) %>% 
    array(., dim = c(D, paramtotal, tt))
  
  return(temp)
}


### EBP functions for joint model
rcppEBPonedraw <-function(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, W_hat, SigmaUV_hat, deltanew_fix, munew_fix, weight, rand_matrix, m)
{
  onedraw = apply(weight, MARGIN = 1, FUN = function(x) sample(1:m, size = 1, prob = x))
  rm(weight)
  rand_one = sapply(1:D,  FUN = function(x) rand_matrix[x, onedraw[x], ]) %>% t()
  rm(rand_matrix)
  rm(onedraw)
  delta_new = deltanew_fix + rep(rand_one[,1:zsize], rep(predictsize, zsize)) %>% matrix(., ncol = zsize)
  mu_new = munew_fix + rep(rand_one[,(zsize+1):(zsize+ysize)], rep(predictsize,ysize)) %>% matrix(ncol = ysize)
  n = dim(deltanew_fix)[1]
  rand_vec = runif(n)
  rm(deltanew_fix)
  rm(munew_fix)
  Yprob_new = multiYprob(P_hat, solve(SigmaZ_hat), delta_new, mu_new, W_hat)
  rcppmulYZdraw(delta_new, mu_new, Yprob_new, t(P_hat) %*% SigmaZ_hat, zsize, ysize, rand_vec) -> temp
  rm(delta_new)
  rm(mu_new)
  rm(Yprob_new)
  
  Z1_draw = rnorm(n, mean = temp[,2+ysize], sd = sqrt(SigmaZ_hat[1]))
  cond_Z2 = temp[,3+ysize] + SigmaZ_hat[1,2]/SigmaZ_hat[1,1]*(Z1_draw - temp[,2+ysize])
  var_Z2 = SigmaZ_hat[2,2] - SigmaZ_hat[1,2]^2/SigmaZ_hat[1,1]
  Z2_draw = rnorm(n,mean = cond_Z2, sd = sqrt(SigmaZ_hat[2,2]))
  
  data.frame(fitting, Y11 = temp[,1], Y21 = temp[,2], Z1 = Z1_draw, 
             Z2 = Z2_draw, Yid = temp[,3] %>% as.factor) %>% 
    rbind(., data.frame(data)) %>% 
    group_by(group) %>% summarise(z1mean = mean(Z1), z2mean = mean(Z2), 
                                  y10mean = mean(1-Y11), y11mean = mean(Y11), 
                                  y20mean = mean(1-Y21), y21mean = mean(Y21), 
                                  z1y10 = sum((1-Y11)*Z1)/sum(1-Y11), 
                                  z1y11 = sum(Y11*Z1)/sum(Y11),
                                  z1y20 = sum((1-Y21)*Z1)/sum(1-Y21), 
                                  z1y21 = sum(Y21*Z1)/sum(Y21),
                                  z2y10 = sum((1-Y11)*Z2)/sum(1-Y11), 
                                  z2y11 = sum(Y11*Z2)/sum(Y11),
                                  z2y20 = sum((1-Y21)*Z2)/sum(1-Y21), 
                                  z2y21 = sum(Y21*Z2)/sum(Y21)) %>% as.matrix() -> mean_draw
  return(mean_draw[,-1])
}


rcppEBPfun <- function(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, W_hat, SigmaUV_hat, m, tt)
{
  
  PSP =  t(P_hat) %*% SigmaZ_hat %*% t(P_hat)/2 + W_hat/2
  probfix = sapply(1:categtotal, FUN = function(x) t(Ychoice[x,]) %*% PSP %*% Ychoice[x,]) 
  
  aux = cbind(1, data$X2)
  delta_fix = aux %*% Bz_hat
  mu_fix = aux %*% By_hat
  Y = cbind(data$Y11, data$Y21)
  Z = cbind(data$Z1, data$Z2)
  A_fix = Y %*% t(P_hat) %*% SigmaZ_hat - Z
  ylabel = data$Yid
  
  deltanew_fix = cbind(1,fitting$X2) %*% Bz_hat
  munew_fix = cbind(1,fitting$X2) %*% By_hat
  
  rand_matrix= rmvnorm(m*D, sigma = SigmaUV_hat) %>% array(., dim = c(D, m, zsize +ysize))
  
  id = rep(1:D, samplesize)
  
  logdens = rcppmulgroupdens(P0 = P_hat, G_hat = solve(SigmaZ_hat), rand_matrix = rand_matrix, delta_fix = delta_fix, mu_fix = mu_fix, 
                             id = id, A_fix =A_fix %>% as.matrix(), ylabel = ylabel, probfix = probfix, m = m, Ychoice = Ychoice, 
                             zsize = zsize, ysize =ysize, categtotal = categtotal, D = D) 
  
  dens = (logdens - median(logdens)) %>% round(3) %>% exp()
  weight = sweep(dens, STATS= rowSums(dens), FUN = "/" , MARGIN = 1)
  
  temp = sapply(1:tt, FUN = function(x) rcppEBPonedraw(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, W_hat, SigmaUV_hat, 
                                                       deltanew_fix, munew_fix, weight, rand_matrix, m)) %>% 
    array(., dim = c(D, paramtotal, tt))
  
  return(temp)
}

#simudata(auxillary, Bz_hat, By_hat,  P_hat, GammaZ_hat, SigmaUV_hat, groupsize, W_hat)
#rcppEBPfun(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, W_hat, SigmaUV_hat, m, tt)
oneboot <- function(data, fitting, Bz0, By0, SigmaZ0, P0, W0, SigmaUV0, m, tt, EMsize, fulledge)
{
  popaux <- rbind(data[,1:2], fitting) %>% arrange(., group)
  popdata0 <- simuldata(popaux$X2 %>% as.vector(), Bz0, By0, P0, solve(SigmaZ0), SigmaUV0, popsize, W0)
  sampledata0 <- left_join(data[,1:2], popdata0)
  
  iniest <- firststep(zsize, ysize,sampledata0,fulledge)
  
  Bz_hat <- iniest$Bz_hat
  By_hat <- iniest$By_hat
  P_hat <- iniest$P_hat
  W_hat <- iniest$W_hat
  #W_hat = matrix(rep(0, 4), 2, 2)
  SigmaZ_hat <- iniest$SigmaZ_hat
  SigmaUV_hat <- iniest$SigmaUV_hat
  
  ###########################
  for(j in 1:100)
    SigmaUV_hat = rcppSigmaUVfun(sampledata0, SigmaZ_hat, P_hat, Bz_hat, By_hat, SigmaUV_hat, W_hat, m)
  
  #EBP for NEW method: (use true data but bootstrap model params)
  rcppEBPfun(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, W_hat, SigmaUV_hat, m, tt) -> temp
  bootmean = apply(temp, MARGIN = c(1,2), FUN = mean) %>% as.vector()
  bootvar = apply(temp, MARGIN = c(1,2), FUN = var) %>% as.vector()
  #EBP for traditional method: 
  # rcppEBPfun(sampledata0, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, W_hat,  SigmaUV_hat, m, tt) -> temp2
  # bootmean2 = apply(temp2, MARGIN = c(1,2), FUN = mean) %>% as.vector()
  # popdata0 %>% group_by(group) %>% 
  #   summarise(z1mean = mean(Z1), z2mean = mean(Z2), 
  #             y10mean = mean(1-Y11), y11mean = mean(Y11), 
  #             y20mean = mean(1-Y21), y21mean = mean(Y21), 
  #             z1y10 = sum((1-Y11)*Z1)/sum(1-Y11), 
  #             z1y11 = sum(Y11*Z1)/sum(Y11),
  #             z1y20 = sum((1-Y21)*Z1)/sum(1-Y21), 
  #             z1y21 = sum(Y21*Z1)/sum(Y21),
  #             z2y10 = sum((1-Y11)*Z2)/sum(1-Y11), 
  #             z2y11 = sum(Y11*Z2)/sum(Y11),
  #             z2y20 = sum((1-Y21)*Z2)/sum(1-Y21), 
  #             z2y21 = sum(Y21*Z2)/sum(Y21)) %>% 
  #   "["(,-1) %>% unlist() %>% as.vector() -> bootpopmean
  # return(list(cbind(bootmean, bootvar, bootpopmean, bootmean2), iniest$edge))
  return(list(cbind(bootmean, bootvar), iniest$edge))
}




oneSimulation <- function(m, t, b, x, EMsize, fulledge)
{
  set.seed(x)
  X = rnorm(D*Pop, sd= 1)  
  truedata = simuldata(X, Bz, By, P, GammaZ, SigmaUV, popsize, W)
  
  mapply(sample_n, split(truedata, truedata$group), samplesize, SIMPLIFY = FALSE) %>%
    bind_rows() -> sampledata
  
  
  iniest <- firststep(zsize, ysize,sampledata, fulledge)
  
  Bz_hat <- iniest$Bz_hat
  By_hat <- iniest$By_hat
  P_hat <- iniest$P_hat
  W_hat <- iniest$W_hat
  #W_hat = matrix(rep(0, 4), 2, 2)
  SigmaZ_hat <- iniest$SigmaZ_hat
  SigmaUV_hat <- iniest$SigmaUV_hat
  
  ###########################
  for(j in 1:100)
    SigmaUV_hat = rcppSigmaUVfun(sampledata, SigmaZ_hat, P_hat, Bz_hat, By_hat, SigmaUV_hat, W_hat, m)
  
  unsampledata = anti_join(truedata, sampledata) %>% "["(,1:2)
  
  #true mean (dim = D * (K + C + K * C))
  truedata %>% group_by(group) %>% summarise(z1mean = mean(Z1), z2mean = mean(Z2), 
                                             y10mean = mean(1-Y11), y11mean = mean(Y11), 
                                             y20mean = mean(1-Y21), y21mean = mean(Y21), 
                                             z1y10 = sum((1-Y11)*Z1)/sum(1-Y11), 
                                             z1y11 = sum(Y11*Z1)/sum(Y11),
                                             z1y20 = sum((1-Y21)*Z1)/sum(1-Y21), 
                                             z1y21 = sum(Y21*Z1)/sum(Y21),
                                             z2y10 = sum((1-Y11)*Z2)/sum(1-Y11), 
                                             z2y11 = sum(Y11*Z2)/sum(Y11),
                                             z2y20 = sum((1-Y21)*Z2)/sum(1-Y21), 
                                             z2y21 = sum(Y21*Z2)/sum(Y21)) -> truemean
  
  
  #only use observed data
  sampledata %>% group_by(group) %>% summarise(z1mean = mean(Z1), z2mean = mean(Z2), 
                                               y10mean = mean(1-Y11), y11mean = mean(Y11), 
                                               y20mean = mean(1-Y21), y21mean = mean(Y21), 
                                               z1y10 = sum((1-Y11)*Z1)/sum(1-Y11), 
                                               z1y11 = sum(Y11*Z1)/sum(Y11),
                                               z1y20 = sum((1-Y21)*Z1)/sum(1-Y21), 
                                               z1y21 = sum(Y21*Z1)/sum(Y21),
                                               z2y10 = sum((1-Y11)*Z2)/sum(1-Y11), 
                                               z2y11 = sum(Y11*Z2)/sum(Y11),
                                               z2y20 = sum((1-Y21)*Z2)/sum(1-Y21), 
                                               z2y21 = sum(Y21*Z2)/sum(Y21)) -> samplemean
  
  #univariate model 
  um_z1 = lmer(Z1 ~ X2 + (1|group), data = sampledata)
  um_z2 = lmer(Z2 ~ X2 + (1|group), data = sampledata)
  um_y11 = glmer(Y11 ~ X2 + (1|group), data = sampledata, family = binomial, nAGQ = 20)
  um_y21 = glmer(Y21 ~ X2 + (1|group), data = sampledata, family = binomial, nAGQ = 20)
  
  zpredict = cbind(predict(um_z1, unsampledata), predict(um_z2, unsampledata))
  
  By1_uni <- fixef(um_y11)%>%unlist()
  By2_uni <- fixef(um_y21)%>%unlist()
  
  V1_uni <- rep(0, D)
  V2_uni <- rep(0, D)
  
  V1_uni[rownames(ranef(um_y11)$group ) %>% as.numeric()] <- ranef(um_y11)$group %>% unlist
  V2_uni[rownames(ranef(um_y21)$group ) %>% as.numeric()] <- ranef(um_y21)$group %>% unlist
  
  SigmaV_hat <- cov(cbind(V1_uni, V2_uni))
  
  By_hat_uni <- cbind(By1_uni, By2_uni)
  EBPfun_uni(sampledata, unsampledata, By_hat_uni, SigmaV_hat, m, 50, zpredict) %>% 
    apply(., MARGIN = c(1,2), FUN = mean) -> unimodelmean
  
  #joint model
  rcppEBPfun(sampledata, unsampledata , Bz_hat, By_hat, SigmaZ_hat, P_hat, W_hat, SigmaUV_hat, m, t)  -> temp
  
  
  jointmean = apply(temp, MARGIN = c(1,2), FUN = mean)
  
  #######Next calculate MSE
  MSE_sample = apply((truemean[,-1] - samplemean[,-1])^2, MARGIN = 2, FUN = mean, na.rm =T)
  MSE_unimodel = apply((truemean[,-1] - unimodelmean)^2, MARGIN = 2, FUN = mean, na.rm =T)
  MSE_jointmodel = apply((truemean[,-1] - jointmean)^2, MARGIN = 2, FUN = mean, na.rm =T)
  #MSE_hiermodel = apply((truemean[,-1] - hiermodelmean)^2, MARGIN = 2, FUN = mean, na.rm =T)
  cbind(MSE_sample, MSE_unimodel, MSE_jointmodel) %>% round(7) %>% "*"(100)
  
  mclapply(1:b, FUN = function(x) oneboot(sampledata, unsampledata, Bz_hat, By_hat, 
                                          SigmaZ_hat, P_hat, W_hat, SigmaUV_hat, m, t, EMsize, fulledge), 
           mc.cores = cores) -> boottemp 
  boot <- NULL
  edge <- list(iniest$edge)
  for(i in 1:length(boottemp))
  {
    boot<- cbind(boot, boottemp[[i]][[1]])
    edge[[i+1]] <- boottemp[[i]][[2]]
  }
  
  result = cbind(truemean[,-1] %>% unlist, samplemean[,-1] %>% unlist, 
                 unimodelmean %>% as.vector, jointmean %>% as.vector, 
                 temp %>% matrix(.,ncol = t), boot)
  return(list(result = result, edge = edge))
  #return(MSEresult)
}
