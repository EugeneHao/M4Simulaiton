relabias <- function(result, D, m, t, b, plot = FALSE, varname, title = NULL)
{
  MSE = NULL
  MSEtrue = NULL
  MSEsample = NULL
  MSEuni = NULL
  Varpart = NULL
  biaspart = NULL
  Varbias = NULL
  truemean = NULL
  estmean = NULL
  id = NULL
  for(i in 1:length(result))
  {
    which(colnames(result[[i]][[1]]) == "") -> a
    if(length(a) > t + 4)
     result[[i]][[1]] <- as.numeric(result[[i]][[1]][,-c(a[-(1:(t+4))])]) %>% matrix(nrow = D*length(varname))
    if(length(dim(result[[i]][[1]])) >1 && is.numeric(result[[i]][[1]]) == TRUE)
      id = c(id, i)
  }
  l = 1
  for(i in id)
  {
    b = (dim(result[[i]][[1]])[2] - 4 - t)/2
    Varpart = cbind(Varpart, result[[i]][[1]][, 5:(4+t)] %>% apply(., MARGIN = 1, FUN = var))
    biaspart = cbind(biaspart, (result[[i]][[1]][, 2*(1:b)+1+t] - result[[i]][[1]][, 4])^2 %>% rowMeans(na.rm = T))
    Varbias = cbind(Varbias, result[[i]][[1]][, 2*(1:b)+2+t] %>% apply(., MARGIN = 1, FUN = mean, na.rm = T) - Varpart[,l])
    #MSE = cbind(MSE, Varpart+biaspart)
    MSEtrue = cbind(MSEtrue, ( result[[i]][[1]][, 1]- result[[i]][[1]][, 4])^2)
    MSEsample = cbind(MSEsample, ( result[[i]][[1]][, 1]- result[[i]][[1]][, 2])^2)
    MSEuni = cbind(MSEuni, ( result[[i]][[1]][, 1]- result[[i]][[1]][, 3])^2)
    truemean = cbind(truemean, result[[i]][[1]][,1])
    estmean = cbind(estmean, result[[i]][[1]][,4])
    l = l+1
  }
  
  MSE =  (Varpart+biaspart-Varbias + abs(Varpart+biaspart-Varbias))/2
  
  #MSE = Varpart + biaspart
  LB = estmean - 1.96 * sqrt(MSE)
  UB = estmean + 1.96 * sqrt(MSE)
  ((LB <= truemean) *(UB >=truemean)) %>% rowSums(na.rm = T) %>% 
    matrix(., nrow = D) %>% colSums(na.rm = T)  -> good
  
  (1 - is.na((LB <= truemean) *(UB >=truemean))) %>% rowSums(na.rm = T) %>% 
    matrix(., nrow = D) %>% colSums(na.rm = T)  -> total
  (good/total) %>% "*"(100) -> CR
  # look at each part: 
  Varpart %>% rowMeans(na.rm = T) %>% matrix(nrow =D) %>% 
    "colnames<-"(c(varname))-> Varpartmean
  biaspart %>% rowMeans(na.rm = T) %>% matrix(nrow =D) %>% 
    "colnames<-"(c(varname)) -> biaspartmean
  Varbias %>% rowMeans(na.rm = T) %>% matrix(nrow =D) %>% 
    "colnames<-"(c(varname))-> Varbiasmean
  MSEtrue %>% rowMeans(na.rm = T) %>% matrix(nrow =D) %>% 
    "colnames<-"(c(varname))-> MSEtruemean
  MSE %>% rowMeans(na.rm = T) %>% matrix(nrow =D) %>% 
    "colnames<-"(c(varname)) -> MSEmean

  MSEtrueavg = MSEtrue %>% rowMeans()
  
  relativebias = (MSEmean - MSEtruemean) %>% colSums(na.rm = T) %>% "/" (MSEtruemean %>% colSums()) %>% as.vector()
  if (plot != FALSE)
  {
    ((MSEmean- MSEtruemean)/MSEtruemean) %>% boxplot(main = " MSE") 
    return(0)
  }
  if(plot == FALSE)
    rbind(relativebias, CR) %>% return()
}



###################
reportMSE <- function(result, D, m, t, b, varname)
{
  MSE = NULL
  MSEtrue = NULL
  MSEsample = NULL
  MSEuni = NULL
  Varpart = NULL
  biaspart = NULL
  Varbias = NULL
  truemean = NULL
  estmean = NULL

  id = NULL
  for(i in 1:length(result))
  {
    which(colnames(result[[i]][[1]]) == "") -> a
    if(length(a) > t + 4)
      result[[i]][[1]] <- as.numeric(result[[i]][[1]][,-c(a[-(1:(t+4))])]) %>% matrix(nrow = D*length(varname))
    if(length(dim(result[[i]][[1]])) >1 && is.numeric(result[[i]][[1]]) == TRUE)
      id = c(id, i)
  }
  
  l = 1
  for(i in id)
  {
    b = (dim(result[[i]][[1]])[2] - 4 - t)/2
    Varpart = cbind(Varpart, result[[i]][[1]][, 5:(4+t)] %>% apply(., MARGIN = 1, FUN = var))
    biaspart = cbind(biaspart, (result[[i]][[1]][, 2*(1:b)+1+t] - result[[i]][[1]][, 4])^2 %>% rowMeans(na.rm = T))
    Varbias = cbind(Varbias, result[[i]][[1]][, 2*(1:b)+2+t] %>% 
                      apply(., MARGIN = 1, FUN = mean, na.rm = T) - Varpart[,l])
    #MSE = cbind(MSE, Varpart+biaspart)
    MSEtrue = cbind(MSEtrue, ( result[[i]][[1]][, 1]- result[[i]][[1]][, 4])^2)
    MSEsample = cbind(MSEsample, ( result[[i]][[1]][, 1]- result[[i]][[1]][, 2])^2)
    MSEuni = cbind(MSEuni, ( result[[i]][[1]][, 1]- result[[i]][[1]][, 3])^2)
    truemean = cbind(truemean, result[[i]][[1]][,1])
    estmean = cbind(estmean, result[[i]][[1]][,4])
    l = l+1
  }
  
  MSE =  (Varpart+biaspart-Varbias + abs(Varpart+biaspart-Varbias))/2
  
  Varpart %>% array(., dim = c(D, length(varname), length(id))) %>% apply(MARGIN = 2, FUN = mean, na.rm = T)-> Varpartmean
  biaspart %>% array(., dim = c(D, length(varname), length(id))) %>% apply(MARGIN = 2, FUN = mean, na.rm = T)-> biaspartmean
  Varbias %>% array(., dim = c(D, length(varname), length(id))) %>% apply(MARGIN = 2, FUN = mean, na.rm = T)-> Varbiasmean
  MSEtrue %>% array(., dim = c(D, length(varname), length(id))) %>% apply(MARGIN = 2, FUN = mean, na.rm = T)-> MSEtruemean
  MSE %>%array(., dim = c(D, length(varname), length(id))) %>% apply(MARGIN = 2, FUN = mean, na.rm = T) -> MSEmean
  MSEsample %>% array(., dim = c(D, length(varname), length(id))) %>% apply(MARGIN = 2, FUN = mean, na.rm = T) -> MSEsamplemean
  MSEuni %>% array(., dim = c(D, length(varname), length(id))) %>% apply(MARGIN = 2, FUN = mean, na.rm = T) -> MSEunimean
  result = cbind(MSEsamplemean, MSEunimean, MSEtruemean, MSEmean) %>% 
    "colnames<-"(c( "sample", "uni", "true", "MSE")) %>% 
    "rownames<-"(varname) %>% "*"(100) %>% round(4)
  return(result)
}


checksize <- function(data)
{
  size = NULL
  for(i in 1:length(data))
  {
      if(length(dim(data[[i]][[1]])) >1)
      size = c(size, i)
  }
  return(length(size))
}

extractPW <- function(data)
{
  temp <- NULL
  for(i in 1:length(data))
  {
    temp <- rbind(temp, c(data[[i]][[2]]$W_hat[2], data[[i]][[2]]$P_hat[1:4], 
                  data[[i]][[2]]$SigmaZ_hat %>% solve %>% "["(2)))
  }
  return(temp)
}

avgMSE <- function(data)
{
  temp = data[[i]][[1]]
  for(i in 2:length(data))
  {
    temp = temp + data[[i]][[1]]
  }
  return(temp/length(data))
}

edgeprob <- function(data)
{
  temp = data[[i]][[2]]$edge
  for(i in 2:length(data))
  {
    temp = temp + data[[i]][[2]]$edge
  }
  temp = temp/length(data)
  return(temp[c(12, 3, 7, 4,8, 2)])
}

bootedgeprob <- function(data)
{
  temp = data[[1]]$edge[[1]]
  for(i in 2:length(data))
  {
    temp = temp + data[[i]]$edge[[1]]
  }
  temp = temp/length(data)
  return(temp[c(12, 3, 7, 5, 8, 2)])
}
