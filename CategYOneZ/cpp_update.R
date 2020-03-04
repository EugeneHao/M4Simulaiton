ratiofuncpp <- ' 
arma::vec x = Rcpp::as<arma::vec>(X);
arma::vec y1 = Rcpp::as<arma::vec>(Y1);
arma::vec y2 = Rcpp::as<arma::vec>(Y2);
arma::vec z = Rcpp::as<arma::vec>(Z);
arma::vec xgroup = Rcpp::as<arma::vec>(group);
double xp = Rcpp::as<double>(p);
double xrand = Rcpp::as<double>(randef);
int xD = Rcpp::as<int>(D);
int xpart = Rcpp::as<int>(part);
arma::vec param = Rcpp::as<arma::vec>(otherparam);
int n = x.n_rows; 
arma::vec numerpart = arma::zeros<arma::vec>(n);
arma::vec linearpart = arma::zeros<arma::vec>(n);
arma::vec result = arma::ones<arma::vec>(xD);
if(xpart == 1)
{
	for(int i = 1; i<=n; i++)
	{
		numerpart(i-1) = (param(0) + param(1) * x(i-1) + xp * z(i-1) + xrand) * y1(i-1);
		linearpart(i-1) = param(0) + param(1) * x(i-1) + xp * z(i-1) + xrand;
	}
}
else
{
	for(int i = 1; i<=n; i++)
	{
		numerpart(i-1) = (param(0) + param(1) * x(i-1) + xp * z(i-1) + xrand) * y2(i-1);	
		linearpart(i-1) = param(0) + param(1) * x(i-1) + xp * z(i-1) + xrand;
	}	
}
for(int i = 1; i<=n; i++)
	result(xgroup(i-1)-1) *= exp(numerpart(i-1))/(1+exp(linearpart(i-1)));
return Rcpp::wrap(result);
'

rcppratiofun <- cxxfunction(signature(X="numeric", Y1="numeric", Y2="numeric",
                                      Z="numeric", group = "int", p="numeric", randef="numeric", D = "int", part = "int",
                                      otherparam="numeric"), ratiofuncpp, plugin ="RcppArmadillo", verbose = TRUE)

#Example: 
#rcppratiofun(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group, 
#              p = 0.5, randef = 0.6, D = 40, part = 1, otherparam = param2)


ratiodifcpp <- ' 
arma::vec x = Rcpp::as<arma::vec>(X);
arma::vec y1 = Rcpp::as<arma::vec>(Y1);
arma::vec y2 = Rcpp::as<arma::vec>(Y2);
arma::vec z = Rcpp::as<arma::vec>(Z);
arma::vec xgroup = Rcpp::as<arma::vec>(group);
double xp = Rcpp::as<double>(p);
double xrand = Rcpp::as<double>(randef);
int xD = Rcpp::as<int>(D);
int xpart = Rcpp::as<int>(part);
arma::vec param = Rcpp::as<arma::vec>(otherparam);
int n = x.n_rows; 
arma::vec numerpart = arma::zeros<arma::vec>(n);
arma::vec linearpart = arma::zeros<arma::vec>(n);
arma::vec Zpart = arma::zeros<arma::vec>(n);
arma::vec Zsum = arma::zeros<arma::vec>(xD);
arma::vec result = arma::ones<arma::vec>(xD);
if(xpart == 1)
{
	for(int i = 1; i<=n; i++)
	{
		numerpart(i-1) = (param(0) + param(1) * x(i-1) + xp * z(i-1) + xrand) * y1(i-1);
		linearpart(i-1) = param(0) + param(1) * x(i-1) + xp * z(i-1) + xrand;
		Zpart(i-1) = z(i-1)*exp(linearpart(i-1))/pow(1+exp(linearpart(i-1)),2) * (2*y1(i-1) - 1);
	}
}
if(xpart == 2)
{
	for(int i = 1; i<=n; i++)
	{
		numerpart(i-1) = (param(0) + param(1) * x(i-1) + xp * z(i-1) + xrand) * y2(i-1);	
		linearpart(i-1) = param(0) + param(1) * x(i-1) + xp * z(i-1) + xrand;
		Zpart(i-1) = z(i-1)*exp(linearpart(i-1))/pow(1+exp(linearpart(i-1)),2) * (2*y2(i-1) - 1);

	}	
}
for(int i = 1; i<=n; i++)
{
	result(xgroup(i-1)-1) *= exp(numerpart(i-1))/(1+exp(linearpart(i-1)));
	Zsum(xgroup(i-1)-1) += Zpart(i-1); 
}
for(int i = 1; i<=xD; i++)
	result(i-1) *= Zsum(i-1);
return Rcpp::wrap(result);
'

rcppratiodif <- cxxfunction(signature(X="numeric", Y1="numeric", Y2="numeric",
                                      Z="numeric", group = "int", p="numeric", randef="numeric", D = "int", part = "int",
                                      otherparam="numeric"), ratiodifcpp, plugin ="RcppArmadillo", verbose = TRUE)

#Example:
#  rcpprationdif(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group, 
#                p = 0.5, randef = 0.6, D = 40, part = 1, otherparam = param2)


rcppYloglike <- function(data, p, sigma, param, part)
{
  sapply(1:length(X_k), FUN = function(x)  
    rcppratiofun(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group, 
                 p, randef = X_k[x]* sigma * sqrt(2),  D=D, part = part,otherparam= param)) -> temp
  (temp %*% W_k/sqrt(pi)) %>% log %>% sum %>% "*"(-1) %>% return()
}


Zlogbasecpp <- ' 
arma::vec x = Rcpp::as<arma::vec>(X);
arma::vec y1 = Rcpp::as<arma::vec>(Y1);
arma::vec y2 = Rcpp::as<arma::vec>(Y2);
arma::vec z = Rcpp::as<arma::vec>(Z);
arma::vec xgroup = Rcpp::as<arma::vec>(group);
arma::vec xP = Rcpp::as<arma::vec>(P);
double xrand = Rcpp::as<double>(randef);
int xD = Rcpp::as<int>(D);
arma::vec xparam = Rcpp::as<arma::vec>(param);
int n = x.n_rows; 
arma::vec A = arma::zeros<arma::vec>(n);
arma::vec result = arma::ones<arma::vec>(xD);
for(int i = 1; i<=n; i++)
{
	A(i-1) = xparam(0)+xparam(1)*x(i-1) + xparam(2)*(y1(i-1)*xP(0)+y2(i-1)*xP(1)) + xrand;
	result(xgroup(i-1)-1) *= R::dnorm4(z(i-1), A(i-1), sqrt(xparam(2)),0);

}
return Rcpp::wrap(result); 
'

rcppZlogbase <- cxxfunction(signature(X="numeric", Y1="numeric", Y2="numeric",
                                      Z="numeric", group = "int", P="numeric", randef="numeric", D = "int", 
                                      param="numeric"), Zlogbasecpp, plugin ="RcppArmadillo", verbose = TRUE)

#Example:
#  rcppZlogbase(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group, 
#	P = P, randef = 0.6, D = 40, param = param1),

Zlogbasedifcpp <- ' 
arma::vec x = Rcpp::as<arma::vec>(X);
arma::vec y1 = Rcpp::as<arma::vec>(Y1);
arma::vec y2 = Rcpp::as<arma::vec>(Y2);
arma::vec z = Rcpp::as<arma::vec>(Z);
arma::vec xgroup = Rcpp::as<arma::vec>(group);
arma::vec xP = Rcpp::as<arma::vec>(P);
double xrand = Rcpp::as<double>(randef);
int xD = Rcpp::as<int>(D);
arma::vec xparam = Rcpp::as<arma::vec>(param);
int xpart = Rcpp::as<int>(part);
int n = x.n_rows; 
arma::vec A = arma::zeros<arma::vec>(n);
arma::vec dif = arma::zeros<arma::vec>(n); 
arma::vec sumdif = arma::zeros<arma::vec>(xD);
arma::vec result = arma::ones<arma::vec>(xD);

if(xpart == 1)
{
	for(int i = 1; i<=n; i++)
	{
		A(i-1) = xparam(0)+xparam(1)*x(i-1) + xparam(2)*(y1(i-1)*xP(0)+y2(i-1)*xP(1)) + xrand;
		dif(i-1) = (z(i-1) - A(i-1))*y1(i-1);
	}
}
else
{
	for(int i = 1; i<=n; i++)
	{
		A(i-1) = xparam(0)+xparam(1)*x(i-1) + xparam(2)*(y1(i-1)*xP(0)+y2(i-1)*xP(1)) + xrand;
		dif(i-1) = (z(i-1) - A(i-1))*y2(i-1);
	}
}

for(int i = 1; i<=n; i++)
{
	result(xgroup(i-1)-1) *= R::dnorm4(z(i-1), A(i-1), sqrt(xparam(2)),0);
	sumdif(xgroup(i-1)-1) += dif(i-1); 
}
for(int i = 1; i<=xD; i++)
	result(i-1) *= sumdif(i-1);
return Rcpp::wrap(result);
'

rcppZlogbasedif <- cxxfunction(signature(X="numeric", Y1="numeric", Y2="numeric",
                                         Z="numeric", group = "int", P="numeric", randef="numeric", D = "int", 
                                         param="numeric", part = "int"), Zlogbasedifcpp, plugin ="RcppArmadillo", verbose = TRUE)


#Example"
# rcppZlogbasedif(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group,
#  P = P, randef = 0.6, D = 40, param = param1, part = 1),

rcppZloglike <- function(data, P, sigma, param)
{
  sapply(1:length(X_k), FUN = function(x)  
    rcppZlogbase(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group, 
                 P = P,randef = sqrt(2) * sigma * X_k[x], D = D, param = param)) -> temp
  (temp %*% W_k/sqrt(pi)) %>% log %>% sum %>% "*"(-1) %>% return()
}


#Example: 
#rcppratiofun(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group, 
#              p = 0.5, randef = 0.6, D = 40, part = 1, otherparam = param2)

##part deriv function in binomial model (neg)
# param is betaY1 or betaY2 
rcppPdevfun <- function(data, p, sigma, param, part)
{
  sapply(1:length(X_k), FUN = function(x) 
    rcppratiofun(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group, 
                 p = p, randef = X_k[x]* sigma * sqrt(2), D = D, part = part, otherparam = param)) -> temp1
  sapply(1:length(X_k), FUN = function(x) 
    rcppratiodif(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group, 
                 p = p, randef = X_k[x]* sigma * sqrt(2), D = D, part = part, otherparam = param)) -> temp2
  ((temp2 %*% W_k) / (temp1 %*% W_k)) %>% sum %>% "*"(-1)%>% return() 
  
}

##part deriv function in continous model (neg)
# return both dev for simple
# sigma = sigmaU (std dev)
# param = [betaZ, sigmaZ]
rcppPdevfun2 <- function(data, P, sigma, param)
{
  Bz <- param[1:2]
  sigmaZ <- param[3]
  sapply(1:length(X_k), FUN = function(x) 
    rcppZlogbasedif(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group, 
                    P = P, randef = X_k[x]* sigma * sqrt(2), D = D, param = param, part = 1)) -> temp1
  sapply(1:length(X_k), FUN = function(x) 
    rcppZlogbasedif(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group, 
                    P = P, randef = X_k[x]* sigma * sqrt(2), D = D, param = param, part = 2)) -> temp2
  sapply(1:length(X_k), FUN = function(x)  
    rcppZlogbase(X = data$X2, Y1 = data$Y1, Y2 = data$Y2, Z = data$Z, group = data$group, 
                 P = P, randef = X_k[x]* sigma * sqrt(2), D = D, param = param)) -> temp
  return(c(((temp1 %*% W_k) / (temp %*% W_k)) %>% sum %>% "*"(-1),
           ((temp2 %*% W_k) / (temp %*% W_k)) %>% sum %>% "*"(-1)))
}

##loglikelihood of whole model (can also used for return other parameters)
# P = [p1, p2]
# Sigma = [sigmaU, sigmaV1, sigmaV2]  (std dev)
# initparam = [bZ, bY1, bY2, sigmaZ]
rcpploglikefun <- function(dataV1, dataV2, data, P, Sigma, initparam, getparam = FALSE)
{
  bZinit = initparam[1:2]
  bY1init = initparam[3:4]
  bY2init = initparam[5:6]
  sigmaZinit = initparam[7]
  o1 = optim(par = c(bZinit, sigmaZinit), fn = rcppZloglike, data = data, 
             P = P, sigma = Sigma[1])
  o2 = optim(par = bY1init, fn = rcppYloglike, data = dataV1, p = P[1], 
             sigma = Sigma[2], part =1)
  o3 = optim(par = bY2init, fn = rcppYloglike, data = dataV2, p = P[2], 
             sigma = Sigma[3], part =2)
  if(getparam == TRUE)
    return(c(o1$par[1:2], o2$par, o3$par, o1$par[3]))
  if(getparam == FALSE)
    return(o1$value + o2$value + o3$value)
}

## derivative of the composite conditional likelihood wrt P
# param = [p1, p2]
# Sigma = [sigmaU, sigmaV1, sigmaV2] (std dev)
rcppfulldevfun <- function(dataV1, dataV2, data, param, Sigma, initparam)
{
  bZinit = initparam[1:2]
  bY1init = initparam[3:4]
  bY2init = initparam[5:6]
  sigmaZinit = initparam[7]
  o1 = optim(par = c(bZinit, sigmaZinit), fn = rcppZloglike, data = data, 
             P = param, sigma = Sigma[1])
  o2 = optim(par = bY1init, fn = rcppYloglike, data = dataV1, p = param[1], 
             sigma = Sigma[2], part =1)
  o3 = optim(par = bY2init, fn = rcppYloglike, data = dataV2, p = param[2], 
             sigma = Sigma[3], part =2)
  param1 = o1$par   #Z   [betaZ, sigmaZ]
  param2 = o2$par   #Y1  betaY1  
  param3 = o3$par   #Y2  betaY2
  devpart1 = c(rcppPdevfun(dataV1, param[1], Sigma[2], param2, part = 1), 
               rcppPdevfun(dataV2, param[2], Sigma[3], param3, part = 2))
  devpart2 = rcppPdevfun2(data, param, Sigma[1], param1)
  return(devpart1 + devpart2)
}

rcppupdatefix <- function(sampledata, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat,maxiter)
{
  data = sampledata
  dataV1 = data[data$Y2 ==0, ]
  dataV2 = data[data$Y1 ==0, ]
  param = P_hat
  Sigma = diag(SigmaUV_hat) %>% sqrt()
  initparam = c(Bz_hat, By_hat %>% as.vector(), SigmaZ_hat)
  optim(par = P_hat, fn = rcpploglikefun, dataV1 = dataV1, dataV2 = dataV2, data = data, Sigma = Sigma, gr = rcppfulldevfun,
        initparam = initparam, control = list(maxit = maxiter))$par -> newP_hat
  rcpploglikefun(dataV1, dataV2, data, newP_hat, Sigma, initparam, getparam = T) ->newother
  return(c(newother, newP_hat))
}
