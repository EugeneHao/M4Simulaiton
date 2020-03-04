mulgroupdenscpp <-'
arma::mat xP0 = Rcpp::as<arma::mat>(P0);
arma::mat xG = Rcpp::as<arma::mat>(G_hat);
arma::cube xarray = Rcpp::as<arma::cube>(rand_matrix);
arma::mat xdelta = Rcpp::as<arma::mat>(delta_fix);
arma::mat xmu = Rcpp::as<arma::mat>(mu_fix);
arma::vec xid = Rcpp::as<arma::vec>(id);
arma::mat xAfix = Rcpp::as<arma::mat>(A_fix);
arma::vec xylabel = Rcpp::as<arma::vec>(ylabel);
arma::vec xprob = Rcpp::as<arma::vec>(probfix);
int xm = Rcpp::as<int>(m);
arma::mat xYchoice = Rcpp::as<arma::mat>(Ychoice);
int zleng = Rcpp::as<int>(zsize);
int yleng = Rcpp::as<int>(ysize);
int xcateg = Rcpp::as<int>(categtotal);
int n = xid.n_rows;
int xD = Rcpp::as<int>(D);
arma::mat temp = arma::zeros<arma::mat>(n,xcateg);
arma::vec temp2 = arma::zeros<arma::vec>(yleng);
arma::vec delta = arma::zeros<arma::vec>(zleng);
double pick;
double sumprob;
double kernel; 
arma::mat result = arma::zeros<arma::mat>(xD, xm);
double check;
for(int i = 1; i<=xm; i++)
{
  for(int j = 1; j<=n; j++)
  {
  for(int a = 1; a <= zleng; a++)
  {
    delta(a-1) = xarray(xid(j-1)-1, i-1, a-1) + xdelta(j-1, a-1);
  }
  for(int l = 1; l<=yleng; l++)
  {
    temp2(l-1) = delta(0) * xP0(0,l-1) + delta(1) * xP0(1,l-1) + xarray(xid(j-1)-1, i-1, zleng+l-1) + xmu(j-1,l-1);
  }
  sumprob = 0;
  for(int k = 1; k <= xcateg; k++) 
  { 
    temp(j-1,k-1) = temp2(0) * xYchoice(k-1, 0) + temp2(1) * xYchoice(k-1,1) + xprob(k-1);
    sumprob += exp(temp(j-1,k-1));
  }
  pick = temp(j-1, xylabel(j-1)-1);
  kernel = -pow(xAfix(j-1,0) + delta(0), 2) * xG(0,0)/2 -pow(xAfix(j-1,1) + delta(1), 2) * xG(1,1)/2 - (xAfix(j-1,0) + delta(0))*(xAfix(j-1,1) + delta(1))*xG(0,1)/2;
  result(xid(j-1)-1, i-1) += pick - log(sumprob) + kernel;
  } 
}
return Rcpp::wrap(result);
'

rcppmulgroupdens <- cxxfunction(signature(P0="numeric", G_hat = "numeric", rand_matrix = "numeric", delta_fix = "numeric", 
                                       mu_fix = "numeric", id = "int", A_fix = "numeric", 
                                       ylabel = "int", probfix = "numeric", m = "int", 
                                       Ychoice = "int",zsize = "int", ysize = "int", categtotal = "int", D = "int"),
                             mulgroupdenscpp,plugin="RcppArmadillo", verbose = TRUE)

#rcppmulgroupdens(P0 = P_hat, G_hat = solve(SigmaZ_hat), rand_matrix = rand_matrix, delta_fix = delta_fix, mu_fix = mu_fix, 
#                 id = id, A_fix =A_fix %>% as.matrix(), ylabel = ylabel, probfix = probfix, m = m, Ychoice = Ychoice, 
#                 zsize = zsize, ysize =ysize, C = C, D = D)

mulSigmaUVweightavgcpp <- ' 
  arma::mat xrand = Rcpp::as<arma::mat>(rand_origin);
  arma::vec xweight = Rcpp::as<arma::vec>(weight);
  int zleng = Rcpp::as<int>(zsize);
  int yleng = Rcpp::as<int>(ysize);
  int n = xweight.n_rows;
  int l = zleng + yleng;
  arma::mat result = arma::zeros<arma::mat>(l,l);
  for(int i = 1; i<=l; i++)
    for(int j = 1; j<=l; j++)
      for(int k = 1; k<=n; k++)
        result(i-1,j-1) += xrand(k-1,i-1) * xrand(k-1,j-1) * xweight(k-1);
  return Rcpp::wrap(result);
  '

rcppmulSigmaUVweightavg <- cxxfunction(signature(rand_origin="numeric", weight = "numeric", zsize = "int", ysize = "int"),
                                    mulSigmaUVweightavgcpp,plugin="RcppArmadillo", verbose = TRUE)

#the last k columns are the conditional mean of Z|Y
#before that there is a line for yid
mulYZdrawcpp <- ' 
  arma::mat xdelta = Rcpp::as<arma::mat>(delta_new);
  arma::mat xmu = Rcpp::as<arma::mat>(mu_new);
  arma::mat xprob = Rcpp::as<arma::mat>(Yprob_new);
  arma::mat xtPSigmaZ = Rcpp::as<arma::mat>(tPSigmaZ);
  int zleng = Rcpp::as<int>(zsize);
  int yleng = Rcpp::as<int>(ysize);
  int n = xmu.n_rows;
  arma::mat Ydraw = arma::zeros<arma::mat>(n,zleng + yleng + 1);
  arma::vec rand = Rcpp::as<arma::vec>(rand_vec);
  double temp = 0; 
  for(int i = 1; i<=n; i++)
  {
    temp = xprob(i-1,0);
    if(rand(i-1) <= temp)
    {
      Ydraw(i-1,yleng) = 1;
      continue;
    }
    temp += xprob(i-1,1);
    if(rand(i-1) <= temp)
    {
      Ydraw(i-1,1) = 1;
      Ydraw(i-1,yleng) = 2;
      continue;
    }  
    temp += xprob(i-1,2);
    if(rand(i-1) <= temp)
    {
      Ydraw(i-1,0) = 1;
      Ydraw(i-1,yleng) = 3;
      continue;
    }  
    temp += xprob(i-1,3);
    if(rand(i-1) <= temp)
    {
      Ydraw(i-1,0) = 1;
      Ydraw(i-1,1) = 1;
      Ydraw(i-1,yleng) = 4;
      continue;
    }  
  }
  for(int i = 1; i<= n; i++)
  {
     Ydraw(i-1,yleng + 1) = xdelta(i-1, 0) + Ydraw(i-1,0) * xtPSigmaZ(0,0) + Ydraw(i-1,1)* xtPSigmaZ(1,0);
     Ydraw(i-1,yleng + 2) = xdelta(i-1, 1) + Ydraw(i-1,0) * xtPSigmaZ(0,1) + Ydraw(i-1,1)* xtPSigmaZ(1,1);
  } 

  return Rcpp::wrap(Ydraw);
'

rcppmulYZdraw <- cxxfunction(signature(delta_new = "numeric", mu_new = "numeric",
                                    Yprob_new = "numeric",  tPSigmaZ = "numeric", 
                                    zsize = "int", ysize = "int", rand_vec = "numeric"), 
                          mulYZdrawcpp, plugin = "RcppArmadillo", verbose = TRUE)