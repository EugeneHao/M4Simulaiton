groupdenscpp <-'
arma::vec xP0 = Rcpp::as<arma::vec>(P0);
double xsigmaZ = Rcpp::as<double>(SigmaZ_hat);
arma::cube xarray = Rcpp::as<arma::cube>(rand_matrix);
arma::vec xdelta = Rcpp::as<arma::vec>(delta_fix);
arma::mat xmu = Rcpp::as<arma::mat>(mu_fix);
arma::vec xid = Rcpp::as<arma::vec>(id);
arma::vec xAfix = Rcpp::as<arma::vec>(A_fix);
arma::vec xylabel = Rcpp::as<arma::vec>(ylabel);
arma::vec xprob = Rcpp::as<arma::vec>(probfix);
int xm = Rcpp::as<int>(m);
arma::mat xYchoice = Rcpp::as<arma::mat>(Ychoice);
int n = xid.n_rows;
int xD = Rcpp::as<int>(D);
arma::mat temp = arma::zeros<arma::mat>(n,3);
arma::vec temp2 = arma::zeros<arma::vec>(2);
double delta;
double pick;
double sumprob;
double kernel; 
arma::mat result = arma::zeros<arma::mat>(xD, xm);
double check;
for(int i = 1; i<=xm; i++)
{
  for(int j = 1; j<=n; j++)
  {
  delta = xarray(xid(j-1)-1, i-1, 0) + xdelta(j-1);
  for(int l = 1; l<=2; l++)
  {
  temp2(l-1) = delta * xP0(l-1) + xarray(xid(j-1)-1, i-1, l) + xmu(j-1,l-1);
  }
  sumprob = 0;
  for(int k = 1; k <=3; k++) 
  { 
  temp(j-1,k-1) = temp2(0) * xYchoice(k-1, 0) + temp2(1) * xYchoice(k-1,1) + xprob(k-1);
  sumprob += exp(temp(j-1,k-1));
  }
  pick = temp(j-1, xylabel(j-1)-1);
  kernel = -pow(xAfix(j-1) + delta, 2)/xsigmaZ/2;
  result(xid(j-1)-1, i-1) += pick - log(sumprob) + kernel;
  } 
}
return Rcpp::wrap(result);
'

rcppgroupdens <- cxxfunction(signature(P0="numeric", SigmaZ_hat = "numeric", rand_matrix = "numeric", delta_fix = "numeric", 
                                   mu_fix = "numeric", id = "int", A_fix = "numeric", 
                                   ylabel = "int", probfix = "numeric", m = "int", 
                                   Ychoice = "int",D = "int"),
                      groupdenscpp,plugin="RcppArmadillo", verbose = TRUE)

SigmaUVweightavgcpp <- ' 
  arma::mat xrand = Rcpp::as<arma::mat>(rand_origin);
  arma::vec xweight = Rcpp::as<arma::vec>(weight);
  int n = xweight.n_rows;
  arma::mat result = arma::zeros<arma::mat>(3,3);
  for(int i = 1; i<=3; i++)
    for(int j = 1; j<=3; j++)
      for(int k = 1; k<=n; k++)
        result(i-1,j-1) += xrand(k-1,i-1) * xrand(k-1,j-1) * xweight(k-1);
  return Rcpp::wrap(result);
  '

rcppSigmaUVweightavg <- cxxfunction(signature(rand_origin="numeric", weight = "numeric"),
                           SigmaUVweightavgcpp,plugin="RcppArmadillo", verbose = TRUE)


YZdrawcpp <- ' 
  arma::vec xdelta = Rcpp::as<arma::vec>(delta_new);
  arma::mat xmu = Rcpp::as<arma::mat>(mu_new);
  arma::mat xprob = Rcpp::as<arma::mat>(Yprob_new);
  arma::vec xP = Rcpp::as<arma::vec>(P_hat);
  double xsimgaZ = Rcpp::as<double>(SigmaZ_hat);
  int n = xmu.n_rows;
  arma::mat Ydraw = arma::zeros<arma::mat>(n,4);
  arma::vec rand = Rcpp::as<arma::vec>(rand_vec);
  for(int i = 1; i<=n; i++)
  {
    if(rand(i-1) < xprob(i-1, 2))
    {
      Ydraw(i-1,1) = 1;
      Ydraw(i-1,2) = 3;
    }
    else if(rand(i-1) < xprob(i-1,1) + xprob(i-1,2) && rand(i-1) >= xprob(i-1, 2))
    {
      Ydraw(i-1,0) = 1;
      Ydraw(i-1,2) = 2; 
    }
    else
      Ydraw(i-1,2) = 1;
    Ydraw(i-1,3) = xdelta(i-1) + (Ydraw(i-1,0) * xP(0) + Ydraw(i-1,1)*xP(1)) * xsimgaZ;
  }
  return Rcpp::wrap(Ydraw);
'

rcppYZdraw <- cxxfunction(signature(delta_new = "numeric", mu_new = "numeric",
  Yprob_new = "numeric", P_hat = "numeric", SigmaZ_hat = "numeric", rand_vec = "numeric"), 
  YZdrawcpp, plugin = "RcppArmadillo", verbose = TRUE)