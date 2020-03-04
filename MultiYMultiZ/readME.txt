The codes under this folder is used for the simulation study which combines two binary variables and two continuous variables.
However, one can change the dimension of discrete variables and continuous variables and all parameters easily to construct different simualtion cases. 

To run the programme correctly, please make sure that the version of R is 3.6.0 or higher. 
The required packages are "dplyr", "tidyr", "mvtnorm", "lme4", "parallel", "Rcpp", "inline", "RcppArmadillo".

Source the other function files before run the main function in runME_multi.R

The simulation raw dataset is an RDS file called multiP0700W03G0full.rds, which is a list, each element contains two parts: 
1. The first part is called result, which is a 560 * 804 matrix. We assume 40 small areas and there are 14 small area parameters in each area. Each row denotes one parameter in one small area. The first 4 columns are true group means, sample group means, group means from univariate approach and group means from multivariate approach. The next 600 columns are the group mean estimates from 600 EBP replicates. Then the last 200 columns are the bootstrap means and bootstrap variance in turn of 100 bootstrap replication study. 
2. The second part is called edge which is a list of 101 4 * 4 matrices. 

The 'relabias' function and 'reportMSE' function are used to analyze the simulation raw data which provides us a table of relative bias (RB) and coverage rate (CR) and a table of MSE results from three approaches and bootstrap procedure. 