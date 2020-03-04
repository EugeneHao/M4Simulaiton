The codes under this folder is used for the simulation study which combines one categorical variables and one continuous variables.

To run the programme correctly, please make sure that the version of R is 3.6.0 or higher. 
The required packages are "dplyr", "tidyr", "mvtnorm", "lme4", "parallel", "Rcpp", "inline", "RcppArmadillo".

Source the other function files before run the main function in runME_oneZoneY.R

The simulation raw dataset is an RDS file called smallcorresult.rds which is a list, each element is a 280 * 804 matrix: 
We assume 40 small areas and there are 7 small area parameters in each area. Each row denotes one parameter in one small area. The first 4 columns are true group means, sample group means, group means from univariate approach and group means from multivariate approach. The next 600 columns are the group mean estimates from 600 EBP replicates. Then the last 200 columns are the bootstrap means and bootstrap variance in turn of 100 bootstrap replication study. 


The 'relabias' function and 'reportMSE' function are used to analyze the simulation raw data which provides us a table of relative bias (RB) and coverage rate (CR) and a table of MSE results from three approaches and bootstrap procedure. 