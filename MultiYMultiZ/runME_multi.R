#This code is used for multi continuous variables and multi categorical variables: 

library('dplyr')
library('tidyr')
library('mvtnorm')
library('lme4')
library('parallel')
library('MASS')
library('Rcpp')
library('inline')
library('RcppArmadillo')
source("base_multi.R")
source("multi_cpp.R")
source("analyse_funs.R")
#source("cpp_update.R")
#source("cpphier.R")

result = NULL   
filename = paste("multiP0700W03G0full", ".rds", sep = "")
saveRDS(result, file = filename)

############################################################
m = 1000              # Monte Carlo size
t = 600               # EBP size
b = 100               # bootstrap size

C = 4                 # number of all combinations of categorical variables
K = 2                 # dimension of Z
L = 2                 # dimension of Y
Cl = c(1,1)           # length of dummy variables 

EMsize = 20           # iteraton times for EM alrorithm
cores = 50            # the number of cores used 
MCsize = 100          # simulation times
seed = 1:100          # set the seed
fulledge = TRUE       # TRUE means not to do edge selection



D = 40                               # number of groups 
Pop = 2000                           # population size 
samplesize = rep(c(10, 30, 50, 70), each = 10) 
popsize = rep(Pop, each = D) 
predictsize = popsize - samplesize


categtotal = prod(Cl+1)
paramtotal = K + C +C*K              # total number of parameters

categtotal = prod(Cl+1)

###########################################################
Bz = matrix(c(1, 0.2, -1, -0.2), nrow = 2)         # fixed coefficients for (1,X) in conditonal model of Z, 
By = matrix(c(-1, 1,
              0.5, -0.5), nrow = 2, byrow = T)     # fixed coefficients for (1,X) in conditonal model of Y , 

P = matrix(c(0.7, 0, 0, -0.7), nrow = 2)           # fixed coefficients for Z in the conditonal model of Y
W = matrix(c(0, 0.3, 0.3, 0), nrow = 2)            # fixed coefficients for Y in the conditonal model of Y

GammaZ = diag(2)                                   # conditional precision matrix of Z

SigmaUV = diag(rep(0.5,C)) %*% 
  matrix(c(1, -0.5, -0.1, 0.1,  
           -0.5, 1, -0.1, 0.1, 
           -0.1,-0.1, 0.25, -0.05, 
           0.1, 0.1, -0.05, 0.25), nrow = C) %*%
  diag(rep(0.5, C))                                #variance of random effects


Ychoice = matrix(c(0,0,
                   0,1,
                   1,0,
                   1,1), nrow = C, byrow = T)      #possible choices for Y (C * (C-L) matrix)
ysize = 2
zsize = 2

#Gauss-Legendre quadrature
X_k = c(-0.1488743389816312, 0.1488743389816312,
        -0.4333953941292472, 0.4333953941292472,
        -0.6794095682990244, 0.6794095682990244,
        -0.8650633666889845, 0.8650633666889845,
        -0.9739065285171717, 0.9739065285171717)

W_k = c(0.2955242247147529, 0.2955242247147529,
        0.2692667193099963, 0.2692667193099963,
        0.2190863625159820, 0.2190863625159820,
        0.1494513491505806, 0.1494513491505806,
        0.0666713443086881, 0.0666713443086881)

####################################################################
# we update the simulation result after every simulation study 
result = list()
for(i in 1: MCsize)
{
  temp <-oneSimulation(m, t, b, seed[i], EMsize, fulledge)
  result[[i]] <- temp
  saveRDS(result, file = filename)
}



########################################################
#analysis part: 
varname <- c("Z1", "Z2", "Y10", "Y11", "Y20", "Y21", "Z1|Y10", "Z1|Y11", "Z1|Y20", "Z1|Y21", 
             "Z2|Y10", "Z2|Y11", "Z2|Y20", "Z2|Y21")

rb = relabias(result, D = D, m = m, t = t, b = b, varname = varname)

# Table of relative bias and coverage rate
rb %>% round(.,4)  %>% "colnames<-"( varname) %>% "rownames<-"(c("RB", "CR"))

# Table of MSEs from direct estimate, univariate model, M4 model and bootstrap estimation 
reportMSE(result, D = D, m = m, t = t, b = b, varname)
