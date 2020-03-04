#This code is used for one continuous variable and one categorical variable: 

library('dplyr')
library('tidyr')
library('mvtnorm')
library('lme4')
library('parallel')
library('MASS')
library('Rcpp')
library('inline')
library('RcppArmadillo')  
source("basefun.R")                            #source three codes 
source("cpp.R")
source("cpp_update.R")
#source("cpphier.R")

#setwd("/home/hao123/GauMultiJoint/turning1126") #/turningparam
result = NULL   
filename = paste("smallcorresult", ".rds", sep = "")     # the name of this file shows that we use a mild correlation 
                                                       # others can change the parameters setting below to do different simulation studies
saveRDS(result, file = filename)

##############################################
m = 1000              # Monte Carlo size
t = 400               # EBP size
b = 100               # bootstrap size
MCsize = 100          # simulation times
seed = 1:100          #set the seed

EMsize = 20           # iteraton times for EM alrorithm
update_all = FALSE    # use EM algorithm to update the fixed coefficient estimates or not 
cores = 50            # the number of cores used 
Ychoice = matrix(c(0,0,
                   1,0,
                   0,1), nrow = 3, byrow = T)    #possible choices for Y (3 categories, so 2 dimension, 3 possible vector)

D = 40                                 # number of groups 
Pop = 2000                             # population size 
samplesize = rep(c(10, 30, 50, 70), each = 10) 
popsize = rep(Pop, each = D) 
fittingsize = popsize - samplesize
numornot = 1                       # use numerical integration approach to find the initial estiamtes
predictsize = rep(Pop, D) - samplesize 

################################################
#Bz = [Bz1, Bz2]
Bz = c(1, 0.2)                      # fixed coefficients for (1,X) in conditonal model of Z, 

By = matrix(c(-1, 0.5,
              0.5, -1), nrow = 2)   # fixed coefficients for (1,X) in conditonal model of Y,
P = c(0.3, -0.3)                    # fixed coefficients for Z in the conditonal model of Y
  
#P = c(0, 0)         #independent 
#P = c(0.7, -0.7)    #strong correlation 

SigmaZ = 1                         #conditonal variance of Z

SigmaUV = diag(rep(0.5,3)) %*% 
  matrix(c(1, -0.1, 0.1,  
           -0.1, 0.25, -0.05, 
           0.1, -0.05, 0.25), nrow = 3) %*%
  diag(rep(0.5, 3))               #variance of random effects

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

# we update the simulation result after every simulation study 
result = list()
for(i in 1: MCsize)
{
    temp <-onestep(m, t, b, seed[i], EMsize, numornot)
       result[[i]] <- temp
    saveRDS(result, file = filename)
}

#mclapply(seed, FUN = function(x) onestep(m, t, b, seed[x], EMsize, numornot), mc.cores = cores) -> result

#saveRDS(result, file = filename)

########################################################
#analysis part: 
varname <- c("Z1", "Y10", "Y11", "Y12", "Z1|Y10", "Z1|Y11", "Z1|Y12")

rb = relabias(result, D = D, m = m, t = t, b = b, varname = varname)

# Table of relative bias and coverage rate
rb %>% round(.,4)  %>% "colnames<-"( varname) %>% "rownames<-"(c("RB", "CR"))

# Table of MSEs from direct estimate, univariate model, M4 model and bootstrap estimation 
reportMSE(result, D = D, m = m, t = t, b = b, varname)
