

##### Output for simulations to compare estimators:

load("StrongCovs.Rdata")

####  Set maxcnt to last iteration number
maxcnt <- 100

biasfun <- function(a,b){ 
	apply(a-b,2,mean, na.rm = TRUE )
}

msehatfun <- function(a,b){ 
	apply( (a-b)^2,2,mean, na.rm = TRUE )
}

 



popestlist <- vector("list",4)
popestlist[[1]] <- zpops[1:maxcnt,]
popestlist[[2]] <- ypops[1:maxcnt,]
popestlist[[3]] <- ratpops[1:maxcnt,]
popestlist[[4]] <- rat0pops[1:maxcnt,]

indestlist <- vector("list",4)
indestlist[[1]] <- zmuEBPINDs[1:maxcnt,]
indestlist[[2]] <- ymuEBPINDs[1:maxcnt,]
indestlist[[3]] <- ratINDs[1:maxcnt,]
indestlist[[4]] <- rat0INDs[1:maxcnt,]

initestlist <- vector("list",4)
initestlist[[1]] <- zmuEBPINITs[1:maxcnt,]
initestlist[[2]] <- ymuEBPINITs[1:maxcnt,]
initestlist[[3]] <- ratEBPINITs[1:maxcnt,]
initestlist[[4]] <- rat0EBPINITs[1:maxcnt,]

msecompare <- c()
for(i in 1:4){
 	msecompare <- cbind(msecompare, 
	#tapply(msehatfun(directestlist[[i]], popestlist[[i]]), Nis, mean),
	tapply(msehatfun(indestlist[[i]], popestlist[[i]]), Nis, mean),
	tapply(msehatfun(initestlist[[i]], popestlist[[i]]), Nis, mean)

	)
}

biascompare <- c()
for(i in 1:4){
 	biascompare <- cbind(biascompare, 
#	tapply(biasfun(directestlist[[i]], popestlist[[i]])^2, Nis, mean),
	tapply(biasfun(indestlist[[i]], popestlist[[i]])^2, Nis, mean),
	tapply(biasfun(initestlist[[i]], popestlist[[i]])^2, Nis, mean)

	)
}

sdcompare <- c()
for(i in 1:4){
 	sdcompare <- cbind(sdcompare, 
#	tapply(biasfun(directestlist[[i]], popestlist[[i]]), Nis, mean),
	tapply(biasfun(indestlist[[i]], popestlist[[i]]), Nis, mean),
	tapply(biasfun(initestlist[[i]], popestlist[[i]]), Nis, mean)

	)
}

library("xtable")

rbind(round(msecompare*100,2),
 round(biascompare/msecompare*100,5)
 )

 xtable(rbind(round(msecompare*100,2),
  round(biascompare/msecompare*100,2)
  ))



tabbar <- round(
	100*rbind(apply( msecompare ,2, mean),
apply(biascompare , 2, mean) /apply( msecompare ,2, mean)
)
,2)

xtable(tabbar, digits = 2)


