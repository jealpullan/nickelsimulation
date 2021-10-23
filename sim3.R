tinput <- data.frame(DOC=2, Ca=2.2, Mg=1,pH=7.3) #input parameters
library(tlpckg)

test_dist <- mlr(input=tinput) #original sensitivity data
sdvals<-read.csv(file.choose(),header=T) # standard errors

testfit<-ssd_fit_dists(test_dist, dists = c("weibull", "lnorm", "llogis", "burrIII3"))
ssd_gof(testfit)

#--------------------------------------------------------------------------#
#Coefficients
#the following code creates simulations of randomness around regression coefficients

#function for coefficient simulations
simMLR <- function (i){ #i 2:DOC, 3:Ca, 4:Mg, 5:pH, 6:DOC.pH, 7:Mg.pH
  newmlr=MLR
    for(j in 1:4){
      newmlr[j,i] <- suppressWarnings(rnorm(1, mean = MLR[j,i], sd = sdvals[j,i]))
    }
  newmlr
}

simMLR <- function (i, j){ #i 2:DOC, 3:Ca, 4:Mg, 5:pH, 6:DOC.pH, 7:Mg.pH #j 1:fish, 2:inv, 3:algae, 4: plants
  newmlr=MLR
    newmlr[j,i] <- suppressWarnings(rnorm(1, mean = MLR[j,i], sd = sdvals[j,i]))
  newmlr
}

#creates multiple models
repFun <- function(N,m){
  resens <- as.data.frame(lapply(1:N, function(i) mlr(tMLR = simMLR(m), input = tinput)))
  colnames(resens)<-rep("Conc", N)
  resens
}

#creates and fits multiple models
repFits<- function(dataConc, fitd){
  N<-ncol(dataConc)
  fitdists <- rep(ssd_fit_dists(dataConc[1], dists=fitd), N)
  for(i in 2:N){
    fitdists[i] <- try(ssd_fit_dists(dataConc[i], dists=fitd)) #try bc burr
  }
  predsout <- data.frame(rep(1:99,N))
  for(i in 1:N){
    predsout[i] <- try(predict(fitdists[[i]])$est)
  }
  predsout <- cbind(rep(1:99,N), stack(predsout))
  colnames(predsout) <- c("percent", "est", "grp")
  predsout
}

#--------------------------------------------------------------------------#
#Sensitivity
#The following code creates simulations of randomness for sensitivity data

sdvz<-test_dist*2 #sensitivity stretched by factor of 2

#simulate sensitivity values
sens<-function(){
  newSens<-test_dist
  for(i in 1:26){
    newSens[i,] <- abs(suppressWarnings(rnorm(1, mean = test_dist[i,], sd = sdvz[i,])))
  }
  newSens
}

#N repeats of simulatio
repSens<-function(N){ ns
  reSensitivity<-as.data.frame(lapply(1:N, function(i) sens()))
  colnames(reSensitivity)<-rep("Conc", N)
  reSensitivity
}

#outputs model predictions
simSens <- function(dataConc, distr){
  N<-ncol(dataConc)
  fitdists <- rep(ssd_fit_dists(dataConc[1], dists=distr), N)
  for(i in 2:N){
    fitdists[i] <- try(ssd_fit_dists(dataConc[i], dists=distr, computable = T))
  }
  preds <- data.frame(1:99)
  predsout <- data.frame(rep(preds,N))
  for(j in 1:N){
    predsout[j] <- try(predict(fitdists[[j]])$est)
  }
  predsout <- cbind(rep(1:99,N), stack(predsout))
  colnames(predsout) <- c("percent", "est", "grp")
  predsout$est <- as.numeric(predsout$est)
  predsout[!(is.na(predsout$est)),]
}

#--------------------------------------------------------------------------#
#Chapter 3

dfDOC<-repFun(1000,2) #1000 sims, changing DOC
burrfitDOC<-repFits(dfDOC, "burrIII3") #runs multiple sims
bfit<-ssd_fit_dists(test_dist, dists="burrIII3") #original fit
bpred<-predict(bfit) #original predictions

ggplot(burrfitDOC, aes(x=est, y=percent, group=grp))+geom_line()+ scale_x_log10() +
  geom_line(data=bpred, aes(x=as.numeric(est), y=percent), colour="red", inherit.aes = F) +
  labs(title="DOC")

dfCa<-repFun(1000,3) #changing Ca
burrfitCa<-repFits(dfCa, "burrIII3") #runs multiple sims

b2<-ggplot(burrfitCa, aes(x=est, y=percent, group=grp))+geom_line()+ scale_x_log10() +
  geom_line(data=bpred, aes(x=as.numeric(est), y=percent), colour="red", inherit.aes = F) +
  labs(title="Ca")

dfMg<-repFun(1000,4) #changing Mg
burrfitMg<-repFits(dfMg, "burrIII3") #runs multiple sims

b3<-ggplot(burrfitMg, aes(x=est, y=percent, group=grp))+geom_line()+ scale_x_log10() +
  geom_line(data=bpred, aes(x=as.numeric(est), y=percent), colour="red", inherit.aes = F) +
  labs(title="Mg")

dfpH<-repFun(1000,5) #changing pH
burrfitpH<-repFits(dfpH, "burrIII3") #runs multiple sims

b4<-ggplot(dd, aes(x=as.numeric(est), y=percent, group=grp))+geom_line()+ scale_x_log10() +
  geom_line(data=bpred, aes(x=as.numeric(est), y=percent), colour="red", inherit.aes = F) +
  labs(title="pH")

dfDOCpH<-repFun(1000,6) #changing DOC.pH
burrfitDOCpH<-repFits(dfDOCpH, "burrIII3") #runs multiple sims

b5<-ggplot(burrfitDOCpH, aes(x=est, y=percent, group=grp))+geom_line()+ scale_x_log10() +
  geom_line(data=bpred, aes(x=as.numeric(est), y=percent), colour="red", inherit.aes = F) +
  labs(title="DOCpH")

dfMgpH<-repFun(1000,7) #changing Mg.pH
burrfitMgpH<-repFits(dfMgpH, "burrIII3") #runs multiple sims

b6<-ggplot(burrfitMgpH, aes(x=est, y=percent, group=grp))+geom_line()+ scale_x_log10() +
  geom_line(data=bpred, aes(x=as.numeric(est), y=percent), colour="red", inherit.aes = F) +
  labs(title="MgpH")

#burr only sensitivity sim
tt<-repSens(1000)

bsfit<-simSens(tt, "burrIII3")

ggplot(bsfit) + geom_line(aes(y = percent, x = as.numeric(est), group = grp)) + scale_x_log10() +
  geom_line(data=bpred, aes(x=as.numeric(est), y=percent), colour="yellow", inherit.aes = F, size=1.2)



#--------------------------------------------------------------------------#
#Chapter 4
#model averaging coefficients
#Only changing DOC coefficient

#gamma+lnorm+llogis+burr
avgdf<-repFun(1000,2)

avgfit<-repFits(avgdf, c("gamma", "lnorm", "llogis", "burrIII3"))

testfit<-ssd_fit_dists(test_dist, dists=c("gamma", "lnorm", "llogis", "burrIII3"))
testpred<-predict(testfit)

ggplot(avgfit, aes(x=est, y=percent, group=grp))+geom_line()+ scale_x_log10()+
  geom_line(data=testpred, aes(x=as.numeric(est), y=percent), colour="red", inherit.aes = F)+
  labs(title="average")

#weibull+llogis+burr
avg2 <- repFits(avgdf, c("weibull", "llogis", "burrIII3"))
testfit2<-ssd_fit_dists(test_dist, dists=c("gamma", "llogis", "burrIII3"))
testpred2<-predict(testfit)

ggplot(avg2, aes(x=as.numeric(est), y=percent, group=grp))+geom_line()+ scale_x_log10()+
  geom_line(data=testpred2, aes(x=as.numeric(est), y=percent), colour="red", inherit.aes = F)

#weibull+burr
avg3 <- repFits(avgdf, c("weibull", "burrIII3"))
testfit3<-ssd_fit_dists(test_dist, dists=c("weibull", "burrIII3"))
testpred3<-predict(testfit)
ggplot(avg3, aes(x=as.numeric(est), y=percent, group=grp))+geom_line()+ scale_x_log10()+
  geom_line(data=testpred3, aes(x=as.numeric(est), y=percent), colour="red", inherit.aes = F)

#model avgeraging using senstivity data
#burr vs weibull+lnorm+llogis+burr

avgSim <- simSens(tt, c("weibull", "lnorm", "llogis", "burrIII3"))

ggplot(avgSim) + geom_line(aes(y = percent, x = est, group = grp)) + scale_x_log10() +
  geom_line(data=testpred, aes(x=as.numeric(est), y=percent), colour="yellow", inherit.aes = F)

#weibull+llogis+burr
avgSim <- simSens(tt, c("weibull", "llogis", "burrIII3"))
avgSim$est<-as.numeric(avgSim$est)

g3<-ggplot(avgSim2) + geom_line(aes(y = percent, x = est, group = grp)) + scale_x_log10() +
  geom_line(data=testpred2, aes(x=as.numeric(est), y=percent), colour="yellow", inherit.aes = F)


#weibull+burr
avgSim3 <- simSens(tt, c("weibull", "burrIII3"))

g4<-ggplot(avgSim3) + geom_line(aes(y = percent, x = est, group = grp)) + scale_x_log10() +
  geom_line(data=testpred3, aes(x=as.numeric(est), y=percent), colour="yellow", inherit.aes = F)


