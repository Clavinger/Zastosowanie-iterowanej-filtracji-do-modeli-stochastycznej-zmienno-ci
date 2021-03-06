#eksperyment symulacyjny dla modelu basic stochastic volatility

#rm(list = ls())
setwd("C:/Users/user/Documents/github/Zastosowanie-iterowanej-filtracji-do-modeli-stochastycznej-zmienno-ci/Dane")


#pakiety
library(forecast)
library(pomp)
library(beepr)
library(doParallel)
library(doSNOW)
library(nimble)
library(FKF)
##########################################################################
##########################################################################
###################################################ustawienia eksperymentu

#liczba symulacji
sym=2

#liczba obserwacji
n=1000

#parametry
mu = -0.5
phi = 0.98
sigma = 0.2


##########################################################################
##########################################################################
#################################################tworzenie macierzy danych

#metoda x parametr x nr symulaci 
parametry.est=array(data=NA, dim=c(4,3,sym+1),dimnames=list(metoda=c("QML","LW","IF","PCMC"),
                                                            parametry=c("mu","phi","sigma_eta"),
                                                            numer_symulacji=c(1:sym,'srednia')))
                                                            
#metoda x parametr  x  nr symulaci x rodzaj bledu 
parametry.blad=array(data=NA,dim=c(4,3,sym+1,5),dimnames=list(metoda=c("QML","LW","IF","PCMC"),
                                                              parametry=c("mu","phi","sigma_eta"),
                                                              numer_symulacji=c(1:sym,'srednia'),
                                                              blad=c("ME","RMSE","MAE",
                                                                  "MPE","MAPE")))
#metoda x  nr symulaci x nr obserwacji 
log.volatility.est=array(data=NA,dim=c(4,sym+1,n),dimnames=list( metoda=c("QML","LW","IF","PCMC"),
                                                                 numer_symulacji=c(1:sym,'srednia'),
                                                                 numer_obserwacji=1:n))
#metoda x  nr symulaci x rodzaj bledu 
log.volatility.blad=array(data=NA,dim=c(4,sym+1,5),dimnames=list(metoda=c("QML","LW","IF","PCMC"),
                                                                 numer_symulacji=c(1:sym,'srednia'),
                                                                 blad=c("ME","RMSE","MAE",
                                                                 "MPE","MAPE")))
#metoda x  nr symulaci  
czas.obliczen=array(data=NA,dim=c(4,sym+1))
rownames(czas.obliczen)<-c("QML","LW","IF","PCMC")  
colnames(czas.obliczen)<-c(1:sym,"srednia")

#nr symulaci x nr obserwacji 
zwroty.sim=array(data=NA,dim=c(sym+1,n))

#nr symulaci x nr obserwacji 
log.volatility.sim=array(data=NA,dim=c(sym+1,n))

#metoda x  nr symulaci  
log.lik=array(data=NA,dim=c(4,sym+1))
rownames(log.lik)<-c("QML","LW","IF","PCMC")
colnames(log.lik)<-c(1:sym,"srednia")

##########################################################################
##########################################################################
##########################################################################
#######################################################implementacja metod
#1) metoda quasi-najwiekszej wiarygodnosci

OUss <- function(mu,phi,sigma){
  Tt <- matrix(phi,ncol=1)
  Zt <- matrix(1,ncol=1)
  ct <- matrix(-1.27,ncol=1)
  dt <- matrix(mu*(1-phi), nrow = 1)
  GGt<- matrix(pi^2/2,nrow = 1,ncol=1)
  HHt<- matrix(sigma^2,nrow=1,ncol=1)
  a0 <-  as.vector(matrix(mu+1.27))
  P0 <- matrix(1,nrow=1,ncol=1)
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
              HHt = HHt))
}

KF.log <- function(theta) {
  sp <- OUss(theta[1], theta[2], theta[3])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt =matrix(y, nrow=1,ncol=n))
  return(-ans$logLik)
}



KF <- function(theta) {
  sp <- OUss(theta[1], theta[2], theta[3])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt =matrix(y, nrow=1,ncol=n))
  return(ans$att[1,])
}


###############################################################################
#2) metoda iterowanej filtracji 


####nazwy
bsv_statenames <- c("H","Y_state")
bsv_rp_names <- c("mu","phi","sigma_eta")
#bsv_ivp_names <- c("H_0")
#bsv_paramnames <- c(bsv_rp_names,bsv_ivp_names)
bsv_paramnames <- c(bsv_rp_names)
bsv_covarnames <- "covaryt"



rproc1 <- "
double omega;
omega = rnorm(0,sigma_eta );
H = mu*(1 - phi) + phi*H + omega;
"

####rownanie procesu pomiaru
rproc2.sim <- "
Y_state = rnorm( 0,exp(H/2) );
"
###do wypelniania danych
rproc2.filt <- "
Y_state = covaryt;
"

###symulacja modelu SVL
bsv_rproc.sim <- paste(rproc1,rproc2.sim)

####filtr czasteczkowy 
bsv_rproc.filt <- paste(rproc1,rproc2.filt)


######inicalizacja
#H=H_0
bsv_initializer <- "
H = rnorm(mu  ,sigma_eta/sqrt((1-phi*phi))) ;
Y_state = rnorm( 0,exp(H/2) );
"
###????y=Y_state;
bsv_rmeasure <- "
y=rnorm( 0,exp(H/2) );
"

####rozk?ad warunkowy zmiennej Y
bsv_dmeasure <- "
lik=dnorm(y,0,exp(H/2),give_log);
"


####przeskalowanie parametr?w 
bsv_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tphi = logit(phi);
"

bsv_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tphi = expit(phi);
"
bsv.model<-pomp(data=data.frame(y=1:n,
                                time=1:n),
                statenames=bsv_statenames,
                paramnames=bsv_paramnames,
                covarnames=bsv_covarnames,
                times="time",
                t0=0,
                covar=data.frame(covaryt=0:n,
                                 time=0:n),
                tcovar="time",
                rmeasure=Csnippet(bsv_rmeasure),
                dmeasure=Csnippet(bsv_dmeasure),
                rprocess=discrete.time.sim(step.fun=Csnippet(bsv_rproc.filt),delta.t=1),
                initializer=Csnippet(bsv_initializer),
                toEstimationScale=Csnippet(bsv_toEstimationScale), 
                fromEstimationScale=Csnippet(bsv_fromEstimationScale)
)

#######################################ustawienia dla pomp

params_test=c(
  mu=mu,
  phi=phi,
  sigma_eta=sigma
)

###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 1

#liczba czasteczek
bsv_Np <-          c(100,1e3,1e3)
bsv_Nmif <-        c(10, 50,150)
bsv_Nreps_eval <-  c(4,  10,  20)
bsv_Nreps_local <- c(4, 10, 10)
bsv_Nreps_global <-c(4, 10, 10)

bsvlist<-list(bsv_Np ,bsv_Nmif,bsv_Nreps_eval,
              bsv_Nreps_local,bsv_Nreps_global )

#parametry do metody mif2
bsv_rw.sd_rp <- 0.02
bsv_rw.sd_ivp <- 0.1
bsv_cooling.fraction.50 <- 0.5
###############################################################################
#3) Liu-West Filter
stochVCode <- nimbleCode({
  
  
  x[1] ~ dnorm(mu*(1-phi) + phi * x0, var = sigmaSquaredInv)
  y[1] ~ dnorm(0, var =  exp(x[1]))
  for(t in 2:T){
    x[t] ~ dnorm(mu*(1-phi) + phi * x[t-1],  var = sigmaSquaredInv)
    y[t] ~ dnorm(0, var = exp(x[t]))
  }
  
  x0 ~ dnorm(0, var = sigmaSquaredInv)
  phi ~ dunif(-1,1)
  sigmaSquaredInv ~ dgamma(5, 20)
  mu ~ dnorm(-0.5, var = 100)      
})

################################################################################
#4) Particle Markov Chain Monte Carlo

hyperparams <- list(min = c(-1,0.9,0), max = c(0,1,1) )

bsv.dprior <- function (params, ..., log) {
  f <- sum(dunif(params, min = hyperparams$min, max = hyperparams$max,
                 log = TRUE))
  if (log) f else exp(f)
}

##########################################################################
##########################################################################
##########################################################################

log.volatility.list  <-pomp::simulate(bsv.model,nsim=sym,seed=123,params=c(mu=mu,phi=phi,sigma_eta=sigma),
        states = TRUE, obs = TRUE)


start_time_total<- Sys.time()
for(s in 1:sym){
  print(s)
  zwroty.sim[s,]= log.volatility.list$obs[1,s,]
  log.volatility.sim[s,]= log.volatility.list$states[1,s,] 
  

  ###############################################################################
  #3) metoda iterowanej filtracji 
  print(paste("numer iteracji: ",s,"(metoda iterowanej filtracji)"))
  bsv.filt<-pomp(data=data.frame(y=zwroty.sim[s,],
                                  time=1:n),
                  statenames=bsv_statenames,
                  paramnames=bsv_paramnames,
                  covarnames=bsv_covarnames,
                  times="time",
                  t0=0,
                  covar=data.frame(covaryt=c(0,zwroty.sim[s,]),
                                   time=0:n),
                  tcovar="time",
                  rmeasure=Csnippet(bsv_rmeasure),
                  dmeasure=Csnippet(bsv_dmeasure),
                  rprocess=discrete.time.sim(step.fun=Csnippet(bsv_rproc.filt),delta.t=1),
                  initializer=Csnippet(bsv_initializer),
                  toEstimationScale=Csnippet(bsv_toEstimationScale), 
                  fromEstimationScale=Csnippet(bsv_fromEstimationScale)
  )
  
  
  
  detectCores()
  cl <- makeCluster(3, type = "SOCK")
  registerDoSNOW(cl)
  
  
  
  
  start_time <- Sys.time()
  t.if.bsv <- system.time({
    if.bsv <- foreach(i=1:bsvlist[[5]][run_level] ,
                      .packages='pomp', .combine=c,.export = "bsvlist", 
                      .options.multicore=list(set.seed=TRUE)) %dopar% try(
                        pomp::mif2(bsv.filt,start=params_test,Np=bsvlist[[1]][run_level] , Nmif=bsvlist[[2]][run_level] ,cooling.type="geometric",
                                   cooling.fraction.50=bsv_cooling.fraction.50,
                                   transform=TRUE,
                                   rw.sd = rw.sd(
                                     mu      = bsv_rw.sd_rp,
                                     phi       = bsv_rw.sd_rp,
                                     sigma_eta = bsv_rw.sd_rp
                                   )
                        )
                        
                      )
    
    L.if.bsv <- foreach(i=1:bsvlist[[5]][run_level] ,.packages='pomp', .export = "bsvlist", 
                        .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                          logmeanexp(
                            replicate(bsvlist[[3]][run_level] ,
                                      logLik(pfilter(bsv.filt,params=coef(if.bsv[[i]]),Np=bsvlist[[1]][run_level]  ))
                            ),
                            se=TRUE)
                        )
    
    H.if.bsv<- foreach(i=1:bsvlist[[5]][run_level] ,.packages='pomp', .export = "bsvlist", 
                       .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                         exp(pfilter(bsv.filt,params=coef(if.bsv[[i]]),Np=bsvlist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                       )
  })
  
  
  stopCluster(cl)
  end_time <- Sys.time()
  czas.obliczen[3,s]=difftime(end_time,start_time, units = "secs")
  
  r.if.bsv <- data.frame(logLik=L.if.bsv[,1],logLik_se=L.if.bsv[,2],t(sapply(if.bsv,coef)))
  parametry.est[3, ,s]=as.vector(r.if.bsv[which.max(r.if.bsv$logLik),3:5],mode="numeric")
  for(i in 1:3) parametry.blad[3,i,s,]=as.vector(accuracy(f= parametry.est[3,i ,s],x=params_test[i]),mode="numeric")
  log.volatility.est[3,s,]=as.vector(pfilter(bsv.filt,params=parametry.est[3,1:3,s],
          Np=bsvlist[[1]][run_level],filter.traj=T)@filter.traj[1,1,2:(n+1)],mode="numeric")
  log.volatility.blad[3,s,]=as.vector(accuracy(f=log.volatility.est[3,s,],x=log.volatility.sim[s,]),mode="numeric")
  log.lik[3,s]=as.numeric(r.if.bsv[which.max(r.if.bsv$logLik),1])
  
  ###############################################################################
  #1) metoda quasi-najwiekszej wiarygodnosci
  print(paste("numer iteracji: ",s,"(metoda quasi-najwiekszej wiarygodnosci)"))
  y=log(zwroty.sim[s,]^2)
  start_time <- Sys.time()
  KF.opt<-optim(c(mu,phi,sigma),KF.log, 
                method="L-BFGS-B",hessian = T,lower=c(-Inf,-1,0),upper=c(Inf,1,Inf))
  end_time <- Sys.time()
  czas.obliczen[1,s]=difftime(end_time,start_time, units = "secs")
  
  
  parametry.est[1, ,s]=KF.opt$par
  for(i in 1:3) parametry.blad[1,i,s,]=as.vector(accuracy(f= parametry.est[1,i ,s],x=params_test[i]),mode="numeric")
  log.volatility.est[1,s,]=(KF(KF.opt$par))
  log.volatility.blad[1,s,]=as.vector(accuracy(f=log.volatility.est[1,s,],x=log.volatility.sim[s,]),mode="numeric")
  log.lik[1,s]=logLik(pfilter(bsv.filt,params=parametry.est[1,1:3,s],
                       Np=bsvlist[[1]][run_level]))
  
  ###############################################################################
  
  #2) Liu and West
  print(paste("numer iteracji: ",s,"(Liu and West)"))
  start_time <- Sys.time()
  stochVolModel <- nimbleModel(code = stochVCode, name ='stochVol',
                               constants = list(T = n),  data = list(y = as.vector(zwroty.sim[s,])),
                               inits = list(mu = mu,  phi = phi,
                                            sigmaSquaredInv = sigma^2))
  CstochVolModel <- compileNimble(stochVolModel)
  stochVolLiuWestFilter <- buildLiuWestFilter(model = stochVolModel, nodes ='x', 
                                              params = c("mu","phi","sigmaSquaredInv"))
  
  CstochVolLiuWestFilter <- compileNimble(stochVolLiuWestFilter,
                                          project = stochVolModel)
  CstochVolLiuWestFilter$run(1000)
  end_time <- Sys.time()
  czas.obliczen[2,s]=difftime(end_time,start_time, units = "secs")
  muSamples <- as.matrix(CstochVolLiuWestFilter$mvEWSamples,'mu')
  phiSamples <- as.matrix(CstochVolLiuWestFilter$mvEWSamples,'phi')
  sigmaSamples <- as.matrix(CstochVolLiuWestFilter$mvEWSamples,'sigmaSquaredInv')

  parametry.est[2, ,s]=c(mean(muSamples),mean(phiSamples),mean(sigmaSamples))
  for(i in 1:3) parametry.blad[2,i,s,]=as.vector(accuracy(f= parametry.est[2,i ,s],x=params_test[i]),mode="numeric")
  log.volatility.est[2,s,]=as.vector(pfilter(bsv.filt,params=parametry.est[2,1:3,s],
                                             Np=bsvlist[[1]][run_level],filter.traj=T)@filter.traj[1,1,2:(n+1)],mode="numeric")
  log.volatility.blad[2,s,]=as.vector(accuracy(f=log.volatility.est[2,s,],x=log.volatility.sim[s,]),mode="numeric")
  log.lik[2,s]=logLik(pfilter(bsv.filt,params=parametry.est[2,1:3,s],
                              Np=bsvlist[[1]][run_level]))
  
  
  ###############################################################################
  
  #4) PMCMC
  print(paste("numer iteracji: ",s,"(PMCMC)"))
  
  hyperparams <- list(min = c(-1,0.9,0), max = c(1,1,1) )
  bsv.dprior <- function (params, ..., log) {
    f <- sum(dunif(params, min = hyperparams$min, max = hyperparams$max,
                   log = TRUE))
    if (log) f else exp(f)
  }
  
  start_time <- Sys.time()
  pmcmc1 <-   pmcmc(pomp(bsv.filt, dprior = bsv.dprior), start = params_test,
                    Nmcmc = 2000, Np = 100, max.fail = Inf,
                    proposal = mvn.diag.rw(c(mu = 0.01, phi = 0.01, sigma_eta = 0.01)))
  
  continue( pmcmc1 ,Nmcmc=3000,proposal=mvn.rw(covmat( pmcmc1 ))) -> pmcmc1 
  end_time <- Sys.time()
  czas.obliczen[4,s]=difftime(end_time,start_time, units = "secs")
  
  parametry.est[4, ,s]=as.vector(coef(pmcmc1 ),mode="numeric")
  for(i in 1:3) parametry.blad[4,i,s,]=as.vector(accuracy(f= parametry.est[4,i ,s],x=params_test[i]),mode="numeric")

  log.volatility.est[4,s,]=as.vector(pfilter(bsv.filt,params=parametry.est[4,1:3,s],
                                             Np=bsvlist[[1]][run_level],filter.traj=T)@filter.traj[1,1,2:(n+1)],mode="numeric")
  
  log.volatility.blad[4,s,]=as.vector(accuracy(f=log.volatility.est[4,s,],x=log.volatility.sim[s,]),mode="numeric")
  log.lik[4,s]=logLik(pmcmc1 )
  
  #koniec petli liczba symulacji
  }
  end_time_total <- Sys.time()
  difftime(end_time_total,start_time_total, units = "mins")
  beep(1)
 
  
  for (i in 1:4){
    parametry.est[i,,sym+1]=apply(parametry.est[i,,1:sym],MARGIN = 1,FUN=mean)
    for(j in 1:5) parametry.blad[i, , sym+1, j]=apply(parametry.blad[i,,1:sym,j],MARGIN = 1,FUN=mean)
    log.volatility.blad[i,sym+1,]=apply(log.volatility.blad[i,1:sym,],MARGIN = 2,FUN=mean)
    log.lik[i,sym+1]=mean(log.lik[i,1:sym])
    czas.obliczen[i,sym+1]=mean(czas.obliczen[i,1:sym])
  }
  #save(parametry.est,parametry.blad,log.volatility.est,log.volatility.blad,czas.obliczen,
   #   zwroty.sim,log.volatility.sim,log.lik,file="bsv_wyniki")
  #load(file="bsv_wyniki")
  
  ##########################################################################
  ##########################################################################
  ##########################################################################
  ######################################################podsumowanie wyników
  
  ##############################################################parametr: mu
  mu
  #oszacowania paramaetrów
  parametry.est[ ,1, ]
  #srednie błedy oszacowan
  parametry.blad[ ,1 , sym+1, ]
  #############################################################parametr: phi
  phi
  #oszacowania paramaetrów
  parametry.est[ ,2, ]
  #srednie błedy oszacowan
  parametry.blad[ ,2 , sym+1, ]
  #######################################################parametr: sigma_eta
  sigma
  #oszacowania paramaetrów
  parametry.est[ ,3, ]
  #srednie błedy oszacowan
  parametry.blad[ ,3 , sym+1, ]  
  
  
  #blad oszaowania procesu ukrytego
  log.volatility.blad[,sym+1,]
  
  #wartosci funkcji wiarygodnosci
  log.lik

  #czas obliczen
  czas.obliczen
  difftime(end_time_total,start_time_total, units = "mins")
  