#model fat stochastic volatility

#rm(list = ls())
#gc()
setwd("C:/Users/user/Documents/github/Zastosowanie-iterowanej-filtracji-do-modeli-stochastycznej-zmienno-ci/Dane")

#dlugosc szeregu czasowego
n=1000
#parametry
mu = -0.5
phi = 0.98
sigma = 0.2
gamma=5

#pakiety
library(pomp)
library(beepr)
library(doParallel)
library(doSNOW)
library(lattice)
library(ggplot2)
library(GGally)
library(gridExtra)
library(reshape2)
library(FKF)

#FSV specyfikacja do symulacji

FSV_statenames <- c("H","G","Y_state")
FSV_rp_names <- c("mu","phi","gamma","sigma_eta")
FSV_paramnames <- FSV_rp_names
FSV_covarnames <- "covaryt"


rproc1 <- "
double lambda,omega;
omega = rnorm(0,sigma_eta);
lambda = rchisq(gamma) ;
H = mu*(1 - phi) + phi*H + omega;
G = lambda/gamma ;
"

####rownanie procesu pomiaru
rproc2.sim <- "
Y_state = rnorm( 0, exp(H/2)/sqrt(G) );
"
###do wypelniania danych
rproc2.filt <- "
Y_state = covaryt;
"

###symulacja modelu SVL
FSV_rproc.sim <- paste(rproc1,rproc2.sim)

####filtr czasteczkowy 
FSV_rproc.filt <- paste(rproc1,rproc2.filt)


######inicalizacja
FSV_initializer <- "
H = rnorm(mu  ,sigma_eta/sqrt((1-phi*phi))) ;
Y_state = rnorm( 0,exp(H/2)/sqrt(G)  );
"
###????
FSV_rmeasure <- "
y=rnorm( 0,exp(H/2)/sqrt(G)  );
"

####rozk?ad warunkowy zmiennej Y
FSV_dmeasure <- "
lik=dnorm(y,0,exp(H/2)/sqrt(G),give_log);
"


####przeskalowanie parametr?w 
FSV_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tphi = logit(phi);
Tgamma = log(gamma);
"

FSV_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tphi = expit(phi);
Tgamma= exp(gamma);
"



####wypelnianie modelu danymi
FSV.model <- pomp(data=data.frame(y=1:n,
                                 time=1:n),
                 statenames=FSV_statenames,
                 paramnames=FSV_paramnames,
                 covarnames=FSV_covarnames,
                 times="time",
                 t0=0,
                 covar=data.frame(covaryt=c(0,1:n),
                                  time=0:n),
                 tcovar="time",
                 rmeasure=Csnippet(FSV_rmeasure),
                 dmeasure=Csnippet(FSV_dmeasure),
                 rprocess=discrete.time.sim(step.fun=Csnippet(FSV_rproc.filt),delta.t=1),
                 initializer=Csnippet(FSV_initializer),
                 toEstimationScale=Csnippet(FSV_toEstimationScale), 
                 fromEstimationScale=Csnippet(FSV_fromEstimationScale)
)

plot(FSV.model)


#######################################ustawienia dla pomp

params_test=c(
  mu=mu,
  phi=phi,
  gamma=gamma,
  sigma_eta=sigma

)

###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 1

#liczba czasteczek
FSV_Np <-          c(100, 1e3, 1e3)
FSV_Nmif <-        c(10,  50,  150)
FSV_Nreps_eval <-  c(4,   10,  20)
FSV_Nreps_local <- c(4,   10,  20)
FSV_Nreps_global <-c(4,   10,  20)

FSVlist<-list(FSV_Np ,FSV_Nmif,FSV_Nreps_eval,
              FSV_Nreps_local,FSV_Nreps_global )

#parametry do metody mif2
FSV_rw.sd_rp <- 0.02
FSV_rw.sd_ivp <- 0.1
FSV_cooling.fraction.50 <- 0.5

##########################################################################
##########################################################################
########################################################symulowanie danych


log.volatility.list  <-simulate(FSV.model,nsim=1,seed=123,params=params_test,
                                states = TRUE, obs = TRUE)
zwroty.sim= log.volatility.list$obs[1,1,]
log.volatility.sim= log.volatility.list$states[1,1,] 
G.sim=log.volatility.list$states[2,1,] 

par(mfrow=c(4,1))
plot(log.volatility.list$obs[1,1,],type='l',ylab="y",main="zwroty") 
plot(log.volatility.list$states[1,1,],type='l',ylab="h",main="log-zmiennosc")
plot(1/sqrt(log.volatility.list$states[2,1,]),type='l',ylab="1/sqrt(G)",main="1/sqrt(G)")
plot(exp(log.volatility.list$states[1,1,]/2)/sqrt(log.volatility.list$states[2,1,]),type='l',ylab="exp(h/2)/sqrt(G)",main="zmiennosc")
par(mfrow=c(1,1))


wykres1= rbind(log.volatility.list$obs[1,1,],
               log.volatility.list$states[1,1,],
               1/sqrt(log.volatility.list$states[2,1,]),
               exp(log.volatility.list$states[1,1,]/2)/sqrt(log.volatility.list$states[2,1,]))
row.names(wykres1) = c("y","H","G1","G2")
wykres1=as.data.frame(t(wykres1))

g1<-ggplot(data = as.data.frame(wykres1), aes(x = 1:n, y = y))  + 
  geom_line(color = "royalblue", size = .5)+ggtitle('zwroty')+labs(x="time")+theme_bw()
g2<-ggplot(data = as.data.frame(wykres1), aes(x = 1:n, y = H))  + 
  geom_line(color = "royalblue", size = .5)+ggtitle('log-zmienność')+labs(x="time",y="h")+theme_bw()
g3<-ggplot(data = as.data.frame(wykres1), aes(x = 1:n, y = G1))  + 
  geom_line(color = "royalblue", size = .5)+ggtitle('1/sqrt(G)')+labs(x="time",y="1/sqrt(G)")+theme_bw()
g4<-ggplot(data = as.data.frame(wykres1), aes(x = 1:n, y = G2))  + 
  geom_line(color = "royalblue", size = .5)+ggtitle('zmienność')+labs(x="time",y="exp(h/2)/sqrt(G)")+theme_bw()

grid.arrange(g1, g2,g3,g4, nrow=4)



##########################################################################
##########################################################################
########################################################Filtrowanie danych

FSV.filt<-pomp(data=data.frame(y=zwroty.sim,
                               time=1:n),
               statenames=FSV_statenames,
               paramnames=FSV_paramnames,
               covarnames=FSV_covarnames,
               times="time",
               t0=0,
               covar=data.frame(covaryt=c(0,zwroty.sim),
                                time=0:n),
               tcovar="time",
               rmeasure=Csnippet(FSV_rmeasure),
               dmeasure=Csnippet(FSV_dmeasure),
               rprocess=discrete.time.sim(step.fun=Csnippet(FSV_rproc.filt),delta.t=1),
               initializer=Csnippet(FSV_initializer),
               toEstimationScale=Csnippet(FSV_toEstimationScale), 
               fromEstimationScale=Csnippet(FSV_fromEstimationScale)
)

#wykres filtracji dla prawdziwych parametr?w
pf1 <- pfilter(FSV.filt,params=params_test,
               Np=1000,filter.traj=T)
plot(pf1)


##########################################################################
##########################################################################
######################################################Estymacja parametrow

FSV_box <- rbind(
  sigma_eta=c(0.1,1),
  phi    = c(0.9,0.99),
  mu = c(-1,01),
  gamma=c(4,6)
)

detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


start_time <- Sys.time()
t.if.FSV <- system.time({
  if.FSV <- foreach(i=1:FSVlist[[5]][run_level] ,
                    .packages='pomp', .combine=c,.export = "FSVlist", 
                    .options.multicore=list(set.seed=TRUE)) %dopar% try(
                      pomp::mif2(FSV.filt,start=apply(FSV_box,1,function(x) runif(1,x[1],x[2])),Np=FSVlist[[1]][run_level] , Nmif=FSVlist[[2]][run_level] ,cooling.type="geometric",
                                 cooling.fraction.50=FSV_cooling.fraction.50,
                                 transform=TRUE,
                                 rw.sd = rw.sd(
                                   mu      = FSV_rw.sd_rp,
                                   phi       = FSV_rw.sd_rp,
                                   sigma_eta = FSV_rw.sd_rp,
                                   gamma     = FSV_rw.sd_rp
                                 )
                      )
                      
                    )
  
  L.if.FSV <- foreach(i=1:FSVlist[[5]][run_level] ,.packages='pomp', .export = "FSVlist", 
                      .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                        logmeanexp(
                          replicate(FSVlist[[3]][run_level] ,
                                    logLik(pfilter(FSV.filt,params=coef(if.FSV[[i]]),Np=FSVlist[[1]][run_level]  ))
                          ),
                          se=TRUE)
                      )
  
  H.if.FSV<- foreach(i=1:FSVlist[[5]][run_level] ,.packages='pomp', .export = "FSVlist", 
                     .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                       exp(pfilter(FSV.filt,params=coef(if.FSV[[i]]),Np=FSVlist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                     )
})


stopCluster(cl)
end_time <- Sys.time()
difftime(end_time,start_time, units = "mins")
plot(if.FSV)


if.FSV.box  <- data.frame(logLik=L.if.FSV[,1],logLik_se=L.if.FSV[,2],t(sapply(if.FSV,coef)))
if.FSV.box [which.max(if.FSV.box$logLik),]
pairs(~logLik+mu+phi+sigma_eta+gamma,data=if.FSV.box )


params_nowe2<- c(
  mu        = as.numeric(if.FSV.box[which.max(if.FSV.box $logLik),'mu']),    
  phi         = as.numeric(if.FSV.box[which.max(if.FSV.box $logLik),'phi']),    
  sigma_eta   = as.numeric(if.FSV.box[which.max(if.FSV.box$logLik),'sigma_eta']),
  gamma       = as.numeric(if.FSV.box[which.max(if.FSV.box$logLik),'gamma'])
)


############################################################################
############################################################################
############################################################################
######################################profile funkcji wiarygodnosci 2 wersja


params_nowe2

#################################################################parametr mu
xx1<-seq(from=-1.5,to=0.5,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.FSV.log<- foreach(i=1:length(xx1) ,.packages='pomp', .export = "FSVlist",.combine=rbind,
                    .options.multicore=list(set.seed=TRUE)) %dopar% {
                      set.seed(87932)
                      logLik(pfilter(FSV.filt,params=c(mu=xx1[i], 
                                                       phi=as.numeric(params_nowe2['phi']),
                                                       gamma=as.numeric(params_nowe2['gamma']),
                                                       sigma_eta=as.numeric(params_nowe2['sigma_eta'])),
                                     Np=FSVlist[[1]][run_level] ))
                    }

stopCluster(cl)
beep(2)

plot(xx1, L.FSV.log, type='l',xlab=expression(mu),ylab="logLik")
#abline(v=params_nowe2['mu_h'],lty=2)
#points(xx1, L.FSV.log)
points(if.FSV.box [,'mu'], if.FSV.box[,'logLik'] ,col='red',pch=16)
p=loess(L.FSV.log~xx1,span=0.5)
lines(xx1,p$fitted,col='blue',lwd=2)


wykres.mu=as.data.frame(t(rbind(xx1,as.vector(L.FSV.log))))
names(wykres.mu)<-c("mu","loglik")
g1<-ggplot(data = wykres.mu, aes(x = mu, y = loglik))  + 
  geom_point(color = "black", size = 1)+
  ggtitle(expression(mu))+labs(x=expression(mu))+geom_smooth()+theme_bw()


#parametr phi
xx2<-seq(from=0.95,to=0.99,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.FSV.log2<- foreach(i=1:length(xx2) ,.packages='pomp', .export = "FSVlist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                       logLik(pfilter(FSV.filt,params=c(mu=as.numeric(params_nowe2['mu']),
                                                        phi=xx2[i],
                                                        gamma=as.numeric(params_nowe2['gamma']),
                                                        sigma_eta=as.numeric(params_nowe2['sigma_eta'])),
                                      Np=FSVlist[[1]][run_level] ))
                     }
stopCluster(cl)
beep(2)

plot(xx2, L.FSV.log2, type='l',xlab=expression(phi),ylab="logLik")
#points(xx2, L.FSV.log2)
points(if.FSV.box[,'phi'], if.FSV.box[,'logLik'] ,col='red',pch=16)
p2=loess(L.FSV.log2~xx2, span=0.5)
lines(xx2,p2$fitted,col='blue',lwd=2)

wykres.phi=as.data.frame(t(rbind(xx2,as.vector(L.FSV.log2))))
names(wykres.phi)<-c("phi","loglik")
g2<-ggplot(data = wykres.phi, aes(x = phi, y = loglik))  + 
  geom_point(color = "black", size = 1)+
  ggtitle(expression(phi))+labs(x=expression(phi))+geom_smooth()+theme_bw()



#parametr sigma_eta
xx3<-seq(from=params_nowe2['sigma_eta']-.1,to=params_nowe2['sigma_eta']+.1,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.FSV.log3<- foreach(i=1:length(xx3) ,.packages='pomp', .export = "FSVlist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                       logLik(pfilter(FSV.filt,params=c(mu=as.numeric(params_nowe2['mu']), 
                                                        gamma=as.numeric(params_nowe2['gamma']),
                                                        phi=as.numeric(params_nowe2['phi']),
                                                        sigma_eta=xx3[i]),
                                      Np=FSVlist[[1]][run_level] ))
                     }

stopCluster(cl)
beep(2)

plot(xx3, L.FSV.log3, type='l',xlab=expression(sigma[eta]),ylab="logLik")
#points(xx3, L.FSV.log3)
points(if.FSV.box[,'sigma_eta'], if.FSV.box[,'logLik'] ,col='red',pch=16)
p3=loess(L.FSV.log3~xx3,span=0.5)
lines(xx3,p3$fitted,col='blue',lwd=2)


wykres.sigma_eta=as.data.frame(t(rbind(xx3,as.vector(L.FSV.log3))))
names(wykres.sigma_eta)<-c("sigma_eta","loglik")
g3<-ggplot(data = wykres.sigma_eta, aes(x = sigma_eta, y = loglik))  + 
  geom_point(color = "black", size = 1)+
  ggtitle(expression(sigma[eta]))+labs(x=expression(sigma[eta]))+geom_smooth()+theme_bw()


#parametr gamma
xx4<-seq(from=params_nowe2['gamma']-2,to=params_nowe2['gamma']+2,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.FSV.log4<- foreach(i=1:length(xx4) ,.packages='pomp', .export = "FSVlist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                       logLik(pfilter(FSV.filt,params=c(mu=as.numeric(params_nowe2['mu']), 
                                                        phi=as.numeric(params_nowe2['phi']),
                                                        sigma_eta=as.numeric(params_nowe2['sigma_eta']),
                                                          gamma=xx4[i]),
                                      Np=FSVlist[[1]][run_level] ))
                     }

stopCluster(cl)
beep(2)

plot(xx4, L.FSV.log4, type='l',xlab=expression(gamma),ylab="logLik")
#points(xx3, L.FSV.log3)
points(if.FSV.box[,'gamma'], if.FSV.box[,'logLik'] ,col='red',pch=16)
p4=loess(L.FSV.log4~xx4,span=0.5)
lines(xx4,p4$fitted,col='blue',lwd=2)


wykres.gamma=as.data.frame(t(rbind(xx4,as.vector(L.FSV.log4))))
names(wykres.gamma)<-c("gamma","loglik")
g4<-ggplot(data = wykres.gamma, aes(x = gamma, y = loglik))  + 
  geom_point(color = "black", size = 1)+
  ggtitle(expression(gamma))+labs(x=expression(gamma))+geom_smooth()+theme_bw()

#################################################################podsumowanie


grid.arrange(g1, g2,g3,g4, ncol=2)

##################################################################################
##################################################################################
########################################estymacja Qauasi-największej wiarygodności



OUss <- function(mu,phi,sigma,gamma){
  Tt <- matrix(phi,ncol=1)
  Zt <- matrix(1,ncol=1)
  ct <- matrix(-1.27-digamma(gamma/2)+log(gamma/2),ncol=1)
  dt <- matrix(mu*(1-phi), nrow = 1)
  GGt<- matrix(pi^2/2+trigamma(gamma/2),nrow = 1,ncol=1)
  HHt<- matrix(sigma^2,nrow=1,ncol=1)
  a0 <-  as.vector(matrix(mu+1.27))
  P0 <- matrix(1,nrow=1,ncol=1)
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
              HHt = HHt))
}

KF.log <- function(theta) {
  sp <- OUss(theta[1], theta[2], theta[3],theta[4])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt =matrix(log(zwroty.sim^2), nrow=1,ncol=n))
  return(-ans$logLik)
}



KF <- function(theta) {
  sp <- OUss(theta[1], theta[2], theta[3])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt =matrix(log(zwroty.sim^2), nrow=1,ncol=n))
  return(ans$att[1,])
}


KF.log(c(mu,phi,sigma,gamma))
KF.opt<-optim(c(mu,phi,sigma,gamma),KF.log, 
              method="L-BFGS-B",hessian = T,lower=c(-Inf,-1,0,0),upper=c(Inf,1,Inf,Inf))
KF.opt$par
KF.opt$value

xx1<-seq(from=-2,to=0,length.out = 100)
p1<-sapply(xx1, function(z)  - KF.log(c(z, KF.opt$par[2],KF.opt$par[3] ,KF.opt$par[4])))
plot(xx1, p1, type='l', xlab=expression(mu))
abline(v=mu)
abline(v=KF.opt$par[1],col='red',lty=2)

xx2<-seq(from=0.95,to=0.999,length.out = 100)
p2<-sapply(xx2, function(z)  - KF.log(c(KF.opt$par[1], z,KF.opt$par[3],KF.opt$par[4])))
plot(xx2, p2, type='l', xlab=expression(phi))
abline(v=phi)
abline(v=KF.opt$par[2],col='red',lty=2)

xx3<-seq(from=0.1,to=0.4,length.out = 100)
p3<-sapply(xx3, function(z)  - KF.log(c(KF.opt$par[1], KF.opt$par[2],z ,KF.opt$par[4])))
plot(xx3, p3, type='l', xlab=expression(sigma[eta]))
abline(v=sigma)
abline(v=KF.opt$par[3],col='red',lty=2)


xx4<-seq(from=1,to=10,length.out = 100)
p4<-sapply(xx4, function(z)  - KF.log(c(KF.opt$par[1], KF.opt$par[2],KF.opt$par[3] ,z)))
plot(xx4, p4, type='l', xlab=expression(gamma))
abline(v=gamma)
abline(v=KF.opt$par[4],col='red',lty=2)


############################################################################
############################################################################
############################################################################
################################################################Liu and West 

##################################nie zrobione
stochVCode <- nimbleCode({
  
  
  x[1] ~ dnorm(mu*(1-phi) + phi * x0, var = 1/sigmaSquaredInv)
  G[1] ~rchisq(gamma)
  y[1] ~ dnorm(0, var =  exp(x[1]))
  
  for(t in 2:T){
    x[t] ~ dnorm(mu*(1-phi) + phi * x[t-1],  var = 1/sigmaSquaredInv)
    y[t] ~ dnorm(0, var = exp(x[t])/G[t])
  }
  
  x0 ~ dnorm(0, var = 1/sigmaSquaredInv)
  phi <- 2 * phiStar - 1
  phiStar ~ dbeta(18, 1)
  sigmaSquaredInv ~ dgamma(5, 20)
  mu ~ dnorm(0, var = 100)      
})

stochVolModel <- nimbleModel(code = stochVCode, name ='stochVol',
                             constants = list(T = length(zwroty)),  data = list(y = as.vector(zwroty)),
                             inits = list(mu = mu,  phiStar = (phi+1)/2,
                                          sigmaSquaredInv = 1/sigma^2))


sigmaSquaredSamples <- 1/sqrt((as.matrix(CstochVolLiuWestFilter$mvEWSamples,'sigmaSquaredInv')))/10
plot(density(sigmaSquaredSamples))
muSamples <- as.matrix(CstochVolLiuWestFilter$mvEWSamples,'mu')
plot(density(muSamples))
phiSamples <- (as.matrix(CstochVolLiuWestFilter$mvEWSamples,'phiStar')+1)/2
plot(density(phiSamples))
round(mean(sigmaSquaredSamples),4)
round(mean(muSamples),4)
round(mean(phiSamples),4)

params_nowe3<- c(
  mu_h        = mean(muSamples) ,  
  phi         = mean(phiSamples),   
  sigma_eta   = mean(sigmaSquaredSamples)
)

logLik(pfilter(bsv.filt,params=params_nowe3 ,Np= 1000))

############################################################################
############################################################################
############################################################################
########################################################################PMCM

hyperparams <- list(min = c(-1,0.9,4,0), max = c(0,1,6,1) )
FSV.dprior <- function (params, ..., log) {
  f <- sum(dunif(params, min = hyperparams$min, max = hyperparams$max,
                 log = TRUE))
  if (log) f else exp(f)
}


pmcmc1 <-   pmcmc(pomp(FSV.filt, dprior = FSV.dprior), start = params_test,
                  Nmcmc = 1000, Np = 100, max.fail = Inf,
                  proposal = mvn.diag.rw(c(mu = 0.01, phi = 0.01, gamma=0.01, sigma_eta = 0.01)))

continue( pmcmc1 ,Nmcmc=2000,proposal=mvn.rw(covmat( pmcmc1 ))) -> pmcmc1 
plot(  pmcmc1 )

plot(  pmcmc1 )
coef(pmcmc1 )
logLik(pmcmc1 )
beep(2)

