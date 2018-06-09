#model  Bandorffa-Nielsena i Sheparda 

#rm(list = ls())
#gc()
setwd("C:/Users/user/Documents/github/Zastosowanie-iterowanej-filtracji-do-modeli-stochastycznej-zmienno-ci/Dane")

##############################################################################
##############################################################################
#################################################################### parametry


############################################
#poczatek okresu badania 
data.poczatkowa='2012-09-04'
#koniec okresu badania
data.koncowa='2016-09-08'
############################################



delta=1

mu=0
beta=0
ksi=.5
omega=.2
alpha=ksi/omega^2
ni=ksi^2/omega^2
lambda=0.01
lambda.f=lambda*ni*delta

prawdziwe.parametry<-c(mu,lambda,ksi,omega)
names(prawdziwe.parametry)<-c("mu","lambda","ksi","omega")


##############################################################################
##############################################################################
##################################################################### pakiety

library(forecast)
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
library(PerformanceAnalytics)

dane=read.csv.zoo('wig_zwroty.csv',sep=',')
 dane=as.xts(dane)
zwroty=dane[paste(data.poczatkowa,"/",data.koncowa,sep=""),]
n=length(zwroty)
x=as.vector(zwroty)
x2=x^2

plot(x,type='l')
plot(x2,type='l')
###########################################################################################
###########################################################################################
############################################################## implementacja filtru Kalmana

OUss <- function(mu, lambda, ksi, omega) {
  Tt <- matrix(c(0, 0, (1-exp(-lambda*delta)), exp(-lambda*delta)), ncol = 2)
  Zt <- matrix(c(0,1/lambda,0, 0), ncol = 2)
  ct <- matrix(c(mu*delta,mu^2*delta^2),ncol=1)
  dt <- matrix(c(ksi*lambda*delta-ksi+ksi*exp(-lambda*delta),ksi-ksi*exp(-lambda*delta)), nrow = 2)
  
  GGt <- matrix(c(ksi*delta, 2*mu*delta^2*ksi,2*mu*delta^2*ksi,
                  4*mu^2*delta^3*ksi+2*ksi^2*delta^2+4*omega^2/(lambda^2)*(exp(-lambda*delta)-1+lambda*delta)), ncol = 2)
  HHt <-2*omega^2* matrix( c( (-3/2-1/2*exp(-2*lambda*delta)+2*exp(-lambda*delta)+lambda*delta),
                              (1-exp(-lambda*delta)-1/2*(1-exp(-2*lambda*delta))),
                              (1-exp(-lambda*delta)-1/2*(1-exp(-2*lambda*delta))),
                              1/2*(1-exp(-2*lambda))),ncol=2)
  
  a0 <- c(ksi,ksi)
  P0 <- 0.1*diag(2)
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
              HHt = HHt))
}



objective <- function(theta) {
  sp <- OUss(theta[1], theta[2], theta[3], theta[4])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = rbind(x,x2))
  return(-ans$logLik)
}

Kalman.filter <- function(theta) {
  sp <- OUss(theta[1], theta[2], theta[3], theta[4])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = rbind(x,x2))
  return(ans$att[1,]/ theta[2])
}



##############################################################################
##############################################################################
##############################################################    Estymacja 



objective(c(mu, lambda, ksi, omega))

start.time <- Sys.time()
model1=tryCatch({optim(c(mu,lambda,alpha,ni),objective, 
                       method="L-BFGS-B", lower=c(-1,0.0001, 0.001,0.001),upper=c(10,1,Inf,Inf))})
end.time <- Sys.time()

difftime(end.time, start.time, units= "secs")
names(model1$par)<-c("mu","lambda","ksi","omega")
round(model1$par,5)

plot(Kalman.filter(model1$par),type='l')

xx1<-seq(from=-1,to=1,length.out = 101)
p1<-sapply(xx1, function(z)  - objective(c(z, model1$par[2],model1$par[3] ,model1$par[4])))
plot(xx1, p1, type='l', xlab=expression(mu))
abline(v=model1$par[1],col='red',lty=2)


wykres.mu=as.data.frame(t(rbind(xx1, p1)))
names(wykres.mu)<-c("mu","loglik")
g1<-ggplot(data = wykres.mu, aes(x = mu, y = loglik))  + 
  ggtitle(expression(mu))+labs(x=expression(mu))+ geom_line(color = "blue", size = 1) +
  geom_segment(mapping=aes(x = model1$par[1], y =loglik[which.min(p1)],xend=model1$par[1],yend=loglik[which.max(p1)] ),color="black", size = 1)


xx2<-seq(from=0.1,to=1.1,length.out = 100)
p2<-sapply(xx2, function(z)  - objective(c(model1$par[1],z,model1$par[3] ,model1$par[4])))
plot(xx2, p2, type='l', xlab=expression(lambda))
abline(v=model1$par[2],col='red',lty=2)


wykres.lambda=as.data.frame(t(rbind(xx2, p2)))
names(wykres.lambda)<-c("lambda","loglik")
g2<-ggplot(data = wykres.lambda, aes(x = lambda, y = loglik))  + 
  ggtitle(expression(lambda))+labs(x=expression(lambda))+ geom_line(color = "blue", size = 1) +
  geom_segment(mapping=aes(x = model1$par[2], y =loglik[which.min(p2)],xend=model1$par[2],yend=loglik[which.max(p2)] ),color="black", size = 1)
 


xx3<-seq(from=0.5,to=1.5,length.out = 100)
p3<-sapply(xx3, function(z)  - objective(c(model1$par[1],model1$par[2], z,model1$par[4])))
plot(xx3, p3, type='l', xlab=expression(xi))
abline(v=model1$par[3],col='red',lty=2)


wykres.xi=as.data.frame(t(rbind(xx3, p3)))
names(wykres.xi)<-c("xi","loglik")
g3<-ggplot(data = wykres.xi, aes(x = xi, y = loglik))  + 
  ggtitle(expression(xi))+labs(x=expression(xi))+ geom_line(color = "blue", size = 1) +
  geom_segment(mapping=aes(x = model1$par[3], y =loglik[which.min(p3)],xend=model1$par[3],yend=loglik[which.max(p3)] ),color="black", size = 1)
 


xx4<-seq(from=0.5,to=1.5,length.out = 100)
p4<-sapply(xx4, function(z)  - objective(c(model1$par[1],model1$par[2], model1$par[3],z)))
plot(xx4, p4, type='l', xlab=expression(omega))
abline(v=model1$par[4],col='red',lty=2)

wykres.omega=as.data.frame(t(rbind(xx4, p4)))
names(wykres.omega)<-c("omega","loglik")
g4<-ggplot(data = wykres.omega, aes(x = omega, y = loglik))  + 
  ggtitle(expression(omega))+labs(x=expression(omega))+ geom_line(color = "blue", size = 1) +
  geom_segment(mapping=aes(x = model1$par[4], y =loglik[which.min(p4)],xend=model1$par[4],yend=loglik[which.max(p4)] ),color="black", size = 1)


grid.arrange(g1, g2,g3,g4, ncol=2)



###########################################################################################
###########################################################################################
########################################################################Iterowana filtracja




rprocess.fun<- function(x, t, params, delta.t, ...){
  
  eta=eta.fun(params['ksi']/params['omega']^2,
              params['ksi']^2/params['omega']^2,
              params['lambda'])
  x2=exp(-params['lambda'])*x["x2"]+eta[1]
  x1=(1-exp(-params['lambda']))*x["x2"]+ eta[2]-eta[1]
  return(c(x1=as.numeric(x1),x2=as.numeric(x2)))
}





rprocess.fun2<- function(x, t, params, delta.t, ...){
  alpha =  params['ksi']/params['omega']^2
  ni    =  params['ksi']^2/params['omega']^2
  lambda =  params['lambda']
  delta=1
  
  ci<-NULL
  sum=0
  i=0
  while(sum<1){
    cip<-rexp(1, rate = lambda*ni*delta)
    ci<-c(ci,cip)
    sum=sum+cip
    i=i+1
  }
  if(length(ci)>1){
    ci<-ci[1:(length(ci)-1)]
    ci<-cumsum(ci)
  }else ci=0
  
  if(sum(ci^2)>0){
    ri<-runif(length(ci))
    temp1=exp(-lambda*delta)*(alpha^(-1))*sum(log(1/ci[1:length(ci)])*exp(lambda*delta*ri[1:length(ci)]))
    temp2=(alpha^(-1))*sum(log(1/ci[1:length(ci)]))
  } else {
    temp1=0
    temp2=0
  }
  eta=c(temp1,temp2)
  
  
  x2=exp(-lambda)*x["x2"]+eta[1]
  x1=(1-exp(-lambda ))*x["x2"]+ eta[2]-eta[1]
  
  
  return(c(x1=as.numeric(x1),x2=as.numeric(x2)))
}


rmeasure.fun <-function (x, t, params,...){
  return(c(y=rnorm(1,params['mu'],sqrt(x['x1']/params['lambda']))))
}


dmeasure.fun <-function (y,x, t, params,log,...){
  return(dnorm(y['y'],params['mu'],sqrt(x['x1']/params['lambda']),log=log))
}


initializer.fun <-function(params, t0,...){
  alpha =  params['ksi']/params['omega']^2
  ni    =  params['ksi']^2/params['omega']^2
  lambda =  params['lambda']
  delta=1
  ci<-NULL
  sum=0
  i=0
  while(sum<1){
    cip<-rexp(1, rate = lambda*ni*delta)
    ci<-c(ci,cip)
    sum=sum+cip
    i=i+1
  }
  if(length(ci)>1){
    ci<-ci[1:(length(ci)-1)]
    ci<-cumsum(ci)
  }else ci=0
  
  if(sum(ci^2)>0){
    ri<-runif(length(ci))
    temp1=exp(-lambda*delta)*(alpha^(-1))*sum(log(1/ci[1:length(ci)])*exp(lambda*delta*ri[1:length(ci)]))
    temp2=(alpha^(-1))*sum(log(1/ci[1:length(ci)]))
  } else {
    temp1=0
    temp2=0
  }
  eta=c(temp1,temp2)
  
  
  
  return(c(x1=as.numeric(params['ksi']+eta[2]-eta[1]),x2=as.numeric(params['ksi'] +eta[1])))
}


logit<-function(p) log(p / (1 - p))
expit<-function(x) exp(x)/(1 + exp(x))

####przeskalowanie parametr?w 
bns_toEstimationScale <- function(params,...){
  mu=as.numeric(params['mu'])
  lambda=as.numeric(log(params['lambda']/(1-params['lambda'])))
  ksi=as.numeric(log(params['ksi']))
  omega=as.numeric(log(params['omega']))
  return(c(mu=mu,lambda=lambda,ksi=ksi,omega=omega))
}



bns_fromEstimationScale <- function(params,...){
  mu=as.numeric(params['mu'])
  lambda=as.numeric(exp(params['lambda'])/(1+exp(params['lambda'])))
  ksi=as.numeric(exp(params['ksi']))
  omega=as.numeric(exp(params['omega']))
  return(c(mu=mu,lambda=lambda,ksi=ksi,omega=omega))
}

####wypelnianie modelu danymi
bns.model <- pomp(data=data.frame(y=x[1:n], time=1:n),
                  times="time",
                  t0=0,
                  rmeasure= rmeasure.fun,
                  dmeasure=dmeasure.fun,
                  rprocess=discrete.time.sim(step.fun =rprocess.fun2, delta.t=1),
                  initializer=initializer.fun,
                  toEstimationScale=bns_toEstimationScale, 
                  fromEstimationScale=bns_fromEstimationScale 
                  
)


plot(bns.model)

params_test=model1$par
round(params_test,5)


###########################################################################################
###########################################################################################
##################################################################################Filtracja

pf1 <- pfilter(bns.model,params=params_test,
               Np=1000,filter.mean=T)
plot(pf1)

plot(Kalman.filter(params_test)[2:n],type='l',col='red')
lines(filter.mean(pf1)[1,2:n]/lambda,col='royalblue')

##########################################################################################
###########################################################################################
##################################################################################Estymacja


###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 3

#liczba czasteczek
bns_Np <-          c(100, 1e3, 1e3)
bns_Nmif <-        c(10,  50,  100)
bns_Nreps_eval <-  c(4,   10,  20)
bns_Nreps_local <- c(4,   10,  20)
bns_Nreps_global <-c(4,   10,  20)

bnslist<-list(bns_Np ,bns_Nmif,bns_Nreps_eval,
              bns_Nreps_local,bns_Nreps_global )

#parametry do metody mif2
bns_rw.sd_rp <- 0.01
bns_rw.sd_ivp <- 0.1
bns_cooling.fraction.50 <- 0.75


bns_box <- rbind(
  mu = c(-.5,.5),
  lambda  = c(0.01,1),
  ksi=c(0.1,1.5),
  omega=c(0.1,1.5)
)


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


start_time <- Sys.time()
t.if.bns <- system.time({
  if.bns <- foreach(i=1:bnslist[[5]][run_level] ,
                    .packages='pomp', .combine=c,.export = "bnslist", 
                    .options.multicore=list(set.seed=TRUE)) %dopar% try(
                      pomp::mif2(bns.model,start=apply(bns_box,1,function(x) runif(1,x[1],x[2])),Np=bnslist[[1]][run_level] , Nmif=bnslist[[2]][run_level] ,cooling.type="geometric",
                                 cooling.fraction.50=bns_cooling.fraction.50,
                                 transform=TRUE,
                                 rw.sd = rw.sd(
                                   mu      = bns_rw.sd_rp,
                                   lambda   = bns_rw.sd_rp,
                                   ksi = bns_rw.sd_rp,
                                   omega    = bns_rw.sd_rp
                                 )
                      )
                    )
  L.if.bns <- foreach(i=1:bnslist[[5]][run_level] ,.packages='pomp',.export ="bnslist", 
                      .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                        logmeanexp(
                          replicate(bnslist[[3]][run_level] ,
                                    logLik(pfilter(bns.model,params=coef(if.bns[[i]]),Np=bnslist[[1]][run_level]  ))
                          ),
                          se=TRUE)
                      )
  
  H.if.bns<- foreach(i=1:bnslist[[5]][run_level] ,.packages='pomp', .export ="bnslist", 
                     .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                       exp(pfilter(bns.model,params=coef(if.bns[[i]]),Np=bnslist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                     )
})


stopCluster(cl)
end_time <- Sys.time()
difftime(end_time,start_time, units = "hours")
plot(if.bns)

if.bns.box  <- data.frame(logLik=L.if.bns[,1],logLik_se=L.if.bns[,2],t(sapply(if.bns,coef)))
if.bns.box [which.max(if.bns.box$logLik),]
round(if.bns.box [which.max(if.bns.box$logLik),],5)



params_new=c(
  mu=if.bns.box [which.max(if.bns.box$logLik),'mu'],
  lambda=if.bns.box [which.max(if.bns.box$logLik),'lambda'],
  ksi=if.bns.box [which.max(if.bns.box$logLik),'lambda'],
  omega=if.bns.box [which.max(if.bns.box$logLik),'omega']
)
round(params_new,5)

wyniki.mif.loglik<-matrix(NaN,ncol=bnslist[[5]][run_level]+1,nrow= bnslist[[2]][run_level])
nazwy=1:bnslist[[5]][run_level]
for(i in 1:bnslist[[5]][run_level]) wyniki.mif.loglik[,i]=conv.rec(if.bns)[[i]][1:bnslist[[2]][run_level] ,1] 
for(i in 1:bnslist[[5]][run_level]) nazwy[i]=paste(i)
wyniki.mif.loglik[,bnslist[[5]][run_level]+1]=1:bnslist[[2]][run_level]
wyniki.mif.loglik=as.data.frame(wyniki.mif.loglik)
names(wyniki.mif.loglik)<-c(nazwy,"nr")
wyniki.mif.loglik.long<-melt(wyniki.mif.loglik,id="nr")
names(wyniki.mif.loglik.long)<-c("Nr_iteracji","IF","value")

g1<-ggplot(data = wyniki.mif.loglik.long, aes(x = Nr_iteracji,y=value,colour=IF,linetype=IF) ) +ggtitle('logLik')+
  geom_line(size = 1)+labs(x="Nr iteracji",y="Loglik")+theme_bw()+
  scale_x_discrete(limits=c(0,25,50,75,100)) +theme(legend.position="none")


wyniki.mif.mu<-matrix(NaN,ncol=bnslist[[5]][run_level]+1,nrow= bnslist[[2]][run_level]+1)
for(i in 1:bnslist[[5]][run_level]) wyniki.mif.mu[,i]=conv.rec(if.bns)[[i]][,3] 
wyniki.mif.mu[,bnslist[[5]][run_level]+1]=0:bnslist[[2]][run_level]
wyniki.mif.mu=as.data.frame(wyniki.mif.mu)
names(wyniki.mif.mu)<-c(nazwy,"nr")
wyniki.mif.mu.long<-melt(wyniki.mif.mu,id="nr")
names(wyniki.mif.mu.long)<-c("Nr_iteracji","IF","value")

g2<-ggplot(data = wyniki.mif.mu.long, aes(x = Nr_iteracji,y=value,colour=IF,linetype=IF) ) +ggtitle(expression(mu))+
  geom_line(size = 1)+labs(x="Nr iteracji",y=expression(mu))+theme_bw()+
  scale_x_discrete(limits=c(0,25,50,75,100))+theme(legend.position="none")



wyniki.mif.lambda<-matrix(NaN,ncol=bnslist[[5]][run_level]+1,nrow= bnslist[[2]][run_level]+1)
for(i in 1:bnslist[[5]][run_level]) wyniki.mif.lambda[,i]=conv.rec(if.bns)[[i]][,4] 
wyniki.mif.lambda[,bnslist[[5]][run_level]+1]=0:bnslist[[2]][run_level]
wyniki.mif.lambda=as.data.frame(wyniki.mif.lambda)
names(wyniki.mif.lambda)<-c(nazwy,"nr")
wyniki.mif.lambda.long<-melt(wyniki.mif.lambda,id="nr")
names(wyniki.mif.lambda.long)<-c("Nr_iteracji","IF","value")

g3<-ggplot(data = wyniki.mif.lambda.long, aes(x = Nr_iteracji,y=value,colour=IF,linetype=IF) ) +ggtitle(expression(lambda))+
  geom_line(size = 1)+labs(x="Nr iteracji",y=expression(lambda))+theme_bw()+
  scale_x_discrete(limits=c(0,25,50,75,100))+theme(legend.position="none")


wyniki.mif.sigma<-matrix(NaN,ncol=bnslist[[5]][run_level]+1,nrow= bnslist[[2]][run_level]+1)
for(i in 1:bnslist[[5]][run_level]) wyniki.mif.sigma[,i]=conv.rec(if.bns)[[i]][,5] 
wyniki.mif.sigma[,bnslist[[5]][run_level]+1]=0:bnslist[[2]][run_level]
wyniki.mif.sigma=as.data.frame(wyniki.mif.sigma)
names(wyniki.mif.sigma)<-c(nazwy,"nr")
wyniki.mif.sigma.long<-melt(wyniki.mif.sigma,id="nr")
names(wyniki.mif.sigma.long)<-c("Nr_iteracji","IF","value")

g5<-ggplot(data = wyniki.mif.sigma.long, aes(x = Nr_iteracji,y=value,colour=IF,linetype=IF) ) +ggtitle(expression(xi))+
  geom_line(size = 1)+labs(x="Nr iteracji",y=expression(xi))+theme_bw()+
  scale_x_discrete(limits=c(0,25,50,75,100))+theme(legend.position="none")



wyniki.mif.sigma<-matrix(NaN,ncol=bnslist[[5]][run_level]+1,nrow= bnslist[[2]][run_level]+1)
for(i in 1:bnslist[[5]][run_level]) wyniki.mif.sigma[,i]=conv.rec(if.bns)[[i]][,6] 
wyniki.mif.sigma[,bnslist[[5]][run_level]+1]=0:bnslist[[2]][run_level]
wyniki.mif.sigma=as.data.frame(wyniki.mif.sigma)
names(wyniki.mif.sigma)<-c(nazwy,"nr")
wyniki.mif.sigma.long<-melt(wyniki.mif.sigma,id="nr")
names(wyniki.mif.sigma.long)<-c("Nr_iteracji","IF","value")

g4<-ggplot(data = wyniki.mif.sigma.long, aes(x = Nr_iteracji,y=value,colour=IF,linetype=IF) ) +ggtitle(expression(omega))+
  geom_line(size = 1)+labs(x="Nr iteracji",y=expression(omega))+theme_bw()+
  scale_x_discrete(limits=c(0,25,50,75,100))+theme(legend.position="none")


grid.arrange( g2,g3,g1, g5,g4,ncol=3)




###########################################################################################
###########################################################################################
##################################################################################Filtracja
pf1 <- pfilter(bns.model,params=params_new,
               Np=1000,filter.mean=T)
pf2<-pfilter(bns.model,params=params_test,
             Np=1000,filter.mean=T)
plot(pf1)

plot(Kalman.filter(params_test)[2:n],type='l',col='red')
lines(filter.mean(pf1)[1,2:n]/params_new[2],col='royalblue')

plot(filter.mean(pf1)[1,2:n]/params_new[2],type='l',col='royalblue')
lines(Kalman.filter(params_test)[2:n], col='red')

 
particle.bns<-filter.mean(pf1)[1,1:n]/params_new[2]
kalman.bns<-filter.mean(pf2)[1,1:n]/params_test[2]
kalman.kalman.bns<-Kalman.filter(params_test)[2:n]


plot(particle.bns,type='l',col='royalblue')
lines(kalman.bns, col='red')




save(x,x2,params_new,params_test,particle.bns,
     kalman.bns,kalman.kalman.bns,file="wyniki_bns_empiryczne")
load(file="wyniki_bns_empiryczne")
