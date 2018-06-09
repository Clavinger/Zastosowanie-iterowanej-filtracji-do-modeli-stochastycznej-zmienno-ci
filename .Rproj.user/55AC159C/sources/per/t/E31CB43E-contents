#model  Bandorffa-Nielsena i Sheparda 

#rm(list = ls())
#gc()
setwd("C:/Users/user/Documents/github/Zastosowanie-iterowanej-filtracji-do-modeli-stochastycznej-zmienno-ci/Dane")

##############################################################################
##############################################################################
#################################################################### parametry


#dlugosc szeregu czasowego
n=1000
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
############################################################## symulacja

z=1:n
eta2=1:n
sigma=1:n

pos<-function(iten){
  ci<-NULL
  sum=0
  i=0
  while(sum<1){
    cip<-rexp(1, rate = iten)
    ci<-c(ci,cip)
    sum=sum+cip
    i=i+1
  }
  if(length(ci)>1){
    ci<-ci[1:(length(ci)-1)]
    ci<-cumsum(ci)
  }else ci=0
  ci
}


for (i in 2:n){
  ci<-pos(lambda.f)
  ri<-runif(length(ci))
  eta2[i]=exp(-lambda*delta)*(alpha^(-1))*sum(log(1/ci[1:length(ci)])*exp(lambda*delta*ri[1:length(ci)]))
}
eta2[eta2==Inf]=0


z[1]=0
sigma[1]=ni/alpha

for (i in 2:n){
  sigma[i]=sigma[(i-1)]*exp(-lambda*delta)+eta2[i]
}

sigma
plot(sigma,type='l')

#set.seed(123) 
u<-rnorm(n)
y=1:n
y2=1:n
lr=1:n
y[1]=0
for(i in 2:n){
  y[i]=y[i-1]+mu*delta+beta*sigma[i-1]+sqrt(sigma[i-1])*u[i-1]*sqrt(delta)
}
plot(y,type='l')

#zwroty
x<-1:(n-1) 
for (i in 2:n) x[i]=y[i]-y[i-1]


#kwadraty zwrot?w
x2<-x^2

plot(x,type='l')
plot(x2,type='l')

##############################################################################
##############################################################################
############################################################## wykres gęstości


den<-density(sigma)
round(mean(sigma),4)
round(sd(sigma),4)
xx<-den$x
yy1<-den$y
yy2<-dgamma(xx,shape=ni ,rate=alpha)
plot(xx,yy1,type='l')
lines(xx,yy2,col='red')
plot(xx,yy2,type='l')

length(yy2)


wykres1.matrix<-rbind(xx,yy1,
                        yy2 )
row.names(wykres1.matrix)<-c("x","Est. jądrowy","Gamma")
wykres1<-as.data.frame(t(wykres1.matrix))
wykres1.long<-melt(wykres1,id="x")
names(wykres1.long)<-c("x","Gęstość","value")

g2<-ggplot(data = wykres1.long, aes(x = x,y=value,colour=Gęstość,linetype=Gęstość) ) +
  ggtitle('(b)')+
  geom_line(size =1)+labs(y="f(x)")+
  scale_colour_manual(values=c("royalblue","tomato"))+
  scale_linetype_manual(values=c( "solid", "solid"))+theme_bw()



wykres3.matrix<-rbind(1:n,sigma)
wykres3<-as.data.frame(t(wykres3.matrix))
names(wykres3)=c("x","sigma")

g1<-ggplot(data = wykres3, aes(x = x,y=sigma))+
  ggtitle('(a)')+
  labs(x="time",y=expression(sigma(t)))+ 
  geom_line(color = "black", size = 1)+theme_bw()


grid.arrange(g1, g2, ncol=2)




##############################################################################
##############################################################################
############################################################## wykres gęstości


objective(c(mu, lambda, ksi, omega))

start.time <- Sys.time()
model1=tryCatch({optim(c(mu,lambda,alpha,ni),objective, 
                       method="L-BFGS-B", lower=c(-1,0.0001, 0.001,0.001),upper=c(10,1,Inf,Inf))})
end.time <- Sys.time()

difftime(end.time, start.time, units= "secs")
names(model1$par)<-c("mu","lambda","ksi","omega")
prawdziwe.parametry
round(model1$par,5)

xx1<-seq(from=-1,to=1,length.out = 101)
p1<-sapply(xx1, function(z)  - objective(c(z, model1$par[2],model1$par[3] ,model1$par[4])))
plot(xx1, p1, type='l', xlab=expression(mu))
abline(v=mu)
abline(v=model1$par[1],col='red',lty=2)


wykres.mu=as.data.frame(t(rbind(xx1, p1)))
names(wykres.mu)<-c("mu","loglik")
g1<-ggplot(data = wykres.mu, aes(x = mu, y = loglik))  + 
  ggtitle(expression(mu))+labs(x=expression(mu))+ geom_line(color = "blue", size = 1) +
  geom_segment(mapping=aes(x = model1$par[1], y =loglik[which.min(p1)],xend=model1$par[1],yend=loglik[which.max(p1)] ),color="black", size = 1)+
  geom_segment(mapping=aes(x =  0, y =loglik[which.min(p1)],xend=0,yend=loglik[which.max(p1)] ),color="red", size = 1)

xx2<-seq(from=0.0001,to=.05,length.out = 100)
p2<-sapply(xx2, function(z)  - objective(c(model1$par[1],z,model1$par[3] ,model1$par[4])))
plot(xx2, p2, type='l', xlab=expression(lambda))
abline(v=lambda)
abline(v=model1$par[2],col='red',lty=2)


wykres.lambda=as.data.frame(t(rbind(xx2, p2)))
names(wykres.lambda)<-c("lambda","loglik")
g2<-ggplot(data = wykres.lambda, aes(x = lambda, y = loglik))  + 
  ggtitle(expression(lambda))+labs(x=expression(lambda))+ geom_line(color = "blue", size = 1) +
  geom_segment(mapping=aes(x = model1$par[2], y =loglik[which.min(p2)],xend=model1$par[2],yend=loglik[which.max(p2)] ),color="black", size = 1)+
  geom_segment(mapping=aes(x =  0.01, y =loglik[which.min(p1)],xend=0.01,yend=loglik[21] ),color="red", size = 1)



xx3<-seq(from=0.1,to=.7,length.out = 100)
p3<-sapply(xx3, function(z)  - objective(c(model1$par[1],model1$par[2], z,model1$par[4])))
plot(xx3, p3, type='l', xlab=expression(xi))
abline(v=ksi)
abline(v=model1$par[3],col='red',lty=2)


wykres.xi=as.data.frame(t(rbind(xx3, p3)))
names(wykres.xi)<-c("xi","loglik")
g3<-ggplot(data = wykres.xi, aes(x = xi, y = loglik))  + 
  ggtitle(expression(xi))+labs(x=expression(xi))+ geom_line(color = "blue", size = 1) +
  geom_segment(mapping=aes(x = model1$par[3], y =loglik[which.min(p3)],xend=model1$par[3],yend=loglik[which.max(p3)] ),color="black", size = 1)+
  geom_segment(mapping=aes(x =  0.5, y =loglik[which.min(p3)],xend=0.5,yend=loglik[67] ),color="red", size = 1)


xx4<-seq(from=0.01,to=.3,length.out = 100)
p4<-sapply(xx4, function(z)  - objective(c(model1$par[1],model1$par[2], model1$par[3],z)))
plot(xx4, p4, type='l', xlab=expression(omega))
abline(v=omega)
abline(v=model1$par[4],col='red',lty=2)

wykres.omega=as.data.frame(t(rbind(xx4, p4)))
names(wykres.omega)<-c("omega","loglik")
g4<-ggplot(data = wykres.omega, aes(x = omega, y = loglik))  + 
  ggtitle(expression(omega))+labs(x=expression(omega))+ geom_line(color = "blue", size = 1) +
  geom_segment(mapping=aes(x = model1$par[4], y =loglik[which.min(p4)],xend=model1$par[4],yend=loglik[which.max(p4)] ),color="black", size = 1)+
  geom_segment(mapping=aes(x =  0.2, y =loglik[which.min(p4)],xend=0.2,yend=loglik[66] ),color="red", size = 1)

grid.arrange(g1, g2,g3,g4, ncol=2)


###########################################################################################
###########################################################################################
########################################################################Iterowana filtracja




eta.fun<-function(alpha,ni,lamba){
  delta=1
  ci<-pos(lambda*ni*delta)
  if(sum(ci^2)>0){
  ri<-runif(length(ci))
  temp1=exp(-lambda*delta)*(alpha^(-1))*sum(log(1/ci[1:length(ci)])*exp(lambda*delta*ri[1:length(ci)]))
  temp2=(alpha^(-1))*sum(log(1/ci[1:length(ci)]))
  } else {
    temp1=0
    temp2=0
  }
  return(c(temp1,temp2))
}

eta.fun(alpha,ni,lamba)



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

params_test=c(
  mu=mu,
  lambda=lambda,
  ksi=ksi,
  omega=omega
  
)



###########################################################################################
###########################################################################################
##################################################################################Filtracja

pf1 <- pfilter(bns.model,params=params_test,
               Np=1000,filter.mean=T)
plot(pf1)
str(filter.mean(pf1))
plot(sigma[2:n],type='l')
lines(Kalman.filter(params_test)[2:n],col='red')
lines(filter.mean(pf1)[1,2:n]/lambda,col='royalblue')

wykres2.matrix<-rbind(sigma[2:n],Kalman.filter(params_test)[2:n],
                      filter.mean(pf1)[1,2:n]/lambda,2:n)
row.names(wykres2.matrix)<-c("Symulacja","Filtr Kalmana","Filtr cząsteczkowy","time")
wykres2<-as.data.frame(t(wykres2.matrix))
wykres2.long<-melt(wykres2,id="time")
names(wykres2.long)<-c("time","Trajektoria","value")




ggplot(data = wykres2.long, aes(x = time,y=value,colour=Trajektoria,linetype=Trajektoria) ) +ggtitle('Estymacja procesu ukrytego')+
  geom_line(size =1)+labs(y=expression(sigma[n]))+scale_colour_manual(values=c("black","royalblue","tomato"))+
  scale_linetype_manual(values=c("solid", "solid", "solid"))


round(accuracy(x=sigma[2:n],f=Kalman.filter(params_test)[2:n]),4)
round(accuracy(x=sigma[2:n],f= filter.mean(pf1)[1,2:n]/lambda),4)



###########################################################################################
###########################################################################################
##################################################################################Estymacja


###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 2

#liczba czasteczek
bns_Np <-          c(100, 1e3, 1e3)
bns_Nmif <-        c(10,  50,  150)
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
  ksi=c(0.1,1),
  lambda    = c(0.001,0.02),
  mu = c(-.5,.5),
  omega=c(0.1,1)
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
difftime(end_time,start_time, units = "mins")
plot(if.bns)

if.bns.box  <- data.frame(logLik=L.if.bns[,1],logLik_se=L.if.bns[,2],t(sapply(if.bns,coef)))
if.bns.box [which.max(if.bns.box$logLik),]
round(if.bns.box [which.max(if.bns.box$logLik),],5)

pairs(~logLik+mu+lambda+ksi+omega,data=if.bns.box )
ggpairs(data=if.bns.box, columns = c("logLik","mu","lambda","ksi","omega"), title = "",  
        axisLabels = "internal", columnLabels = c("logLik","mu","lambda","ksi","omega"),
        upper = list(continuous = "points"), diag=list(continuous = "blankDiag"),
        lower = list(continuous = "cor")) 



params_new=c(
  mu=if.bns.box [which.max(if.bns.box$logLik),'mu'],
  lambda=if.bns.box [which.max(if.bns.box$logLik),'lambda'],
  ksi=if.bns.box [which.max(if.bns.box$logLik),'lambda'],
  omega=if.bns.box [which.max(if.bns.box$logLik),'omega']
)


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
 
save(x,x2,params_new,file="wyniki_bns")
load(file="wyniki_bns")
