#eksperyment symulacyjny dla modelu basic stochastic volatility

#rm(list = ls())
#gc()
setwd("C:/Users/user/Documents/github/Zastosowanie-iterowanej-filtracji-do-modeli-stochastycznej-zmienno-ci/Dane")


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
##########################################################################
##########################################################################
################################################################ustawienia

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
######################################################implementacja modelu


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

##########################################################################
##########################################################################
########################################################symulowanie danych

log.volatility.list  <-simulate(bsv.model,nsim=1,seed=123,params=c(mu=mu,phi=phi,sigma_eta=sigma),
                                      states = TRUE, obs = TRUE)
zwroty.sim= log.volatility.list$obs[1,1,]
log.volatility.sim= log.volatility.list$states[1,1,] 

par(mfrow=c(3,1))
plot(log.volatility.list$obs[1,1,],type='l',ylab="y",main="zwroty") 
plot(log.volatility.list$states[1,1,],type='l',ylab="h",main="log-zmiennosc")
plot(exp(log.volatility.list$states[1,1,]/2),type='l',ylab="exp(h/2)",main="zmiennosc")
par(mfrow=c(1,1))


wykres1= rbind(log.volatility.list$obs[1,1,],
              log.volatility.list$states[1,1,],
              exp(log.volatility.list$states[1,1,]/2))
row.names(wykres1) = c("y","H","exp(H/2)")
wykres1=as.data.frame(t(wykres1))

g1<-ggplot(data = as.data.frame(wykres1), aes(x = 1:n, y = y))  + 
  geom_line(color = "royalblue", size = .5)+ggtitle('zwroty')+labs(x="time")+theme_bw()
g2<-ggplot(data = as.data.frame(wykres1), aes(x = 1:n, y = H))  + 
  geom_line(color = "royalblue", size = .5)+ggtitle('log-zmienność')+labs(x="time")+theme_bw()
g3<-ggplot(data = as.data.frame(wykres1), aes(x = 1:n, y = exp(H/2)))  + 
  geom_line(color = "royalblue", size = .5)+ggtitle('zmienność')+labs(x="time")+theme_bw()
grid.arrange(g1, g2,g3, nrow=3)



##########################################################################
##########################################################################
########################################################Filtrowanie danych



bsv.filt<-pomp(data=data.frame(y=zwroty.sim,
                               time=1:n),
               statenames=bsv_statenames,
               paramnames=bsv_paramnames,
               covarnames=bsv_covarnames,
               times="time",
               t0=0,
               covar=data.frame(covaryt=c(0,zwroty.sim),
                                time=0:n),
               tcovar="time",
               rmeasure=Csnippet(bsv_rmeasure),
               dmeasure=Csnippet(bsv_dmeasure),
               rprocess=discrete.time.sim(step.fun=Csnippet(bsv_rproc.filt),delta.t=1),
               initializer=Csnippet(bsv_initializer),
               toEstimationScale=Csnippet(bsv_toEstimationScale), 
               fromEstimationScale=Csnippet(bsv_fromEstimationScale)
)

#wykres filtracji dla prawdziwych parametr?w
pf1 <- pfilter(bsv.filt,params=params_test,
               Np=1000,filter.traj=T)
plot(pf1)
coef(pf1)


wykres2.matrix<-rbind(pf1@data,
               pf1@eff.sample.size,
               pf1@cond.loglik,
               pf1@filter.traj[1, 1,2:1001],
               exp(pf1@filter.traj[1, 1,2:1001]/2),
               log.volatility.list$states[1,1,],
               exp(log.volatility.list$states[1,1,]/2),
               1:n)
row.names(wykres2.matrix)<-c("y","eff.sample.size","cond.loglik",
                     "filtered_H", "filtered_exp","H","exp","time")
wykres2<-as.data.frame(t(wykres2.matrix))
wykres2a<-as.data.frame(t(wykres2.matrix[c(6,4,8),]))
names(wykres2a)<-c("Symulacja","Filtracja", "time")  
wykres2a.long<-melt(wykres2a,id="time")
names(wykres2a.long)<-c("time","Rodzaj","value")

wykres2b<-as.data.frame(t(wykres2.matrix[c(7,5,8),]))
names(wykres2b)<-c("Symulacja","Filtracja", "time")  
wykres2b.long<-melt(wykres2b,id="time")
names(wykres2b.long)<-c("time","Rodzaj","value")


g1<-ggplot(data = wykres2, aes(x = 1:n, y = y))  + 
  geom_line(color = "royalblue", size = .5)+ggtitle('zwroty')+labs(x="time")+theme_bw()


g2<-ggplot(data = wykres2, aes(x = time))  + 
  geom_line(aes(y = H), colour="royalblue",size = .5)+
  geom_line(aes(y = filtered_H), colour="red",size = .5,linetype = "dashed")+
  ggtitle('log-zmienność')+labs(x="time")+theme_bw()

#druga wersja
g2a<-ggplot(data = wykres2a.long, aes(x = time,y=value,colour=Rodzaj,linetype=Rodzaj) ) +ggtitle('log-zmienność')+
  geom_line(size = .5)+labs(y="H")+theme_bw()+scale_colour_manual(values=c("royalblue","tomato"))+
  scale_linetype_manual(values=c("solid", "dashed"))


g3<-ggplot(data = wykres2, aes(x = 1:n ))  + 
  geom_line( aes(y = exp), colour="royalblue",size = .5)+
  geom_line( aes(y = filtered_exp), colour="red",size = .5,linetype = "dashed")+
  ggtitle( 'zmienność')+labs(x="time",y="exp(H/2)")+theme_bw()
#druga wersja
g3a<-ggplot(data = wykres2b.long, aes(x = time,y=value,colour=Rodzaj,linetype=Rodzaj) ) +ggtitle('zmienność')+
  geom_line(size = .5)+labs(y="exp(H/2)")+theme_bw()+scale_colour_manual(values=c("royalblue","tomato"))+
  scale_linetype_manual(values=c("solid", "dashed"))


g4<-ggplot(data = wykres2, aes(x = 1:n, y = eff.sample.size))  + 
  geom_line(color = "royalblue", size = .5)+ggtitle('eff.sample.size')+labs(x="time")+theme_bw()


g5<-ggplot(data = wykres2, aes(x = 1:n, y = cond.loglik))  + 
  geom_line(color = "royalblue", size = .5)+ggtitle('cond.loglik')+labs(x="time")+theme_bw()

grid.arrange(g1, g4,g2,g5,g3, nrow=3,ncol=2)

grid.arrange(g1, g2a,g4, g3a,g5, nrow=3,ncol=2)

##########################################################################
##########################################################################
######################################################Estymacja parametrow

bsv_box <- rbind(
  sigma_eta=c(0.1,1),
  phi    = c(0.9,0.99),
  mu = c(-1,01)
)

detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


start_time <- Sys.time()
t.if.bsv <- system.time({
  if.bsv <- foreach(i=1:bsvlist[[5]][run_level] ,
                    .packages='pomp', .combine=c,.export = "bsvlist", 
                    .options.multicore=list(set.seed=TRUE)) %dopar% try(
                      pomp::mif2(bsv.filt,start=apply(bsv_box,1,function(x) runif(1,x[1],x[2])),Np=bsvlist[[1]][run_level] , Nmif=bsvlist[[2]][run_level] ,cooling.type="geometric",
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
difftime(end_time,start_time, units = "mins")
plot(if.bsv)


if.bsv.box  <- data.frame(logLik=L.if.bsv[,1],logLik_se=L.if.bsv[,2],t(sapply(if.bsv,coef)))
if.bsv.box [which.max(if.bsv.box$logLik),]
pairs(~logLik+mu+phi+sigma_eta,data=if.bsv.box )


params_nowe2<- c(
  mu        = as.numeric(if.bsv.box[which.max(if.bsv.box $logLik),'mu']),    
  phi         = as.numeric(if.bsv.box[which.max(if.bsv.box $logLik),'phi']),    
  sigma_eta   = as.numeric(if.bsv.box[which.max(if.bsv.box$logLik),'sigma_eta']) 
)

ggpairs(data=if.bsv.box, columns = c("logLik","mu","phi","sigma_eta"), title = "",  
        axisLabels = "internal", columnLabels = c("logLik","mu","phi","sigma_eta"),
        upper = list(continuous = "points"), diag=list(continuous = "blankDiag"),
        lower = list(continuous = "cor"))+theme_bw()

     
wyniki.mif.loglik<-matrix(NaN,ncol=bsvlist[[5]][run_level]+1,nrow= bsvlist[[2]][run_level])
nazwy=1:bsvlist[[5]][run_level]
for(i in 1:bsvlist[[5]][run_level]) wyniki.mif.loglik[,i]=conv.rec(if.bsv)[[i]][1:bsvlist[[2]][run_level] ,1] 
for(i in 1:bsvlist[[5]][run_level]) nazwy[i]=paste(i)
wyniki.mif.loglik[,bsvlist[[5]][run_level]+1]=1:bsvlist[[2]][run_level]
wyniki.mif.loglik=as.data.frame(wyniki.mif.loglik)
names(wyniki.mif.loglik)<-c(nazwy,"nr")
wyniki.mif.loglik.long<-melt(wyniki.mif.loglik,id="nr")
names(wyniki.mif.loglik.long)<-c("Nr_iteracji","IF","value")

g1<-ggplot(data = wyniki.mif.loglik.long, aes(x = Nr_iteracji,y=value,colour=IF,linetype=IF) ) +ggtitle('logLik')+
  geom_line(size = 1)+labs(x="Nr iteracji",y="Loglik")+theme_bw()+
  scale_x_discrete(limits=1:bsvlist[[2]][run_level])


wyniki.mif.mu<-matrix(NaN,ncol=bsvlist[[5]][run_level]+1,nrow= bsvlist[[2]][run_level]+1)
for(i in 1:bsvlist[[5]][run_level]) wyniki.mif.mu[,i]=conv.rec(if.bsv)[[i]][,5] 
wyniki.mif.mu[,bsvlist[[5]][run_level]+1]=0:bsvlist[[2]][run_level]
wyniki.mif.mu=as.data.frame(wyniki.mif.mu)
names(wyniki.mif.mu)<-c(nazwy,"nr")
wyniki.mif.mu.long<-melt(wyniki.mif.mu,id="nr")
names(wyniki.mif.mu.long)<-c("Nr_iteracji","IF","value")

g2<-ggplot(data = wyniki.mif.mu.long, aes(x = Nr_iteracji,y=value,colour=IF,linetype=IF) ) +ggtitle(expression(mu))+
  geom_line(size = 1)+labs(x="Nr iteracji",y=expression(mu))+theme_bw()+
  scale_x_discrete(limits=0:bsvlist[[2]][run_level])



wyniki.mif.phi<-matrix(NaN,ncol=bsvlist[[5]][run_level]+1,nrow= bsvlist[[2]][run_level]+1)
for(i in 1:bsvlist[[5]][run_level]) wyniki.mif.phi[,i]=conv.rec(if.bsv)[[i]][,4] 
wyniki.mif.phi[,bsvlist[[5]][run_level]+1]=0:bsvlist[[2]][run_level]
wyniki.mif.phi=as.data.frame(wyniki.mif.phi)
names(wyniki.mif.phi)<-c(nazwy,"nr")
wyniki.mif.phi.long<-melt(wyniki.mif.phi,id="nr")
names(wyniki.mif.phi.long)<-c("Nr_iteracji","IF","value")

g3<-ggplot(data = wyniki.mif.phi.long, aes(x = Nr_iteracji,y=value,colour=IF,linetype=IF) ) +ggtitle(expression(phi))+
  geom_line(size = 1)+labs(x="Nr iteracji",y=expression(phi))+theme_bw()+
  scale_x_discrete(limits=0:bsvlist[[2]][run_level])



wyniki.mif.sigma<-matrix(NaN,ncol=bsvlist[[5]][run_level]+1,nrow= bsvlist[[2]][run_level]+1)
for(i in 1:bsvlist[[5]][run_level]) wyniki.mif.sigma[,i]=conv.rec(if.bsv)[[i]][,3] 
wyniki.mif.sigma[,bsvlist[[5]][run_level]+1]=0:bsvlist[[2]][run_level]
wyniki.mif.sigma=as.data.frame(wyniki.mif.sigma)
names(wyniki.mif.sigma)<-c(nazwy,"nr")
wyniki.mif.sigma.long<-melt(wyniki.mif.sigma,id="nr")
names(wyniki.mif.sigma.long)<-c("Nr_iteracji","IF","value")

g4<-ggplot(data = wyniki.mif.sigma.long, aes(x = Nr_iteracji,y=value,colour=IF,linetype=IF) ) +ggtitle(expression(sigma[eta]))+
  geom_line(size = 1)+labs(x="Nr iteracji",y=expression(sigma[eta]))+theme_bw()+
  scale_x_discrete(limits=0:bsvlist[[2]][run_level])


grid.arrange(g1, g2,g3,g4,  nrow=2,ncol=2)


##########################################################################
##########################################################################
###########################################################Profile funkcji 

library(magrittr)
library(plyr)



###############################################################parametr: mu

#szybka wersja
profileDesign(
  mu=seq(from=-1,to=0,length=21),
  lower=c(phi=phi,sigma_eta=sigma),upper=c(phi=phi,sigma_eta=sigma),
  nprof=2
) -> pd
str(pd)

#dokladna wersja
#profileDesign(
#  mu=seq(from=-1,to=0,length=21),
#  lower=c(phi=0.95,sigma_eta=0.01),upper=c(phi=0.99,sigma_eta=0.5),
#  nprof=10
#) -> pd
#str(pd)

pairs(~mu+phi+sigma_eta,data=pd)


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


foreach (p=iter(pd,"row"),
         .combine=rbind,
         .errorhandling="remove",
         .packages=c("pomp","magrittr","reshape2","plyr"),
         .export="bsv.filt",.inorder=FALSE
) %dopar%
{
  bsv.filt %>% 
    mif2(start=unlist(p),Nmif=1,Np=1000,transform=TRUE,
         cooling.fraction.50=0.8,cooling.type="geometric",
         rw.sd= rw.sd(phi = 0.02,
                sigma_eta = 0.02)) %>%
    mif2() -> mf
  
  pf <- replicate(5,pfilter(mf,Np=1000))  ## independent particle filters
  ll <- sapply(pf,logLik)
  ll <- logmeanexp(ll,se=TRUE)
 nfail <- sapply(pf,getElement,"nfail")  ## number of filtering failures
  
  data.frame(as.list(coef(mf)),
             loglik = ll[1],
             loglik.se = ll[2],
             nfail.min = min(nfail),
             nfail.max = max(nfail))
} %>% arrange(mu,-loglik) -> mu_prof

stopCluster(cl)





pairs(~loglik+mu+phi+sigma_eta,data=mu_prof,subset=loglik>max(loglik)-10)

library(plyr)
mu_prof %>% 
  mutate(mu=signif(mu,8)) %>%
  ddply(~mu,subset,loglik==max(loglik)) %>%
  ggplot(aes(x=mu,y=loglik))+
  geom_point()+geom_smooth()+
  theme_bw()->g1




###############################################################parametr: phi

#szybka wersja
profileDesign(
  phi=seq(from=0.95,to=0.99,length=20),
  lower=c(mu=mu,sigma_eta=sigma),upper=c(mu=mu,sigma_eta=sigma),
  nprof=2
) -> pd

#dokladniejsza wersja
#profileDesign(
#  phi=seq(from=0.9,to=0.99,length=20),
#  lower=c(mu=-1,sigma_eta=0.01),upper=c(mu=0,sigma_eta=0.5),
#  nprof=10
#) -> pd

str(pd)

pairs(~phi+mu+sigma_eta,data=pd)


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


foreach (p=iter(pd,"row"),
         .combine=rbind,
         .errorhandling="remove",
         .packages=c("pomp","magrittr","reshape2","plyr"),
         .export="bsv.filt",.inorder=FALSE
) %dopar%
{
  bsv.filt %>% 
    mif2(start=unlist(p),Nmif=1,Np=1000,transform=TRUE,
         cooling.fraction.50=0.8,cooling.type="geometric",
         rw.sd= rw.sd(phi = 0.02,
                      mu = 0.02)) %>%
    mif2() -> mf
  
  pf <- replicate(5,pfilter(mf,Np=1000))  ## independent particle filters
  ll <- sapply(pf,logLik)
  ll <- logmeanexp(ll,se=TRUE)
  nfail <- sapply(pf,getElement,"nfail")  ## number of filtering failures
  
  data.frame(as.list(coef(mf)),
             loglik = ll[1],
             loglik.se = ll[2],
             nfail.min = min(nfail),
             nfail.max = max(nfail))
} %>% arrange(phi,-loglik) -> phi_prof

stopCluster(cl)





pairs(~loglik+mu+phi+sigma_eta,data=phi_prof,subset=loglik>max(loglik)-10)

library(plyr)
phi_prof %>% 
  mutate(phi=signif(phi,8)) %>%
  ddply(~phi,subset,loglik==max(loglik)) %>%
  ggplot(aes(x=phi,y=loglik))+
  geom_point()+geom_smooth()+
  theme_bw()->g2





########################################################parametr: sigma_eta

#szybka wersja 
profileDesign(
  sigma_eta=seq(from=0.01,to=0.5,length=20),
  lower=c(mu=-1,phi=0.95),upper=c(mu=0,phi=0.99),
  nprof=2
) -> pd

#dokladna wersja
profileDesign(
  sigma_eta=seq(from=0.01,to=0.5,length=20),
  lower=c(mu=mu,phi=phi),upper=c(mu=mi,phi=phi),
  nprof=10
) -> pd

str(pd)

pairs(~mu+phi+sigma_eta,data=pd)


detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


foreach (p=iter(pd,"row"),
         .combine=rbind,
         .errorhandling="remove",
         .packages=c("pomp","magrittr","reshape2","plyr"),
         .export="bsv.filt",.inorder=FALSE
) %dopar%
{
  bsv.filt %>% 
    mif2(start=unlist(p),Nmif=1,Np=1000,transform=TRUE,
         cooling.fraction.50=0.8,cooling.type="geometric",
         rw.sd= rw.sd(phi = 0.02,
                      mu = 0.02)) %>%
    mif2() -> mf
  
  pf <- replicate(5,pfilter(mf,Np=1000))  ## independent particle filters
  ll <- sapply(pf,logLik)
  ll <- logmeanexp(ll,se=TRUE)
  nfail <- sapply(pf,getElement,"nfail")  ## number of filtering failures
  
  data.frame(as.list(coef(mf)),
             loglik = ll[1],
             loglik.se = ll[2],
             nfail.min = min(nfail),
             nfail.max = max(nfail))
} %>% arrange(sigma_eta,-loglik) -> sigma_eta_prof

stopCluster(cl)





pairs(~loglik+mu+phi+sigma_eta,data=sigma_eta_prof,subset=loglik>max(loglik)-10)

library(plyr)
sigma_eta_prof %>% 
  mutate(sigma_eta=signif(sigma_eta,8)) %>%
  ddply(~sigma_eta,subset,loglik==max(loglik)) %>%
  ggplot(aes(x=sigma_eta,y=loglik))+
  geom_point()+geom_smooth()+
  theme_bw()->g3

#################################################################podsumowanie


grid.arrange(g1, g2,g3, ncol=3)

############################################################################
############################################################################
############################################################################
######################################profile funkcji wiarygodnosci 2 wersja


params_nowe2=params_test

#################################################################parametr mu
xx1<-seq(from=-1.5,to=0.5,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.bsv.log<- foreach(i=1:length(xx1) ,.packages='pomp', .export = "bsvlist",.combine=rbind,
                    .options.multicore=list(set.seed=TRUE)) %dopar% {
                      set.seed(87932)
                      logLik(pfilter(bsv.filt,params=c(mu=xx1[i], phi=as.numeric(params_nowe2['phi']),sigma_eta=as.numeric(params_nowe2['sigma_eta'])),
                                     Np=bsvlist[[1]][run_level] ))
                    }

stopCluster(cl)
beep(2)

plot(xx1, L.bsv.log, type='l',xlab=expression(mu),ylab="logLik")
#abline(v=params_nowe2['mu_h'],lty=2)
#points(xx1, L.bsv.log)
points(if.bsv.box [,'mu'], if.bsv.box[,'logLik'] ,col='red',pch=16)
p=loess(L.bsv.log~xx1,span=0.5)
lines(xx1,p$fitted,col='blue',lwd=2)


wykres.mu=as.data.frame(t(rbind(xx1,as.vector(L.bsv.log))))
names(wykres.mu)<-c("mu","loglik")
g1<-ggplot(data = wykres.mu, aes(x = mu, y = loglik))  + 
  geom_point(color = "black", size = 1)+
  ggtitle(expression(mu))+labs(x=expression(mu))+geom_smooth()+theme_bw()


#parametr phi
xx2<-seq(from=0.95,to=0.99,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.bsv.log2<- foreach(i=1:length(xx2) ,.packages='pomp', .export = "bsvlist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                       logLik(pfilter(bsv.filt,params=c(mu=as.numeric(params_nowe2['mu']), phi=xx2[i],sigma_eta=as.numeric(params_nowe2['sigma_eta'])),
                                      Np=bsvlist[[1]][run_level] ))
                     }
stopCluster(cl)
beep(2)

plot(xx2, L.bsv.log2, type='l',xlab=expression(phi),ylab="logLik")
#points(xx2, L.bsv.log2)
points(if.bsv.box[,'phi'], if.bsv.box[,'logLik'] ,col='red',pch=16)
p2=loess(L.bsv.log2~xx2, span=0.5)
lines(xx2,p2$fitted,col='blue',lwd=2)

wykres.phi=as.data.frame(t(rbind(xx2,as.vector(L.bsv.log2))))
names(wykres.phi)<-c("phi","loglik")
g2<-ggplot(data = wykres.phi, aes(x = phi, y = loglik))  + 
  geom_point(color = "black", size = 1)+
  ggtitle(expression(phi))+labs(x=expression(phi))+geom_smooth()+theme_bw()



#parametr sigma_eta
xx3<-seq(from=params_nowe2['sigma_eta']-.1,to=params_nowe2['sigma_eta']+.1,length.out = 100)
detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)
L.bsv.log3<- foreach(i=1:length(xx3) ,.packages='pomp', .export = "bsvlist",.combine=rbind,
                     .options.multicore=list(set.seed=TRUE)) %dopar% {
                       set.seed(87932)
                       logLik(pfilter(bsv.filt,params=c(mu=as.numeric(params_nowe2['mu']), phi=as.numeric(params_nowe2['phi']),sigma_eta=xx3[i]),
                                      Np=bsvlist[[1]][run_level] ))
                     }

stopCluster(cl)
beep(2)

plot(xx3, L.bsv.log3, type='l',xlab=expression(sigma[eta]),ylab="logLik")
#points(xx3, L.bsv.log3)
points(if.bsv.box[,'sigma_eta'], if.bsv.box[,'logLik'] ,col='red',pch=16)
p3=loess(L.bsv.log3~xx3,span=0.5)
lines(xx3,p3$fitted,col='blue',lwd=2)


wykres.sigma_eta=as.data.frame(t(rbind(xx3,as.vector(L.bsv.log3))))
names(wykres.sigma_eta)<-c("sigma_eta","loglik")
g3<-ggplot(data = wykres.sigma_eta, aes(x = sigma_eta, y = loglik))  + 
  geom_point(color = "black", size = 1)+
  ggtitle(expression(sigma[eta]))+labs(x=expression(sigma[eta]))+geom_smooth()+theme_bw()


#################################################################podsumowanie


grid.arrange(g1, g2,g3, ncol=3)
