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

plot(fsv.model)


#######################################ustawienia dla pomp

params_test=c(
  mu=mu,
  phi=phi,
  sigma_eta=sigma,
  gamma=gamma
)

###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 1

#liczba czasteczek
FSV_Np <-          c(100,1e3,1e3)
FSV_Nmif <-        c(10, 50,150)
FSV_Nreps_eval <-  c(4,  10,  20)
FSV_Nreps_local <- c(4, 10, 10)
FSV_Nreps_global <-c(4, 10, 10)

FSVlist<-list(bsv_Np ,bsv_Nmif,bsv_Nreps_eval,
              bsv_Nreps_local,bsv_Nreps_global )

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
               statenames=bsv_statenames,
               paramnames=bsv_paramnames,
               covarnames=bsv_covarnames,
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
