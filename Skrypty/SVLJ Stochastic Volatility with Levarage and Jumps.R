#model  stochastic volatility with levarage and jumps

#rm(list = ls())
#gc()
setwd("C:/Users/user/Documents/github/Zastosowanie-iterowanej-filtracji-do-modeli-stochastycznej-zmienno-ci/Dane")

#dlugosc szeregu czasowego
n=1000
#parametry
mu = -0.5
phi = 0.98
sigma = 0.2
rho=-0.3
lambda=0.01
sigma_J=2

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




####nazwy
svlj_statenames <- c("H","Y_state","J")
svlj_rp_names <- c("mu","phi","sigma_eta","rho","lambda","sigma_J")
svlj_paramnames <- c(svlj_rp_names)
svlj_covarnames <- "covaryt"



rproc1 <- "
double betal,omega,pt,v,s;
pt=dnorm(Y_state,0,sqrt(exp(H)+sigma_J*sigma_J),0)*lambda/( dnorm(Y_state,0,sqrt(exp(H)+sigma_J*sigma_J),0)*lambda +dnorm(Y_state,0,exp(H/2),0)*(1-lambda));
v=Y_state*exp(-H/2)/(exp(H)+sigma_J*sigma_J);
s=sqrt(sigma_J*sigma_J/(exp(H)+sigma_J*sigma_J));
J=rbinom(1,pt);
betal=(1-J)*Y_state*exp(-H/2)+J*rnorm(v,s);
omega = rnorm(0,sigma_eta * sqrt(1-rho*rho));
H = mu*(1 - phi) + phi*H + betal * rho + omega;
"

####rownanie procesu pomiaru
rproc2.sim <- "
Y_state = rnorm( 0,exp(H/2) )+rbinom(1,lambda)*rnorm(0,sigma_J);
"
###do wypelniania danych
rproc2.filt <- "
Y_state = covaryt;
"

###symulacja modelu SVL
svlj_rproc.sim <- paste(rproc1,rproc2.sim)

####filtr czasteczkowy 
svlj_rproc.filt <- paste(rproc1,rproc2.filt)


######inicalizacja
svlj_initializer <- "
H = rnorm(mu,sigma_eta/sqrt((1-phi*phi))) ;
Y_state = rnorm( 0,exp(H/2) )+rbinom(1,lambda)*rnorm(0,sigma_J);
"
###????
svlj_rmeasure <- "
y=rnorm( 0,exp(H/2) )+rbinom(1,lambda)*rnorm(0,sigma_J);
"

####rozk?ad warunkowy zmiennej Y
svlj_dmeasure <- "
lik=dnorm(y,0,exp(H/2),give_log)*(1-lambda)+dnorm(y,0,sqrt(exp(H)+sigma_J*sigma_J),give_log)*lambda;
"


####przeskalowanie parametr?w 
svlj_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tsigma_J=log(sigma_J);
Tphi = logit(phi);
Tlambda = logit(lambda);
Trho = 0.5*log((rho+1)/(1-rho));
"

svlj_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tsigma_J=exp(sigma_J);
Tlambda = expit(lambda);
Tphi = expit(phi);
Trho= (exp(2*rho)-1)/(exp(2*rho)+1);
"


####wypelnianie modelu danymi
svlj.model <- pomp(data=data.frame(y=1:n, time=1:n),
                  statenames=svlj_statenames,
                  paramnames=svlj_paramnames,
                  covarnames=svlj_covarnames,
                  times="time",
                  t0=0,
                  covar=data.frame(covaryt=c(0:n),
                                   time=0:n),
                  tcovar="time",
                  rmeasure=Csnippet(svlj_rmeasure),
                  dmeasure=Csnippet(svlj_dmeasure),
                  rprocess=discrete.time.sim(step.fun=Csnippet(svlj_rproc.sim),delta.t=1),
                  initializer=Csnippet(svlj_initializer),
                  toEstimationScale=Csnippet(svlj_toEstimationScale), 
                  fromEstimationScale=Csnippet(svlj_fromEstimationScale)
)


plot(svlj.model)




##testowa wersja parametr?w
params_test <- c(
  mu = mu,       
  phi =phi,  
  sigma_eta = sigma,
  rho=rho,
  lambda=lambda,
  sigma_J= sigma_J
)



#######################################ustawienia dla pomp


###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 1
#liczba czasteczek
svl_Np <-          c(100, 1e3, 1e3)
svl_Nmif <-        c(10,  50,  150)
svl_Nreps_eval <-  c(4,   10,  20)
svl_Nreps_local <- c(4,   10,  20)
svl_Nreps_global <-c(4,   10,  20)

svllist<-list(svl_Np ,svl_Nmif,svl_Nreps_eval,
              svl_Nreps_local,svl_Nreps_global )

#parametry do metody mif2
svl_rw.sd_rp <- 0.01
svl_rw.sd_ivp <- 0.1
svl_cooling.fraction.50 <- 0.5

##########################################################################
##########################################################################
########################################################symulowanie danych


log.volatility.list  <-simulate(svlj.model,nsim=1,seed=12,params=params_test,
                                states = TRUE, obs = TRUE)
zwroty.sim= log.volatility.list$obs[1,1,]
log.volatility.sim= log.volatility.list$states[1,1,] 

par(mfrow=c(3,1))
plot(log.volatility.list$obs[1,1,],type='l',ylab="y",main="zwroty") 
plot(log.volatility.list$states[1,1,],type='l',ylab="h",main="log-zmiennosc")
plot(log.volatility.list$states[3,1,],type='l',ylab="J",main="Skoki")
par(mfrow=c(1,1))




##########################################################################
##########################################################################
########################################################Filtrowanie danych




svlj.filt<-pomp(data=data.frame(y=zwroty.sim,
                                time=1:n),
                statenames=svlj_statenames,
                paramnames=svlj_paramnames,
                covarnames=svlj_covarnames,
                times="time",
                t0=0,
                covar=data.frame(covaryt=c(0,zwroty.sim),
                                 time=0:n),
                tcovar="time",
                rmeasure=Csnippet(svlj_rmeasure),
                dmeasure=Csnippet(svlj_dmeasure),
                rprocess=discrete.time.sim(step.fun=Csnippet(svlj_rproc.filt),delta.t=1),
                initializer=Csnippet(svlj_initializer),
                toEstimationScale=Csnippet(svlj_toEstimationScale), 
                fromEstimationScale=Csnippet(svlj_fromEstimationScale)
)

#wykres filtracji dla prawdziwych parametr?w
pf1 <- pfilter(svlj.filt,params=params_test,
               Np=1000,filter.traj=T)
plot(pf1)






##########################################################################
##########################################################################
######################################################Estymacja parametrow

svlj_box <- rbind(
  sigma_eta=c(0.1,1),
  phi    = c(0.9,0.99),
  mu = c(-1,01),
  sigma_eta=c(0.1,4),
  lambda=c(0.001,0.1), 
  rho=c(-0.5,0)
)

detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


start_time <- Sys.time()
t.if.svlj <- system.time({
  if.svlj <- foreach(i=1:svljlist[[5]][run_level] ,
                     .packages='pomp', .combine=c,.export = "svljlist", 
                     .options.multicore=list(set.seed=TRUE)) %dopar% try(
                       pomp::mif2(svlj.filt,start=apply(svlj_box,1,function(x) runif(1,x[1],x[2])),Np=svljlist[[1]][run_level] , Nmif=svljlist[[2]][run_level] ,cooling.type="geometric",
                                  cooling.fraction.50=svlj_cooling.fraction.50,
                                  transform=TRUE,
                                  rw.sd = rw.sd(
                                    mu      = svlj_rw.sd_rp,
                                    phi       = svlj_rw.sd_rp,
                                    sigma_eta = svlj_rw.sd_rp,
                                    sigma_J = svlj_rw.sd_rp,
                                    lambda = svlj_rw.sd_rp,
                                    rho    = svlj_rw.sd_rp
                                  )
                       )
                       
                     )
  
  L.if.svlj <- foreach(i=1:svljlist[[5]][run_level] ,.packages='pomp', .export = "svljlist", 
                       .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                         logmeanexp(
                           replicate(svljlist[[3]][run_level] ,
                                     logLik(pfilter(svlj.filt,params=coef(if.svlj[[i]]),Np=svljlist[[1]][run_level]  ))
                           ),
                           se=TRUE)
                       )
  
  H.if.svlj<- foreach(i=1:svljlist[[5]][run_level] ,.packages='pomp', .export = "svljlist", 
                      .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                        exp(pfilter(svlj.filt,params=coef(if.svlj[[i]]),Np=svljlist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                      )
})


stopCluster(cl)
end_time <- Sys.time()
difftime(end_time,start_time, units = "mins")
plot(if.svlj)


if.svlj.box  <- data.frame(logLik=L.if.svlj[,1],logLik_se=L.if.svlj[,2],t(sapply(if.svlj,coef)))
if.svlj.box [which.max(if.svlj.box$logLik),]
pairs(~logLik+mu+phi+sigma_eta+rho,data=if.svlj.box )




