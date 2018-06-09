#model  stochastic volatility with levarage

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
svl_statenames <- c("H","Y_state")
svl_rp_names <- c("mu","phi","sigma_eta","rho")
svl_paramnames <- c(svl_rp_names)
svl_covarnames <- "covaryt"



rproc1 <- "
double beta,omega;
omega = rnorm(0,sigma_eta * sqrt(1-rho*rho));
beta = Y_state * sigma_eta ;
H = mu*(1 - phi) + phi*H + beta * rho * exp(-H/2) + omega;
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
svl_rproc.sim <- paste(rproc1,rproc2.sim)

####filtr czasteczkowy 
svl_rproc.filt <- paste(rproc1,rproc2.filt)


######inicalizacja
svl_initializer <- "
H = rnorm(mu,sigma_eta/sqrt((1-phi*phi))) ;
Y_state = rnorm( 0,exp(H/2) );
"
###????
svl_rmeasure <- "
y=rnorm(0,exp(H/2));
"

####rozk?ad warunkowy zmiennej Y
svl_dmeasure <- "
lik=dnorm(y,0,exp(H/2),give_log);
"


####przeskalowanie parametr?w 
svl_toEstimationScale <- "
Tsigma_eta = log(sigma_eta);
Tphi = logit(phi);
Trho = 0.5*log((rho+1)/(1-rho));
"

svl_fromEstimationScale <- "
Tsigma_eta = exp(sigma_eta);
Tphi = expit(phi);
Trho= (exp(2*rho)-1)/(exp(2*rho)+1);
"


####wypelnianie modelu danymi
svl.model <- pomp(data=data.frame(y=1:n, time=1:n),
                 statenames=svl_statenames,
                 paramnames=svl_paramnames,
                 covarnames=svl_covarnames,
                 times="time",
                 t0=0,
                 covar=data.frame(covaryt=c(0:n),
                                  time=0:n),
                 tcovar="time",
                 rmeasure=Csnippet(svl_rmeasure),
                 dmeasure=Csnippet(svl_dmeasure),
                 rprocess=discrete.time.sim(step.fun=Csnippet(svl_rproc.sim),delta.t=1),
                 initializer=Csnippet(svl_initializer),
                 toEstimationScale=Csnippet(svl_toEstimationScale), 
                 fromEstimationScale=Csnippet(svl_fromEstimationScale)
)


plot(svl.model)




##testowa wersja parametr?w
params_test <- c(
  mu = mu,       
  phi =phi,  
  sigma_eta = sigma,
  rho=rho
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


log.volatility.list  <-simulate(svl.model,nsim=1,seed=123,params=params_test,
                                states = TRUE, obs = TRUE)
zwroty.sim= log.volatility.list$obs[1,1,]
log.volatility.sim= log.volatility.list$states[1,1,] 

par(mfrow=c(3,1))
plot(log.volatility.list$obs[1,1,],type='l',ylab="y",main="zwroty") 
plot(log.volatility.list$states[1,1,],type='l',ylab="h",main="log-zmiennosc")
plot(exp(log.volatility.list$states[1,1,]/2),type='l',ylab="exp(h/2)",main="zmiennosc")
par(mfrow=c(1,1))


##########################################################################
##########################################################################
########################################################Filtrowanie danych


svl.filt<-pomp(data=data.frame(y=zwroty.sim,
                               time=1:n),
               statenames=svl_statenames,
               paramnames=svl_paramnames,
               covarnames=svl_covarnames,
               times="time",
               t0=0,
               covar=data.frame(covaryt=c(0,zwroty.sim),
                                time=0:n),
               tcovar="time",
               rmeasure=Csnippet(svl_rmeasure),
               dmeasure=Csnippet(svl_dmeasure),
               rprocess=discrete.time.sim(step.fun=Csnippet(svl_rproc.filt),delta.t=1),
               initializer=Csnippet(svl_initializer),
               toEstimationScale=Csnippet(svl_toEstimationScale), 
               fromEstimationScale=Csnippet(svl_fromEstimationScale)
)

#wykres filtracji dla prawdziwych parametr?w
pf1 <- pfilter(svl.filt,params=params_test,
               Np=1000,filter.traj=T)
plot(pf1)


##########################################################################
##########################################################################
######################################################Estymacja parametrow

svl_box <- rbind(
  sigma_eta=c(0.1,1),
  phi    = c(0.9,0.99),
  mu = c(-1,01),
  rho=c(-0.5,0)
)

detectCores()
cl <- makeCluster(3, type = "SOCK")
registerDoSNOW(cl)


start_time <- Sys.time()
t.if.svl <- system.time({
  if.svl <- foreach(i=1:svllist[[5]][run_level] ,
                    .packages='pomp', .combine=c,.export = "svllist", 
                    .options.multicore=list(set.seed=TRUE)) %dopar% try(
                      pomp::mif2(svl.filt,start=apply(svl_box,1,function(x) runif(1,x[1],x[2])),Np=svllist[[1]][run_level] , Nmif=svllist[[2]][run_level] ,cooling.type="geometric",
                                 cooling.fraction.50=svl_cooling.fraction.50,
                                 transform=TRUE,
                                 rw.sd = rw.sd(
                                   mu      = svl_rw.sd_rp,
                                   phi       = svl_rw.sd_rp,
                                   sigma_eta = svl_rw.sd_rp,
                                   rho    = svl_rw.sd_rp
                                 )
                      )
                      
                    )
  
  L.if.svl <- foreach(i=1:svllist[[5]][run_level] ,.packages='pomp', .export = "svllist", 
                      .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                        logmeanexp(
                          replicate(svllist[[3]][run_level] ,
                                    logLik(pfilter(svl.filt,params=coef(if.svl[[i]]),Np=svllist[[1]][run_level]  ))
                          ),
                          se=TRUE)
                      )
  
  H.if.svl<- foreach(i=1:svllist[[5]][run_level] ,.packages='pomp', .export = "svllist", 
                     .combine=cbind,.options.multicore=list(set.seed=TRUE)) %dopar% try(
                       exp(pfilter(svl.filt,params=coef(if.svl[[i]]),Np=svllist[[1]][run_level],pred.mean=TRUE)@pred.mean[1,])
                     )
})


stopCluster(cl)
end_time <- Sys.time()
difftime(end_time,start_time, units = "mins")
plot(if.svl)


if.svl.box  <- data.frame(logLik=L.if.svl[,1],logLik_se=L.if.svl[,2],t(sapply(if.svl,coef)))
if.svl.box [which.max(if.svl.box$logLik),]
pairs(~logLik+mu+phi+sigma_eta+rho,data=if.svl.box )


params_nowe2<- c(
  mu        = as.numeric(if.svl.box[which.max(if.svl.box $logLik),'mu']),    
  phi         = as.numeric(if.svl.box[which.max(if.svl.box $logLik),'phi']),    
  sigma_eta   = as.numeric(if.svl.box[which.max(if.svl.box$logLik),'sigma_eta']),
  rho       = as.numeric(if.svl.box[which.max(if.svl.box$logLik),'rho'])
)

############################################################################
########################################################################PMCM

hyperparams <- list(min = c(-1,0.9,0,-0.5), max = c(0,1,1,0) )
svl.dprior <- function (params, ..., log) {
  f <- sum(dunif(params, min = hyperparams$min, max = hyperparams$max,
                 log = TRUE))
  if (log) f else exp(f)
}


pmcmc1 <-   pmcmc(pomp(svl.filt, dprior = svl.dprior), start = params_test,
                  Nmcmc = 1000, Np = 100, max.fail = Inf,
                  proposal = mvn.diag.rw(c(mu = 0.01, phi = 0.01,sigma_eta = 0.01, rho=0.01)))

continue( pmcmc1 ,Nmcmc=2000,proposal=mvn.rw(covmat( pmcmc1 ))) -> pmcmc1 
plot(  pmcmc1 )
coef(pmcmc1 )
logLik(pmcmc1 )
beep(1)
