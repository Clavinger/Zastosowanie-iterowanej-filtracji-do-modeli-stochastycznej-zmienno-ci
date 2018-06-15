#model  Bandorffa-Nielsena i Sheparda napisany w C

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
#library(ggplot2)
#library(GGally)
#library(gridExtra)
#library(reshape2)
library(FKF)


##############################################################################
##############################################################################
############################################################## symulacja

z=1:n
sigma=1:n


###########################################################################################
###########################################################################################
############################################################## symulacja 2 wersja 
#C:\Users\user\Documents\R\win-library\3.4\pomp\include
#https://github.com/kingaa/pomp/issues/61


eta.vec<-function(lambda,ksi,omega){
  alpha=ksi/omega^2
  ni=ksi^2/omega^2
  
  ci=0
  licznik=0
  eta1=0
  eta2=0
  while(ci<1){
    licznik=licznik+1
    ci<-ci+rexp(1, rate = lambda*ni)
    ri=runif(1)
    if(ci>1) break
    eta1=eta1+log(1/ci)*exp(lambda*ri)
    eta2=eta2+log(1/ci)
  }
  if(licznik==1){
    eta1=0
    eta2=0
  } else{
    eta1=1/alpha*exp(-lambda)*eta1
    eta2=1/alpha*eta2
  }
  return(c(eta1,eta2))
}



start_time <- Sys.time()
z[1]=0
sigma[1]=ksi

for (i in 2:n){
  eta=eta.vec(lambda,ksi,omega)
  sigma[i]=sigma[(i-1)]*exp(-lambda*delta)+eta[1]
  z[i]=z[i-1]+eta[2]
}
end_time <- Sys.time()
difftime(end_time,start_time, units = "sec")

plot(z,type='l')
plot(sigma,type='l')


set.seed(123) 
u<-rnorm(n)
y=1:n
y2=1:n
lr=1:n
y[1]=0
for(i in 2:n){
  y[i]=y[i-1]+mu*delta + beta*sigma[i-1] + sqrt(sigma[i-1])*u[i-1]*sqrt(delta)
}
plot(y,type='l')

#zwroty
x<-1:(n-1) 
for (i in 2:n) x[i]=y[i]-y[i-1]


#kwadraty zwrot?w
x2<-x^2

plot(sigma,type='l')
plot(x,type='l')
plot(x2,type='l')

###########################################################################################
###########################################################################################
########################################################################Iterowana filtracja
?pomp

#bsn specyfikacja do symulacji

bsn_statenames <- c("sigma","sigma_n","Y_state")
bsn_rp_names <- c("mu","lambda","xi","omega")
bsn_paramnames <- bsn_rp_names
bsn_covarnames <- "covaryt"


 

rproc1 <- "
double *eta=etafun(lambda,xi,omega);
sigma_n=(1-exp(-lambda))*sigma + eta[1]-eta[0];
sigma=exp(-lambda)*sigma + eta[0];
free(eta);
"

####rownanie procesu pomiaru
rproc2.sim <- "
Y_state =  rnorm( mu,sqrt(sigma_n/lambda) );
"
###do wypelniania danych
rproc2.filt <- "
Y_state = covaryt;
"

###symulacja modelu SVL
bsn_rproc.sim <- paste(rproc1,rproc2.sim)

####filtr czasteczkowy 
bsn_rproc.filt <- paste(rproc1,rproc2.filt)


######inicalizacja
bsn_initializer <- "
  double alpha=xi/(omega*omega);
  double ni=xi*xi/(omega*omega);
sigma=rgamma(ni,1/alpha);
sigma_n=sigma*lambda;
Y_state = rnorm( mu,sqrt(sigma_n/lambda) );
"

###????
bsn_rmeasure <- "
y=rnorm( mu,sqrt(sigma_n/lambda) );
"

####rozk?ad warunkowy zmiennej Y
bsn_dmeasure <- "
lik=dnorm(y,mu,sqrt(sigma_n/lambda) ,give_log);
"


####przeskalowanie parametr?w 
bsn_toEstimationScale <- "
Tlambda = log(lambda);
Txi = log(xi);
Tomega= log(omega);
"

bsn_fromEstimationScale <- "
Tlambda = exp(lambda);
Txi = exp(xi);
Tomega= exp(omega);
"



####wypelnianie modelu danymi
bsn.model.C <- pomp(data=data.frame(y=x,
                                  time=1:n),
                  statenames=bsn_statenames,
                  paramnames=bsn_paramnames,
                  covarnames=bsn_covarnames,
                  times="time",
                  t0=0,
                  covar=data.frame(covaryt=c(0,x),
                                   time=0:n),
                  tcovar="time",
                  rmeasure=Csnippet(bsn_rmeasure),
                  dmeasure=Csnippet(bsn_dmeasure),
                  rprocess=discrete.time.sim(step.fun=Csnippet(bsn_rproc.sim),delta.t=1),
                  initializer=Csnippet(bsn_initializer),
                  toEstimationScale=Csnippet(bsn_toEstimationScale), 
                  fromEstimationScale=Csnippet(bsn_fromEstimationScale)
)

plot(bsn.model.C)


params_test=c(
  mu=mu,
  lambda=lambda,
  xi=ksi,
  omega=omega
)



#wykres filtracji dla prawdziwych parametr?w
start_time <- Sys.time()
pf1 <- pfilter(bsn.model.C,params=params_test,
               Np=1000,filter.traj=T)
end_time <- Sys.time()
difftime(end_time,start_time, units = "sec")

plot(pf1)



symulacja <-simulate(bsn.model.C ,nsim=1,seed=123,params=params_test,
                                states = TRUE, obs = TRUE)


plot(symulacja $ states[1,1,],type='l')

