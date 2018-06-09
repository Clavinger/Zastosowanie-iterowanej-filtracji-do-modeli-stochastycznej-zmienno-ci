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
library(ggplot2)
library(GGally)
library(gridExtra)
library(reshape2)
library(FKF)


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

start_time <- Sys.time()
z[1]=0
sigma[1]=ni/alpha

for (i in 2:n){
  sigma[i]=sigma[(i-1)]*exp(-lambda*delta)+eta2[i]
}
end_time <- Sys.time()
difftime(end_time,start_time, units = "sec")


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
double lambda,omega;
omega = rnorm(0,sigma_eta);
lambda = rchisq(gamma) ;
H = mu*(1 - phi) + phi*H + omega;
G = lambda/gamma ;
"

####rownanie procesu pomiaru
rproc2.sim <- "
Y_state =  rnorm( mu,sigma_n/lambda );
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
sigma_n=xi;
sigma=xi;
Y_state = rnorm( mu,sigma_n/lambda );
"

###????
bsn_rmeasure <- "
y=rnorm( mu,sigma_n/lambda);
"

####rozk?ad warunkowy zmiennej Y
bsn_dmeasure <- "
lik=dnorm(y,mu,sigma_n/lambda,give_log);
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
bsn.model.C <- pomp(data=data.frame(y=1:n,
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

