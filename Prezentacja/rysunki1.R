library(pomp)
library(latticeExtra)
#rm(list = ls())
setwd("C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Dane")


#dlugosc szeregu czasowego
n=1000
#parametry
mu = -0.5
phi = 0.98
sigma = 0.2

load(file ='dane do przykladu2')
plot(sim1.sim)

zwroty= sim1.sim@data

####nazwy
bsv_statenames <- c("H","Y_state")
bsv_rp_names <- c("mu_h","phi","sigma_eta")
#bsv_ivp_names <- c("H_0")
#bsv_paramnames <- c(bsv_rp_names,bsv_ivp_names)
bsv_paramnames <- c(bsv_rp_names)
bsv_covarnames <- "covaryt"



rproc1 <- "
double omega;
omega = rnorm(0,sigma_eta );
H = mu_h*(1 - phi) + phi*H + omega;
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
H = rnorm(mu_h  ,sigma_eta/sqrt((1-phi*phi))) ;
Y_state = rnorm( 0,exp(H/2) );
"
###????
bsv_rmeasure <- "
y=Y_state;
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


####wypelnianie modelu danymi
bsv.filt <- pomp(data=data.frame(y=as.vector(zwroty),
                                 time=1:length(zwroty)),
                 statenames=bsv_statenames,
                 paramnames=bsv_paramnames,
                 covarnames=bsv_covarnames,
                 times="time",
                 t0=0,
                 covar=data.frame(covaryt=c(0,as.vector(zwroty)),
                                  time=0:length(zwroty)),
                 tcovar="time",
                 rmeasure=Csnippet(bsv_rmeasure),
                 dmeasure=Csnippet(bsv_dmeasure),
                 rprocess=discrete.time.sim(step.fun=Csnippet(bsv_rproc.filt),delta.t=1),
                 initializer=Csnippet(bsv_initializer),
                 toEstimationScale=Csnippet(bsv_toEstimationScale), 
                 fromEstimationScale=Csnippet(bsv_fromEstimationScale)
)


plot(bsv.filt)

##testowa wersja parametr?w
params_test <- c(
  mu_h = mu,       
  phi = phi,     
  sigma_eta =sigma
)


#wykres filtracji dla prawdziwych parametr?w
pf1 <- pfilter(bsv.filt,params=params_test,
               Np=10,save.states=TRUE)




plot(pf1)

dane=matrix(0,nrow=10,ncol=20)
srednia=rep(0,20)
for(i in 1:20){

srednia[i]=mean(dane[,i])
}
plot(1:20,seq(-2,1,length.out = 20),pch=1,col="white",xlab="t",ylab="h")
for(i in 1:10) {
  lines(1:20,dane[i, ],col='red')
points(1:20,dane[i, ],col='red')
}
lines(sim1.sim@states[1,1:20],col='blue',lwd=3)
points(sim1.sim@states[1,1:20],col='blue',pch=1)
lines(1:20,srednia,col="darkgreen",lwd=3)
points(1:20,srednia,col='darkgreen',pch=1)

legend(x="topleft",legend=c("true","mean","particles"),
       col=c("blue","darkgreen","red"),lty=1,lwd=2)



#####################################################


dane=matrix(0,nrow=10,ncol=10)
srednia=rep(0,10)
for(i in 1:10){
  dane[,i]=pf1@saved.states[[i]][1, ]
  srednia[i]=mean(dane[,i])
}
plot(1:20,seq(-2,1,length.out = 20),pch=1,col="white",xlab="t",ylab="h")
for(i in 1:10) {
  lines(1:10,dane[i, ],col='red')
  points(1:10,dane[i, ],col='red')
}
lines(sim1.sim@states[1,1:10],col='blue',lwd=3)
points(sim1.sim@states[1,1:10],col='blue',pch=1)
lines(1:10,srednia,col="darkgreen",lwd=3)
points(1:10,srednia,col='darkgreen',pch=1)
legend(x="topleft",legend=c("true","mean","particles"),
       col=c("blue","darkgreen","red"),lty=1,lwd=2)


########################################################
########################################################



#wykres filtracji dla prawdziwych parametr?w
set.seed(1)
pf1 <- pfilter(bsv.filt,params=params_test,
               Np=10,save.states=TRUE)
set.seed(1)
pf2 <- pfilter(bsv.filt,params=params_test,
               Np=100,save.states=TRUE)
set.seed(1)
pf3 <- pfilter(bsv.filt,params=params_test,
               Np=1000,save.states=TRUE)


vals <-   data.frame(
p1<-pf1@saved.states[[20]][1, ],
p2<-pf2@saved.states[[20]][1, ],
)

ecdfplot(~ p1 + p2, data=vals, auto.key=list(space='right'))




ecdf1 <- ecdf(pf1@saved.states[[20]][1, ])
ecdf2 <- ecdf(pf2@saved.states[[20]][1, ])
ecdf3 <- ecdf(pf3@saved.states[[20]][1, ])
plot(ecdf3, verticals=TRUE, do.points=FALSE, col='blue',,main="",lwd=2)
plot(ecdf2, verticals=TRUE, do.points=FALSE, add=TRUE, col='orange',lwd=2)
plot(ecdf1, verticals=TRUE, do.points=FALSE, add=TRUE,col='red',lwd=2)
legend(x="topleft",legend=c("J=10","J=100","J=1000"),
       col=c("red","orange","blue"),lty=1,lwd=2)
