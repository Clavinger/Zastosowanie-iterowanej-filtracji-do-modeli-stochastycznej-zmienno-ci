dane2<-read.csv('wig_d.csv', header = T)
dane3<-xts(dane2[, -1], order.by=as.POSIXct(dane2$Data))
dane<-xts(dane3$Zamkniecie)


zwroty=x


par(mfrow=c(2,1))
#wykres WIG20 dzienne ceny zamkniecia
plot(dane,major.ticks="years", main="WIG dzienne ceny zamkniecia",
     minor.ticks=NULL,grid.ticks.on = "years")

#wykres WIG20 dzienne logarytmiczne zwroty
plot(zwroty,major.ticks="years",main="WIG dzienne logarytmiczne zwroty",
     minor.ticks=NULL,grid.ticks.on = "years")
par(mfrow=c(1,1))



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
H = rnorm(mu_h,sigma_eta/sqrt((1-phi*phi))) ;
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



##--------Likelihood maximization using randomized starting values--------


bsv_box <- rbind(
  sigma_eta=c(0.001,1),
  phi    = c(0.9,1),
  mu_h = c(-1,1)
)

###trzy szybkosci filtru: 1 -szybki, 2 -sredni, 3 - wolny
run_level <- 2

#liczba czasteczek
bsv_Np <-          c(100,1e3,2e3)
bsv_Nmif <-        c(10, 50,200)
bsv_Nreps_eval <-  c(4,  10,  20)
bsv_Nreps_local <- c(4, 10, 20)
bsv_Nreps_global <-c(4, 10, 100)

bsvlist<-list(bsv_Np ,bsv_Nmif,bsv_Nreps_eval,
              bsv_Nreps_local,bsv_Nreps_global )

#parametry do metody mif2
bsv_rw.sd_rp <- 0.02
bsv_rw.sd_ivp <- 0.1
bsv_cooling.fraction.50 <- 0.5




 

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
                                   mu_h      = bsv_rw.sd_rp,
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
beep(1)

plot(if.bsv  )

#save(if.bsv.box,L.bsv.box,H.bsv.box, file="bsv_box_eval.rda")


r.box <- data.frame(logLik= L.if.bsv [,1],logLik_se= L.if.bsv  [,2],t(sapply(if.bsv,coef)))
summary(r.box$logLik,digits=5)
r.box [which.max(r.box $logLik),]

#rysunki procesu zmiennosci

params_nowe2<- c(
  mu_h        = r.box [which.max(r.box $logLik),'mu_h'],    
  phi         = r.box [which.max(r.box $logLik),'phi'],   
  sigma_eta   = r.box [which.max(r.box $logLik),'sigma_eta'] 
)

pf3 <- pfilter(bsv.filt,params=params_nowe2,
               Np=bsvlist[[1]][run_level],filter.traj=T)

par(mfrow=c(2,1))
plot(zwroty,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
sigma=as.xts(exp(pf2@filter.traj[1,1,2:(dim(pf2@filter.traj)[3])]/2),order.by=index(zwroty))
plot(sigma,minor.ticks=NULL,grid.ticks.on = "years",major.ticks="years")
par(mfrow=c(1,1))

particle.bvs<- exp(pf3@filter.traj[1,1,])
plot(particle.bns,type='l',col='royalblue')
lines(kalman.bns, col='red')
lines(particle.bvs,col='green')


wykres2.matrix<-rbind(particle.bvs[2:(n+1)],kalman.bns,1:n)
row.names(wykres2.matrix)<-c("bsv","Kalman","time")
wykres2<-as.data.frame(t(wykres2.matrix))
wykres2.long<-melt(wykres2,id="time")
names(wykres2.long)<-c("time","Rodzaj","value")

g1<-ggplot(data = wykres2.long, aes(x = time,y=value,colour=Rodzaj,linetype=Rodzaj) ) +ggtitle('(c)')+
  geom_line(size = .8)+labs(y="zmiennność")+theme_bw()+scale_colour_manual(values=c("royalblue","orange"))+
  scale_linetype_manual(values=c("solid", "dashed"))

wykres2.matrix<-rbind(particle.bvs[2:(n+1)],particle.bns[2:(n+1)],1:n)
row.names(wykres2.matrix)<-c("bsv","particle","time")
wykres2<-as.data.frame(t(wykres2.matrix))
wykres2.long<-melt(wykres2,id="time")
names(wykres2.long)<-c("time","Rodzaj","value")

g2<-ggplot(data = wykres2.long, aes(x = time,y=value,colour=Rodzaj,linetype=Rodzaj) ) +ggtitle('(b)')+
  geom_line(size = .8)+labs(y="zmiennność")+theme_bw()+scale_colour_manual(values=c("royalblue","darkgreen"))+
  scale_linetype_manual(values=c("solid", "dotdash"))

wykres2.matrix<-rbind(kalman.bns,particle.bns[2:(n+1)],1:n)
row.names(wykres2.matrix)<-c("Kalman","particle","time")
wykres2<-as.data.frame(t(wykres2.matrix))
wykres2.long<-melt(wykres2,id="time")
names(wykres2.long)<-c("time","Rodzaj","value")

g3<-ggplot(data = wykres2.long, aes(x = time,y=value,colour=Rodzaj,linetype=Rodzaj) ) +ggtitle('(a)')+
  geom_line(size = .8)+labs(y="zmiennność")+theme_bw()+scale_colour_manual(values=c("orange","darkgreen"))+
  scale_linetype_manual(values=c( "dashed","dotdash"))

grid.arrange( g3,g2,g1,ncol=1)

save(x,x2,params_new,params_test,particle.bns,
     kalman.bns,kalman.kalman.bns,file="wyniki_bns_empiryczne")
load(file="wyniki_bns_empiryczne")


 accuracy(f=particle.bvs,x=particle.bns)
 accuracy(f=particle.bvs,x=kalman.bns)
 accuracy(f=particle.bns,x=kalman.bns)