
par(mfrow=c(1,2))
xx<-seq(-0.999,.999,length.out = 1000)
p1<-sapply(xx, function(z)   atanh(z))
p2<-sapply(xx, function(z)   z)
plot(xx,p1,type='l',col='red',xlab="x",ylab="y")
lines(xx,p2,col='blue')
legend(x='topleft',legend = c("f(x)=atanh(x)","h(x)=x"),
       lty=c(1,1),col=c('red','blue'))
 
xx<-seq(-4.999,4.999,length.out = 1000)
p1<-sapply(xx, function(z)   tanh(z))
p2<-sapply(xx, function(z)   z)
plot(xx,p1,type='l',col='red',xlab="x",ylab="y")
lines(xx,p2,col='blue')
legend(x='topleft',legend = c("g(x)=tanh(x)","h(x)=x"),
       lty=c(1,1),col=c('red','blue'))
par(mfrow=c(1,1))
log(0.5)/log(.9)
exp(5)
log(0.13)
