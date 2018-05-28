xx<-seq(-4,4,length.out = 1000)
y<-dnorm(xx)
par(mar=c(2,2,2,2))
plot(xx,y,type='l',xlab="x",ylab="f(x)")
proba=rnorm(100)
y.proba=dnorm(proba)

licznik=0
suma=0
for(i in 1:length(proba)) 
  {
  segments(x0=proba[i], y0=0, x1 = proba[i], y1 = y.proba[i])
  if( (proba[i]>-1) & (proba[i]<1)) {
    licznik=licznik +1
    suma=suma+y.proba[i]
    }
}
licznik/length(proba)
suma/sum(y.proba)
pnorm(1)-pnorm(-1)


abline(v=-1,col='red',lty=2,lwd=2)
abline(v=1,col='red',lty=2,lwd=2)
