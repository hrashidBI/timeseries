

### Assignment 1, Summer 2022
### Md. Humayun Rashid
### ID- 2021 6011
### Dataset Name: UKDriverDeaths


library(dlm) ##install.packages("dlm")
library(TTR) ##install.packages("TTR")
library(tseries) ##install.packages("tseries")
library(lawstat) ##install.packages("lawstat")


par(mfrow=c(1,1))


help(UKDriverDeaths)

UKDriverDeaths

dt=log(UKDriverDeaths)
ar(dt)

par(mfrow=c(1,1))

plot(dt,main="TS Plot")

plot(diff(dt),main="Diff(dt)")

plot(diff(diff(dt)),main="diff(diff(dt))")

dev.copy2eps(file="UKDriverDeaths.eps")

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE)) 


par(mfrow=c(1,1))

plot(dt,main="TS Plot")


lines(SMA(dt, 12),lty="dotdash")


legend("bottomright", col = c("black","black"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("log(UKDriverDeaths)", "12 simple moving average"))


acf(dt,main="ACF")
pacf(dt,main="PACF")

dev.copy2eps(file="acf-pacf.eps")


par(mfrow=c(1,1))
slope0 <- mean(diff(dt))
freq=10

init_par=log(var(dt))

mod1 <- function(u) {
  trend <- dlmModPoly(order=2,
                      #dV = 1e-7,
                      #dV = exp(u[4]),
                      dW = c(exp(u[2]),exp(u[3])),
                      dV = exp(u[1]),
                      #dW = 0,
                      #m0 = c(dt[1])
                      m0 = c(dt[1],slope0),
                      C0 = (1e-7)*diag(2)
  )
  mod <- dlmModSeas(f=freq,
                    #dV = 0,
                    dV = exp(u[4]),
                    dW = c(exp(u[5]),rep(0,freq-2)),
                    #dW = c(0,0)
                    m0 = rep(0, freq - 1),
                    C0 = 1e+07 * diag(nrow = freq - 1)
  )
  return(trend+mod)
}

init_par=log(var(dt))

outMLE1 <- dlmMLE(dt,dlmMLE(dt,c(init_par,init_par,init_par,init_par,init_par), mod1)$par, mod1)
outMLE1$par

## [1]  -6.836786  -4.397117 -25.512625  -6.836786 -25.605640


dlm1 <- mod1(outMLE1$par)

dlm1$FF
dlm1$V

res1 <- dlmFilter(dt,dlm1) #Filtering

ress1<- dlmSmooth(dt,dlm1) #Smoothing

par(mfrow=c(1,1))

plot(dt,xlab="",ylab="value",type="o",col=("green"))
lines(res1$f,col=2,lty="dotdash")  



legend("bottomleft", col = c("blue","blue"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("log(UKDriverDeaths)", "filtering value"))

par(mfrow=c(1,1))
plot(dt,xlab="",ylab="value",type="o",col=("red"))


ts.plot(ress1$s[-1,1],col=c(1),lty=c("solid"))
ts.plot(ress1$s[-1,3],col=c(1),lty=c("solid"))


par(mfrow=c(1,1))
ar(dt, AIC=T)$aic
slope0=mean(diff(dt))


mod2 <- function(u) {
  trend <- dlmModPoly(order=2,
                      dV = 0.000001,
                      dW = exp(u[1:2]),
                      m0 = c(dt[1],slope0),
                      C0 = (1e-7)*diag(2)
  )
  ar<- dlmModARMA(
    ar=c(u[4:5],rep(0,9),u[6]),
    #ar=c(u[4:5]),
    #ma=c(u[1]),
    sigma2 = exp(u[3]),
    #C0 = 10^(-10)*diag(11)
  )
  
  return(trend+ar)
}



outMLE2 <- dlmMLE(dt,dlmMLE(dt,
                              c(init_par,init_par,
                                log(var(ar(dt, AIC=F, order.max=13)$resid[-c(1:13)])),
                                ar(dt, AIC=T, order.max=13)$ar[1:2],
                                ar(dt, AIC=T, order.max=13)$ar[13]
                              ),
                              mod2)$par,mod2)

outMLE2$par

##[1]  -7.339188309 -25.595884248  -5.029440304   0.218193379
##[5]  -0.002124301   0.675004101

dlm2 <- mod2(outMLE2$par)
dlm2$FF

res2 <- dlmFilter(dt,dlm2)

par(mfrow=c(1,1)) 

plot(dt,xlab="",ylab="value",type="o",col=("green"))
lines(res2$f,lty="dotdash")

ress2<- dlmSmooth(dt,dlm2)



par(mfrow=c(1,1), mar=c(1.5,4,0.5,0.5) + 0.1, cex=0.6)


plot.ts(dropFirst(ress2$s)[,1],ylab="level+trend",ann=F,yax.flip=TRUE)

legend("bottomright",col=c(1),lty="solid",legend=c("level+trend"))


plot.ts(dropFirst(ress2$s)[,3],main="",ann=F,yax.flip=TRUE, type="o")

legend("bottomright",col=c(1),lty="solid",legend=c("Seasonal"))


par(mfrow=c(1, 1))

zansa1=residuals(res1)$res


qqnorm(zansa1)
qqline(zansa1)

acf(zansa1)
summary(zansa1)


shapiro.test(zansa1)



ks.test(zansa1, "pnorm", mean=mean(zansa1), sd=sqrt(var(zansa1)))

jarque.bera.test(zansa1)

rjb.test(zansa1)

Box.test(zansa1,type= "Ljung-Box")


Box.test(zansa1,type= "Box-Pierce")



sum(abs(zansa1))/length(zansa1) ##MAD
## [1] 0.7577441


sum(zansa1^2)/length(zansa1) ##MSE
## [1] 0.9531191


sum(abs(zansa1)/dt)/length(zansa1) ##MAPE

##[1] 0.1021601


sum((dt-res1$f)^2)/sum((dt[-1]-dt[-length(dt)])^2) ##U
##[1] 1.45821


par(mfrow=c(1, 1))
zansa2=residuals(res2)$res
qqnorm(zansa2, main="")
qqline(zansa2)

acf(zansa2)
shapiro.test(zansa2)


ks.test(zansa2, "pnorm", mean=mean(zansa2), sd=sqrt(var(zansa2)))

jarque.bera.test(zansa2)

rjb.test(zansa2)


Box.test(zansa2,type= "Ljung-Box")
Box.test(zansa2,type= "Box-Pierce")

mean(abs(res1$f-dt)) ##MAD
## [1] 0.1173216
mean(abs(res2$f-dt))
##[1] 0.07598289

mean((res1$f-dt)^2) ###MSE
## [1] 0.02360453

mean((res2$f-dt)^2) 
## [1] 0.009086485



mean(abs(res1$f-dt)/dt) ## MAPE
# [1] 0.01581873

mean(abs(res2$f-dt)/dt)
#[1] 0.01028233


#U
sqrt(sum((res1$f-dt)[-(1:5)]^2)/
       sum(diff(dt[-(1:4)])^2))
##[1] 1.204158
sqrt(sum((res2$f-dt)[-(1:5)]^2)/
       sum(diff(dt[-(1:4)])^2))
##[1] 0.7438305



par(mfrow=c(1,1),mar=c(2,2,2,0.5))



fore1 <- dlmForecast(res1, nAhead=12)


ciTheory1 <- (outer(sapply(fore1$Q, FUN=function(x) sqrt(diag(x))), qnorm(c(0.25,0.75))) +
                as.vector(t(fore1$f)))

plot(dt,xlim=c(1970,1980),type="o",col="green",main="model1",ylab="")

lines(cbind(ciTheory1,fore1$f[,1])[,1], lty="dotdash")

lines(cbind(ciTheory1,fore1$f[,1])[,2], lty="dotdash")

lines(cbind(ciTheory1,fore1$f[,1])[,3], lty="solid")


legend("bottomleft", col = c("black","black"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("one-step-ahead forecast","theoretical bounds"))

fore2 <- dlmForecast(res2, nAhead=12)


ciTheory2 <- (outer(sapply(fore2$Q, FUN=function(x) sqrt(diag(x))), qnorm(c(0.25,0.75))) +
                as.vector(t(fore2$f)))

plot(dt,xlim=c(1970,1980),type="o",col="red",main="model2")




lines(cbind(ciTheory2,fore2$f[,1])[,1], lty="dotdash")
lines(cbind(ciTheory2,fore2$f[,1])[,2], lty="dotdash")
lines(cbind(ciTheory2,fore2$f[,1])[,3], lty="solid")



legend("bottomleft", col = c("black","black"),
       lty = c("solid", "dotdash"), pch = c(NA,NA), bty = "n",
       legend = c("one-step-ahead forecast","95% Confidence Interval"))
dev.copy2eps(file="UKDriverDeathsforecast.eps",width=10)





