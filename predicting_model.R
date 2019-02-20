library(readxl)
batman = read_xlsx("batman.xlsx", col_names = FALSE)
names(batman) = c("Receipts", "Date", "Day")
batman = ts(batman[,1])
batmanT1 = log(batman)
batmanT2 = sqrt(batman)
batmanT3 = batman^{-1/2}
par(mfrow = c(2,2))
plot.ts(batman, ylab = 'batman', main = 'Plot of batman')
plot.ts(batmanT1, ylab = 'ln(batman)', main = 'Plot of ln(batman)')
plot.ts(batmanT2, ylab = 'sqrt(batman)', main = 'Plot of sqrt(batman)')
plot.ts(batmanT3, ylab = '1/sqrt(batman)', main = 'Plot of 1/sqrt(batman)')

par(mfrow = c(2,2))
hist(batman)
hist(batmanT1)
hist(batmanT2)
hist(batmanT3)

time= 1:119
Model0 = lm(batman~time)
plot(Model0$residuals)
Model1 = lm(batmanT1~time)
plot(Model1$residuals)
Model2 = lm(batmanT2~time)
plot(Model2$residuals)
Model3 = lm(batmanT3~time)
plot(Model3$residuals)

t = 1:119
model02 = lm(batmanT1~poly(t,2))
summary(model02)
par(mfrow = c(2,2))
plot(model02$fitted.values)
plot(model02$residuals)
acf(model02$residuals)
pacf(model02$residuals)
AIC(model02)

model03 = lm(batmanT1 ~poly(t,3))
summary(model03)
par(mfrow = c(2,2))
plot(model03$fitted.values)
plot(model03$residuals)
acf(model03$residuals)
pacf(model03$residuals)
AIC(model03)

model04 = lm(batmanT1~poly(t,4))
summary(model04)
par(mfrow = c(2,2))
plot(model04$fitted.values)
plot(model04$residuals)
acf(model04$residuals)
pacf(model04$residuals)
AIC(model04)

trndseas=function(y,seas=1,lam=1,degtrnd=0){
  
  # requires the R-package 'pracma'
  
  # fits  a trend plus seasonal for the "best" Box-Cox 
  # transformation.
  
  # input: y, observed series; seas, seasons
  
  # input: lam, the grid of Box-Cox transformations (lambda values)
  
  # input: degtrnd, degree of the polynomial trend, if
  # degtrnd=0, then the fitted trend is constant.
  
  # output:  coef, regression coefficients - the
  # first degtrnd+1 values for the trend part and the
  # rest associated with the seasonals
  
  # output: fit, fitted y-values; res, residuals,
  
  # output: trend, fitted trend; season, fitted seasonals
  
  # output: rsq, r-square values for different lambda in the
  
  # output: lamopt, the value of lambda (among those supplied 
  # in the vector lam) at which r-square is maximum.
  
  m=length(lam)
  n=length(y)
  
  # Part of design matrix for estimating trend
  if(degtrnd>0) {
    tm=seq(1/n,1,by=1/n) #use normalized time
    x1=poly(tm,degree=degtrnd,raw=TRUE) #generate the matrix of x
    x1=cbind(rep(1,n),x1)
  } else {
    x1=as.matrix(rep(1,n),ncol=1) #x matrix with only intercept
  }
  
  # Part of design matrix for estimating seasonality
  x2=NULL
  if(seas>1){
    sn=rep(1:seas,length.out=n)
    x2=factor(sn,levels=unique(sn),ordered=TRUE)
    x2=model.matrix(~x2-1) #matrix without the intercept
    m2=ncol(x2)
    m21=m2-1
    x2=x2[,1:m21]-matrix(rep(x2[,m2],m21),ncol=m21,nrow=nrow(x2),byrow=FALSE) #include the negative ones
  }
  
  x=cbind(x1,x2)  # design matrix
  
  xx=t(x)%*%x
  rsq=rep(1,m) #empty vector of length m for computation of r squared of each transformation fit
  m1=ncol(x1)     #degtrnd+1
  m11=m1+1
  mx=ncol(x)      # degtrnd+1+seas-1
  
  for(i in 1:m) { #m is the length of lambda, do Boxcox transformation
    if (lam[i]==0) {
      yt=log(y)
    } else {
      yt=y^lam[i]
    }
    xy=t(x)%*%yt
    coef=solve(xx,xy)
    fit=x%*%coef #this is the trend and seasonal fit
    res=yt-fit
    ssto=(n-1)*var(yt)
    sse=t(res)%*%res
    rsq[i]=1-sse/ssto
  }
  
  ii=which.max(rsq) 
  lamopt=lam[ii]   #choose lambda optimal according to r squared
  if (lamopt==0) {
    yt=log(y)
  } else {
    yt=y^lamopt
  } #optimal transformation done ; yt
  xy=t(x)%*%yt
  coef=solve(xx,xy)
  fit=x%*%coef
  trnd=x1%*%coef[1:m1] #extract the trend part
  season=NULL
  if(seas>1){
    season=c(coef[m11:mx],-sum(coef[m11:mx]))
    #season = x2%*%coef[m11:mx]
  }
  res=yt-fit
  
  result=list(coef=coef,fit=fit,trend=trnd,res=res,season=season,rsq=rsq,lamopt=lamopt)
  return(result)
}
lam = seq(-1,1,by=0.05)
ff = trndseas(batmanT1,seas = 5,lam = 1,degtrnd = 2)
rsq = ff$rsq
rsq
attributes(ff)
m.fit = ff$trend
ff$season
n = length(batmanT1)
s.fit = rep(ff$season,length.out=n)
smooth.fit = ff$fit
par(mfrow=c(2,2))
plot.ts(batmanT1)
plot.ts(m.fit, main='Estimated Trend')
plot.ts(s.fit,main='Estimated Seasonal Component')
plot.ts(batmanT1,main='Estimated Smooth Part')
points(smooth.fit,type='l',col='red')
months = 1:5
plot.ts(batmanT1)
plot(months, ff$season, type='l', ylab = 'Seasonal', main = 'Seasonals for log(batman)')
plot.ts(ff$res,type = 'l', main = "Estimated rough")

par(mfrow=c(1,1))
x = batmanT1-m.fit-s.fit
acf(x)
pacf(x)
hist(x)
qqnorm(x)
qqline(x)

fitAR0 = arima(x,order=c(0,0,0))
fitAR1 = arima(x,order=c(1,0,0))
fitAR2 = arima(x,order=c(2,0,0))
fitAR3 = arima(x,order=c(3,0,0))
fitAR4 = arima(x,order=c(4,0,0))
fitAR5 = arima(x,order=c(5,0,0))
fitAR6 = arima(x,order=c(6,0,0))
aicc = function(model){
  n = model$nobs
  p = length(model$coef)
  aicc = model$aic + 2*p*(p+1)/(n-p-1)
  return(aicc)
}
aiccAR0 = aicc(fitAR0)
aiccAR1 = aicc(fitAR1)
aiccAR2 = aicc(fitAR2)
aiccAR3 = aicc(fitAR3)
aiccAR4 = aicc(fitAR4)
aiccAR5 = aicc(fitAR5)
aiccAR6 = aicc(fitAR6)
AICC = c(aiccAR0,aiccAR1,aiccAR2,aiccAR3,aiccAR4,aiccAR5,aiccAR6)
AICC
par(mfrow=c(1,1))
plot.ts(fitAR6$residuals)
acf(fitAR6$residuals)
pacf(fitAR6$residuals)
Box.test(fitAR6$residuals,lag=10,'Ljung-Box')

t2 = 1:112
Yt = batmanT1[1:112]
modelT2 = lm(Yt~poly(t2,2))
summary(modelT2)
par(mfrow = c(2,2))
plot(modelT2$fitted.values)
plot(modelT2$residuals)
acf(modelT2$residuals)
pacf(modelT2$residuals)
AIC(modelT2)

lam = seq(-1,1,by=0.05)
ff1 = trndseas(Yt,seas = 5,lam = 1,degtrnd = 2)
rsq1 = ff1$rsq
rsq1
attributes(ff1)
m.fit1 = ff1$trend
season1 = ff1$season
season1
n1 = length(Yt)
s.fit1 = rep(season1,length.out=n1)
smooth.fit1 = ff1$fit
par(mfrow=c(2,2))
plot.ts(Yt)
plot.ts(m.fit1, main='Estimated Trend on first 112 observations')
plot.ts(s.fit1,main='Estimated Seasonal Component on first 112 observations')
plot.ts(Yt,main='Estimated Smooth Part on first 112 observations')
points(smooth.fit1,type='l',col='red')

plot.ts(Yt)
plot.ts(m.fit1, main='Estimated Trend on first 112 observations')
month = 1:5
plot(month, season1, type='l', ylab = 'Seasonal', main = 'Seasonals on first 112 observation')
plot.ts(ff1$res,type = 'l', main = "Estimated rough on first 112 observation")

par(mfrow=c(2,2))
x1 = Yt-m.fit1-s.fit1
acf(x1)
pacf(x1)
hist(x1)
qqnorm(x1)
qqline(x1)
fitTAR6 = arima(x1,order=c(6,0,0))
aiccTAR6 = fitTAR6$aic + 2*7*(7+1)/(n1-7-1)
aiccTAR6
par(mfrow=c(1,1))
plot.ts(fitTAR6$residuals)
acf(fitTAR6$residuals)
pacf(fitTAR6$residuals)
Box.test(fitTAR6$residuals,lag=10,'Ljung-Box')

library(Hmisc)
trend7Days = approxExtrap(m.fit1[1:112],m.fit1[1:105],xout=m.fit1[106:112], method="linear")[1]
trend7Days
ff2 = trndseas(Yt[106:112],seas = 5,lam = 1,degtrnd = 2)
ff2$season
h = 7
deg = 2
coef = ff1$coef[1:(deg+1)]
time1 = (n1+(1:h))/n1
predmat = matrix(rep(time1,deg)^rep(1:deg,each=h),nrow=h,byrow=FALSE)
predmat = cbind(rep(1,h),predmat)
predmat
m.fc = predmat %*% coef
s.fc = rep(ff1$season,length.out=n1+h)
s.fc = s.fc[-(1:n1)]
s.fc
fcast = predict(fitTAR6,n.ahead=h)
x.fc = fcast$pred
x.fc
y.fc = m.fc + s.fc + x.fc
y.fc
plot.ts(Yt,xlim=c(0,n1+h))
points(x=n1+1:h, y=y.fc, col='purple',type='b',pch=19)
