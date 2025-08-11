source('fnexpectile.R')
source('fn.R')
library(fGarch)
library(quantreg)
library(expectreg)
library(fields)

b0 = 1
b1 = 0.9
b2 = 1
omega0 = 2*pi*0.09
omega1 = 2*pi*0.12
omegac = 2*pi*0.25
r = 0.6   
n = 200
tau = 0.9 
ncores = 2

kep_smo7 = 0
kep_smo9 = 0
qp_smo = 0
pg_smo = 0

kep_res7 = 0
kep_res9 = 0
qp_res = 0
pg_res = 0
k = 2000

set.seed(1)
for (i in 1:k){
  xt =  as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=n))
  at = b0 + b1*cos(omega0*(1:n)) + b2*sin(omega1*(1:n))
  yt = xt*at
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  
  ###xt
  kep_no7 = kep(xt, freq, tau = 0.7, k=1.5 ,n.cores=1)
  kep_smo7 = kep_smo7 + smooth.spline(kep_no7, cv=T)$y/k
  kep_no9 = kep(xt, freq, tau = 0.9, k=1.5 ,n.cores=1)
  kep_smo9 = kep_smo9 + smooth.spline(kep_no9, cv=T)$y/k

  pg_no = spec.pgram(xt, plot=F)$spec
  pg_smo = pg_smo + smooth.spline(pg_no, cv=T)$y/k 
  
  qp_no = lap.spec.new2(xt,freq,0.5,intercept=T,type=1,n.cores=1)
  qp_smo = qp_smo + smooth.spline(qp_no, cv=T)$y/k
  
  ####yt
  kep7 = kep(yt, freq, tau = 0.7, k=1.5 ,n.cores=1)
  kep_res7 = kep_res7 + kep7/k
  kep9 = kep(yt, freq, tau = 0.9, k=1.5 ,n.cores=1)
  kep_res9 = kep_res9 + kep9/k
  
  pg = spec.pgram(yt, plot=F)$spec
  pg_res = pg_res + pg/k
  
  qp = lap.spec.new2(yt,freq,0.5,intercept=T,type=1,n.cores=1)
  qp_res = qp_res + qp/k

  print(i)
}

kep_smo7 = kep_smo7/sum(kep_smo7)
kep_smo9 = kep_smo9/sum(kep_smo9)
qp_smo = qp_smo/sum(qp_smo)
pg_smo = pg_smo/sum(pg_smo)

kep_res7 = kep_res7/sum(kep_res7)
kep_res9 = kep_res9/sum(kep_res9)
qp_res = qp_res/sum(qp_res)
pg_res = pg_res/sum(pg_res)

ts.plot(kep_smo7)
ts.plot(kep_smo9)
ts.plot(qp_smo)
ts.plot(pg_smo)

ts.plot(kep_res7)
ts.plot(kep_res9)
ts.plot(qp_res)
ts.plot(pg_res)

#image
tau=seq(from=5,to=95,2)/100
omegac = 2*pi*0.25
k = 2000
set.seed(1)
kep_ima11 = 0
# kep_ima15 = 0
kep_ima19 = 0
for (i in 1:k){
  xt =as.numeric(arima.sim(list(order=c(2,0,0),ar=c(2*r*cos(omegac),-r^2),sd=1),n=n))
  at = b0 + b1*cos(omega0*(1:n)) + b2*sin(omega1*(1:n))
  yt = xt*at
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  
  # kep_no_ima11 = kep(xt, freq, tau, k=1.1 ,n.cores=2)
  # kep_ima11 = kep_ima11 + kep_no_ima11/k
  kep_no_ima15 = kep(xt, freq, tau, k=1.5 ,n.cores=2)
  kep_ima15 = kep_ima15 + kep_no_ima15/k
  # kep_no_ima19 = kep(xt, freq, tau, k=1.9 ,n.cores=2)
  # kep_ima19 = kep_ima19 + kep_no_ima19/k
  print(i)
}

for(i in 1:length(tau)){
  # kep_ima11[,i] = smooth.spline(kep_ima11[,i],cv=T)$y 
  kep_ima15[,i] = smooth.spline(kep_ima15[,i],cv=T)$y
  # kep_ima19[,i] = smooth.spline(kep_ima19[,i],cv=T)$y 
  
}
######plot
par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
set.panel(1,3)

plot((0:99)/200, pg_smo/sum(pg_smo) ,type='l',xaxt="n",yaxt="n",xlab='',ylab='',col ='black',main='(a)',ylim=c(0.001,0.023))
lines((0:98)/198, qp_smo/sum(qp_smo), col='blue')
lines((0:98)/198, kep_smo7/sum(kep_smo7), col='red')
lines((0:98)/198, kep_smo9/sum(kep_smo9), col='green')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq",ylab=expression('Periodograms'), line=2, cex.lab=1.2)
grid()
legend("topright", c('PG','LP',expression('KEP (' ~ alpha == 0.7 ~')'),expression('KEP (' ~ alpha == 0.9 ~')')),inset = 0.01,
       lty=c(1, 1),  col=c( "black", "blue",'red','green'),
       seg.len=2,bty = "n")

# plot((0:98)/198, kep_smo7/sum(kep_smo7) ,type='l',xaxt="n",yaxt="n",xlab='',ylab='',col ='black')
# lines((0:98)/198, kep_smo9/sum(kep_smo9), col='blue')
# axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
# title(xlab="Freq",ylab=expression('Periodograms'), line=2, cex.lab=1.2)
# grid()
# legend("bottom", c(expression('KEP (' ~ alpha == 0.7 ~')'),expression('KEP (' ~ alpha == 0.9 ~')')),inset = 0.01,
#        lty=c(1, 1),  col=c( "black", "blue"),
#        seg.len=2,bty = "n")

image((0:98)/198, y = seq(from=5,to=95,2)/100, z=spec.normalize( kep_ima11),xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0,0.5),col = tim.colors(), main='(b)')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq",ylab="Expectile ( k = 1.5 )", line=2, cex.lab=1.2)

## out

plot((0:98)/198, kep_res7/sum(kep_res7) ,type='l',xaxt="n",yaxt="n",xlab='',ylab='',col ='red',main='(c)',ylim=c(0.003,0.055))
lines((0:98)/198, kep_res9/sum(kep_res9), col='green')
lines((0:98)/198, pg_res[2:100]/sum(pg_res[1:99]), col='black')
lines((0:98)/198, qp_res/sum(qp_res), col='blue')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq",ylab=expression('Periodograms'), line=2, cex.lab=1.2)
grid()
legend("topright", c('PG','LP',expression('KEP (' ~ alpha == 0.7 ~')'),expression('KEP (' ~ alpha == 0.9 ~')')),inset = 0.01,
       lty=c(1, 1),  col=c( "black", "blue",'red','green'),
       seg.len=2,bty = "n")



#####MIXTURE
s = 1000
n = 200
tau=seq(from=5,to=95,2)/100
mix_kep = matrix(0,100-1,46)
mix_pg = rep(0,99)
mix_qp = rep(0,99)
set.seed(1)
for(i in 1:s){   ##GARCH
  data1 = mymix(200)
  n <- length(data1)
  freq <- c(0:(n-1))/n
  freq <- freq[freq>0 & freq <0.5]
  q1 <- kep(data1, freq, tau,k =1.5, n.cores=1)
  q1 = spec.normalize(q1)
  mix_kep = mix_kep + q1/s
  
  mix_pg = mix_pg + smooth.spline(pgram(data1))$y/s
  mix_qp = mix_qp + smooth.spline(lap.spec.new2(data1,freq,0.5,intercept=T,type=1,n.cores=1),cv=T)$y/s
  print(i)}

for(i in 1:46){
  mix_kep[,i] = smooth.spline(mix_kep[,i], cv=T)$y
}

set.panel(1,2)
par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
image((0:98)/198, y = seq(from=5,to=95,2)/100, z=spec.normalize(mix_kep),xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0,0.5),col = tim.colors(), main='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq",ylab="Expectile ( k = 1.5 )", line=2, cex.lab=1.2) 

plot((1:99)/198, mix_pg/sum(mix_pg),type='l',xlab='',ylab='',col='blue',xaxt='n',yaxt='n',ylim=c(0,0.041))
lines((1:99)/198, mix_kep[,4]/sum(mix_kep[,4]), type='l',col='red')
lines((1:99)/198, mix_kep[,44]/sum(mix_kep[,44]), type='l',col='black')
title(xlab="Freq",ylab="Periodograms", line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
legend(x = 'topright', 
       legend = c('PG',expression('KEP (' ~ alpha == 0.1 ~ ", " ~ k == 1.5~')'),expression('KEP (' ~ alpha == 0.9 ~ ", " ~ k == 1.5~')')),
       lty = 1,
       col = c( 'blue','red','black') ,inset = 0.01,bty = "n")
