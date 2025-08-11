source('fnexpectile.R')
source('fn.R')
library(fGarch)
library(quantreg)
library(expectreg)
library(fields)
tau=seq(from=5,to=95,2)/100


data = read.csv('sp.csv')
data = data$Close[14540:22103]
data = diff(log(data))

n <- length(data)
freq <- c(0:(n - 1)) / n 
freq <- freq[freq>0 & freq <0.025]
ep = expectile_peri_ls(data, freq, tau=(5:95)/100, n.cores=2)
ep = spec.normalize(ep) 
kep11 = kep(data, freq, tau, k=1.1 ,n.cores=4)
kep11 = spec.normalize(kep11)
kep15 = kep(data, freq, tau, k=1.5 ,n.cores=4)
kep15 = spec.normalize(kep15)
kep19 = kep(data, freq, tau, k=1.9 ,n.cores=4)
kep19 = spec.normalize(kep19)
pg = (1/length(data))*(abs(fft(data)[2:190]))^2


# pgsmo = smooth.spline(pg,spar=0.8)$y
# pgsmo = pgsmo/sum(pgsmo)
# ep = spec.normalize(ep)
# kep2 = spec.normalize(kep2)
# kep15 = spec.normalize(kep15)

# epsmo = ep
# kep2smo = ep
# kep15smo = ep
# for(i in 1:91){
#   epsmo[,i] = smooth.spline(ep[,i],spar=0.8)$y
#   kep2smo[,i] = smooth.spline(kep2[,i],spar=0.8)$y
#   kep15smo[,i] = smooth.spline(kep15[,i],spar=0.8)$y
# }


#1
layout_matrix <- matrix(c(1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11), nrow = 4, byrow = TRUE)
layout(mat = layout_matrix, widths = c(1, 1, 1))
par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
plot((1:7563)/252+1986, data ,type='l',xaxt="n",yaxt="n",xlab='',ylab='',ylim=c(-0.1,0.113), xlim = c(1986,2016))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Year",'Daily log return', line=2, cex.lab=1.2)
grid()
lines((1:7563)/252+1986, 0.014+0.007*sin((1:7563)/410 -0.5),col='blue')
lines((1:7563)/252+1986, -0.014-0.007*sin((1:7563)/410 -0.5),col='blue')
legend("top", c('Daily log returns','10-year cycle'),inset = 0.01,
       lty=c(1, 1),  col=c( "black", "blue"),
       seg.len=2)


plot((0:188)/7560,pg/sum(pg) ,type='l',xaxt="n",yaxt="n",xlab='',ylab='',ylim=c(0,0.16))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per day)",ylab="PG", line=2, cex.lab=1.2)
grid()
text(16/3000, 0.5 ,'No 10-year cycle')

# text(8/3000, 0.12 ,'10-year cycle')

#2
plot((0:188)/7560, (kep11[,44]/sum(kep11[,44]))[1:189] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='',ylim=c(0,0.16))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('KEP (' ~ alpha == 0.9 ~ ", " ~ k == 1.1~')'), line=2, cex.lab=1.2)
text(16/3000, 0.12 ,'10-year cycle')
grid()

plot((0:188)/7560, (kep15[,44]/sum(kep15[,44]))[1:189] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='',ylim=c(0,0.16))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('KEP (' ~ alpha == 0.9 ~ ", " ~ k == 1.5~')'), line=2, cex.lab=1.2)
text(16/3000, 0.12 ,'10-year cycle')
grid()

plot((0:188)/7560, (kep19[,44]/sum(kep19[,44]))[1:189] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='',ylim=c(0,0.16))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('KEP (' ~alpha == 0.9 ~ ", " ~ k == 1.9~')'), line=2, cex.lab=1.2)
text(16/3000, 0.12 ,'10-year cycle')
grid()

#####
image((0:188)/7560, tau , spec.normalize(kep11)[1:189,] ,xaxt="n",yaxt="n",xlab='',ylab='', col = tim.colors())
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('Expectile ( k = 1.1 )'), line=2, cex.lab=1.2)

image((0:188)/7560, tau , spec.normalize(kep15)[1:189,] ,xaxt="n",yaxt="n",xlab='',ylab='',col = tim.colors())
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('Expectile ( k = 1.5 )'), line=2, cex.lab=1.2)

image((0:188)/7560, tau , spec.normalize(kep19)[1:189,] ,xaxt="n",yaxt="n",xlab='', ylab='',col = tim.colors())
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('Expectile ( k = 1.9 )'), line=2, cex.lab=1.2)

######
plot(tau, kep11[3,] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Expectile",ylab=expression('KEP (' ~ alpha == 0.9 ~ ", " ~ k == 1.1~')'), line=2, cex.lab=1.2)
grid()

plot(tau, kep15[3,],type='l',xaxt="n",yaxt="n",xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Expectile",ylab=expression('KEP (' ~ alpha == 0.9 ~ ", " ~ k == 1.5~')'), line=2, cex.lab=1.2)
grid()

plot(tau, kep19[3,] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Expectile",ylab=expression('KEP (' ~ alpha == 0.9 ~ ", " ~ k == 1.9~')'), line=2, cex.lab=1.2)
grid()


#########
smoothed19 = kep19
 for(i in 1:46){
   smoothed19[,i] = smooth.spline(kep19[,i], cv=T)$y
 }
load('garch_sdf.rdata')
par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
set.panel(1,2)
image((0:188)/7560, tau ,smoothed19 ,xaxt="n",yaxt="n",xlab='',ylab='',xlim=c(0,0.025),col = tim.colors())
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per day)",ylab="Expectile ( k = 1.9 )", line=2, cex.lab=1.2)

image(1:199/398, y = (5:95)/100, z=garch_sdf,xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(0,0.5),col = tim.colors())
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq",ylab="Expectile ( k = 1.9 )", line=2, cex.lab=1.2)

