source('fnexpectile.R')
source('fn.R')
library(fGarch)
library(quantreg)
library(expectreg)
library(fields)

###############eeg
tau=seq(from=5,to=95,2)/100
library(R.matlab)
data = readMat('Patient-3_seizure-2_Q5-Q6_17-21.mat')
data = data$tmpdata
data = data[1,1000:4000]

n <- length(data)
freq <- c(0:(n-1))/n
freq <- freq[freq>0 & freq <0.1]
ep = expectile_peri_ls(data, freq, tau, n.cores=2)
ep = spec.normalize(ep)
pg = (1/n)*(abs(fft(data)[2:301]))^2
pg=pg/sum(pg)

kep11 = kep(data, freq, tau, k=1.1 ,n.cores=4)
kep11 = spec.normalize(kep11)
kep15 = kep(data, freq, tau, k=1.5 ,n.cores=4)
kep15 = spec.normalize(kep15)
kep19 = kep(data, freq, tau, k=1.9 ,n.cores=4)
kep19 = spec.normalize(kep19)

#1
layout_matrix <- matrix(c(1, 1, 2, 3, 4, 5,6,7,8,9,10,11), nrow = 4, byrow = TRUE)
layout(mat = layout_matrix, widths = c(1, 1, 1))

par(mar=c(4.1,4,1.5,1.5),mgp=c(3,0.5,0)) 
plot((1:3001)/1000, data ,type='l',xaxt="n",yaxt="n",xlab='',ylab='',ylim=c(-100,510))
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Time (Second)",ylab=bquote('Voltage: '~mu~'V'), line=2, cex.lab=1.2)
grid()
lines((1:3001)/1000, rep(generalized_expectile(data,0.1,1.1),3001),col='blue',lty=1)
lines((1:3001)/1000, rep(generalized_expectile(data,0.1,1.5),3001),col='blue',lty=2)
lines((1:3001)/1000, rep(generalized_expectile(data,0.1,1.9),3001),col='blue',lty=3)

lines((1:3001)/1000, rep(generalized_expectile(data,0.9,1.1),3001),col='red',lty=1)
lines((1:3001)/1000, rep(generalized_expectile(data,0.9,1.5),3001),col='red',lty=2)
lines((1:3001)/1000, rep(generalized_expectile(data,0.9,1.9),3001),col='red',lty=3)

legend("top", 
       legend = c(expression(alpha == 0.1 ~ ", " ~ k == 1.1 ~'   '),
                  expression(alpha == 0.9 ~ ", " ~ k == 1.1),
                  expression(alpha == 0.1 ~ ", " ~ k == 1.5~'   '),
                  expression(alpha == 0.9 ~ ", " ~ k == 1.5),
                  expression(alpha == 0.1 ~ ", " ~ k == 1.9~'   '),
                  expression(alpha == 0.9 ~ ", " ~ k == 1.9)),
       col = c("red", "blue", "red", "blue", "red", "blue"),
       lty = c(1, 1, 2, 2, 3, 3), 
       lwd = 1,
       ncol = 3, bty = "n")

lines(c(1.4,1.85),c(-40,-40),lty=2)
lines(c(1.4,1.85),c(90,90),lty=2)
lines(c(1.4,1.4),c(-40,90),lty=2)
lines(c(1.85,1.85),c(-40,90),lty=2)
text(1.625,-80,'Main spike',cex=1.3)

lines(c(1.45,1.65),c(160,160),lty=2)
lines(c(1.45,1.65),c(380,380),lty=2)
lines(c(1.45,1.45),c(160,380),lty=2)
lines(c(1.65,1.65),c(160,380),lty=2)
text(1.3,350,'Bursts')

plot((0:149)/3000,(pg/sum(pg))[1:150] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab="PG", line=2, cex.lab=1.2)
grid()
text(30/3000, 0.3 ,'Main spikes')
text(75/3000, 0.05 ,'Bursting')

######
plot((0:149)/3000, (kep11[,44]/sum(kep11[,44]))[1:150] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('KEP (' ~ alpha == 0.9 ~ ", " ~ k == 1.1~')'), line=2, cex.lab=1.2)
text(30/3000, 0.11 ,'Main spikes',cex=1)
text(75/3000, 0.04 ,'Bursting',cex=1)
grid()

plot((0:149)/3000, (kep15[,44]/sum(kep15[,44]))[1:150] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('KEP (' ~ alpha == 0.9 ~ ", " ~ k == 1.5~')'), line=2, cex.lab=1.2)
text(30/3000, 0.15 ,'Main spikes',cex=1)
text(75/3000, 0.05 ,'Bursting',cex=1)
grid()

plot((0:149)/3000, (kep19[,44]/sum(kep19[,44]))[1:150] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('KEP (' ~ alpha == 0.9 ~ ", " ~ k == 1.9~')'), line=2, cex.lab=1.2)
text(30/3000, 0.155 ,'Main spikes',cex=1)
text(75/3000, 0.05 ,'Bursting',cex=1)
grid()

plot((0:149)/3000, (kep11[,4]/sum(kep11[,4]))[1:150] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('KEP (' ~ alpha == 0.1 ~ ", " ~ k == 1.1~')'), line=2, cex.lab=1.2)
text(30/3000, 0.09*3 ,'Main spikes',cex=1)
text(75/3000, 0.038 ,'No bursting',cex=1)
grid()

plot((0:149)/3000, (kep15[,4]/sum(kep15[,4]))[1:150] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('KEP (' ~ alpha == 0.1 ~ ", " ~ k == 1.5~')'), line=2, cex.lab=1.2)
text(30/3000, 0.15*2 ,'Main spikes',cex=1)
text(75/3000, 0.05 ,'No bursting',cex=1)
grid()

plot((0:149)/3000, (kep19[,4]/sum(kep19[,4]))[1:150] ,type='l',xaxt="n",yaxt="n",xlab='',ylab='')
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('KEP (' ~ alpha == 0.1 ~ ", " ~ k == 1.9~')'), line=2, cex.lab=1.2)
text(30/3000, 0.178*2 ,'Main spikes',cex=1)
text(75/3000, 0.05 ,'No bursting',cex=1)
grid()

image((0:149)/3000, tau , spec.normalize(kep11)[1:150,]^0.7 ,xaxt="n",yaxt="n",xlab='',ylab='', col = tim.colors())
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('Expectile ( k = 1.1 )'), line=2, cex.lab=1.2)

image((0:149)/3000, tau , spec.normalize(kep15)[1:150,]^0.7 ,xaxt="n",yaxt="n",xlab='',ylab='',col = tim.colors())
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('Expectile ( k = 1.5 )'), line=2, cex.lab=1.2)

image((0:149)/3000, tau , spec.normalize(kep19)[1:150,]^0.7 ,xaxt="n",yaxt="n",xlab='', ylab='',col = tim.colors())
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
title(xlab="Freq (cycles per ms)",ylab=expression('Expectile ( k = 1.9 )'), line=2, cex.lab=1.2)