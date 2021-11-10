require(R.matlab)
rm(list = ls())
ecgpath='/data/p_02186/TMS_ECG2/analyses/EEGmotor/EKGdata/'
setwd(ecgpath)
eventpath='/data/p_02186/TMS_ECG2/analyses/EEGmotor/FDIdata/events/'
bigdata=c()
subids=c('VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11',
         'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23',
         'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37')
s=1
b=1
ind=1
fs=5000
# I deleted kubios fil

for (s in 1:length(subids)){
  subdata=c()
  Rname=paste('ecgdata', subids[s], '_', b, '.mat', sep='')
  ecgname=paste(ecgpath, 'filtecg', subids[s],'_', 'motor', b, '.txt', sep='')
  
  if (! file.exists(ecgname)) next
  data=c()
  data= readMat(paste(eventpath, 'motor_event', subids[s], '.mat' ,sep=''))
  data=data[[1]]
  data=t(data)
  ECG=read.table(ecgname)
  
  R_peaks= readMat(Rname)
  R_peaks= R_peaks[[1]]
  #data=na.omit(data)
  data=cbind(matrix(data=NA,nrow=nrow(data),ncol=2), data, matrix(data=NA,nrow=nrow(data),ncol=5))
  for (ind in 1:nrow(data)) { 
    # the time of the peak just before stimulus onset
    if (is.na(data[ind,3])) next
    pos=max(which(R_peaks < data[ind,3]))
    ecgpos1=R_peaks[pos]*fs 
    RRint=(R_peaks[pos+1]-R_peaks[pos])*fs #in data points
    ecgpos2=ecgpos1+RRint-300   
    twave1=ECG[(ecgpos1-600):ecgpos2+600,]
    twave=ECG[(ecgpos1+700):ecgpos2,]
    tmaxpos=which.max(twave[1:((RRint-700)/3),2]) 
    twave2=twave[tmaxpos:dim(twave)[1],]
    
    dp=0.12*fs 
    if (dp>dim(twave2)[1]) {
      xm=which(diff(twave2[,2])==min(diff(twave2[,2]))) 
    } else {
      xm=which(diff(twave2[1:dp,2])==min(diff(twave2[1:dp,2]))) 
    }
    xm=xm[1]
    ym=twave2[xm,2]
    xr=500+xm 
    xseq=xm:xr
    yseq=twave2[xm:xr,2]
    
    #write a function find the end of twave
    trapez_area <- function(xm, ym, xseq, yseq, xr) {
      a <- numeric()
      for (i in seq_along(xseq)){
        a[i] <- 0.5 * (ym - yseq[i]) * ((2*xr) - xseq[i] - xm)
      }
      x_tend <- which.max(a)+xm-1
      return(x_tend)
    }
    tend=trapez_area(xm, ym, xseq, yseq, xr)
    
    
    nextR=(R_peaks[pos+1]*fs)-10
    syslen=twave2[tend,1]-R_peaks[pos]
    diastole_dp=(nextR-syslen*fs):(nextR-1)
    diastole_ecg=ECG[diastole_dp,]
    diastole_ecg[1,1]
    
    # par(mfrow=c(1,2))
    # # #To check whether twave end detection worked well
    # plot(twave1,col='black',xlab='time(ms)', ylab='electrical potential(mV)')
    # points(twave[tmaxpos,1],twave[tmaxpos,2],col='magenta',pch='+',cex=2)
    # points(twave2[tend,1],twave2[tend,2],col='green',pch='+',cex=2)
    # points(diastole_ecg, col='blue')
    # 
    # plot(twave2,col='black', xlab='time(ms)', ylab='electrical potential(mV)')
    # title(paste('subject ',s,'_block ',b,'_trial ', ind, sep=''),line=-2, outer=TRUE)
    # points(twave2[xm,1],twave2[xm,2],col='blue',pch='+',cex=2)
    # points(twave2[xr,1],twave2[xr,2],col='blue',pch='+',cex=2)
    # points(twave2[tend,1],twave2[tend,2],col='green',pch='+',cex=2)
    # points(twave[tmaxpos,1],twave[tmaxpos,2],col='magenta',pch='+',cex=2)
    
    # 
    
    # the difference between the onset and previous R peak
    diff2peak=data[ind,3] - R_peaks[pos]
    # relative position of onset in cardiac cycle
    stim_degree=360 * diff2peak/(R_peaks[pos+1] - R_peaks[pos])
    
    data[ind,1] = s
    data[ind,2] = b
    data[ind,4] = data[ind,3]<twave2[tend,1] # whether in systole
    data[ind,5] = data[ind,3]>diastole_ecg[1,1] # whether in diastole
    data[ind,6] = diff2peak
    data[ind,7] = stim_degree
    data[ind,8] = twave2[tend,1]-R_peaks[pos]
  }
  
  bigdata=rbind(bigdata,data)
}
bigdata=data.frame(bigdata)
names(bigdata)=c('subject','block','tms','systole', 'diastole', 'diff2peak', 'stim_degree', 'syslength')
write.csv2(bigdata,file="/data/p_02186/TMS_ECG2/analyses/EEGmotor/FDIdata/sys_dys_motor_onset.csv")
#####
require(R.matlab)
rm(list=ls())
ecgdata = read.csv2("/data/p_02186/TMS_ECG2/analyses/EEGmotor/FDIdata/sys_dys_motor_onset.csv")
writeMat(con="/data/p_02186/TMS_ECG2/analyses/EEGmotor/FDIdata/sys_dys_motor_onset.mat", x=as.matrix(ecgdata))

