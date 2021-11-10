require(R.matlab)
rm(list = ls())
ecgpath='/data/p_02186/TMS_ECG2/analyses/EEGsham/EKGsham/'
setwd(ecgpath)
fdipath='/data/p_02186/TMS_ECG2/analyses/EEGsham/FDIdata/'
bigdata=c()
subids=c('VP01', 'VP02', 'VP03', 'VP04', 'VP05','VP06', 'VP07', 'VP08', 'VP09', 'VP10', 'VP11',
         'VP12', 'VP13', 'VP14', 'VP15', 'VP16', 'VP17', 'VP18', 'VP19','VP20', 'VP21', 'VP22', 'VP23',
         'VP24', 'VP25', 'VP26', 'VP28', 'VP29', 'VP30','VP31', 'VP32', 'VP33', 'VP34', 'VP35', 'VP36', 'VP37')
s=1
b=1
ind=74
fs=5000
# I deleted kubios fil

for (s in 1:length(subids)){
  subdata=c()
  for (b in 1:4){
    Rname=paste('ecgdata', subids[s], '_', b, '.mat', sep='')
    ecgname=paste(ecgpath, 'filtecg', subids[s],'_', 'tms', b, '.txt', sep='')
    
    if (! file.exists(ecgname)) next
    data=c()
    #data= readMat(paste(fdipath, 'data_all', subids[s], '_', b, '.mat' ,sep='')) #since it is sham no meps,
    data= readMat(paste(fdipath, 'data_tms_nomep', subids[s], '_', b, '.mat' ,sep='')) #since it is sham no meps,
    data=data[[1]]
    
    ECG=read.table(ecgname)
    
    R_peaks= readMat(Rname)
    R_peaks= R_peaks[[1]]
    #data=na.omit(data)
    data=cbind(matrix(data=NA,nrow=nrow(data),ncol=3), data, matrix(data=NA,nrow=nrow(data),ncol=5))
    for (ind in 1:nrow(data)) { 
      # the time of the peak just before stimulus onset
      data[ind,1] = s
      data[ind,2] = b
      data[ind,3] = ind
      if (is.na(data[ind,4])) next
      pos=max(which(R_peaks < data[ind,4]))
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
      diff2peak=data[ind,4] - R_peaks[pos]
      # relative position of onset in cardiac cycle
      stim_degree=360 * diff2peak/(R_peaks[pos+1] - R_peaks[pos])
      
      
      data[ind,6] = data[ind,4]<twave2[tend,1] # whether in systole
      data[ind,7] = data[ind,4]>diastole_ecg[1,1] # whether in diastole
      data[ind,8] = diff2peak
      data[ind,9] = stim_degree
      data[ind,10] = twave2[tend,1]-R_peaks[pos]
    }
    subdata=rbind(subdata,data)
  }
  bigdata=rbind(bigdata,subdata)
}
bigdata=data.frame(bigdata)
names(bigdata)=c('subject','block','trial','tms','mep','systole', 'diastole', 'diff2peak', 'stim_degree', 'syslength')
#write.csv2(bigdata,file="sys_dys_sham.csv", row.names = F)
write.csv2(bigdata,file="sys_dys_sham_nomep.csv", row.names = F)
#############################

# cancel the trials where Twave detection did not work -- systole detection
bigdata = read.csv2("/data/p_02186/TMS_ECG2/analyses/EEGsham/EKGsham/sys_dys_sham.csv")
rejmin=matrix(data=NA, nrow =length(unique(bigdata$subject)), ncol= 3)
rejmax=matrix(data=NA, nrow =length(unique(bigdata$subject)), ncol= 3)
meansys=matrix(data=NA,nrow=36,ncol=2)

for (s in 1:36){
  
  min_th=mean(bigdata$syslength[bigdata$subject==s],na.rm=T) -4*sd(bigdata$syslength[bigdata$subject==s],na.rm=T)  
  meansys[s,1]=min_th
  
  max_th=mean(bigdata$syslength[bigdata$subject==s],na.rm=T) +4*sd(bigdata$syslength[bigdata$subject==s],na.rm=T)  
  meansys[s,2]=max_th
}

minbigrej=c()
maxbigrej=c()
for (s in 1:36){
  minsubrej=c()  
  maxsubrej=c()
  if (sum(bigdata$subject[bigdata$syslength<meansys[s,1] & bigdata$subject==s ],na.rm=T)>0) {
    minsubrej=cbind(bigdata$subject[bigdata$syslength<meansys[s,1] & bigdata$subject==s ], bigdata$block[bigdata$syslength<meansys[s,1] & bigdata$subject==s],bigdata$X[bigdata$syslength<meansys[s,1] & bigdata$subject==s])
    minbigrej=rbind(minbigrej, minsubrej)
  }
  
  if (sum(bigdata$subject[bigdata$syslength>meansys[s,2] & bigdata$subject==s ],na.rm=T)) {
    maxsubrej=cbind(bigdata$subject[bigdata$syslength>meansys[s,2] & bigdata$subject==s ], bigdata$block[bigdata$syslength>meansys[s,2] & bigdata$subject==s],bigdata$X[bigdata$syslength>meansys[s,2] & bigdata$subject==s])
    maxbigrej=rbind(maxbigrej, maxsubrej)
  }
}

bigdata_cor=bigdata
rejs=rbind(minbigrej, maxbigrej)
rejs=na.omit(rejs)
for (r in 1:nrow(rejs)) {
  bigdata_cor=bigdata_cor[!(bigdata_cor$subject==rejs[r,1] & bigdata_cor$block==rejs[r,2] & bigdata_cor$X==rejs[r,3]),]
}
hist(bigdata_cor$syslength)
#write.csv2(bigdata_cor,file='/data/p_02186/TMS_ECG2/analyses/EEGsham/EKGsham/sys_dys_sham_final.csv', row.names = F)
#I did not exclude any trials since detection seems worked well even for these trials, manual check!


##
require(R.matlab)
rm(list=ls())
# ecgdata = read.csv2("/data/p_02186/TMS_ECG2/analyses/EEGsham/EKGsham/sys_dys_sham.csv")
# writeMat(con="/data/p_02186/TMS_ECG2/analyses/EEGsham/EKGsham/ecgdata.mat", x=as.matrix(ecgdata))

ecgdata = read.csv2("/data/p_02186/TMS_ECG2/analyses/EEGsham/EKGsham/sys_dys_sham_nomep.csv")
writeMat(con="/data/p_02186/TMS_ECG2/analyses/EEGsham/EKGsham/ecgdata_nomep.mat", x=as.matrix(ecgdata))
