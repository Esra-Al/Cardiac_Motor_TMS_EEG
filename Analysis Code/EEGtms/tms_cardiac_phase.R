require(R.matlab)
rm(list = ls())
ecgpath='/data/p_02186/TMS_ECG2/analyses/EKGdata/'
setwd(ecgpath)
fdipath='/data/p_02186/TMS_ECG2/analyses/EEGtms/FDIdata/'
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
    data= readMat(paste(fdipath, 'data', subids[s], '_', b, '.mat' ,sep=''))
    data=data[[1]]
    
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
      data[ind,5] = data[ind,3]<twave2[tend,1] # whether in systole
      data[ind,6] = data[ind,3]>diastole_ecg[1,1] # whether in diastole
      data[ind,7] = diff2peak
      data[ind,8] = stim_degree
      data[ind,9] = twave2[tend,1]-R_peaks[pos]
    }
    subdata=rbind(subdata,data)
  }
  bigdata=rbind(bigdata,subdata)
}
bigdata=data.frame(bigdata)
names(bigdata)=c('subject','block','tms','mep','systole', 'diastole', 'diff2peak', 'stim_degree', 'syslength')
write.csv2(bigdata,file="sys_dys_data.csv")
#############################

# cancel the trials where Twave detection did not work -- systole detection
bigdata = read.csv2("/data/p_02186/TMS_ECG2/analyses/EKGdata/sys_dys_data.csv")
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
#write.csv2(bigdata_cor,file= paste('/data/p_02186/TMS_ECG2/analyses/EKGdata/sys_dys_data_final.csv', sep=''), row.names = F)
#I did not exclude any trials since detection seems worked well even for these trials, manual check!



#######################
rm(list=ls())
ecgdata = read.csv2("/data/p_02186/TMS_ECG2/analyses/EKGdata/sys_dys_data.csv")
writeMat(con="/data/p_02186/TMS_ECG2/analyses/EKGdata/det_diff2peak/ecgdata.mat", x=as.matrix(ecgdata))

################
rm(list=ls())
ecgdata = read.csv2("/data/p_02186/TMS_ECG2/analyses/EKGdata/sys_dys_data.csv")
rrint =  aggregate(RRinterval ~ subject, ecgdata, FUN=mean)
mean(rrint$RRinterval*1000)*0.32

systole = aggregate(syslength ~ subject, ecgdata, FUN=mean)
mean(systole$syslength*1000)
sd(systole$syslength*1000)

mep_amp = aggregate(mep ~ subject, ecgdata, FUN=mean)
mean(mep_amp$mep)
sd(mep_amp$mep)

require(reshape)
library(svglite)
#bigdata=subset(bigdata, subject!=1)
ecgdata = na.omit(ecgdata)
systole = ecgdata[ecgdata$systole==1,]
diastole = ecgdata[ecgdata$diastole==1,]
subnum = length(unique(ecgdata$subject))
subjects = unique(ecgdata$subject)

mep_sys = matrix(data=NA,nrow=subnum,ncol=1)
mep_dys = matrix(data=NA,nrow=subnum,ncol=1)

for (s in 1:subnum) {
  mep_sys[s] = mean(systole$mep[systole$subject==subjects[s]])
  mep_dys[s] = mean(diastole$mep[diastole$subject==subjects[s]])
}

mean(mep_sys, na.rm=T)
mean(mep_dys, na.rm=T)
t.test(mep_sys, mep_dys, paired = T)
wilcox.test(mep_sys,mep_dys,paired=TRUE)
mep_box=data.frame(heart = c(rep('systole',subnum),rep('diastole',subnum)),
                    mep = c(mep_sys, mep_dys))

mep_box=na.omit(mep_box)
mep_box$heart=as.character(mep_box$heart)
mep_box$heart=factor(mep_box$heart, levels=unique(mep_box$heart))
library(ggplot2)
#png('mep_tms_boxplot.png',width=2000, height=1400,res=300)
figcor=ggplot(mep_box, aes(x=heart, y=mep, color=heart)) +
  #stat_compare_means(aes(group = heart),  label = "p.signif", size=6, label.y = c(64,82), method = 't.test', paired=TRUE, tip.length = c(0.1, 0.03)) +
  geom_boxplot(alpha=0, outlier.size=NA, position=position_dodge(0.8), outlier.shape = NA) +
  geom_point(alpha=0.5, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0.1), fill="white", shape=21) +
  scale_color_manual(name = "heart", values = c("red", "blue"))+
  labs(x=' ', y="mep ", fill=' ') +
  #scale_y_continuous(breaks=c(0, 1000, 2000), limits = c(0,2000))+
   theme_bw() +
  theme(axis.text.x= element_blank(),
        axis.text.y=  element_blank(),
        axis.title.y= element_blank(),
        #legend.text=element_text(size=15),  legend.key.size=unit(2,"line"), legend.position=c(.87, .2),
       # legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave('mep_tms_bin_boxplot.svg',figcor,width=5,height=3)

library(ggpubr)
ggqqplot(mep_box$mep)
shapiro.test(mep_box$mep)


sum(mep_sys>mep_dys)
#######

require(R.matlab)
tepdiff=readMat('/data/p_02186/TMS_ECG2/analyses/EEGtms/merge/cardiacphase/oldcleaning/sysdys_tep.mat')
tepdiff=tepdiff$sysdys.tep
mepdiff=t(mep_sys-mep_dys)
shapiro.test(mepdiff)
cormat=cor.test(tepdiff, mepdiff, method = c("pearson"))

