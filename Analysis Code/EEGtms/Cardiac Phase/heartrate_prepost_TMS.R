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
    data=cbind(matrix(data=NA,nrow=nrow(data),ncol=2), data, matrix(data=NA,nrow=nrow(data),ncol=10))
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
      if (pos>2) {
        data[ind,10]=(R_peaks[pos-1] - R_peaks[pos-2]) #pre2
      }
      data[ind,11]=(R_peaks[pos] - R_peaks[pos-1]) #pre1
      data[ind,12]=(R_peaks[pos+1] - R_peaks[pos]) #current
      data[ind,13]=(R_peaks[pos+2] - R_peaks[pos+1]) #post1
      data[ind,14]=(R_peaks[pos+3] - R_peaks[pos+2]) #post2
    }
    subdata=rbind(subdata,data)
  }
  bigdata=rbind(bigdata,subdata)
}
bigdata=data.frame(bigdata)
names(bigdata)=c('subject','block','tms','mep','systole', 'diastole', 'diff2peak', 'stim_degree', 'syslength',  'RRpre2','RRpre1','RRcur','RRpost1','RRpost2')
write.csv2(bigdata,file="heartrate_prepost_TMS.csv")
##########

rm(list = ls())
library(svglite)
library(ggplot2)
data=read.csv2("/data/p_02186/TMS_ECG2/analyses/EKGdata/heartrate_prepost_TMS.csv")
data=subset(data, systole==1 | diastole==1)

RRsub=aggregate(RRcur ~ subject, data, FUN=mean)
RRpre2=aggregate(RRpre2 ~ subject*systole, data, FUN=mean)
RRpre1=aggregate(RRpre1 ~ subject*systole, data, FUN=mean)
RRcur=aggregate(RRcur ~ subject*systole, data, FUN=mean)

t.test(RRcur$RRcur[RRcur$systole==0], RRcur$RRcur[RRcur$systole==1], paired = T)
RRpost1=aggregate(RRpost1 ~ subject*systole, data, FUN=mean)
RRpost2=aggregate(RRpost2 ~ subject*systole, data, FUN=mean)


hrdata= data.frame(subject=factor(), interval=factor(), time=factor(),  hr=numeric())
#hrdata=rbind(hrdata, data.frame(subject=RRpre2$subject, interval= RRpre2$systole, time=rep(1, length(RRpre2$subject)),hr=RRpre2$RRpre2))
hrdata=rbind(hrdata, data.frame(subject=RRpre1$subject, interval= RRpre1$systole, time=rep(2, length(RRpre2$subject)),hr=RRpre1$RRpre1))
hrdata=rbind(hrdata, data.frame(subject=RRcur$subject, interval= RRcur$systole, time=rep(3, length(RRpre2$subject)),hr=RRcur$RRcur))
hrdata=rbind(hrdata, data.frame(subject=RRpost1$subject, interval= RRpost1$systole, time=rep(4, length(RRpre2$subject)),hr=RRpost1$RRpost1))
#hrdata=rbind(hrdata, data.frame(subject=RRpost2$subject, interval= RRpost2$systole, time=rep(5, length(RRpre2$subject)),hr=RRpost2$RRpost2))


hrdata$interval=as.factor(hrdata$interval)
hrdata$time=as.factor(hrdata$time)
hrdata$subject=as.factor(hrdata$subject)

## Create three functions used to calculate within-subjects variation. These fragments of code are credited to Winston Chang and were used and copied here under the CC0 license (http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#understanding-within-subjects-error-bars)
# Function 1: normDataWithin
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL, na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  # Remove this subject mean column
  data$subjMean <- NULL
  return(data)
}
# Function 2: summarySEwithin
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL, idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}
# Function 3: summarySE
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}
## Calculate within-subjects variation
d <- summarySEwithin(data=hrdata, measurevar = 'hr', betweenvars=NULL, withinvars=c('time','interval'), idvar='subject', na.rm=FALSE, conf.interval=.95, .drop=TRUE)

# Change into milliseconds 
d$hr <- d$hr * 1000
d$ci <- d$ci * 1000
# Prepare descriptions
#labs <- c("S-2", "S-1", "Stimulus", "S+1", "S+2", "S+3")
labs <- c( "S-1", "Stimulus", "S+1", "S+2")

# Save interval as factor
d$interval <- factor(d$interval)

pd <- position_dodge(0)
fig1 <- ggplot(d, aes(x=time, y=hr, group=interval,colour=interval)) + 
  geom_ribbon(aes(ymin=hr-ci, ymax=hr+ci, fill = factor(interval)), position=pd, show.legend = F, alpha = 0.3, colour = NA) +
  scale_fill_manual(values=c("blue3", "firebrick3")) +
  geom_line(position=pd, size = 0.7)+
  geom_point(position=pd)+theme_classic()+
  labs(x = " ", y = "Interbeat interval (ms)") +
  scale_x_discrete(labels=labs) +
  scale_y_continuous(limits=c(960,980)) +
  scale_color_manual(values=c("blue3", "firebrick3")) 
print(fig1)
ggsave('/data/p_02186/TMS_ECG2/analyses/EKGdata/heartrate_prepost_TMS.svg',fig1)

library(afex)
## Perform two-way repeated measures ANOVA
# Run to see ucorrected degrees of freedom
fit_all <- aov_ez("subject","hr", hrdata, within=c("time", "interval"), return = "nice")
fit_all 
# Perform test with corrected degrees of freedom
fit_all <- aov_ez("subject","hr", hrdata, within=c("time", "interval"))
fit_all # to see corrected degrees of freedom 
summary(fit_all) # see epsilon values
# Calculate confidence intervals
ref <- emmeans::emmeans(fit_all, specs = c("interval", "time"))
ref 
# Post-hoc Bonferroni-corrected paired t tests. 
res=emmeans::contrast(ref,method="pairwise",adjust="bonferroni")
summary(res)




#load ezANOVA package
library("ez") 
hrdata$time=as.factor(hrdata$time)
hrdata$interval=as.factor(hrdata$interval)
# call to anova wrapper
output_anova = ezANOVA(
  data = hrdata  # dataframe containing all relevant variables, see below
  , dv = hr     # dependent variable: P50amp/N150amp/...
  , wid = subject     # array with the subject number within the dataframe
  , within= .(time, interval) # specify the names of within-subject factors
  , type =3        # sum-of-squares-type, should be '3' in your case, but check function help to be sure
  , detailed = TRUE   # some output options
  , return_aov = TRUE
)
print(output_anova)

t.test(RRpre1$RRpre1[RRpre1$systole ==0], RRpre1$RRpre1[RRpre1$systole == 1], paired = T)
t.test(RRcur$RRcur[RRcur$systole ==0], RRcur$RRcur[RRcur$systole == 1], paired = T)
t.test(RRpost1$RRpost1[RRpost1$systole ==0], RRpost1$RRpost1[RRpost1$systole == 1], paired = T)

t.test(RRpre1$RRpre1[RRpre1$systole ==1], RRcur$RRcur[RRcur$systole ==1], paired = T)
t.test(RRcur$RRcur[RRcur$systole ==1], RRpost1$RRpost1[RRpost1$systole ==1], paired = T)

t.test(RRpre1$RRpre1[RRpre1$systole ==0], RRcur$RRcur[RRcur$systole ==0], paired = T)
t.test(RRcur$RRcur[RRcur$systole ==0], RRpost1$RRpost1[RRpost1$systole ==0], paired = T)

t.test(RRpre1$RRpre1[RRpre1$systole ==0], RRpost1$RRpost1[RRpost1$systole ==0], paired = T)
####

