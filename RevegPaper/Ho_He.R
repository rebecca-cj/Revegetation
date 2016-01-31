## average population level He 
# for microcarpa individuals excluding all outliers (RCV, RDK, RML, FDK, FML, NGR)
# using loci selected through population filtering using ALL SITES including site with <15 samples (see 'Round4_Analysis.txt' for more information)

## He ##
file<-dir(path="He_data/", pattern='*.He')
summary<-data.frame(row.names=c('mean','sd','SE'))
for(j in 1:length(file))
{
  tmp<-read.table(paste("He_data/",file[j],sep=""),header=T,row.names=NULL)
  obj<-c(mean(tmp[,7]),sd(tmp[,7]),(sd(tmp[,7])/sqrt(length(tmp[,7]))))
  summary<-data.frame(summary,obj)
  colnames(summary)[j]<-sub('.He','',(paste(file[j])))
}

summary_sorted<-data.frame(summary$VWGn,summary$VWGd,summary$FCV,summary$RCV,
                           summary$VBUn,summary$VBUd,summary$FPS,summary$RPS,
                           summary$VMCn,summary$VMCd,summary$FHH,summary$RHH,
                           summary$VDBn,summary$VDBd,summary$FDK,summary$RDK,
                           summary$NTRn,summary$NTRd,summary$FRJ,summary$RRJ,
                           summary$NCRn,summary$NCRd,summary$FDN,summary$RDN,
                           summary$NGRn,summary$NGRd,summary$FML,summary$RML)
rownames(summary_sorted)<-rownames(summary)
colnames(summary_sorted)<-sub('summary.','',colnames(summary_sorted))

He<-as.matrix(summary_sorted)

# add row of numbers for plotting
He<-as.matrix(rbind(He,c(seq(1,4),seq(6,9),seq(11,14),seq(16,19),seq(21,24),seq(26,29),seq(31,34))))

obj=He
png('MeanSE-He.png',width=12,height=6,units="in",res=300)
plot(obj[4,],obj[1,],pch=c(rep(c(16,1,6,3),7)),cex=1.25,axes=F,ylim=c(0.07,0.085),xlab=NA,ylab="He",cex.lab=1,cex.axis=1)
for(i in 1:28){arrows(obj[4,i],(obj[1,i]+obj[3,i]),obj[4,i],(obj[1,i]-obj[3,i]),angle=90,length=0.1,code=3,lwd=0.8)}
axis(side=2,tck=-0.02,las=2,cex.axis=0.8)
axis(side=1,labels=F,at=c(1,6,11,16,21,26,31,35),tck=-0.02,line=-0.75)
mtext(c("B. Marsh","Horsham","Bendigo","Dookie","The Rock","Crowther","Grenfell"),1,at=c(3,8,13,18,23,28,33),las=0,cex=1)
legend(27,0.085,c("Natural (normal sampling)","Natural (dense sampling)","Fragment","Revegetation"),pch=c(16,1,6,3),col=c(rep("black",4)),cex=0.6,pt.cex=0.8)
dev.off()


## average population level Ho 
# for microcarpa individuals excluding all outliers (RCV, RDK, RML, FDK, FML, NGR)
# using loci selected through population filtering using ALL SITES including site with <15 samples (see 'Round4_Analysis.txt' for more information)

# libraries
library(car)
library(agricolae)
library(lme4)
library(lsmeans)

## Ho ##
file<-dir(path="Ho_data/", pattern='*.het')
summary<-data.frame(row.names=c('mean','sd','SE'))
for(j in 1:length(file))
{
  tmp<-read.table(paste("Ho_data/",file[j],sep=""),header=T,row.names=NULL)
  obj<-c(mean(tmp[,7]),sd(tmp[,7]),(sd(tmp[,7])/sqrt(length(tmp[,7]))))
  summary<-data.frame(summary,obj)
  colnames(summary)[j]<-sub('.het','',(paste(file[j])))
}

summary_sorted<-data.frame(summary$VWGn,summary$VWGd,summary$FCV,summary$RCV,
                           summary$VBUn,summary$VBUd,summary$FPS,summary$RPS,
                           summary$VMCn,summary$VMCd,summary$FHH,summary$RHH,
                           summary$VDBn,summary$VDBd,summary$FDK,summary$RDK,
                           summary$NTRn,summary$NTRd,summary$FRJ,summary$RRJ,
                           summary$NCRn,summary$NCRd,summary$FDN,summary$RDN,
                           summary$NGRn,summary$NGRd,summary$FML,summary$RML)
rownames(summary_sorted)<-rownames(summary)
colnames(summary_sorted)<-sub('summary.','',colnames(summary_sorted))

Ho<-as.matrix(summary_sorted)

# add row of numbers for plotting
Ho<-as.matrix(rbind(Ho,c(seq(1,4),seq(6,9),seq(11,14),seq(16,19),seq(21,24),seq(26,29),seq(31,34))))

png('MeanSE-Ho.png',width=12,height=6,units="in",res=300)
plot(Ho[4,],Ho[1,],pch=c(rep(c(16,1,6,3),7)),cex=1.25,axes=F,ylim=c(0.07,0.085),xlab=NA,ylab="Avg Indv Observed Heterozygosity",cex.lab=1,cex.axis=1)
for(i in 1:28){arrows(Ho[4,i],(Ho[1,i]+Ho[3,i]),Ho[4,i],(Ho[1,i]-Ho[3,i]),angle=90,length=0.1,code=3,lwd=0.8)}
axis(side=2,tck=-0.02,las=2,cex.axis=0.8)
axis(side=1,labels=F,at=c(1,6,11,16,21,26,31,35),tck=-0.02,line=-0.75)
#for(i in c(1:8,13:16,25:28)){text(obj[4,i],(obj[1,i]+(obj[3,i]+0.0004)),significance[5,i],cex=0.7,)}
mtext(c("B. Marsh","Horsham","Bendigo","Dookie","The Rock","Crowther","Grenfell"),1,at=c(3,8,13,18,23,28,33),las=0,cex=1)
text(1,0.085,"(a)")
legend(25,0.085,c("Natural (normal sampling)","Natural (dense sampling)","Fragment","Revegetation"),pch=c(16,1,6,3),col=c(rep("black",4)),cex=0.6,pt.cex=0.8)
dev.off()

# Combine He and Ho plots

par(mfrow=c(2,1))
plot(Ho[4,],Ho[1,],pch=c(rep(c(16,1,6,3),7)),cex=1.25,axes=F,ylim=c(0.07,0.085),xlab=NA,ylab="Avg Indv Observed Heterozygosity",cex.lab=1,cex.axis=1)
for(i in 1:28){arrows(Ho[4,i],(Ho[1,i]+Ho[3,i]),Ho[4,i],(Ho[1,i]-Ho[3,i]),angle=90,length=0.1,code=3,lwd=0.8)}
axis(side=2,tck=-0.02,las=2,cex.axis=0.8)
axis(side=1,labels=F,at=c(1,6,11,16,21,26,31,35),tck=-0.02,line=-0.75)
mtext(c("B. Marsh","Horsham","Bendigo","Dookie","The Rock","Crowther","Grenfell"),1,at=c(3,8,13,18,23,28,33),las=0,cex=1)
text(1,0.085,"(a)")
legend(25,0.085,c("Natural (normal sampling)","Natural (dense sampling)","Fragment","Revegetation"),pch=c(16,1,6,3),col=c(rep("black",4)),cex=0.6,pt.cex=0.8)
plot(He[4,],He[1,],pch=c(rep(c(16,1,6,3),7)),cex=1.25,axes=F,ylim=c(0.07,0.085),xlab=NA,ylab="Expected heterozygosity",cex.lab=1,cex.axis=1)
for(i in 1:28){arrows(He[4,i],(He[1,i]+He[3,i]),He[4,i],(He[1,i]-He[3,i]),angle=90,length=0.1,code=3,lwd=0.8)}
axis(side=2,tck=-0.02,las=2,cex.axis=0.8)
axis(side=1,labels=F,at=c(1,6,11,16,21,26,31,35),tck=-0.02,line=-0.75)
for(i in c(1:8,13:16,25:28)){text(He[4,i],(He[1,i]+(He[3,i]+0.0005)),He.sig.plotting[5,i],cex=0.7,)}
mtext(c("B. Marsh","Horsham","Bendigo","Dookie","The Rock","Crowther","Grenfell"),1,at=c(3,8,13,18,23,28,33),las=0,cex=1)
text(1,0.085,"(b)")
legend(27,0.085,c("Natural (normal sampling)","Natural (dense sampling)","Fragment","Revegetation"),pch=c(16,1,6,3),col=c(rep("black",4)),cex=0.6,pt.cex=0.8)

# alternative plot

png('Ho_He_alt_MeanSE_a.png',width=12,height=9,units="in",res=300)
par(fig=c(0,1,0.4,1))
plot(He[4,],He[1,],pch=c(rep(c(16,1,6,3),7)),cex=1.25,axes=F,ylim=c(0.07,0.085),xlab=NA,ylab="Expected heterozygosity",cex.lab=1.1,cex.axis=0.8)
for(i in 1:28){arrows(He[4,i],(He[1,i]+He[3,i]),He[4,i],(He[1,i]-He[3,i]),angle=90,length=0.1,code=3,lwd=0.8)}
axis(side=2,tck=-0.02,las=2,cex.axis=0.8)
axis(side=1,labels=F,at=c(1,6,11,16,21,26,31,35),tck=-0.02,line=-0.75)
text(1,0.084,"(a)",cex=1.5)
legend(25,0.085,c("Large remnant (normal)","Large remnant (dense)","Small remnant","Revegetation"),pch=c(16,1,6,3),col=c(rep("black",4)),cex=1,pt.cex=1.2)
par(fig=c(0,1,0,0.6),new=T)
plot(Ho[4,],Ho[1,],pch=c(rep(c(16,1,6,3),7)),cex=1.25,axes=F,ylim=c(0.07,0.085),xlab=NA,ylab="Avg Indv Observed Heterozygosity",cex.lab=1.1,cex.axis=0.8)
for(i in 1:28){arrows(Ho[4,i],(Ho[1,i]+Ho[3,i]),Ho[4,i],(Ho[1,i]-Ho[3,i]),angle=90,length=0.1,code=3,lwd=0.8)}
axis(side=2,tck=-0.02,las=2,cex.axis=0.8)
axis(side=1,labels=F,at=c(1,6,11,16,21,26,31,35),tck=-0.02,line=-0.75)
mtext(c("B. Marsh","Horsham","Bendigo","Dookie","The Rock","Crowther","Grenfell"),1,at=c(3,8,13,18,23,28,33),las=0,cex=1.1)
text(1,0.084,"(b)",cex=1.5)
mtext("Location",side=1,line=2,cex=1.1)
dev.off()

