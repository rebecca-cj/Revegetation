## Calculate and plot average genetic distance within a site
# Matrix of individual genetic distances (1 - Identity-By-State) generated in Plink

# libraries

library(gdata)

# create list containing population sample info
# each .pop file contains list of samples for that site

popinfo<-list()
files<-dir(pattern='*.pop')
for(i in 1:length(files)){
  input<-read.table(files[i],stringsAsFactors = F)
  obj<-as.vector(input[,1])
  name<-sub('.pop','',files[[i]])
  popinfo[[paste(name)]]<- obj
}

# read distance matrices
NFR.dist<-as.matrix(read.table('NFR_plink.mdist'))
NFR.names<-read.table('NFR_plink.mdist.id',stringsAsFactors =F)
rownames(NFR.dist)<-NFR.names[,2]
colnames(NFR.dist)<-NFR.names[,2]

# For within site comparisons (samples to themselves = "triangular matrix")
# calculate distance between individuals within a site based on using "triangular matrix" worth of data
  # removes duplicate comparisons found in square matrix when comparing samples to themselves (affects SE calculation)

distance.matrix=NFR.dist
output<-data.frame(stringsAsFactors = F,row.names=c('Site','Mean','SD','SE','count'))
for(p in names(popinfo)){
  total<-vector()
  count=0
  pop<-popinfo[[p]]
  for(i in 1:((length(pop))-1)){
    sample=pop[i]
    x<-distance.matrix[sample,]
    remaining<-pop[-(1:i)]
    for(j in 1:length(remaining)){
      compare=remaining[j]
      total<-c(total, x[[compare]])
      count<-count + 1
    }
  }
  popdist<-c(p,mean(total),sd(total),(sd(total)/sqrt(length(total))),count)
  output<-data.frame(stringsAsFactors = F,output,popdist)
}

within.site.dist<-data.frame(output[2:5,1:28],stringsAsFactors = F)
colnames(within.site.dist)<-output[1,1:28]

# order sites in data frame for "correct" plotting order

within.site.dist.sorted<-data.frame(stringsAsFactors = F,
                                    within.site.dist$VWGn,within.site.dist$VWGd,within.site.dist$FCV,within.site.dist$RCV,
                                    within.site.dist$VBUn,within.site.dist$VBUd,within.site.dist$FPS,within.site.dist$RPS,
                                    within.site.dist$VMCn,within.site.dist$VMCd,within.site.dist$FHH,within.site.dist$RHH,
                                    within.site.dist$VDBn,within.site.dist$VDBd,within.site.dist$FDK,within.site.dist$RDK,
                                    within.site.dist$NTRn,within.site.dist$NTRd,within.site.dist$FRJ,within.site.dist$RRJ,
                                    within.site.dist$NCRn,within.site.dist$NCRd,within.site.dist$FDN,within.site.dist$RDN,
                                    within.site.dist$NGRn,within.site.dist$NGRd,within.site.dist$FML,within.site.dist$RML)
rownames(within.site.dist.sorted)<-rownames(within.site.dist)
colnames(within.site.dist.sorted)<-sub('within.site.dist.','',colnames(within.site.dist.sorted))

# add row of x axis values for easier plotting

within.site.dist.sorted<-rbind(within.site.dist.sorted,x=c(seq(1,4),seq(6,9),seq(11,14),seq(16,19),seq(21,24),seq(26,29),seq(31,34)))


# plot within site distance (1-IBS)

obj=data.matrix(within.site.dist.sorted)
png('GD.png',width=12,height=6,units="in",res=300)
plot(obj[5,],obj[1,],axes=F,xlab=NA,ylim=c(0.06,0.075),ylab="Avg distance between individuals (1-IBS)",pch=c(rep(c(16,1,6,3),7)))
for(i in 1:28){arrows(obj[5,i],(obj[1,i]+obj[3,i]),obj[5,i],(obj[1,i]-obj[3,i]),angle=90,length=0.1,code=3)}
axis(side=2,tck=-0.02,las=2,cex.axis=0.8)
axis(side=1,labels=F,at=c(1,6,11,16,21,26,31,35),tck=-0.02,line=-0.75)
mtext(c("B. Marsh","Horsham","Bendigo","Dookie","The Rock","Crowther","Grenfell"),1,at=c(3,8,13,18,23,28,33),las=0,cex=1)
legend(26,0.075,c("Natural (normal sampling)","Natural (dense sampling)","Fragment","Revegetation"),pch=c(16,1,6,3),col=c(rep("black",4)),cex=0.8,pt.cex=1)
mtext("Location",side=1,line=2,cex=1.1)
dev.off()

