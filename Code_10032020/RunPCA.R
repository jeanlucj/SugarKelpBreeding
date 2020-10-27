#Run PCA on the P1P2 SNPs data

##!!!##rownames(PCA$Scores)[which(rownames(PCA$Scores)=="SL18-LD-13-Mg-3")]<-"SL18-LD-13-MG3"


#GPSNP<-read.csv("GPsnp_mainGenome_NA0.8_P1P2P3.csv",sep=",",header=T)

rm(list=ls())
load("GPsAmat_NA0.8.Rdata")
  ls()
  dim(GPSNP)
  GPSNP[1:4,1:6]  
  GPSNP[1:4,(ncol(GPSNP)-4):ncol(GPSNP)]
GPSNP<-imputedFromAmat

#rownames(GPSNP)[which(rownames(GPSNP)=="SL18-LD-13-Mg-3")]<-"SL18-LD-13-RepMG3"


pop<-"P1P2P3"
PCA<-PCA_fcn(GPSNP,pop = pop)

xlim=c(-400,600)
ylim=c(-400,600)

PCA2<-PCA_lim_fcn(xlim=xlim,ylim=ylim)

#
Link<-read.csv("3_plates_sequenced_names_Update.csv",sep=",",header=T)
  head(Link)
  
Link$Parent<-as.character(Link$Parent)  
Link$Entry<-Link$Name_InGeno #!!
Link$Program<-Link$Region #!!


PCA3<-PCA_Link_fcn(Link,M.scores=PCA$Scores,pop=pop)


#PCA4 part by hand
# PCA.link<-merge(Link,as.matrix(PCA$Scores),by.x="Entry",by.y="row.names",all.y=TRUE) ### Link file has the "Entry" and "Program
#   PCA.link[1:4,1:11]
#   
# #PCA.link<-PCA3$PCA.link[-which(PCA3$PCA.link$Entry=="SL18-LD-13-RepMG3"),]

# 
# PCA.link$Program<-factor(PCA.link$Program)



color<-function (x){
  if (x=="GOM"){
    return("red")
  }else if(x=="SNE"){
    return("black")
  }else if(x=="Check"){
    return("blue")
  }
}



#sapply(PCA.link$Program, color)


### Code change to ggplot, so that color and shapes separate by groups
#
PCA.link<-PCA3$PCA.link
Prop<-PCA$Portion

dev.set(3) 	
tiff(file=paste(pop,"_PCA_x,y,lim_with_color_Shape.tiff",sep=""),width=1200,height=800,units="px",pointsize=12,res=300)
par(cex=0.8)

xlim=c(-400,800)
ylim=c(-400,800)

###### Still need to fix this
plot(PCA.link$PC1,PCA.link$PC2,main=paste("PCA for",pop,sep=" "),xlim=xlim,ylim=ylim,xlab=paste("PC1 (",Prop[1],"%)",sep=""),ylab=paste("PC2 (",Prop[2],"%)",sep=""),col=sapply(PCA.link$Region, color),pch=as.numeric(PCA.link$Location))
legend(610,650,legend=levels(PCA.link$Program),pch=unique(as.numeric(PCA.link$Location)),col="black",cex=0.4)
######

#plot(PCA.link$PC1,PCA.link$PC2,main=paste("PCA for",pop,sep=" "),xlim=xlim,ylim=ylim,xlab=paste("PC1 (",Prop[1],"%)",sep=""),ylab=paste("PC2 (",Prop[2],"%)",sep=""),col=as.numeric(PCA.link$Parent),pch=as.numeric(PCA.link$Location))

dev.off()






PCA4<-PCA_color_fcn(pop=pop,Prop=Prop,PCA.link=PCA.link,xlim=xlim,ylim=ylim,color=color,xlegend=400,ylegend=500,cex=2)



