#### MUST NAME your first PCA_fcn object "PCA", in order for the other ones to work
#### 1. PCA_fcn()
#### G.num: the nxm matrix in -1,0,1 form, with rownames being the varieties
#### pop: defines which population you are working with
#### It writes out the PCA scores
#### It returns 1) Portion, 2) Scores and 3) Pop as the % of variance explained by PC as well as the PCscores
#### 
#### 2. PCA_lim_fcn() plots PCA again with the user defined xlim and ylim
####
#### 3. PCA_color_fcn() plots PCA with colors identified for each Group, 
#### 4. The Link file: with 3 columns, the first two has to be "Entry" and "Program"


#### 1.
PCA_fcn<-function (G.num, pop) {

G.num[G.num==1]<-2
G.num[G.num==0]<-1
G.num[G.num==-1]<-0

M<- apply(G.num,2,as.numeric)  
rownames(M)<-rownames(G.num)

#M.pca <- prcomp(M, retx=TRUE, center=TRUE, scale.=TRUE)
M.pca<-prcomp(M,retx=TRUE,center=TRUE)
M.pca.scores<-M.pca$x

# The order of the M.pca.scores is the same as Var order in M and in G.num

Prop<-M.pca$sdev^2/sum(M.pca$sdev^2)
Prop<-100*round(Prop,digits=4)

write.csv(M.pca.scores[,c(1:30)],paste(pop,"_PCA_scores.csv",sep="")) # Save only the first 30 PCs scores 
	# M.pca.scores[1:5,1:6]

### Plot without Link file
dev.set(1) 		### Plot PC1 VS PC2
tiff(file=paste(pop,"_PCA_testFig.tiff",sep=""),width=1200,height=800,units="px",pointsize=12,res=150)
par(cex=0.5)
plot(M.pca.scores[,1],M.pca.scores[,2],pch=1,main=paste("PCA for",pop,sep=" "),xlab=paste("PC1 (",Prop[1],"%)",sep=""),ylab=paste("PC2 (",Prop[2],"%)",sep=""))
dev.off()
 
 return(list(Portion=Prop, Scores=M.pca.scores, Pop=pop))	
}


#### 2.Re-run with xlim and ylim
PCA_lim_fcn<-function (xlim,ylim){
		Prop<-PCA$Portion
		M.scores<-PCA$Scores
		pop<-PCA$Pop
		
dev.set(2) 		### Plot PCA
tiff(file=paste(pop,"_PCA_x,y,lim.tiff",sep=""),width=1200,height=800,units="px",pointsize=12,res=150)
par(cex=0.5)

plot(M.scores[,1],M.scores[,2],pch=1,main=paste("PCA for",pop,sep=" "),xlim=xlim,ylim=ylim,xlab=paste("PC1 (",Prop[1],"%)",sep=""),ylab=paste("PC2 (",Prop[2],"%)",sep=""))
dev.off()
	
}


#### 3.

PCA_Link_fcn<-function(Link,M.scores,pop){
PCA.link<-merge(Link,M.scores,by.x="Entry",by.y="row.names",all.y=TRUE) ### Link file has the "Entry" and "Program

	#dim(PCA.link); PCA.link[1:5,1:7] 
	
write.csv(PCA.link,paste(pop,"_PCA_scores_link.csv",sep=""))

return(list(PCA.link=PCA.link))
}


#### 4.
PCA_color_fcn<-function(pop,Prop,PCA.link,xlim,ylim,color,xlegend,ylegend,cex) {
		
dev.set(3) 	
tiff(file=paste(pop,"_PCA_x,y,lim_with_color.tiff",sep=""),width=1200,height=800,units="px",pointsize=12,res=150)
par(cex=0.5)

plot(PCA.link$PC1,PCA.link$PC2,pch=1,main=paste("PCA for",pop,sep=" "),xlim=xlim,ylim=ylim,xlab=paste("PC1 (",Prop[1],"%)",sep=""),ylab=paste("PC2 (",Prop[2],"%)",sep=""),col=sapply(PCA.link$Program, color))


legend(xlegend,ylegend,legend=levels(PCA.link$Program),pch=1,col=unique(sapply(PCA.link$Program,color)),cex=cex)

dev.off()

}


save(PCA_fcn,PCA_lim_fcn,PCA_Link_fcn,PCA_color_fcn,file="PCA_function_Update.Rdata")

