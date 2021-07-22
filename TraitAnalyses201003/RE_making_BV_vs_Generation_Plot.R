rm(list=ls())
library(here)
here()
datafdr<-here("TraitAnalyses201003/data/")
#wd<-here("TraitAnalyses201003/UpdateAsh/")
wd<-here("TraitAnalyses201003/")

yr<-"Three"
diphap<-"MixedPloidy"
allBLUPs<-read.csv(paste0(datafdr,"allBLUPs_Plots+Individuals_withSGP_950_AddfndrsMrkData_",diphap,"_",yr,"_0715_2021.csv"),sep=",",header=TRUE,row.names=1)
  head(allBLUPs)
  dim(allBLUPs)
colnames(allBLUPs)<-c("DWpM","WWpM","pDW","BD","BL","BmWid","BTh","SL","SDia")  

### 1. fndr SPs, do not have FG/MG in it
sampledSP<-which(!(grepl("FG",rownames(allBLUPs))|grepl("MG",rownames(allBLUPs))))  #104
Fndr_BLUPs<-allBLUPs[sampledSP,]

### 2. All the Crosses we made that has "x" in the name, has FG/MG in it, BUT no "UCONN_S" in it
### Note, this are always UCONN-S crossing with a fdr progenySP. So no need to separate
 fndrSP<-which(grepl("x",rownames(allBLUPs))& !grepl("UCONN-S",rownames(allBLUPs)))
 UCONNSP<-which(grepl("x",rownames(allBLUPs))& grepl("UCONN-S",rownames(allBLUPs)))
	str(progenySP)
	str(UCONNSP)
 progenySP<-c(fndrSP,UCONNSP)	
SPCross_BLUPs<-allBLUPs[progenySP,]
	str(SPCross_BLUPs)                                    #329
#rownames(allBLUPs) is the same as that in u, and hMat

### 3. Find out SProg_GPs (UCONN_GPs)    # has FG/MG in it
	#### !!!!This is not correct, Because UCONN-S would also be in a cross that has used UCONN GP
SProg_GPs<-which(grepl("UCONN-S",rownames(allBLUPs)) & (!grepl("x",rownames(allBLUPs))))
SP_GPs_BLUPs<-allBLUPs[Sprog_GPs,]  
	str(SP_GPs_BLUPs)   # 78
 
#### 4. GPs from fnders (not SP_GPs, not SP_Cross, has the FG/MG in it)
releaseGP<-which(!(rownames(allBLUPs)%in%c(rownames(SP_GPs_BLUPs),rownames(SPCross_BLUPs))) & (grepl("FG",rownames(allBLUPs))|grepl("MG",rownames(allBLUPs))))   #439

releaseGP_BLUPs<-allBLUPs[releaseGP,]
# progenySP <- which(nchar(rownames(allBLUPs)) > 15)   #4. the Cross name is larger than 15 characters

fndNames <- rownames(allBLUPs)[sampledSP]
hasDesc <- sapply(fndNames, function(n) grep(n, rownames(allBLUPs)[progenySP])) # grep the fndr pattern in progenySP string
nDesc <- sapply(hasDesc, length)  #sapply: function applied to each element of a vector
  nDesc


traitname<-"DWpM" # !!!!!!!!!!

# Different colors for fndr, GP, progeny SP crosses
# colSPGP <- c(rep("Founder SPs", length(sampledSP)), rep("GPs from Founder", length(releaseGP)), rep("Farm SPs", length(progenySP)),rep("GPs from Farm SPs",nrow(SProg_GPs)))

colSPGP <- c(rep("Founder SPs", length(sampledSP)), rep("GPs from Founder", length(releaseGP)), rep("SPs progeny", length(fndrSP)),rep("SPs progeny used Farm GPs", length(UCONNSP)),rep("GPs from Farm SPs",nrow(SProg_GPs)))

colSPGP[which(nDesc == 0)] <- "Founder SPs without progeny"

OrderedBLUPs<-as.data.frame(rbind(allBLUPs[sampledSP,],allBLUPs[releaseGP,],allBLUPs[fndrSP,],allBLUPs[UCONNSP,],SProg_GPs))
  dim(OrderedBLUPs)
  str(colSPGP)
  sum(rownames(allBLUPs)%in%rownames(OrderedBLUPs))
  
OrderedBLUPs$colSPGP<-colSPGP  
OrderedBLUPs$Trait<-OrderedBLUPs[,traitname]  
  (dim(OrderedBLUPs))

write.csv(OrderedBLUPs,paste0("OrderedBLUPs_3yrs_",nrow(OrderedBLUPs),"Individuals.csv"))
##### Need to add and find out which SPs is from 2019 or 2020 ???????????/
  #OrderedBLUPs$Year<-expss::vlookup(OrderedBLUPS)
#### FIGURE 3 !!!!!!! Plot out each generation of BLUPs
pdf(paste0(datafdr,traitname,"vs Individual",yr,"_PhotoScore23_with_",diphap,"_Figure3_0718_2021.pdf"))
#plot(OrderedBLUPs[,colnames(OrderedBLUPs)%in%traitname], pch=16, col=colSPGP, xlab="Individual number", ylab=traitname) # !!The second col is DWpM

ggplot(data=OrderedBLUPs,mapping=aes(x=1:950,y=Trait,color=colSPGP))+
         geom_point()+
    scale_color_manual(name="Category",
                       values=c("Founder SPs"="black","Founder SPs without progeny"="grey","GPs from Founder"="dark green","SPs progeny"="dark red","SPs progeny used Farm GPs"="darksalmon","GPs from Farm SPs"="red"))+
  labs(x="Individual",y="Dry Weight per Meter (kg/m)")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank())
  

dev.off()  