# Within and Between Locations
setwd("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/SugarKelpBreeding/TraitAnalyses200820/")
rm(list=ls())
#install.packages("splitstackshape")

#### Subset for Site, Sex, clean or green
selection<-function(site,sex,clean,green,data){
  #if (green=="x"){
  #  set<-subset(data,Site==site & Gender==sex & Green==green & is.na(Clean))	
  #} else{
    set<-subset(data,Site==site & Gender==sex & Clean==clean & is.na(Green))		
  #}
  ### Calculate the Cross.num per Loc and per SP
  num.Loc<-aggregate(set$Cross.num,by=list(set$Loc),FUN=sum,na.rm=TRUE)  ### Cross# per location
  num.SP<-aggregate(set$Cross.num,by=list(set$SP),FUN=sum,na.rm=TRUE)    ### Cross# per SP
  return(list(set=set,num.Loc=num.Loc,num.SP=num.SP))
}


####
summary<-read.csv("Active_culture_biomass_assesment_for_Mao_100120.csv",sep=",",header=T,na.strings=c("","NA"))  ##!! 
### Use na.string to replace empty cells to NAs

summary2<-summary[,colnames(summary)%in%c("GametophyteID","Location","Site","Sex","X..Clean","X..Green","X..Crosses")]
  head(summary2)
  str(summary2)
colnames(summary2)<-c("GPID","Site","Gender","Loc","Cross.num","Clean","Green")

# summary2$Gender<-as.character(summary2$Gender)
# summary2$Gender[which(summary2$Gender=="Male")]="M"
# summary2$Gender[which(summary2$Gender=="Female")]="F"
  summary2[140,]  ###### This is what we have so far
  
summary2$Cross.num<-as.numeric(as.character(summary2$Cross.num))  #### Make numeric !!!! otherwise factor changes the cross num
# There is no NAs in the summary2
#summary2$Cross.num[is.na(summary2$Cross.num)] <- 0     #### Make NA to be 0

GPName<-str_split_fixed(summary2$GPID,"-",4)
summary2$SP<-paste0(GPName[,1],"-",GPName[,2],"-",GPName[,3])  #



### Within GOM Site, Clean; Green  !!!!!
site<-"GOM"


### Within SNE Site, Clean; Green  !!!!!
site<-"SNE"




setC.F<-selection(site=site,sex="F",clean="x",green="NA",summary2)
setC.M<-selection(site=site,sex="M",clean="x",green="NA",summary2)

# setG.F<-selection(site=site,sex="F",clean="NA",green="x",summary2)
# setG.M<-selection(site=site,sex="M",clean="NA",green="x",summary2)

lsSite<-list(setC.F,setC.M)

setwd("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/SugarKelpBreeding/TraitAnalyses200820/SNE_NewCrosses")

################## Summarize this part for manual checking 10 and 10 crosses
##### Sum Cross.num per Location
i=1
#### Cross.num summed up by Loc Clean(C)
num.Loc.C<-merge(lsSite[[i]]$num.Loc,lsSite[[i+1]]$num.Loc,by.x="Group.1",by.y="Group.1",all=TRUE)
colnames(num.Loc.C)<-c(paste("C_",site,"_Loc",sep=""),"Cross.num_F","Cross.num_M")

# #### Cross.num summed up by Loc Green(G)
# num.Loc.G<-merge(lsSite[[i+2]]$num.Loc,lsSite[[i+3]]$num.Loc,by.x="Group.1",by.y="Group.1",all=TRUE)
# colnames(num.Loc.G)<-c(paste("C_",site,"_Loc",sep=""),"Cross.num_F","Cross.num_M")

# Site.num<-merge(num.Loc.C,num.Loc.G,by.x=paste("C_",site,"_Loc",sep=""),by.y=paste("C_",site,"_Loc",sep=""),all=TRUE)
# name1<-gsub(".x",".Clean",colnames(Site.num))  ### data x is the Clean one: num.Loc.C
# name2<-gsub(".y",".Green",name1)       ### data y is the Green one: num.Loc.G
# name2
# colnames(Site.num)<-name2
write.csv(num.Loc.C,paste(site,"_Crossnum_per_Loc_1005.csv",sep=""))

##### Cross.num summed up by SP Clean(C)
i=1
#### Cross.num summed up by SP Clean(C)
num.C<-merge(lsSite[[i]]$num.SP,lsSite[[i+1]]$num.SP,by.x="Group.1",by.y="Group.1",all=TRUE)
colnames(num.C)<-c(paste("SP_",site,sep=""),"Cross.num_F","Cross.num_M") ### !!!

# #### Cross.num summed up by SP Green(G)
# num.G<-merge(lsSite[[i+2]]$num.SP,lsSite[[i+3]]$num.SP,by.x="Group.1",by.y="Group.1",all=TRUE)
# colnames(num.G)<-c(paste("SP_",site,sep=""),"Cross.num_F","Cross.num_M")

# Site.num<-merge(num.C,num.G,by.x=paste("SP_",site,sep=""),by.y=paste("SP_",site,sep=""),all=TRUE)
# name1<-gsub(".x",".Clean",colnames(Site.num))  ### data x is the Clean one: num.Loc.C
# name2<-gsub(".y",".Green",name1)       ### data y is the Green one: num.Loc.G
# name2
# colnames(Site.num)<-name2
# Site.num
  num.C
write.csv(num.C,paste(site,"_Crossnum_per_SP_1005.csv",sep="")) ### !!!

colSums(num.Loc.C[,2:3],na.rm=TRUE)
#colSums(num.Loc.G[,2:3],na.rm=TRUE)

####################################

### From same Loc #10 (CC, JS, NC, NL, OI--- GOM. BL, FI, FW, PI, TI --- SNE)
### From same SP # 10	
### Assign random numbers to the F and then to M
### After (somewhat manually)remove the ones that are 10, 10, 10
### Randomly assign across locations


library("splitstackshape")
Random<-function(site,sex,clean,green,df){	
  ### Expand based on how many Cross.num that GP has
  new_df<-expandRows(df, "Cross.num") 
  new_df$Random<-runif(1:nrow(new_df))
  
  ### Randomly assign a number to the GPs
  dataY<-new_df
  nlines<-nrow(dataY)
  sets<-rep(1:nlines,1)               ### Make it # nlines !!!
  sample.num<-sets[order(runif(nrow(dataY)))]
  dataY$Random<-sample.num
  
  write.csv(dataY,file=paste(site,"_",sex,"_Clean",clean,"_Green",green,"_Cross.num_with_Ran.csv",sep="")) 
  return(dataY)
}


setC.F.Ran<-Random(site=site,sex="F",clean="x",green="NA",df=lsSite[[i]]$set)
setC.M.Ran<-Random(site=site,sex="M",clean="x",green="NA",df=lsSite[[i+1]]$set)

# setG.F.Ran<-Random(site=site,sex="F",clean="NA",green="x",df=lsSite[[i+2]]$set)
# setG.M.Ran<-Random(site=site,sex="M",clean="NA",green="x",df=lsSite[[i+3]]$set)

dim(setC.F.Ran) # 110 8
dim(setC.M.Ran) # 208 8

# dim(setG.F.Ran) # 165     #/35       #165      #35
# dim(setG.M.Ran) # 62      #/44       #62       #44



# setwd("/Users/maohuang/Desktop/Kelp/Crossing_Design/RMed")
# rm(list=ls())
# 
# ##### !!!
# site<-"GOM"
# 
# ##### !!!
# site<-"SNE"
# 
# ### GOM Trial 2
# SetC.F<-read.csv(paste(site,"_F_RMed_SameRMed.csv",sep=""),sep=",",header=T,row.names=1)
# SetC.F.Ran<-SetC.F[!SetC.F$RM=="RM",]
# 
# ### SNE, GOM Trial 1
# SetC.F<-read.csv(paste(site,"_F_RMed.csv",sep=""),sep=",",header=T,row.names=1)
# SetC.F.Ran<-SetC.F  ## SNE set only
# 
# dim(SetC.F)
# head(SetC.F)
# dim(SetC.F.Ran)
# head(SetC.F.Ran)
# 
# ### GOM Trial 2
# #SetC.M<-read.csv(paste(site,"_M_RMed_SameRMed.csv",sep=""),sep=",",header=T,row.names=1)
# 
# ### SNE, GOM Trial 1
# SetC.M<-read.csv(paste(site,"_M_RMed.csv",sep=""),sep=",",header=T,row.names=1)
# 
# dim(SetC.M)
# SetC.M.Ran<-SetC.M[!SetC.M$RM=="RM",]
# dim(SetC.M.Ran)
# head(SetC.M.Ran)



#### Randome sampling within locations
setC.F.Ran<-droplevels(setC.F.Ran)
setC.M.Ran<-droplevels(setC.M.Ran)

AllwithinCrosses<-NULL
for (i in levels(setC.M.Ran$Loc)){
  F.Ran<-setC.F.Ran[setC.F.Ran$Loc==i,]
  F.Ran<-F.Ran[order(F.Ran$Random),]
  M.Ran<-setC.M.Ran[setC.M.Ran$Loc==i,]
  M.Ran<-M.Ran[order(M.Ran$Random),]
  withinCross<-cbind(as.character(F.Ran$GPID)[1],F.Ran$Random[1],as.character(M.Ran$GPID)[1],M.Ran$Random[1],as.character(M.Ran$Loc)[1])
  AllwithinCrosses<-rbind(AllwithinCrosses,withinCross)
}

AllwithinCrosses<-as.data.frame(AllwithinCrosses)
colnames(AllwithinCrosses)<-c("FG","FGRan","MG","MGRan","Loc")
for (col in c("FGRan","MGRan")){
  AllwithinCrosses[,col]<-as.numeric(as.character(AllwithinCrosses[,col]))
}

str(AllwithinCrosses)

AllwithinCrosses$FGuniq<-paste0(AllwithinCrosses$FG,"_",AllwithinCrosses$FGRan) # GP plus its randome number
AllwithinCrosses$MGuniq<-paste0(AllwithinCrosses$MG,"_",AllwithinCrosses$MGRan)

setC.F.Ran$UniqID<-paste0(setC.F.Ran$GPID,"_",setC.F.Ran$Random)
setC.M.Ran$UniqID<-paste0(setC.M.Ran$GPID,"_",setC.M.Ran$Random)

setC.F.Ran$RM<-ifelse(setC.F.Ran$UniqID%in%AllwithinCrosses$FGuniq,"RM",0)  
setC.M.Ran$RM<-ifelse(setC.M.Ran$UniqID%in%AllwithinCrosses$MGuniq,"RM",0) 

write.csv(AllwithinCrosses,paste0(site,"_WithinLocCrosses_1005.csv"))

#########################################
#Try to remove the lines that are already being used. Reupload these four files
#setC.F.Ran
#setC.M.Ran
#setG.F.Ran
#setG.M.Ran

######################################### Set Random

SetC.F.Ran<-setC.F.Ran[!setC.F.Ran$RM=="RM",]
SetC.M.Ran<-setC.M.Ran[!setC.M.Ran$RM=="RM",]


#Run20<-read.csv("TO_pick_20DiversityCrosses_StillHasBiomass_toRunBetweenLoc.csv",sep=",",header=T)


##### Do the between crosses matching

Set.small<-setC.F.Ran ### GOM C!!!   # <-SetC.F.Ran
Set.large<-setC.M.Ran ### GOM C!!!   # <-SetC.M.Ran
culture<-"clean" ### !!!

# ######################################
# ######################################
# Set.small<-setG.M.Ran ### GOM G!!!
# Set.large<-setG.F.Ran ### GOM G!!!
# culture<-"green"

Loc.num<-unique(Set.small$Loc) 
Set.small$Uniq<-paste("X",1:nrow(Set.small),sep="")
Set.large$Uniq<-paste("Y",1:nrow(Set.large),sep="")

tmp.s<-NULL
tmp.L<-NULL

Set.sample<-Set.large ### Run ONLY ONCE !!
Match.s<-NULL
Match.L<-NULL

for (i in 1:length(Loc.num)){
  
  set.s<-Set.small[Set.small$Loc%in%Loc.num[i],]
  set.L<-Set.sample[!Set.sample$Loc%in%Loc.num[i],]
  
  set.s<-set.s[order(set.s$Random),]  ### Order based on Random assigned value
  set.L<-set.L[order(set.L$Random),]  ### Order based on Random assigned value
  
  Match.s<-c(Match.s,as.character(set.s$GPID))               ### Get the set.s IDs at location i
  Match.L<-c(Match.L,as.character(set.L[1:nrow(set.s),]$GPID))   ### Match up the # of rows to it
  
  ### remove the ones already sampled
  tmp.s<-c(tmp.s,set.s$Uniq)   
  tmp.L<-c(tmp.L,set.L[1:nrow(set.s),]$Uniq)   ### the matching set.s rows from set.L
  leaveout<-c(tmp.L) 	
  Set.sample<-Set.sample[!Set.sample$Uniq%in%leaveout,]	
  
}

# Match.s.table<-subset(Set.small,GPID%in%Match.s)
# Match.s.table2<-Match.s.table[order(match(Match.s.table$GPID,Match.s)),]
# identical(Match.s.table2$GPID,Match.s)

Match.all<-matrix(nrow=nrow(Set.small),ncol=3)
Match.all[,1]<-Match.s
Match.all[,2]<-Match.L
Match.all[,3]<-c(paste(Match.s," x ",Match.L,sep=""))

str(Match.s)
str(Match.L)
str(Match.all)

write.csv(Match.all,paste(site,"_",culture,"_Match.all_",length(tmp.L),"crosses_AfterSameRMed_1005.csv",sep=""))
write.csv(Set.sample,paste(site,"_",culture,"_LeftOver.L_",nrow(Set.sample),"GPs_AfterSameRMed_1005.csv",sep=""))






#### Merge the repeat crosses to biomass

RepeatCross<-read.csv("RepeatCrosses_DwPM_SNE.csv",sep=",",header=TRUE) ## SNE!!
  dim(summary2)  
  head(summary2)
  
for (col in c("FG","MG")){
  RepeatCross[,col]<-as.character(RepeatCross[,col])
}  
summary2$GPID<-as.character(summary2$GPID)
  

#install.packages("expss")  ## VLOOKUP Vlooup vlooup
library(expss)  
RepeatCross$FGbio<-vlookup(RepeatCross$FG,dict=summary2,result_column="Cross.num",lookup_column="GPID")
RepeatCross$MGbio<-vlookup(RepeatCross$MG,dict=summary2,result_column="Cross.num",lookup_column="GPID")

write.csv(RepeatCross,paste0(site,"_RepeatCross_withCrossnum_1005.csv"))

##### Estimate and calculate statistics on # crosses being made per loation F/M
#Match.s.used<-Set.small[Set.small$GPID%in%Match.all[,1],] ### Not necessarily in the order of crossing designs
#Match.L.used<-Set.large[Set.large$GPID%in%Match.all[,2],]

#Match.s.sum<-table(Match.s.used$Loc)
#Match.L.sum<-table(Match.L.used$Loc)

#Match.both.sum<-
##### Need to write: D1: Ran one Loc M (CLEAN: GOM, M is limit; SNE, F is limit), 
##### then D2: Ran the other Loc, Match what's to the D1, remove that line from D2
##### Move to the next Loc


########## Run the 20 diversity
Run20<-read.csv("TO_pick_20DiversityCrosses_StillHasBiomass_toRunBetweenLoc.csv",sep=",",header=T)
head(Run20)

Run20F<-Run20[!is.na(Run20$FbioLeft),1:2]
colnames(Run20F)<-c("GPID","Cross.num")

Run20M<-Run20[!is.na(Run20$MbioLeft),3:4]
colnames(Run20M)<-c("GPID","Cross.num")

#setC.F.Ran<-Random(site=site,sex="F",clean="x",green="NA",df=Run20F)
#setC.M.Ran<-Random(site=site,sex="M",clean="x",green="NA",df=Run20M)

  df<-Run20F
  df<-Run20M
#library("splitstackshape")
#Random<-function(site,sex,clean,green,df){	
  ### Expand based on how many Cross.num that GP has
  new_df<-expandRows(df, "Cross.num") 
  new_df$Random<-runif(1:nrow(new_df))
  
  ### Randomly assign a number to the GPs
  dataY<-new_df
  nlines<-nrow(dataY)
  sets<-rep(1:nlines,1)               ### Make it # nlines !!!
  sample.num<-sets[order(runif(nrow(dataY)))]
  dataY$Random<-sample.num
  
  write.csv(dataY,file=paste("Cross.num_with_Run20M.csv",sep="")) 
 # return(dataY)
#}

setC.F.Ran<-dataY
setC.F.Ran$Loc<-str_split_fixed(string=setC.F.Ran$GPID,"-",4)[,2]

setC.M.Ran<-dataY
setC.M.Ran$Loc<-str_split_fixed(string=setC.M.Ran$GPID,"-",4)[,2]
