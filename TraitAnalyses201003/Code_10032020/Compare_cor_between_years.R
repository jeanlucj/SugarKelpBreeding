## evaluate cor between 2019_2020 data VS 2021 data
## Phenotypic Selection cor

WD<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003/"

datafdr<-paste0(WD,"data/")
pop<-"Yr1920to21"
load(paste0(datafdr,"dBLUPs_3_yrs_cal_WithinYr_",pop,".Rdata"))
## Note, the year.x and year.y is not correct for Yr209 nor 2020
yrs3<-WithinYr_Both_dBLUPs
  dim(yrs3)  #341 
X1920<-yrs3$Crosses.x[!yrs3$Year.x==2021]
X21<-yrs3$Crosses.x[yrs3$Year.x==2021]

comX<-X21[X21%in%X1920]
  
X1920_comX<-yrs3[yrs3$Crosses.x%in%comX &!(yrs3$Year.x==2021),]
X21_comX<-yrs3[yrs3$Crosses.x%in%comX & (yrs3$Year.x==2021),]
  identical(X1920_comX$Crosses.x,X21_comX$Crosses.x)
X21_comX<-X21_comX[order(match(X21_comX$Crosses.x,X1920_comX$Crosses.x)),] 
  identical(X1920_comX$Crosses.x,X21_comX$Crosses.x)
  
traits<-c("wetWgtPerM", "percDryWgt", "dryWgtPerM", "densityBlades","bladeLength","bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter")
cor<-NULL
for (j in 1:length(traits)){
  cor<-c(cor,cor(X1920_comX[,traits[j]],X21_comX[,traits[j]],use="complete"))
}
  cor
  names(cor)<-traits
  # wetWgtPerM     percDryWgt     dryWgtPerM  densityBlades    bladeLength  bladeMaxWidth bladeThickness 
  # 0.006236214   -0.080638326    0.004768809    0.389343450    0.041159072    0.627555522    0.065313101 
  # stipeLength  stipeDiameter 
  # 0.006479333    0.340854886 
  
### Within Yr dBLUPs, get the  cor within each year
  
pop<-"Yr19to21"
load(paste0(datafdr,"dBLUPs_3_yrs_cal_WithinYr_",pop,".Rdata"))
  dim(WithinYr_Both_dBLUPs)  #346
  
  
traits<-c("dryWgtPerM","wetWgtPerM", "percDryWgt","densityBlades","bladeLength","bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter")
newColName<-c("DWpM","WWpM","pDW","DB","BL","BmWid","Bth","SL","SDia")

for (col in traits){
  colnames(WithinYr_Both_dBLUPs)[col]
}


years<-c(2019,2020,2021)
cor<-matrix(nrow=3,ncol=length(traits))

for (i in 1:length(years)){
  yr<-years[i]
  data<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==yr,traits]
    print(dim(data))
    head(data)
## RM NAs if there are any
  data<-na.omit(data)
    print(dim(data))
  for (j in 1:length(traits)){
   colnames(data)[colnames(data)==traits[j]]<-newColName[j]
  }
    dev.set(i)
    pdf(paste0(datafdr,"CorrPlot_SPplots",yr,"_PhotoScore23_dBLUPs.pdf"))
    corrplot::corrplot.mixed(cor(data), diag="n", tl.cex=0.6, tl.col=1)
    dev.off()
}  


### Hist of BV for three years' data
rm(list=ls())
TraitWant<-"DWpM"

WD<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003/"
datafdr<-paste0(WD,"data/")

forPlot<-read.csv(paste0(datafdr,"allBLUPs_3yrs_all_scenarios_sort_",TraitWant,".csv"))

years<-c(2019,2020,2021,"Two","Three")
  #329 invidiuals of SPs, all these will be having an estimated SPs
#for (i in 1:length(years)){
i=5
  data<-forPlot[,colnames(forPlot)%in%c("X",paste0(TraitWant,"_",years[i]))]
    print(dim(data))
    print(head(data))
  data$Trait<-data[,paste0(TraitWant,"_",years[i])]  
  library(ggplot2)
  
histo<-ggplot(data=data,aes(x=Trait))+
  geom_histogram(binwidth=0.01)+
  labs(x="Trait GEBV",y="Plot Count")+
  theme_bw()+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"))
#}


load(paste0(datafdr,"dataNHpi_withChk_3Yrs_PhotoScore123_07152021.rdata"))  
  ls(dataNHpi3yrs_C)

data2019<-droplevels(dataNHpi3yrs_C[dataNHpi3yrs_C$Year==2019 &(!dataNHpi3yrs_C$PhotoScore==1),])
  dim(data2019) #139
data2020<-droplevels(dataNHpi3yrs_C[dataNHpi3yrs_C$Year==2020 &(!dataNHpi3yrs_C$PhotoScore==1),])
  dim(data2020) #144
data2021<-droplevels(dataNHpi3yrs_C[dataNHpi3yrs_C$Year==2021 &(!dataNHpi3yrs_C$PhotoScore==1),])
  dim(data2021) #95
  
  head(data2019)
  head(data2020)
  head(data2021)
Cross2019<-unique(data2019$Crosses) 
Cross2020<-unique(data2020$Crosses)
Cross2021<-unique(data2021$Crosses)
  str(Cross2019)   #124
  str(Cross2020)  #129
  str(Cross2021)  #87

  head(data)
  dim(data)
data$Cross2019<-data$X%in%Cross2019
data$Cross2020<-data$X%in%Cross2020
data$Cross2021<-data$X%in%Cross2021

segment_data1 = data.frame(
  x = c(data$Trait[data$Cross2019==TRUE]),
  xend = c(data$Trait[data$Cross2019==TRUE]), 
  y = rep(0,length(data$Trait[data$Cross2019==TRUE])),
  yend = rep(1, length(data$Trait[data$Cross2019==TRUE]))
)

segment_data2 = data.frame(
  x = c(data$Trait[data$Cross2020==TRUE]),
  xend = c(data$Trait[data$Cross2020==TRUE]), 
  y = rep(0,length(data$Trait[data$Cross2020==TRUE])),
  yend = rep(2, length(data$Trait[data$Cross2020==TRUE]))
)

segment_data3 = data.frame(
  x = c(data$Trait[data$Cross2021==TRUE]),
  xend = c(data$Trait[data$Cross2021==TRUE]), 
  y = rep(0,length(data$Trait[data$Cross2021==TRUE])),
  yend = rep(3, length(data$Trait[data$Cross2021==TRUE]))
)

  dim(segment_data3)
  head(segment_data)
histo +
  geom_segment(data=segment_data1,aes(x = x, y = y, xend = xend, yend = yend,color="green"),linetype="dashed",show.legend=TRUE)+
geom_segment(data=segment_data2,aes(x = x, y = y, xend = xend, yend = yend,color="red"),linetype="dashed",show.legend=TRUE)+
  geom_segment(data=segment_data3,aes(x = x, y = y, xend = xend, yend = yend,color="blue"),linetype="dashed",show.legend=TRUE)+
  scale_color_manual(name = "",
                     values = c( "green","red","blue"),
                     labels = c("Crosses in Yr2019", "Crosses in Yr2020","Crosses in Yr2021"))


data[data$Cross2021=="TRUE"&data$Cross2020=="TRUE",]

CrossRepeat<-data[data$Trait>0.2,]
  dim(CrossRepeat)
write.csv(CrossRepeat,paste0(datafdr,"Crosses_Worth_Repeating_3_yearsBVanalysis.csv")) 

CrossRepeatRaw<-dataNHpi3yrs_C[dataNHpi3yrs_C$Crosses%in%CrossRepeat$X,]

yrs3<-yrs3[order(yrs3$Year.x,-yrs3$dryWgtPerM),]

for (Yr in c(2019,2020,2021)){
  yrs3[yrs3$Year.x==Yr,]$order<-1:sum(yrs3$Year.x==Yr)
}

write.csv(yrs3,paste0(datafdr,"yrs3_ordered.csv"))
low10<-c("SL18-CC-8-FG1xSL18-JS-9-MG3","SL18-JS-9-FG3xSL18-NL-2-MG1","SL18-JS-16-FG3xSL18-JS-9-MG3","SL18-JS-9-FG3xSL18-JS-5-MG2","SL18-CC-8-FG1xSL18-NL-2-MG1","SL18-CC-8-FG1xSL18-NC-8-MG2","SL18-JS-16-FG3xSL18-NC-8-MG2")

yrs3[yrs3$Crosses.x%in%low10,]
