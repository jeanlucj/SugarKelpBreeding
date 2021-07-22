
WD<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003/"
#WD<-"/local/workdir/mh865/GCA_SCA/OneTime192021/"  # run in terminal

datafdr<-paste0(WD,"data/")

  pop<-"Yr19to21"
#save(WithinYr_Both_dBLUPs,file=paste0(datafdr,"dBLUPs_3_yrs_cal_WithinYr_",pop,".Rdata"))
 load(paste0(datafdr,"dBLUPs_3_yrs_cal_WithinYr_",pop,".Rdata"))
##### 
Y<-WithinYr_Both_dBLUPs
Y$Year<-as.factor(WithinYr_Both_dBLUPs$Year.x)
Y$Crosses<-WithinYr_Both_dBLUPs$Crosses.x
library(ggplot2)
plot<-ggplot(data=Y,aes(dryWgtPerM,Year))+
  geom_point(aes(color=as.factor(Year)))+ 
  geom_line(aes(group=as.factor(Crosses)))
print(plot)
identify(x=Y$dryWgtPerM,y=Y$Year.x)
Y[Y$dryWgtPerM==max(Y[Y$Year==2021,]$dryWgtPerM),]

 # SA18-CB-10-FG6xSL18-NC-2-MG1
load(paste0(datafdr,"dataNHpi_withChk_3Yrs_PhotoScore123_07152021.rdata"))   ## Plot
  ls()
  dataNHpi3yrs_C[dataNHpi3yrs_C$Crosses=="SA18-CB-10-FG6xSL18-NC-2-MG1",]
#####
