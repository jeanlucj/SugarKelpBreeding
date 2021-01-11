## Testing lines and block effects
library(rrBLUP)
rm(list=ls())

#load("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/dataNHpi_withChk_3_sets_PhotoScore23.rdata")   ## Plot
#load("dataNHim_withChk_3_sets_PhotoScore0123.rdata")  ## Indi


dataAll <-read.csv("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/dataNHpi_ManuallyAddGrowth_9.29.20_Ireland_Data_AshFDWpM.csv",sep=",",header=TRUE) ##!!!

dataAll$AshFDwPM<-(dataAll$wetWgtPerM*dataAll$percDryWgt/100)*(1-(dataAll$Ash/100))
dataAll$popChk <- ifelse(substr(dataAll$plotNo, 1, 1) == "Z", substr(dataAll$plotNo, 1, 2), "ES")  # Checks VS ES

# exptlSP <- as.character(dataAll$popChk) == "ES"
# dataAll$entry <- as.character(dataAll$popChk)
# dataAll$entry[exptlSP] <- as.character(dataAll$plotNo[exptlSP])
# dataAll$group <- as.factor(ifelse(exptlSP, 1, 0))

ls()
head(dataAll)
dataAll[1:10,]

########################### A. Plot level data input, estimate pedNH
dataNHpi<-dataAll[dataAll$Region=="GOM",] ## !!! RM SNE, 530 rows

dataNHpi <- dataNHpi[order(dataNHpi$plotNo),]  #### Plots in alphabetic order
dim(dataNHpi)
str(dataNHpi)    
dataNHpi<-dataNHpi[!dataNHpi$crossID=="Buffer",]  ## !!! RMed Buffer lines

dataNHpi<-dataNHpi[!dataNHpi$PhotoScore==0,] ## !!! RM PhotoScore=0,  447 rows
#dataNHpi<-dataNHpi[dataNHpi$PhotoScore>1,] ## !!! RM PhotoScore<=1
dim(dataNHpi)
colnames(dataNHpi)
dataNHpi19<-dataNHpi[dataNHpi$Year==2019,] ## 2019 data   !!!!!!
dataNHpi19_C<-dataNHpi19
dataNHpi19_C<-dataNHpi19_C[order(dataNHpi19_C$plotNo),]  ## Order plotNo alphabetically


dataNHpi<-dataNHpi19_C  ##!!!!!

library(lme4)
exptlSP <- as.character(dataNHpi$popChk) == "ES"
dataNHpi$entry <- as.character(dataNHpi$popChk)
dataNHpi$entry[exptlSP] <- as.character(dataNHpi$plotNo[exptlSP])
dataNHpi$group <- as.factor(ifelse(exptlSP, 1, 0))
print("Wet Weight Per Plot")

library(lmerTest) # can give p-values

 for (col in c( "line", "block","popChk","group","entry")) 
  dataNHpi[,col] <- factor(dataNHpi[,col])

fitAug <- lmer(dryWgtPerM ~ popChk + line*block + (1|entry:group), data=dataNHpi)
print(aov <- anova(fitAug))
#> anova(fitAug)
#Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq   Mean Sq NumDF   DenDF F value Pr(>F)
#popChk     0.008459 0.0042294     2  10.812  0.4603 0.6429
#line       0.029086 0.0096952     3 143.641  1.0552 0.3702
#block      0.049054 0.0061317     8 139.696  0.6674 0.7195
#line:block 0.168385 0.0070161    24 129.949  0.7636 0.7751

fitAug <- lmer(wetWgtPlot ~ popChk + line*block + (1|entry:group), data=dataNHpi) ##!!!popChk and group is identical
print(aov <- anova(fitAug))
#Type III Analysis of Variance Table with Satterthwaite's method
#            Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
#popChk      0.7268 0.36341     2  11.322  0.2536 0.7803
#line        3.6786 1.22620     3 144.801  0.8558 0.4657
#block       6.4044 0.80055     8 140.540  0.5587 0.8101
#line:block 29.4904 1.22877    24 129.486  0.8576 0.6575


dataNHim <- read.csv("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Indi_2019_2020_Compile_07202020_Edit_CorrectplotNo.csv",sep=",",header=TRUE)
dataNHim<-dataNHim[dataNHim$Region=="GOM",] ## !!! RM SNE, 530 rows
  dim(dataNHim)  #5046

dataNHim<-dataNHim[!dataNHim$crossID=="Buffer",]
  dim(dataNHim) #4996

dataNHim19_C<-dataNHim[dataNHim$Year==2019,] 
dataNHim19_C<-dataNHim19_C[order(dataNHim19_C$plotNo),]  ### Plots in Alphabetic order; Here the plot did not RM plots with photoScore < 2,3
  dim(dataNHim19_C)  # 2800 x 12
  tail(dataNHim19_C)
  colnames(dataNHim19_C)
  
data<-dataNHim19_C  ### !!!!!!
  data$popChk <- ifelse(substr(data$plotNo, 1, 1) == "Z", substr(data$plotNo, 1, 2), "ES") # if it contains Z, be it, otherwise ES
  
#### vlookup the line and block based on plotNo
  library(expss)
  data$line<-vlookup(data$plotNo,dict=dataNHpi,result_column="line",lookup_column="plotNo")  
  data$block<-vlookup(data$plotNo,dict=dataNHpi,result_column="block",lookup_column="plotNo")
  
  exptlSP <- as.character(data$popChk) == "ES"
  data$entry <- as.character(data$popChk)
  data$entry[exptlSP] <- as.character(data$plotNo[exptlSP])
  data$group <- as.factor(ifelse(exptlSP, 1, 0))

  library(lmerTest) # can give p-values
  
  for (col in c( "line", "block","popChk","group","entry")) 
    data[,col] <- factor(data[,col])
  
  for (col in c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")){
    data[,col]<-as.numeric(data[,col])
  }
  
  fitAug <- lmer(stipeLength ~ popChk + line*block + (1|entry:group), data=data)
  print(aov <- anova(fitAug))
#Type III Analysis of Variance Table with Satterthwaite's method
#  Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
#  popChk       3128    1564     2  68.043  0.1147 0.891784   
 # line       199434   66478     3 152.679  4.8770 0.002882 **
  #  block      184648   23081     8 155.523  1.6933 0.103985   
  #line:block 318458   13269    24 161.595  0.9735 0.503850   
  
  fitAug <- lmer(stipeDiameter ~ popChk + line*block + (1|entry:group), data=data)
  print(aov <- anova(fitAug))
#Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
#popChk        228     114     2  98.077  0.0111    0.9890    
#line       321655  107218     3 164.219 10.4246 2.573e-06 ***
#block       59707    7463     8 167.541  0.7256    0.6687    
#line:block 347478   14478    24 184.603  1.4077    0.1078 
  
  fitAug <- lmer(bladeLength ~ popChk + line*block + (1|entry:group), data=data)
  print(aov <- anova(fitAug))
#Type III Analysis of Variance Table with Satterthwaite's method
#  Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
#  popChk      137531   68766     2  18.723  0.4597 0.6384
#  line        273300   91100     3 169.147  0.6091 0.6100
#  block      1986050  248256     8 172.561  1.6598 0.1114
#  line:block 3427438  142810    24 165.719  0.9548 0.5284  
  
  fitAug <- lmer(bladeMaxWidth ~ popChk + line*block + (1|entry:group), data=data)
  print(aov <- anova(fitAug))
#Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
#popChk      17544  8772.0     2  55.924  0.8402 0.43700  
#line        28523  9507.6     3 161.124  0.9107 0.43726  
#block       27021  3377.7     8 163.844  0.3235 0.95617  
#line:block 468098 19504.1    24 166.228  1.8682 0.01214 *  
  
  fitAug <- lmer(bladeThickness ~ popChk + line*block + (1|entry:group), data=data)
  print(aov <- anova(fitAug))

#Type III Analysis of Variance Table with Satterthwaite's method
#  Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
#  popChk        28.5   14.23     2  84.309  0.0457 0.9553370    
 # line        5789.4 1929.79     3 148.506  6.2008 0.0005357 ***
  #  block       5300.6  662.57     8 151.574  2.1290 0.0363016 *  
   # line:block 13136.3  547.35    24 164.859  1.7587 0.0213345 * 
    