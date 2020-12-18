setwd("/Users/maohuang/Desktop/Other ppls project/Stagonospora/")

ls()
load("GAPIT_StagTP464lines_AA,AC,CC.Rdata")


write.csv(myG,"myG_Oct01.csv",row.names=FALSE)
write.csv(myKI,"myKI_Oct01.csv",row.names=FALSE)
write.csv(myPC,"myPC_Oct01.csv",row.names=FALSE)
write.csv(myY,"myY_Oct01.csv",row.names=FALSE)