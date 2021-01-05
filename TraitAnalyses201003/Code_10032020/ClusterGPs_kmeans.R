#Clustering GPS

library(magrittr)
library(ClusterR)
hMat <- readRDS(file="hMat_analyzeNH.rds")
gpMat <- hMat[196:312, 196:312] ## !!
gpEig <- eigen(gpMat, symmetric=T)
clustDat <- gpEig$vectors %*% diag(sqrt(gpEig$values))
tst <- KMeans_rcpp(clustDat, 8, num_init=100)
plot(gpEig$vectors[,1:2], col=tst$clusters, pch=16)
xlim <- range(gpEig$vectors[,1]); ylim <- range(gpEig$vectors[,2])
op <- par(mfrow=c(2, 2))
for (clust in 1:4){
  plot(gpEig$vectors[tst$clusters == clust,1:2], pch=16, xlim=xlim, ylim=ylim)
}
for (clust in 5:8){
  plot(gpEig$vectors[tst$clusters == clust,1:2], pch=16, xlim=xlim, ylim=ylim)
}
par(op)
