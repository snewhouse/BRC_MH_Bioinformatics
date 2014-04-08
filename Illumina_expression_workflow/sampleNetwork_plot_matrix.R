sampleNetwork_plot_matrix <- function(exprsDat, sampleDat, colBy=c("chip","group"), thrs=2 ) {
  
  gp_col <- sampleDat[,colBy]
  cat(" setting up data for qc plots","\r","\n")
  ## expression matrix and IAC
  cat(" expression matrix and IAC","\r","\n")
  datExprs <- exprsDat
  IAC <- cor(datExprs)
  IAC_d <- 1-IAC
  samle_names <- colnames(datExprs)
  IAC=cor(datExprs, method="p",use="p")
  diag(IAC)=0
  A.IAC=((1+IAC)/2)^2  ## ADJACENCY MATRIX
  cat(" fundamentalNetworkConcepts","\r","\n")
  FNC=fundamentalNetworkConcepts(A.IAC) ## WGCNA
  K2=FNC$ScaledConnectivity
  Z.K=(K2-mean(K2))/sd(K2)
  Z.C=(FNC$ClusterCoef-mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef)
  rho <- signif(cor.test(Z.K,Z.C,method="s")$estimate,2)
  rho_pvalue <- signif(cor.test(Z.K,Z.C,method="s")$p.value,2)
  Z.K_outliers <- Z.K < -thrs
  Z.K_score <- Z.K[Z.K_outliers==TRUE]
  Z.K_outliers <- colnames(datExprs)[Z.K_outliers==TRUE]
  # set colours
  cat(" colorvec [",paste(colBy),"]","\r","\n")
  
  colorvec <- labels2colors(as.character(gp_col))
  
  mean_IAC <- mean(IAC[upper.tri(IAC)])
  
  ## samplenetwork
  local( 
{colLab <<- function(n,treeorder) { 
  if(is.leaf(n)) {
    a <- attributes(n)
    i <<- i+1
    attr(n, "nodePar") <-   c(a$nodePar, list(lab.col = colorvec[treeorder][i], lab.font = i%%3))
  } 
  n 
} 
i <- 0
})
cat(" begin SampleNetwork plots","\r","\n")
## Cluster for pics
cluster1 <- hclust(as.dist(1-A.IAC),method="average")
cluster1order <- cluster1$order
cluster2 <- as.dendrogram(cluster1,hang=0.1)
cluster3 <- dendrapply(cluster2,colLab,cluster1order)
## PLOTS
## cluster IAC
par(mfrow=c(2,2))
par(mar=c(5,6,4,2))
plot(cluster3,nodePar=list(lab.cex=1,pch=NA),
     main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),
     xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
mtext(paste("distance: 1 - ISA ",sep=""),cex=0.8,line=0.2)
## Connectivity
par(mar=c(5,5,4,2))
plot(Z.K,main="Connectivity", ylab="Z.K",xaxt="n",xlab="Sample",type="n",cex.main=1.8,cex.lab=1.4)
text(Z.K,labels=samle_names,cex=0.8,col=colorvec)
abline(h=-2)
abline(h=-3)
## ClusterCoef
par(mar=c(5,5,4,2))
plot(Z.C,main="ClusterCoef", ylab="Z.C",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4,type="n")
text(Z.C,labels=samle_names,cex=0.8,col=colorvec)
abline(h=-2)
abline(h=-3)
## Connectivity vs ClusterCoef
par(mar=c(5,5,4,2))
plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=colorvec,cex.main=1.8,cex.lab=1.4)
abline(lm(Z.C~Z.K),col="black",lwd=2)
mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
abline(v=-2,lty=2,col="grey")
abline(h=-2,lty=2,col="grey")
res <- list(Z.K_outliers=Z.K_outliers, Z.K_score=Z.K_score)
# return(Z.K_outliers, return(Z.K_score)
return(res)
}