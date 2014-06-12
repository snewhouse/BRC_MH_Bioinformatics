library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(affy)
library(cluster)
library(impute)
library(WGCNA)
allowWGCNAThreads(nThreads=6)
library(gplots)
library(limma)
library(vsn)
library(MBCB)
library(lumiHumanIDMapping)
library(scatterplot3d)
library(relaimpo)
library(plyr)
library(ggplot2)
library(gdata)
library(sva)
library(pamr)
library(glmnet)
library(snow)
library(parallel)
library(foreach)
library(doParallel)
library(caret)

####
def.par <- par(no.readonly = TRUE)

dev.off()
par(def.par)


############################
## subset eset
############################

removeSamples_eset_lumi <- function(eset, sampleRemove) {
sample <- sampleNames(eset)
samples_to_remove <- sampleRemove
samples_to_keep <- (sample %in% samples_to_remove)==FALSE
sel_samp_names <- sampleNames(eset)[samples_to_keep]
eset <- eset[,samples_to_keep]
ControlData <- getControlData(eset)
ControlData  <- ControlData[,c("controlType","ProbeID",sel_samp_names)]
eset <- addControlData2lumi(ControlData , eset)
return(eset)
}

############################
## shuffle_cols
############################
shuffle_cols <- function(data) {
data_var <- names(data)
data_var_shuffled <- sample(data_var,size=length(data_var),replace=FALSE)
return(data_var_shuffled)
}
############################
## shuffle_rows
############################
shuffle_rows <- function(data) {
data_var <- rownames(data)
data_var_shuffled <- sample(data_var,size=length(data_var),replace=FALSE)
return(data_var_shuffled)
}

############################
## data_summary_plots
############################
data_summary_plots <- function(data,results_prefix) {
##library(car)
##library(psych)
##library(mi)
if(require("car")){
    print("car is loaded correctly")
} else {
    print("trying to install car")
    install.packages("car")
    if(require(car)){
        print("car installed and loaded")
    } else {
        stop("could not install car")
    }
}
##
if(require("psych")){
    print("psych is loaded correctly")
} else {
    print("trying to install psych")
    install.packages("psych")
    if(require(psych)){
        print("psych installed and loaded")
    } else {
        stop("could not install psych")
    }
}
if(require("mi")){
    print("mi is loaded correctly")
} else {
    print("trying to install mi")
    install.packages("mi")
    if(require(mi)){
        print("mi installed and loaded")
    } else {
        stop("could not install mi")
    }
}
### get data

## res fil names ##
plotfile <- paste(results_prefix,".data_summary_plots.pdf",sep="")
resfile_num <- paste(results_prefix,".numerical_data_summary.csv",sep="")
resfile_cat <- paste(results_prefix,".categorical_data_summary.csv",sep="")

## save default, for resetting...c(5, 4, 4, 2) + 0.1 ## c(bottom, left, top, right)
def.par <- par(no.readonly = TRUE)

## missing.data.analysis
cat(" missing.data.analysis ","\r","\n")
missing_TF <- is.na(data)
missing_counts <- rowSums(missing_TF)
n_sample_with_missing_data <- sum(missing_counts>=1)
n_sample_with_complete_data <- sum(missing_counts==0)
data_with_missing_var <- data[missing_counts>=1,]
data_with_complete_var <- data[missing_counts==0,]
cat(" missing.data.plot ","\r","\n")
cat(" missing.pattern.plot ","\r","\n")
pdf(file=paste(results_prefix,".missing.pattern.plot.pdf",sep=""),width=11,height=8)
par(mar=c(4.1, 12.1, 1.1, 2.1))
missing.pattern.plot ( data, y.order = FALSE, x.order = FALSE,obs.col="white", mis.col="black" )
par(def.par)
dev.off()

## data classes
data_class <- sapply(data ,class)
class_list <- unique(data_class)
cat(" The following data classes are observed [",unique(data_class),"]","\r","\n")

## get numerical data and summarize
sel_numerical <- sapply(data ,class) %in% c("numeric","integer")
num_data <- data[,sel_numerical]
## write to csv
sink(resfile_num)
cat("VARNAME,N,NMISS,PROP_MISS,MEAN,SD,VAR,MEDIAN,IQR,MAD,MIN,MAX,SHAPIRO_WILKS_W,SHAPIRO_WILKS_P,SKEW,KURTOSIS","\n")
for(myvar in names(num_data) ){
X <- num_data[,myvar]
name_var <- paste(myvar,sep="")
mean_var <- signif(mean(X,na.rm=TRUE),3)
sd_var <- signif(sd(X,na.rm=TRUE),3)
Var_var <- signif(var(X,na.rm=TRUE),3)
med_var <- median(X, na.rm = FALSE)
iqr_var <- IQR(X, na.rm = FALSE, type = 7)
mad_var <- mad(X, center = median(X), constant = 1.4826, na.rm = FALSE,low = FALSE, high = FALSE)
range_var_min <- range(X, na.rm = FALSE)[1]
range_var_max <- range(X, na.rm = FALSE)[2]
normtest_sw_p <- signif(shapiro.test(X)$p.value,3)
normtest_sw_w <- signif(shapiro.test(X)$statistic,3)
skew_var <- skew(X, na.rm = TRUE)
kurtosi_var <- kurtosi(X, na.rm = TRUE)
nmiss <- sum(is.na(X))
nsamp <- length(X)
pecentage_missing <- round(nmiss/nsamp,3)
res <- paste(name_var,nsamp,nmiss,pecentage_missing,mean_var,sd_var,Var_var,med_var,iqr_var,mad_var,range_var_min,range_var_max,normtest_sw_w,normtest_sw_p,skew_var,kurtosi_var,sep=",")
cat(res,"\n")
}
sink()

## get categorical data and make tables of counts
sel_categorical <- sapply(data ,class) %in% c("character","factor")
categorica_data <- as.data.frame(data[,sel_categorical])
colnames(categorica_data) <- names(data)[sel_categorical]
sink(resfile_cat)
for(myvar in names(categorica_data) ){
X <- categorica_data[,myvar]
name_var <- paste(myvar,sep="")
table_var <- table(X)
nmiss <- sum(is.na(X))
nsamp <- length(X)
pecentage_missing <- round(nmiss/nsamp,3)
###observed_var <- paste(unique(X),collapse=",")
observed_var <- paste(X,collapse=",")
SINGLUARITY <- min(table_var)==1; ## any var seen only once??
min_count <- min(table_var)
cat("VARNAME,",name_var,"\n",sep="")
cat("N,",nsamp,"\n",sep="")
cat("NMISS,",nmiss,"\n",sep="")
cat("PROP_MISS,",pecentage_missing,"\n",sep="")
cat("MIN_COUNT,",min_count,"\n",sep="")
cat("SINGLUARITY,",SINGLUARITY,"\n",sep="")
cat("OBSERVED_TERMS,",observed_var,"\n",sep="")
cat("COUNTS,",paste(table_var,collapse=",",sep=""),"\n","\n",sep="")
}
sink()

#################
## Start plots ##
#################
pdf(file=plotfile,width=11,height=8)

par(mar=c(4.1, 12.1, 1.1, 2.1))
###missing.pattern.plot(data, y.order = FALSE, x.order = FALSE,obs.col="white", mis.col="black" )
par(def.par)

for( class_type in class_list ) {

cat(" doing ",class_type,"\r","\n")

	      if(class_type =="character") {
					new_data <- data[,data_class==class_type]
					    for(myvar in names(new_data) ){
					    var_table <- sort(table(new_data[,myvar]),decreasing=FALSE)
					    par(mar=c(5.1, 8, 4, 2.1))
					    bp <- barplot(var_table,las=2,main=paste(myvar), horiz=TRUE,col="light blue")
					    text(0,bp,round(var_table, 1),cex=1,pos=4)
						    }
	      }
	      else if(class_type=="factor") {
					      new_data <- data[,data_class==class_type]
					      for(myvar in names(new_data) ){
					      var_table <- sort(table(new_data[,myvar]),decreasing=FALSE)
					      par(mar=c(5.1, 8, 4, 2.1))
					      bp <- barplot(var_table,las=2,main=paste(myvar), horiz=TRUE,col="light blue")
					      text(0,bp,round(var_table, 1),cex=1,pos=4)
	      	            	      }
	      }
	      else if(class_type=="numeric") {
							  new_data <- data[,data_class==class_type]
							  for(myvar in names(new_data) ){
							  nf <- layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,3))
							  par(mar=c(4.1, 4.1, 1.1, 2.1))
							  X <- as.numeric(new_data[,myvar])
							  boxplot(X, horizontal=TRUE,  outline=TRUE,main=paste(myvar))
							  hist(X,xlab=paste(myvar),breaks=50,main="",prob=TRUE)
							  lines(density(X),col="blue")
							  meanvar <- signif(mean(X,na.rm=TRUE),3)
							  sdvar <- signif(sd(X,na.rm=TRUE),3)
							  normtest_sw <- signif(shapiro.test(X)$p.value,3)
							  nmiss <- sum(is.na(X))
							  pecentage_missing <- round(nmiss/length(X),3)
							  mtext(paste("\nMean:[",meanvar,"]\nSD:[",sdvar,"]\nnormP:[",normtest_sw,"]\n NSAMPLE:[",length(X),"]\nNMISS:[",nmiss,"]\n %Missing:[",pecentage_missing,"]",sep=" "), side = 3, adj=1, padj=1)
							  par(def.par) # 1=bottom, 2=left, 3=top, 4=right
							  qqPlot(X,main=paste(myvar),pch=20,ylab=paste(myvar),col="blue")
							  							  }
      	      }
      	      else if(class_type=="integer") {
							  new_data <- data[,data_class==class_type]
							  new_data <- apply(new_data,2,as.numeric)
							  for(myvar in names(new_data) ){
							  nf <- layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,3))
							  par(mar=c(4.1, 4.1, 1.1, 2.1))
							  X <- as.numeric(new_data[,myvar])
							  boxplot(X, horizontal=TRUE,  outline=TRUE,main=paste(myvar))
							  hist(X,xlab=paste(myvar),breaks=50,main="",prob=TRUE)
							  lines(density(X),col="blue")
							  meanvar <- signif(mean(X,na.rm=TRUE),3)
							  sdvar <- signif(sd(X,na.rm=TRUE),3)
							  normtest_sw <- signif(shapiro.test(X)$p.value,3)
							   nmiss <- sum(is.na(X))
							  pecentage_missing <- round(nmiss/length(X),3)
							  mtext(paste("\nMean:[",meanvar,"]\nSD:[",sdvar,"]\nnormP:[",normtest_sw,"]\n NSAMPLE:[",length(X),"]\nNMISS:[",nmiss,"]\n %Missing:[",pecentage_missing,"]",sep=" "), side = 3, adj=1, padj=1)
							  par(def.par)
							  qqPlot(X,main=paste(myvar),pch=20,ylab=paste(myvar),col="blue")
							  }
	      }
	      else if(class_type=="logical") {
					      new_data <- data[,data_class==class_type]
					      for(myvar in names(new_data) ){
					      var_table <- sort(table(new_data[,myvar]),decreasing=FALSE)
					      par(mar=c(5.1, 8, 4, 2.1))
					      bp <-  barplot(var_table,las=2,main=paste(myvar), horiz=TRUE,col="light blue")
					      text(0,bp,round(var_table, 1),cex=1,pos=4)
					      par(def.par)
	      	      }
	      }
	  }
	  dev.off()
	  	  dev.off()
	  	  	  dev.off()
	  	  	  	  dev.off()
	  	  	  	  dev.off()
	  }

############################
## negBeadOutlierRepMean
#############################
negBeadOutlierRepMean <- function(x) {
		z_out_samp <- abs(  as.numeric( scale(x) )  ) > 2
		mean_pop <- mean(x[z_out_samp==FALSE])
		sd_pop <- sd(x[z_out_samp==FALSE])
		new_x <- ifelse( abs(  as.numeric( scale(x) )  ) > 2, mean_pop, x )
		return(new_x)
}
##############
## quantfun ##
##############
quantfun <- function(x,probs=c(seq(0,1,0.25))) {
# changed from function(x,probs=c(seq(0,1,0.20))) to function(x,probs=c(seq(0,1,0.25)))
as.integer(cut(x, quantile(x, probs ), include.lowest=TRUE))
}
##############
## var_gene ##
##############
zero_var_probe <- function(gx_matrix) {
gx <- gx_matrix
var_gx <- apply(gx,1,var)
zero_var <- var_gx==0
return(zero_var)
}
mean_probe <- function(gx_matrix) {
gx <- gx_matrix
mean_gx <- apply(gx,1,mean)
return(mean_gx)
}
sd_probe <- function(gx_matrix) {
gx <- gx_matrix
sd_gx <- apply(gx,1,sd)
return(sd_gx)
}
var_probe <- function(gx_matrix) {
gx <- gx_matrix
var_gx <- apply(gx,1,var)
return(var_gx)
}
max_probe <- function(gx_matrix) {
gx <- gx_matrix
max_gx <- apply(gx,1,max)
return(max_gx)
}
min_probe <- function(gx_matrix) {
gx <- gx_matrix
min_gx <- apply(gx,1,min)
return(min_gx)
}
has_var_probe <- function(gx_matrix) {
gx <- gx_matrix
var_gx <- apply(gx,1,var)
has_var <- var_gx!=0
return(has_var)
}
has_var_probe2 <- function(gx_matrix) {
gx <- gx_matrix
var_gx <- apply(gx,1,var)
has_var <- var_gx >= quantile(var_gx, probs="0.1" );
return(has_var)
}
############################
## write_expression_files ##
############################
write_expression_files <- function(eset, outfile) {

cat(" Writing probe exprs matrix [", paste(outfile,".exprs_matrix.txt",sep="")  ,"]","\r","\n")
gx <- exprs(eset)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(outfile,".exprs_matrix.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing probe se.exprs matrix [", paste(outfile,".se.exprs_matrix.txt",sep="")  ,"]","\r","\n")
gx <- se.exprs(eset)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(outfile,".se.exprs_matrix.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing probe detection matrix [", paste(outfile,".detection_matrix.txt",sep="")  ,"]","\r","\n")
gx <- detection(eset)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(outfile,".detection_matrix.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing probe beadNum matrix [", paste(outfile,".beadNum_matrix.txt",sep="")  ,"]","\r","\n")
gx <- beadNum(eset)
gx <- as.data.frame(gx)
gx$nuID <- rownames(gx)
gx <- merge(fData(eset),gx,by.x="nuID",by.y="nuID",sort=FALSE)
write.table(gx , file=paste(outfile,".beadNum_matrix.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing probe PCA matrix [", paste(outfile,".pca_matrix.txt",sep="")  ,"]","\r","\n")
gx <- exprs(eset)
gx <- gx[zero_var_probe(gx)==FALSE,]
pca_gx <- prcomp(t(gx))$x
pca_gx <- pca_gx[,1:20]
pca_gx <- cbind(rownames(pca_gx),pca_gx)
grep_pc <- grep("PC",colnames(pca_gx))
colnames(pca_gx) <-   c("sampleID",colnames(pca_gx)[grep_pc])
write.table(pca_gx  , file=paste(outfile,".pca_matrix.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing pData slot of eset and adding PCA data to [", paste(outfile,".pData.txt",sep="")  ,"]","\r","\n")
pgx <- pData(eset)
pgx <- merge(pgx, pca_gx, by.x="sampleID",by.y="sampleID",sort=FALSE,all.x=TRUE)
write.table(pgx  , file=paste(outfile,".pData.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

cat(" Writing fData slot of eset [", paste(outfile,".fData.txt",sep="")  ,"]","\r","\n")
fgx <- fData(eset)
sel_probe <- paste(fgx$nuID,sep="")
gx <- exprs(eset)
gx <- gx[sel_probe,]
fgx$mean_probe <- mean_probe(gx);
fgx$sd_probe <- sd_probe(gx);
fgx$var_probe <- var_probe(gx);
fgx$min_probe <- min_probe(gx);
fgx$max_probe <- max_probe(gx);
write.table(fgx  , file=paste(outfile,".fData.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)
rm("gx","pca_gx","pgx","fgx")
}

##############
## qc_plots ##
##############


gx_qc_plots_lumi <- function(eset, outfile ,do_pca=TRUE ) {
  outpdf <- paste(outfile,".qc_plots.pdf",sep="");
  cat(" startin qc plots","\r","\n")
  cat(" saving all plots to ",outpdf,"\r","\n")
  ## get pheno data
  pheno <- pData(eset)
  ## basic colours
  chip_col <- labels2colors( as.character(pData(eset)$Sentrix.Barcode))
  group_col <- labels2colors( as.character(pData(eset)$GROUPS))
  pheno_col <- labels2colors( as.character(pData(eset)$PHENOTYPE))
  gender_col <- labels2colors( as.character(pData(eset)$SEX))
  tissue_col <- labels2colors( as.character(pData(eset)$TISSUE))
  ## batch pheno data
  sel_tech <- grep("tech",names(pheno))
  batch_pheno <- pheno[,sel_tech]
  ## id what is char/fav versus numerical
  sel_batch <- sapply(batch_pheno ,class) %in% c("character","factor")
  sel_num_batch <- sapply(batch_pheno ,class) %in% c("numeric")
  ## get names
  batch_var_names <- names( batch_pheno[sel_batch])
  batch_var_names_numeric <- names(batch_pheno[sel_num_batch])
  ## quantfun numeric data
  numeric_tech <- as.data.frame(batch_pheno[,sel_num_batch])
  colnames(numeric_tech) <- batch_var_names_numeric
  quant_numeric <- apply( numeric_tech,2,quantfun)
  ## colours
  batch_col <- apply(batch_pheno[,sel_batch],2,labels2colors) ## colours
  quant_numeric_col <- apply(quant_numeric,2,numbers2colors)
  datColors <- cbind(chip_col,group_col,pheno_col,gender_col,tissue_col,batch_col,quant_numeric_col)
  ## expression matrix and IAC
  gx <- t(exprs(eset));
  datExprs <- exprs(eset)
  IAC <- cor(datExprs)
  IAC_d <- 1-IAC
  #
  #
  samle_names <- sampleNames(eset)
  IAC=cor(datExprs, method="p",use="p")
  diag(IAC)=0
  A.IAC=((1+IAC)/2)^2  ## ADJACENCY MATRIX
  FNC=fundamentalNetworkConcepts(A.IAC) ## WGCNA
  K2=FNC$ScaledConnectivity
  Z.K=(K2-mean(K2))/sd(K2)
  Z.C=(FNC$ClusterCoef-mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef)
  Z.MAR=(FNC$MAR-mean(FNC$MAR))/sd(FNC$MAR)
  rho <- signif(cor.test(Z.K,Z.C,method="s")$estimate,2)
  rho_pvalue <- signif(cor.test(Z.K,Z.C,method="s")$p.value,2)
  colorvec <- labels2colors(as.character(pData(eset)$Sentrix.Barcode))
  mean_IAC <- mean(IAC[upper.tri(IAC)])
  ## flashClust
  dist_exprs <- dist( t(datExprs), method="e" )
  sampleTree <- flashClust( dist_exprs, method = "average");
  ## PLOTS
  if(do_pca==TRUE) {
    pca_raw <- prcomp(gx)$x;
    pdf(outpdf, width=16.5,height=11.7)
    ## Standard plots
    plot(eset, what='boxplot', col=chip_col )
    plot(eset, what='density' )
    plot(eset, what='cv'  )
    ## PCA plots
    scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by chip ",color="black", pch=21,bg=chip_col)
    scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by Group ",color="black", pch=21,bg=group_col)
    scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by Phenotype ",color="black", pch=21,bg=pheno_col)
    scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by Gender ",color="black", pch=21,bg=gender_col)
    scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main="3D Scatterplot coloured by Gender ",color="black", pch=21,bg=tissue_col)
    ## loop PCA and plot by tech var
    for(tech_var in batch_var_names) {
      tech_var_name <- paste(tech_var,sep="")
      tech_var_col <- labels2colors( as.character(batch_pheno[,tech_var_name]) )
      scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main=paste("3D Scatterplot coloured by ",tech_var_name,sep=""),color="black",pch=21,bg=tech_var_col)
    }
    for(tech_var in batch_var_names_numeric) {
      tech_var_name <- paste(tech_var,sep="")
      tech_var_col <- numbers2colors( quant_numeric[,tech_var_name])
      scatterplot3d(pca_raw[,"PC1"],pca_raw[,"PC2"],pca_raw[,"PC3"], main=paste("3D Scatterplot coloured by ",tech_var_name,sep=""),color="black",bg=tech_var_col,pch=21)
    }
    ## other plots
    #tmp_iac <- outlierSamples(datExprs,thresh=iac_sd_thrs, showplots=TRUE)
    plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,dendroLabels = pData(eset)$Sample.ID, main="Sample dendrogram and trait heatmap")
    heatmap.2(IAC , trace="none", ColSideColors=labels2colors( as.character(pData(eset)$PHENOTYPE)), RowSideColors=labels2colors( as.character(pData(eset)$GROUPS)), col="bluered", main="IAC",hclustfun=function(x) hclust(x,method='complete'), distfun= function(x) dist(x,method='euclidean') )
    #
    # heatmap A.IAC
    heatmap.2(A.IAC , trace="none", ColSideColors=labels2colors( as.character(pData(eset)$PHENOTYPE)), RowSideColors=labels2colors( as.character(pData(eset)$GROUPS)), col="bluered", main="A.IAC",hclustfun=function(x) hclust(x,method='complete'), distfun= function(x) dist(x,method='euclidean') )
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
    ## Cluster for pics
    cluster1 <- hclust(as.dist(1-A.IAC),method="average")
    cluster1order <- cluster1$order
    cluster2 <- as.dendrogram(cluster1,hang=0.1)
    cluster3 <- dendrapply(cluster2,colLab,cluster1order)
    ## PLOTS
    ## cluster IAC
    par(mfrow=c(2,2))
    par(mar=c(5,6,4,2))
    plot(cluster3,nodePar=list(lab.cex=1,pch=NA),main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
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
    dev.off()
  } else {
    pdf(outpdf, width=16.5,height=11.7)
    plot(eset, what='boxplot', col=chip_col )
    plot(eset, what='density' )
    plot(eset, what='cv'  )
    #tmp_iac <- outlierSamples(datExprs,thresh=iac_sd_thrs, showplots=TRUE)
    plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,dendroLabels = pData(eset)$Sample.ID, main="Sample dendrogram and trait heatmap")
    heatmap.2(IAC , trace="none", ColSideColors=labels2colors( as.character(pData(eset)$PHENOTYPE)), RowSideColors=labels2colors( as.character(pData(eset)$GROUPS)), col="bluered", main="IAC",hclustfun=function(x) hclust(x,method='complete'), distfun= function(x) dist(x,method='euclidean') )
    dev.off()
  }
}





#################
## qc_plots_2 ##
##  no 3D scatter plot, no save pdf
##############

# basic_qc_plot_lumi
basic_qc_plot_lumi <- function(eset) {
  
  ## flashClust
  cat(" Running flashClust","\r","\n")
  chip_col <- labels2colors( as.character(pData(eset)$Sentrix.Barcode))
  datExprs <- exprs(eset)
  dist_exprs <- dist( t(datExprs), method="e" )
  sampleTree <- flashClust( dist_exprs, method = "average");
  ## Standard plots
  cat(" beging plotting boxplot","\r","\n")
  plot(eset, what='boxplot', col=chip_col )
#  dev.off()
 # par(def.par)
  cat(" beging plotting outlier","\r","\n")
  plot(eset, what='outlier'  )
#  dev.off()
 # par(def.par)
  cat(" beging plotting sampleTree <- flashClust( dist_exprs, method = \"average\")","\r","\n")
  plot(sampleTree)
#  dev.off()
  #par(def.par)
  cat(" beging plotting density","\r","\n")
  plot(eset, what='density' )
#  dev.off()
  #par(def.par)
  cat(" beging plotting cv","\r","\n")
  plot(eset, what='cv'  )
#  dev.off()
  #par(def.par)
    
}

## pca_plot_lumi
pca_plot_lumi <- function(eset, probe_list) {
  
  cat(" setting up data for qc plots","\r","\n")
  ## get pheno data
  cat(" get pheno data","\r","\n")
  sel_probes <- probe_list
  pheno <- pData(eset)
  ## basic colours
  cat(" basic colours","\r","\n")
  chip_col <- labels2colors( as.character(pData(eset)$Sentrix.Barcode))
  group_col <- labels2colors( as.character(pData(eset)$GROUPS))
  pheno_col <- labels2colors( as.character(pData(eset)$PHENOTYPE))
  gender_col <- labels2colors( as.character(pData(eset)$SEX))
  tissue_col <- labels2colors( as.character(pData(eset)$TISSUE))
  ## batch pheno data
  cat(" batch pheno data","\r","\n")
  sel_tech <- grep("tech",names(pheno))
  batch_pheno <- pheno[,sel_tech]
  ## id what is char/fav versus numerical
  sel_batch <- sapply(batch_pheno ,class) %in% c("character","factor")
  sel_num_batch <- sapply(batch_pheno ,class) %in% c("numeric")
  ## get names
  batch_var_names <- names( batch_pheno[sel_batch])
  batch_var_names_numeric <- names(batch_pheno[sel_num_batch])
  ## quantfun numeric data
  numeric_tech <- as.data.frame(batch_pheno[,sel_num_batch])
  colnames(numeric_tech) <- batch_var_names_numeric
  quant_numeric <- apply( numeric_tech,2,quantfun)
  ## colours
  batch_col <- apply(batch_pheno[,sel_batch],2,labels2colors) ## colours
  quant_numeric_col <- apply(quant_numeric,2,numbers2colors)
  datColors <- cbind(chip_col,group_col,pheno_col,gender_col,tissue_col,batch_col,quant_numeric_col)
  #
  cat(" begin PCA plots","\r","\n")
  cat(" calculating PCs using prcomp()","\r","\n")
  def.par <- par(no.readonly = TRUE)
  gx <- t(exprs(eset[sel_probes,]));
  pca_gx <- prcomp(gx);
  pca_summary <- summary(pca_gx)
  pca_importance_var_exp <- summary(pca_gx)$importance[2,]
  pca_importance_var_exp_cum <- summary(pca_gx)$importance[3,]
  cat(" plotting Proportion of Variance Explained","\r","\n") 
  par(mfrow=c(1,2))
  plot(pca_importance_var_exp[1:20],ylab="PCA Proportion of Variance Explained", type="b",col="blue")
  plot(pca_importance_var_exp_cum[1:20],ylab="PCA Cumulative Proportion of Variance Explained", ylim=c(0,1.1),type="b",col="blue");abline(h=0.90);abline(h=1.00)
  par(def.par)
  #dev.off()
  # PCA MATRIX
  pca_raw <- pca_gx$x;
  ## PCA plots
  plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=" PCA plot coloured by chip ",col="black", pch=21,bg=chip_col)
  plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=" PCA plot coloured by Group ",col="black", pch=21,bg=group_col)
  plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=" PCA plot coloured by Phenotype ",col="black", pch=21,bg=pheno_col)
  plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=" PCA plot coloured by Gender ",col="black", pch=21,bg=gender_col)
  plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=" PCA plot coloured by Tissue ",col="black", pch=21,bg=tissue_col)
  ## loop PCA and plot by tech var
  for(tech_var in batch_var_names) {
    cat(" begin looping through batch variable PCA plots ",tech_var,"\r","\n")
    tech_var_name <- paste(tech_var,sep="")
    tech_var_col <- labels2colors( as.character(batch_pheno[,tech_var_name]) )
    plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=paste(" PCA plot coloured by ",tech_var_name,sep=""),col="black",pch=21,bg=tech_var_col)
  }
  for(tech_var in batch_var_names_numeric) {
    cat(" begin looping through batch variable PCA plots ",tech_var,"\r","\n")
    tech_var_name <- paste(tech_var,sep="")
    tech_var_col <- numbers2colors( quant_numeric[,tech_var_name])
    plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=paste(" PCA plot coloured by ",tech_var_name,sep=""),col="black",bg=tech_var_col,pch=21)
  }
}

# coloured_dendrogram_lumi

coloured_dendrogram_lumi <- function(eset){
  ## get pheno data
  cat(" get pheno data","\r","\n")
  pheno <- pData(eset)
  ## basic colours
  cat(" basic colours","\r","\n")
  chip_col <- labels2colors( as.character(pData(eset)$Sentrix.Barcode))
  group_col <- labels2colors( as.character(pData(eset)$GROUPS))
  pheno_col <- labels2colors( as.character(pData(eset)$PHENOTYPE))
  gender_col <- labels2colors( as.character(pData(eset)$SEX))
  tissue_col <- labels2colors( as.character(pData(eset)$TISSUE))
  ## batch pheno data
  cat(" batch pheno data","\r","\n")
  sel_tech <- grep("tech",names(pheno))
  batch_pheno <- pheno[,sel_tech]
  ## id what is char/fav versus numerical
  sel_batch <- sapply(batch_pheno ,class) %in% c("character","factor","Date")
  sel_num_batch <- sapply(batch_pheno ,class) %in% c("numeric")
  ## get names
  batch_var_names <- names( batch_pheno[sel_batch])
  batch_var_names_numeric <- names(batch_pheno[sel_num_batch])
  ## quantfun numeric data
  numeric_tech <- as.data.frame(batch_pheno[,sel_num_batch])
  colnames(numeric_tech) <- batch_var_names_numeric
  quant_numeric <- apply( numeric_tech,2,quantfun)
  ## colours
  cat(" batch colours","\r","\n")
  batch_col <- apply(batch_pheno[,sel_batch],2,labels2colors) ## colours
  quant_numeric_col <- apply(quant_numeric,2,numbers2colors)
  datColors <- cbind(chip_col,group_col,pheno_col,gender_col,tissue_col,batch_col,quant_numeric_col)
  # sampleTree
  cat(" datExprs and sampleTree","\r","\n")
  datExprs <- exprs(eset)
  dist_exprs <- dist( t(datExprs), method="e" )
  sampleTree <- flashClust( dist_exprs, method = "average");
  # plotDendroAndColors
  cat(" plotDendroAndColors","\r","\n")
  #def.par <- par(no.readonly = TRUE)
  par(mar=c(5.1, 20, 4, 2.1)) #bottom, left, top and right margins
  plotDendroAndColors(sampleTree,
                      groupLabels=names(datColors),
                      colors=datColors,
                      dendroLabels = pData(eset)$Sample.ID,marAll = c(1, 10, 3, 1), 
                      main="Sample dendrogram and trait heatmap")
  #par(def.par)
  
}


# heatmap_plot_lumi_eset_raw
heatmap_plot_lumi <- function(eset){
  
  cat(" setting up data for qc plots","\r","\n")
  ## get pheno data
  cat(" get pheno data","\r","\n")
  pheno <- pData(eset)
  ## basic colours
  cat(" basic colours","\r","\n")
  chip_col <- labels2colors( as.character(pData(eset)$Sentrix.Barcode))
  group_col <- labels2colors( as.character(pData(eset)$GROUPS))
  pheno_col <- labels2colors( as.character(pData(eset)$PHENOTYPE))
  gender_col <- labels2colors( as.character(pData(eset)$SEX))
  tissue_col <- labels2colors( as.character(pData(eset)$TISSUE))
  ## batch pheno data
  cat(" batch pheno data","\r","\n")
  sel_tech <- grep("tech",names(pheno))
  batch_pheno <- pheno[,sel_tech]
  ## id what is char/fav versus numerical
  sel_batch <- sapply(batch_pheno ,class) %in% c("character","factor")
  sel_num_batch <- sapply(batch_pheno ,class) %in% c("numeric")
  ## get names
  batch_var_names <- names( batch_pheno[sel_batch])
  batch_var_names_numeric <- names(batch_pheno[sel_num_batch])
  ## quantfun numeric data
  numeric_tech <- as.data.frame(batch_pheno[,sel_num_batch])
  colnames(numeric_tech) <- batch_var_names_numeric
  quant_numeric <- apply( numeric_tech,2,quantfun)
  ## colours
  batch_col <- apply(batch_pheno[,sel_batch],2,labels2colors) ## colours
  quant_numeric_col <- apply(quant_numeric,2,numbers2colors)
  datColors <- cbind(chip_col,group_col,pheno_col,gender_col,tissue_col,batch_col,quant_numeric_col)
  ## expression matrix and IAC
  cat(" expression matrix and IAC","\r","\n")
  gx <- t(exprs(eset));
  datExprs <- exprs(eset)
  IAC <- cor(datExprs)
  IAC_d <- 1-IAC
  samle_names <- sampleNames(eset)
  IAC=cor(datExprs, method="p",use="p")
  diag(IAC)=0
  A.IAC=((1+IAC)/2)^2  ## ADJACENCY MATRIX
  
  #  heatmap IAC
  heatmap.2(IAC , trace="none", 
            ColSideColors=labels2colors( as.character(pData(eset)$PHENOTYPE)), 
            RowSideColors=labels2colors( as.character(pData(eset)$GROUPS)), 
            col="bluered", main="IAC",
            hclustfun=function(x) hclust(x,method='complete'), 
            distfun= function(x) dist(x,method='euclidean') )
  
  # heatmap A.IAC
  heatmap.2(A.IAC , trace="none", 
            ColSideColors=labels2colors( as.character(pData(eset)$PHENOTYPE)),
            RowSideColors=labels2colors( as.character(pData(eset)$GROUPS)), 
            col="bluered", main="A.IAC",
            hclustfun=function(x) hclust(x,method='complete'), 
            distfun= function(x) dist(x,method='euclidean') )
  
}


# sampleNetwork_plot_all_lumi

sampleNetwork_plot_all_lumi <- function(eset, colBy=c("chip","group") ) {
  
  gp_col <- colBy
  cat(" setting up data for qc plots","\r","\n")
  ## expression matrix and IAC
  cat(" expression matrix and IAC","\r","\n")
  datExprs <- exprs(eset)
  IAC <- cor(datExprs)
  IAC_d <- 1-IAC
  samle_names <- sampleNames(eset)
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
  Z.K_outliers <- Z.K < -3
  Z.K_score <- Z.K[Z.K_outliers==TRUE]
  Z.K_outliers <- colnames(datExprs)[Z.K_outliers==TRUE]
# set colours
  cat(" colorvec [",paste(gp_col),"]","\r","\n")
if(gp_col=="chip") { colorvec <- labels2colors(as.character(pData(eset)$Sentrix.Barcode)) }
if(gp_col=="group") { colorvec <- labels2colors(as.character(pData(eset)$GROUPS)) }
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

# gx_qc_plots_lumi_2
gx_qc_plots_lumi_2 <- function(eset, do_pca=TRUE ) {
  cat(" starting qc plots","\r","\n")
  cat(" setting up data for qc plots","\r","\n")
  ## get pheno data
  cat(" get pheno data","\r","\n")
  pheno <- pData(eset)
  ## basic colours
  cat(" basic colours","\r","\n")
  chip_col <- labels2colors( as.character(pData(eset)$Sentrix.Barcode))
  group_col <- labels2colors( as.character(pData(eset)$GROUPS))
  pheno_col <- labels2colors( as.character(pData(eset)$PHENOTYPE))
  gender_col <- labels2colors( as.character(pData(eset)$SEX))
  tissue_col <- labels2colors( as.character(pData(eset)$TISSUE))
  ## batch pheno data
  cat(" batch pheno data","\r","\n")
  sel_tech <- grep("tech",names(pheno))
  batch_pheno <- pheno[,sel_tech]
  ## id what is char/fav versus numerical
  sel_batch <- sapply(batch_pheno ,class) %in% c("character","factor")
  sel_num_batch <- sapply(batch_pheno ,class) %in% c("numeric")
  ## get names
  batch_var_names <- names( batch_pheno[sel_batch])
  batch_var_names_numeric <- names(batch_pheno[sel_num_batch])
  ## quantfun numeric data
  numeric_tech <- as.data.frame(batch_pheno[,sel_num_batch])
  colnames(numeric_tech) <- batch_var_names_numeric
  quant_numeric <- apply( numeric_tech,2,quantfun)
  ## colours
  batch_col <- apply(batch_pheno[,sel_batch],2,labels2colors) ## colours
  quant_numeric_col <- apply(quant_numeric,2,numbers2colors)
  datColors <- cbind(chip_col,group_col,pheno_col,gender_col,tissue_col,batch_col,quant_numeric_col)
  ## expression matrix and IAC
  cat(" expression matrix and IAC","\r","\n")
  gx <- t(exprs(eset));
  datExprs <- exprs(eset)
  IAC <- cor(datExprs)
  IAC_d <- 1-IAC
  samle_names <- sampleNames(eset)
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
  colorvec <- labels2colors(as.character(pData(eset)$Sentrix.Barcode))
  mean_IAC <- mean(IAC[upper.tri(IAC)])
  ## flashClust
  cat(" flashClust","\r","\n")
  dist_exprs <- dist( t(datExprs), method="e" )
  sampleTree <- flashClust( dist_exprs, method = "average");
  ## Standard plots
  cat(" beging plotting boxplot","\r","\n")
  plot(eset, what='boxplot', col=chip_col )
  cat(" beging plotting density","\r","\n")
  plot(eset, what='density' )
  cat(" beging plotting cv","\r","\n")
  plot(eset, what='cv'  )
  cat(" beging plotting outlier","\r","\n")
  plot(eset, what='outlier'  )
  cat(" beging plotting sampleTree","\r","\n")
  plot(sampleTree)
  if(do_pca==TRUE) {
    cat(" begin PCA plots","\r","\n")
    pca_raw <- prcomp(gx)$x;
    ## PCA plots
    plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=" PCA plot coloured by chip ",col="black", pch=21,bg=chip_col)
    plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=" PCA plot  coloured by Group ",col="black", pch=21,bg=group_col)
    plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=" PCA plot  coloured by Phenotype ",col="black", pch=21,bg=pheno_col)
    plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=" PCA plot  coloured by Gender ",col="black", pch=21,bg=gender_col)
    plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=" PCA plot  coloured by Tissue ",col="black", pch=21,bg=tissue_col)
    ## loop PCA and plot by tech var
    for(tech_var in batch_var_names) {
      cat(" begin looping through batch variable PCA plots ",tech_var,"\r","\n")
      tech_var_name <- paste(tech_var,sep="")
      tech_var_col <- labels2colors( as.character(batch_pheno[,tech_var_name]) )
      plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=paste(" PCA plot coloured by ",tech_var_name,sep=""),col="black",pch=21,bg=tech_var_col)
    }
    for(tech_var in batch_var_names_numeric) {
      cat(" begin looping through batch variable PCA plots ",tech_var,"\r","\n")
      tech_var_name <- paste(tech_var,sep="")
      tech_var_col <- numbers2colors( quant_numeric[,tech_var_name])
      plot(pca_raw[,"PC1"],pca_raw[,"PC2"], main=paste(" PCA plot coloured by ",tech_var_name,sep=""),col="black",bg=tech_var_col,pch=21)
    }
  }
    ## other plots
    #tmp_iac <- outlierSamples(datExprs,thresh=iac_sd_thrs, showplots=TRUE)
    cat(" begin plotDendroAndColors and  heatmap.2 plots","\r","\n")
    plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,dendroLabels = pData(eset)$Sample.ID, main="Sample dendrogram and trait heatmap")
    heatmap.2(IAC , trace="none", ColSideColors=labels2colors( as.character(pData(eset)$PHENOTYPE)), RowSideColors=labels2colors( as.character(pData(eset)$GROUPS)), col="bluered", main="IAC",hclustfun=function(x) hclust(x,method='complete'), distfun= function(x) dist(x,method='euclidean') )
    # heatmap A.IAC
    heatmap.2(A.IAC , trace="none", ColSideColors=labels2colors( as.character(pData(eset)$PHENOTYPE)), RowSideColors=labels2colors( as.character(pData(eset)$GROUPS)), col="bluered", main="A.IAC",hclustfun=function(x) hclust(x,method='complete'), distfun= function(x) dist(x,method='euclidean') )
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
    plot(cluster3,nodePar=list(lab.cex=1,pch=NA),main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
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
  }





####################
## bgcor_mbcb ######
####################
bgcor_mbcb <- function(eset, outfile) {
  
cat(" Start background correction ", "\r","\n")

sig <- paste(outfile,".exprsSignal.txt",sep="")
neg <- paste(outfile,".negative.control.txt",sep="")

cat(" get negativeControl bead data ","\r","\n")

negativeControl <- getControlData(eset)
negativeControl <- subset(negativeControl, negativeControl$controlType=="NEGATIVE")
negativeControl_orig <- negativeControl[,c(3:dim(negativeControl)[2])]

#############################################################
## removing outlier neg beads and replace with mean
#############################################################

negBeadOutlierRepMean2 <- function(x) {
    z_out_samp <- abs(  as.numeric( scale(x) )  ) > 2
    mean_pop <- mean(x[z_out_samp==FALSE])
    sd_pop <- sd(x[z_out_samp==FALSE])
    new_x <- ifelse( abs(  as.numeric( scale(x) )  ) > 2, mean_pop, x )
    return(new_x)
  }

negativeControl <- apply(negativeControl_orig ,2, negBeadOutlierRepMean2)

neg_max <- apply(negativeControl,2,max)
neg_sd <- apply(negativeControl,2,sd)
neg_mean <- apply(negativeControl,2,mean)

#############################################################
## save negativeControl bead expression data
#############################################################
write.table(negativeControl, file=neg, sep="\t");

#############################################################
## get expression data
#############################################################
cat(" get expression bead data ","\r","\n")

expressionSignal <- exprs(eset)

#############################################################
## save expression data
#############################################################
write.table(expressionSignal, file=sig, sep="\t");

#############################################################
## set data for mbcb
#############################################################
cat(" set data for mbcb ","\r","\n")

data <- mbcb.parseFile(sig, neg);

signal <- data$sig;
negCon <- data$con;

#############################################################
## run mbcb
#############################################################
cat(" run background correct using mbcb.correct(method=\"MLE\")  ","\r","\n")

gx_mbcb <- mbcb.correct(signal,negCon,npBool=FALSE,mleBool=TRUE, isRawBead=FALSE)

###############
## save mbcb ##
###############
cat(" mbcb complete ","\r","\n")
cat(" saveing raw mbcb matrix to ",paste(outfile,".mbcb.correct.output.RData",sep=""), "\r","\n")
save(gx_mbcb , file=paste(outfile,".mbcb.correct.output.RData",sep=""))

################################
## select mbcb method results ##
################################
cat(" Selected mbcb method: ",mbcb_method,"\r","\n")
##if(mbcb_method=="NP") { gx_mbcb <- gx_mbcb$NP }
if(mbcb_method=="MLE") { gx_mbcb <- gx_mbcb$MLE }

#############################################################
## replace names with original as R adds "X" to numbers #####
#############################################################
cat(" replace names with original sampleNames(eset_raw), as R adds X to numbers ","\r","\n")
colnames(gx_mbcb) <- sampleNames(eset)
gx_mbcb[1:10,1:10]
#############################################################
## make new eset for bk corrected data
#############################################################
cat(" Creating Background Corrected Data set: eset_bg  ","\r","\n")
eset_bg <- eset
#############################################################
## replace old exprs data with new mbcb bk corrected data
#############################################################
cat(" replace old exprs data with new mbcb Background corrected data ","\r","\n")
exprs(eset_bg) <- as.matrix(gx_mbcb)
sampleNames(eset_bg) <- sampleNames(eset)
sampleNames(eset_bg)
## return 
cat(" returning new mbcb Background corrected data: eset_bg ","\r","\n")
return(eset_bg)
}




##################
##lumi.N <- lumiExpresso(eset_raw, variance.stabilize=TRUE,  varianceStabilize.param=list(method='vst'),  normalize.param=list(method='rsn'), bg.correct=FALSE)


# Austin Hilliard, White lab UCLA, Sep 2009
#
# group of functions for removing outlier probes and samples in microarray data
#
# removeOutlierProbes removes probes outside of specified stdev range from mean, runs on probesets (rows) or samples (cols) of input data 
#
# removeOutlierProbesIterate runs removeOutlierProbes iteratively until no outliers remain
# 
# removeTooManyNAs looks for probesets (rows) and samples (cols) with more than a specified number of missing values and removes them
# 
# outlierSamples computes inter-sample correlations and performs hierarchical clustering to find sample outliers as taught by Mike Oldham (formerly of the Geschwind lab), 
# written prior to MO giving me beta version of his RemoveOutliers function.
### MO's function RemoveOutliers function does this iteratively, as well as running ComBat, with user interactivity
#
# outlierSamplesIterate runs outlierSamples iteratively until user quits or chooses not to remove any more samples
#
# preProc runs all of the above and saves their output lists, including each iteration of processed data
#
#########################################################################

### libraries
library(lattice)
library(Biobase)
library(affy)
library(limma)
library(vsn)
library(preprocessCore)
library(lumi)
library(annotate)
library(lumiHumanAll.db)
library(affy)
library(cluster)
library(impute)
library(WGCNA)
library(gplots)
library(limma)
library(MBCB)
library(lumiHumanIDMapping)
library(scatterplot3d)
########################################

########################################

# utility functions

cv = function(x) {
  # compute coefficient of variation
  stdev = sd(as.numeric(x), na.rm=TRUE);
  avg = mean(as.numeric(x), na.rm=TRUE);
  cv = stdev/avg;
  return(cv)
}

# run vsn
runVSN = function(data, plot=T) {
  if (!is.matrix(data)) {data = as.matrix(data)};
  dataSet = new("ExpressionSet", exprs = data);
  dataSetVSN = justvsn(dataSet);
  dataVSN = exprs(dataSetVSN);
  dataVSN = as.data.frame(dataVSN);
  
  if (plot) {par(mfrow=c(1,2)); meanSdPlot(data); meanSdPlot(dataSetVSN)};
  
  return(dataVSN);
  
}
##########################################################################

removeOutlierProbes = function(data, deviate=3, rowORcol=1) {
  
  # deviate is number of stdevs away from mean probe must surpass to be outlier.
  # 	could be assessing distribution of single probe across samples,
  #	or distribution of all probes on single sample
  #
  # rowORcol chooses whether computation is over probesets(1) or samples(2)
  # 	if 1, looking for single probe outliers within probeset across samples
  #	if 2, looking for outliers within samples
  
  if(rowORcol==1){label=rownames(data)}; 
  if(rowORcol==2){label=names(data)};
  
  data_dim = dim(data);
  
  # set up data frame to hold mean, sd, cv, and cut-offs for each probeset/sample
  data_stats = data.frame(mean=numeric(data_dim[rowORcol]),
                          sd=numeric(data_dim[rowORcol]),
                          cv=numeric(data_dim[rowORcol]),
                          up_thresh=numeric(data_dim[rowORcol]), 
                          low_thresh=numeric(data_dim[rowORcol]), 
                          row.names=label 
  );
  
  # compute mean, sd, cv, and cut-offs for each probeset/sample
  data_stats[,1] = apply(data, rowORcol, mean, na.rm=TRUE);
  data_stats[,2] = apply(data, rowORcol, sd, na.rm=TRUE);
  data_stats[,3] = data_stats[,2] / data_stats[,1];
  data_stats[,4] = data_stats[,1] + (deviate * data_stats[,2]);
  data_stats[,5] = data_stats[,1] - (deviate * data_stats[,2]);  
  
  # get indices and values of outliers
  outlier_positions = which(data > data_stats[,4] | data < data_stats[,5], arr.ind=T); 
  outlier_vals = data[outlier_positions];
  outlier_vals = cbind(outlier_positions,outlier_vals);
  
  # get number of outliers and compute percentage of total probes they represent
  total_outliers = dim(outlier_positions)[1];
  total_probes = dim(data)[1] * dim(data)[2];
  percent_outliers = 100 * total_outliers / total_probes;
  
  # make copy of data and insert NAs for outliers
  dataClean = data;
  dataClean[data > data_stats[,4] | data < data_stats[,5]] = NA; 
  
  out = list(data_stats=data_stats,
             deviate=deviate,
             outlier_positions=outlier_positions,
             outlier_vals=outlier_vals,
             total_probes=total_probes,
             total_outliers=total_outliers,
             percent_outliers=percent_outliers,
             dataClean=dataClean
  );
  return(out);
  
}

########################################################

removeOutlierProbesIterate = function(data, deviate=3, rowORcol=1) {
  
  if(rowORcol==1){cat('Removing probes >',deviate,'stdevs from probeset mean...\n')}; 
  if(rowORcol==2){cat('Removing probes >',deviate,'stdevs from sample mean...\n')};
  
  single_round_outliers = 1;
  outlier_running_count = 0;
  iteration = 0;
  total_probes = dim(data)[1] * dim(data)[2];
  outlier_positions = list();
  outlier_vals = list(); 
  outliers_by_round = c();
  
  while (single_round_outliers > 0) {
    
    out = removeOutlierProbes(data, deviate, rowORcol);
    single_round_outliers = out$total_outliers;
    outlier_running_count = outlier_running_count + out$total_outliers;
    data = out$dataClean;
    iteration = iteration + 1;
    outlier_positions[[iteration]] = out$outlier_positions;
    outlier_vals[[iteration]] = out$outlier_vals; 
    outliers_by_round[iteration] = single_round_outliers;
    cat('Round',iteration, 'outliers:', single_round_outliers, '\n');
    
  };
  
  names(outlier_positions) = paste('round', c(1:iteration), sep='');
  names(outlier_vals) = paste('round', c(1:iteration), sep='');
  names(outliers_by_round) = paste('round', c(1:iteration), sep='');
  percent_outliers = 100 * outlier_running_count / total_probes;
  dataClean = out$dataClean;
  
  cat('\n');
  cat('Total outliers:', outlier_running_count, '\n');
  cat('Percentage of probes that were outliers:', percent_outliers, '\n');
  
  output = list(dataClean=dataClean, 
                outlier_positions=outlier_positions,
                outlier_vals=outlier_vals,
                outliers_by_round=outliers_by_round,
                total_outliers=outlier_running_count,
                total_probes=total_probes,
                percent_outliers=percent_outliers, 
                sd_cutoff=deviate
  );
  return(output);
  
}

#######################################################

removeTooManyNAs = function (data, probe_thresh=NULL, sample_thresh=NULL) {
  
  if(is.numeric(probe_thresh)){
    probe_thresh=probe_thresh;
  } else {
    probe_thresh = floor(ncol(data)/2);
  }
  cat('\nremoving probes with >', probe_thresh, ' missing measurements\n', sep='');
  
  if(is.numeric(sample_thresh)){
    sample_thresh=sample_thresh;
  } else {
    sample_thresh = floor(nrow(data)/2);
  }
  cat('removing samples with >', sample_thresh, ' missing measurements\n', sep='');
  
  countNAs = apply(is.na(data), 1, sum);
  probes_over_thresh = (1:dim(data)[1])[countNAs > probe_thresh];
  cat('\n');
  cat(length(probes_over_thresh),'probes removed \n');
  
  countNAsSamps = apply(is.na(data), 2, sum);
  samples_over_thresh = (1:dim(data)[2])[countNAsSamps > sample_thresh]; 
  cat(length(samples_over_thresh),'samples removed \n');
  
  dataClean = data;
  if (length(probes_over_thresh) > 0) {dataClean = data[-probes_over_thresh, ] };
  if (length(samples_over_thresh) > 0) {dataClean = data[ , -samples_over_thresh] }; 	
  names(probes_over_thresh) = names(countNAs)[probes_over_thresh];
  names(samples_over_thresh) = names(countNAsSamps)[samples_over_thresh];
  
  output = list(dataClean=dataClean,
                countNAs=countNAs, 
                countNAsSamps=countNAsSamps,
                probes_over_thresh=probes_over_thresh,
                samples_over_thresh=samples_over_thresh,
                probe_thresh=probe_thresh,  
                sample_thresh=sample_thresh 
  );
  
  return(output);
  
}

##########################################################

outlierSamples = function(data, thresh=2, showplots=T) {
  
  if (thresh < 0) {thresh = -thresh};
  
  IAC=cor(data,method="p",use="complete.obs");
  clust=hclust(as.dist(1-IAC),method="average");
  
  meanIAC=mean(IAC[upper.tri(IAC)]);
  meanIACdiag=mean(IAC);
  samplemeanIAC=apply(IAC,2,mean);
  sdCorr=sd(samplemeanIAC);    
  numbersd=(samplemeanIAC-meanIACdiag)/sdCorr;
  
  cat('\n');
  cat('Looking for outlier samples (>',thresh,'stdevs from meanIACdiag)...\n');
  cat('meanIAC =',meanIAC,'\n');
  cat('meanIACdiag =',meanIACdiag,'\n');
  cat('\n') 
  
  over_thresh = numbersd < -thresh | numbersd > thresh; 
  samples_to_remove = numbersd[over_thresh];
  #formatted = data.frame(z.IAC=samples_to_remove);
  
  dataClean = data[, !over_thresh]; 
  
  cat('All samples z.IAC: \n');
  print(numbersd);
  cat('\n\n');
  if (length(samples_to_remove)!=0) {
    cat('Possible outliers: \n')
    #print(formatted);
    print(samples_to_remove);
    cat('\n');
  } else {
    cat('No samples >',thresh,'stdevs from meanIACdiag \n');
  };
  
  output = list(dataClean=dataClean,
                IAC=IAC,
                meanIAC=meanIAC,
                meanIACdiag=meanIACdiag,
                samplemeanIAC=samplemeanIAC,
                numbersd=numbersd,
                clust=clust,
                samples_to_remove=samples_to_remove
  );
  
  if (showplots) {
    par(mfrow=c(1,2)); 
    plot(clust); 
    plot(numbersd, type='n'); 
    text(numbersd, labels=names(numbersd), cex=0.75);
  };
  
  return(output);
  
}

##########################################################

outlierSamplesIterate = function (data, IACthresh=2, showplots=T) {
  
  if (IACthresh < 0) {IACthresh = -IACthresh};
  samples_removed = c();
  temp = data;
  
  while (TRUE) {
    
    out = outlierSamples(temp,as.numeric(IACthresh),showplots);
    to_remove = out$samples_to_remove; 
    if (length(to_remove) < 1) {break}; 
    
    answer_raw = readline(prompt='List samples (0 if none) to remove with single spaces in between (no commas): ');
    if (answer_raw==0) {cat('You didn\'t remove any samples!!! \n'); break};
    answer = strsplit(answer_raw, ' ');
    answer = answer[[1]]; 
    
    temp = temp[, -match(answer, names(temp))];
    samples_removed = c(samples_removed, to_remove[names(to_remove)==answer]); 
    cat('Sample(s)', answer_raw, 'removed \n');
    
  };
  
  cat('\n');
  cat('Any more suspicious samples to remove?');
  choose = menu(c('Yes','No'));
  
  if (choose==1) {
    answer_raw = readline(prompt='List samples to remove with single spaces in between (no commas): ');
    answer = strsplit(answer_raw, ' ');
    answer = answer[[1]]; 
    temp = temp[, -match(answer, names(temp))];
    samples_removed = c(samples_removed, out$numbersd[answer]); 
    cat('Sample(s)', answer_raw, 'removed \n');
  } else {
    cat('Ok... finished\n')
  };
  
  output = list(dataClean=temp, 
                IAC=out$IAC, 
                meanIACdiag=out$meanIACdiag,
                samplemeanIAC=out$samplemeanIAC,
                numbersd=out$numbersd,
                samples_removed=samples_removed
  );
  
  return(output);
  
}

##########################################################


preProcess = function (datIN,
                       removeOutlierProbes=T, deviate=3, rowORcol=1,
                       removeTooManyNAs=T, probe_thresh=NULL, sample_thresh=NULL,
                       removeOutlierSamples=T, IACthresh=2, showplots=T,
                       Qnorm=T,
                       vsn=F, vsnPlot=F,
                       CVsort=F) {
  
  # check input, if ok, assign input to 'temp'. if not, quit 
  if (is.data.frame(datIN) || is.matrix(datIN)) {
    temp = datIN;
    cat('Input data has',nrow(temp),'rows (genes) and',ncol(temp),'columns (samples) \n\n');
  } else {
    stop('Input data must be in the form of a data frame or matrix!');
  };
  
  # if removeOutlierProbes=T, run on 'temp' and save processed data in 'temp' for further processing
  # if removeOutlierProbes=F, skip. 'temp' continues to hold input data.
  # 	set outlierProbesOUTPUT to NULL
  if (removeOutlierProbes) {
    outlierProbesOUTPUT = removeOutlierProbesIterate(temp, deviate, rowORcol);
    temp = outlierProbesOUTPUT$dataClean;
    cat('Processed data available in output as $data_removedOutlierProbes \n');
  } else {
    cat('Skipping removal of outlier probes, checking for probesets and samples with too much missing data...\n');
    outlierProbesOUTPUT = NULL;
  };
  
  # if removeTooManyNAs=T, run on 'temp' and save processed data in 'temp' for further processing	
  # if removeTooManyNAs=F, skip. 'temp' holds output of removeOutlierProbes or input data.
  # 	set checkMissingDataOUTPUT to NULL
  if (removeTooManyNAs) {
    checkMissingDataOUTPUT = removeTooManyNAs(temp, probe_thresh, sample_thresh);
    temp = checkMissingDataOUTPUT$dataClean;		
    cat('...Processed data ($data_checkedMissingData) has',nrow(temp),'rows and',ncol(temp),'columns \n');
  } else {
    cat('Skipping removal of probesets and samples with too much missing data, checking for outlier samples...\n');
    checkMissingDataOUTPUT = NULL;	
  };
  
  # if removeOutlierSamples=T, run on 'temp' and save processed data in 'temp' for further processing
  # if removeOutlierSamples=F, skip. 'temp' holds output of removeTooManyNAs, removeOutlierProbes, or input 
  # 	set outlierSamplesOUTPUT to NULL
  if (removeOutlierSamples) {
    outlierSamplesOUTPUT = outlierSamplesIterate(temp, IACthresh, showplots);
    temp = outlierSamplesOUTPUT$dataClean;
    cat('...Processed data ($data_removedOutlierSamples) now has',nrow(temp),'rows and',ncol(temp),'columns \n\n');
  } else {
    cat('Skipping removal of outlier samples...\n');
    outlierSamplesOUTPUT = NULL;
  };
  
  # if Qnorm=T, run on 'temp', save output to 'tempQnorm' for function output as 'data_Qnorm'.
  # 	set 'temp' equal to 'tempQnorm' for further processing.
  # 	ask if user wants to re-check for outlier samples but don't run iterative version as it is just a check
  # if Qnorm=F, skip and set data_Qnorm to NULL.
  # 	'temp' holds output of removeOutlierSamples, removeTooManyNAs, removeOutlierProbes, or input  
  if (Qnorm) {
    cat('..........\n');
    cat('Performing quantile normalization...\n')
    data_Qnorm = as.data.frame(normalize.quantiles(as.matrix(temp)));
    names(data_Qnorm) = names(temp); rownames(data_Qnorm) = rownames(temp);
    temp = data_Qnorm;
    cat('Normalized data available in output as $data_Qnorm \n\n');
    cat('Re-check for sample outliers?\n')
    answer = menu(c('Yes','No'));
    if (answer==1) {
      postNormOutlierSamples = outlierSamples(data_Qnorm, IACthresh, showplots);
      cat('\n');
      cat('Don\'t remove suspicious samples after normalizing!!!\n');
      cat('Instead, re-run and remove right before normalizing \n\n');
    } else {
      cat('You may want to re-check for outlier samples \n\n');
    };
  } else {
    cat('Skipping quantile normalization...\n\n');
    data_Qnorm = NULL;
  };
  
  # if vsn=T, run on 'temp', save output to 'tempVSN' for function output as 'data_VSN'.
  # 	set 'temp' equal to 'tempVSN' for further processing
  # if vsn=F, skip and set data_VSN to NULL.
  # 	'temp' holds output of Qnorm, removeOutlierSamples, removeTooManyNAs, removeOutlierProbes, or input
  if (vsn) {
    cat('Performing variance stabilization...\n\n');
    data_VSN = runVSN(temp, vsnPlot);
    temp = data_VSN;
    cat('Variance stabilized data available in output as $data_VSN \n\n');
  } else {
    cat ('Skipping variance stabilization...\n\n');
    data_VSN = NULL;
  };
  
  # if CVsort=T, compute CVs of rows in 'temp' and sort rows of 'temp' by CVs
  # 	save sorted 'temp' as 'tempCVsort' for function output as 'data_Sorted'
  # 	'temp' holds output of vsn, Qnorm, removeOutlierSamples, removeTooManyNAs, removeOutlierProbes, or input
  # if CVsort=F, skip and set 'data_Sorted' to NULL
  if (CVsort) {
    cat('Sorting probes by CV...\n')
    data_CVs = apply(temp, 1, cv);
    data_Sorted = temp[order(data_CVs, decreasing=T), ];
    cat('Sorted data available in output as $data_Sorted, CVs in $data_CVs \n\n');
  } else {
    cat('Skipping CV sort...\n\n');
    data_Sorted = NULL;
    data_CVs = NULL;
  };
  
  cat('Creating list for output...\n');
  output = list(outlierProbesOUTPUT=outlierProbesOUTPUT,
                checkMissingDataOUTPUT=checkMissingDataOUTPUT,
                outlierSamplesOUTPUT=outlierSamplesOUTPUT,
                data_removedOutlierProbes=outlierProbesOUTPUT$dataClean,
                data_checkedMissingData=checkMissingDataOUTPUT$dataClean,
                data_removedOutlierSamples=outlierSamplesOUTPUT$dataClean,
                data_Qnorm=data_Qnorm,
                data_VSN=data_VSN,
                data_CVs=data_CVs,
                data_Sorted=data_Sorted
  );
  
  cat('Write any of the output to .csv file?\n');
  choose = menu(c('Yes','No'));
  if (choose==1) {
    choices = names(output)[4:10];
    cat('Choose one...\n');
    choose = 3 + menu(choices=choices);
    file = readline('Enter a name for the output file, including .csv extension... ');
    write_out = output[[choose]];
    write.csv(write_out, file);
  } else {
    cat('You didn\'t write any of the output to a file\n');
  };
  
  cat('All done!\n');
  return(output); 
  
}

