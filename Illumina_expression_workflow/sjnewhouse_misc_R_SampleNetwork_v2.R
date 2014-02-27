###################
## sampleNetwork ##
###################


phenotype_data <- pData(eset) 
gene_expression_matrix <- exprs(eset)
outfile <- "test"
groupsLabel="GROUPS"

basic_sampleNetwork <- function(gene_expression_matrix , phenotype_data, groupsLabel="GROUPS", outfile, sd_thrs=2 ) {
    
    #' get list of groups
    group_list <- unique(phenotype_data[,groupsLabel])
    n_groups <- length(group_list)
    cat(" Groups [",group_list,"]","\r","\n")    
    
    #' make results matrix  
    res_names <- c("Group","mean_IAC","n_outliers","min_Z.K","rho","rho_pvalue","Z.K_outliers")
    res <- matrix(ncol=length(res_names),nrow=1000, dimnames=list(1:1000,res_names))

    #' set row = 1 to start   
    row <- 1
    
    #' loop throgh groups and calc network concepts
  for(mgroup in group_list) {
    
      cat(" doing group [",mgroup,"]","\r","\n")
      
      #' subset gene expression matrix to samples in group X
      
      gx <- gene_expression_matrix[,phenotype_data$GROUPS==mgroup]
      
      #' subset phenotype data to samples in group X       
      sample_info <- phenotype_data[phenotype_data$GROUPS==mgroup,]
      samle_names <- rownames(sample_info)
      
      #' set as character
      groups <- as.character(sample_info$GROUPS)
      pheno <- as.character(sample_info$PHENOTYPE)
      gpcolors <- labels2colors(groups)
      
      ## IAC, ADJACENCY & fundamentalNetworkConcepts, Z.K, Z.C, Z.MAR
      cat(" Calculating fundamentalNetworkConcepts Metrics ","\r","\n")
      
      IAC=cor(gx,method="p",use="p")
      diag(IAC)=0
      A.IAC=((1+IAC)/2)^2  ## ADJACENCY MATRIX
      
      FNC <- fundamentalNetworkConcepts(A.IAC) ## WGCNA
      
      K2 <- round(as.numeric(FNC$ScaledConnectivity),2)
      Z.K <- round(as.numeric((K2-mean(K2))/sd(K2)),2)
      Z.C <- round( as.numeric((FNC$ClusterCoef-mean(FNC$ClusterCoef))/sd(FNC$ClusterCoef)),2)
      Z.MAR <- round(as.numeric((FNC$MAR-mean(FNC$MAR))/sd(FNC$MAR)),2)
      rho <- as.numeric(signif(cor.test(Z.K,Z.C,method="s")$estimate,2))
      rho_pvalue <- as.numeric(signif(cor.test(Z.K,Z.C,method="s")$p.value,2))
      
      ## OUTLIERS
      Z.K_outliers <- Z.K < -sd_thrs
      Z.K_outliers <- colnames(gx)[Z.K_outliers==TRUE]
      n_outliers <- as.numeric(length(Z.K_outliers))
      
      mean_IAC <- as.numeric(round(as.numeric(mean(IAC[upper.tri(IAC)])),2))
      min_Z.K <- as.numeric(round(as.numeric(min(Z.K)),2))
      
      cat(" Number of Z.K outliers [", n_outliers,"]","\r","\n")
      cat(" mean_IAC [", mean_IAC,"]","\r","\n")
      
      ## Data frame of .SampleNetwork_Stats.txt
      cat(" Making Data fram of fundamentalNetworkConcepts Metrics ","\r","\n")
      FNC_DF <- as.data.frame(FNC)
      FNC_DF$Z.K <- round(as.numeric( (K2-mean(K2))/sd(K2) ),2)
      FNC_DF$Z.C <- round(as.numeric((FNC_DF$ClusterCoef-mean(FNC_DF$ClusterCoef))/sd(FNC_DF$ClusterCoef)),2)
      FNC_DF$Z.MAR <- round(as.numeric((FNC_DF$MAR-mean(FNC_DF$MAR))/sd(FNC_DF$MAR)),2)
      FNC_DF$mean_IAC <- round(as.numeric(mean(IAC[upper.tri(IAC)])),2)
      FNC_DF$Mean_Connectivity <- round(as.numeric(mean(FNC$Connectivity)),2)
      FNC_DF$Mean_ScaledConnectivity <- round(as.numeric(mean(FNC$ScaledConnectivity)),2)
      FNC_DF$Mean_ClusterCoef <- round(as.numeric(mean(FNC$ClusterCoef)),2)
      FNC_DF$Mean_MAR <- round(as.numeric(mean(FNC$MAR)),2)
      FNC_DF$Decentralization <- round(as.numeric(1-FNC_DF$Centralization),2)
      FNC_DF$Homogeneity <- round(as.numeric(1-FNC_DF$Heterogeneity),2)
      FNC_DF$rho <- round(as.numeric(signif(cor.test(FNC_DF$Z.K,FNC_DF$Z.C,method="s")$estimate)),2)
      FNC_DF$rho_pvalue <- as.numeric(signif(cor.test(FNC_DF$Z.K,FNC_DF$Z.C,method="s")$p.value))
      FNC_DF <- cbind(rownames(FNC_DF),FNC_DF)
      colnames(FNC_DF) <- c("Sample.ID",names(FNC_DF[-1]))
     
      ## write data
      cat(" Saving Data fram of fundamentalNetworkConcepts Metrics [",paste(outfile,".SampleNetwork_Stats.txt",sep=""),"]","\r","\n")
      write.table(FNC_DF,file=paste(outfile,".group.",mgroup,".SampleNetwork_Stats.txt",sep=""),sep="\t",row.name=FALSE,quote=FALSE)
      
      ## write data for IAC
      cat(" Saving Data fram of fundamentalNetworkConcepts Z.K outliers [",paste(outfile,".SampleNetwork_Stats_Z.K_outliers.txt",sep=""),"]","\r","\n")
      write.table(FNC_DF[Z.K_outliers,],file=paste(outfile,".group.",mgroup,".SampleNetwork_Stats_Z.K_outliers.txt",sep=""),sep="\t",row.name=FALSE,quote=FALSE)
      
      ## set colours by chip of GROUPS
      colorvec <- labels2colors(as.character(sample_info$Sentrix.Barcode))
      ## plots of fundamentalNetworkConcepts
      local({
        colLab <<- function(n,treeorder) {
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
      cat(" Plotting SampleNetwork Metrics [",paste(outfile,".group.",mgroup,".SampleNetwork.qc.pdf",sep=""),"]","\r","\n")
      meanIAC <- mean(IAC[upper.tri(IAC)],na.rm=T)
      cluster1 <- hclust(as.dist(1-A.IAC),method="average")
      cluster1order <- cluster1$order
      cluster2 <- as.dendrogram(cluster1,hang=0.1)
      cluster3 <- dendrapply(cluster2,colLab,cluster1order)
      ## PLOTS
      ## cluster IAC
      pdf(file=paste(outfile,".group.",mgroup,".SampleNetwork.qc.pdf",sep=""),width=11,height=8)
      par(mfrow=c(2,2))
      par(mar=c(5,6,4,2))
      plot(cluster3,nodePar=list(lab.cex=1,pch=NA),main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
      mtext(paste("distance: 1 - ISA ",sep=""),cex=0.8,line=0.2)
      ## Connectivity
      par(mar=c(5,5,4,2))
      plot(Z.K,main="Connectivity", ylab="Z.K",xaxt="n",xlab="Sample",type="n",cex.main=1.8,cex.lab=1.4)
      text(Z.K,labels=samle_names,cex=0.8,col=colorvec)
      abline(h=-2);	abline(h=-3)
      ## ClusterCoef
      par(mar=c(5,5,4,2))
      plot(Z.C,main="ClusterCoef", ylab="Z.C",xaxt="n",xlab="Sample",cex.main=1.8,cex.lab=1.4,type="n")
      text(Z.C,labels=samle_names,cex=0.8,col=colorvec)
      abline(h=-2);	abline(h=-3)
      ## Connectivity vs ClusterCoef
      par(mar=c(5,5,4,2))
      plot(Z.K,Z.C,main="Connectivity vs ClusterCoef",xlab="Z.K",ylab="Z.C",col=gpcolors,cex.main=1.8,cex.lab=1.4)
      abline(lm(Z.C~Z.K),col="black",lwd=2)
      mtext(paste("rho = ",signif(cor.test(Z.K,Z.C,method="s")$estimate,2)," p = ",signif(cor.test(Z.K,Z.C,method="s")$p.value,2),sep=""),cex=0.8,line=0.2)
      abline(v=-2,lty=2,col="grey");abline(h=-2,lty=2,col="grey")
      dev.off()
      
    ## RESULTS
 
      ##
      Z.K_outliers <- paste(Z.K_outliers, collapse=";")
      
      res[row,] <-  c(
                      Group=mgroup,
                      mean_IAC=as.numeric(round(mean_IAC,2)),
                      n_outliers=as.numeric(n_outliers),
                      min_Z.K=as.numeric(round(min_Z.K,2)),
                      rho=as.numeric(round(rho,2)),
                      rho_pvalue=rho_pvalue,
                      Z.K_outliers=Z.K_outliers);
            
      row <- row + 1
  }
    
    ## FINAL RESULTS
    res_out <- as.data.frame(res)
    res_out <- subset(res_out, res_out$mean_IAC!="NA")
    res_out[,c("mean_IAC","n_outliers","min_Z.K","rho","rho_pvalue")] <- sapply(res_out[,c("mean_IAC","n_outliers","min_Z.K","rho","rho_pvalue")],as.numeric )
    save(res_out,file=paste(outfile,".group.",mgroup,".SampleNetwork.qc.RData",sep=""))
    return(res_out)
  }

#################################
##  basic_sampleNetworkIterate ##
#################################
##```{r basic_sampleNetworkIterate, eval=FALSE}
## function(gene_expression_matrix , phenotype_data, groupsLabel="GROUPS", outfile, sd_thrs=2 ) 
basic_sampleNetworkIterate <- function(gene_expression_matrix , phenotype_data, groupsLabel="GROUPS", outfile, sd_thrs=2) {

  
    ## make results matrix basic_sampleNetworkIterate_summary
    basic_sampleNetworkIterate_summary_names <- c("Group","round","nSamp","nOutlier","mean_IAC","min_Z.K","rho","rho_pvalue","Z.K_outliers")
    basic_sampleNetworkIterate_summary <- matrix(ncol=length(basic_sampleNetworkIterate_summary_names),nrow=1000, dimnames=list(1:1000,basic_sampleNetworkIterate_summary_names))
    
    ## get group names
    group_list <- unique(phenotype_data[,groupsLabel])
    
    cat(" Number of Groups [",length(group_list),"]\nGroup Names [",group_list,"]","\r","\n")
    ## set outlier running count to 0
    outlier_running_count <- 0;
    ## iteration count
    iteration <- 1;
    ## sd threshold
    sd_thrs <- sd_thrs
    ## outlier count
    iac_outlier_samples <- c();
    ## n groups
    n_groups <- length(group_list)
    
    ## loop through groups and run basic_sampleNetwork
    
    for(mgroup in group_list)  {
      
      ## subset eset to group samples
  
      group_samples <- rownames(phenotype_data)[phenotype_data$GROUPS==mgroup]
      
      cat(" Subset gene_expression_matrix to group [",mgroup,"]","\r","\n")
      
      gx <- gene_expression_matrix
      gx_group <- gene_expression_matrix[,group_samples]
      phenotype_data_group <- phenotype_data[group_samples,]
     
      cat(" doing group [",mgroup,"]","\r","\n")
      
      ## run basic_sampleNetwork
      
      out <- basic_sampleNetwork(gx_group ,phenotype_data_group, groupsLabel="GROUPS",outfile=paste(outfile,".group.",mgroup,".round.",iteration,sep="" ))
      
      ## get stats from basic_sampleNetwork
      n_samp <-  length(colnames(gx_group))
      mgroup <- out$Group
      mean_IAC <- round(as.numeric(out$mean_IAC),2);
      min_Z.K <- round(as.numeric(out$min_Z.K),2);
      rho_pvalue <- as.numeric(signif(out$rho_pvalue,2));
      rho <- round(as.numeric(out$rho),2)
      
      ## Outlier samples
      Z.K_outliers <- out$Z.K_outliers
      
      Z.K_outliers_list <- paste(out$Z.K_outliers)
      
      Z.K_outliers_list <- strsplit(Z.K_outliers_list,";")[[1]]
      
      ## update outlier_running_count
      outlier_running_count <- as.numeric(outlier_running_count + length(Z.K_outliers_list));
      
      ## make list of outlier samples
      iac_outlier_samples <- c(Z.K_outliers_list, iac_outlier_samples)
      
      ## store in basic_sampleNetworkIterate_summary results matrix
      #basic_sampleNetworkIterate_summary_names <- c("Group","round","nSamp","nOutlier","mean_IAC","min_Z.K","rho","rho_pvalue","Z.K_outliers")
      res <- cbind(
        Group=mgroup,
        iteration=as.numeric(iteration),
        n_samp=as.numeric(n_samp),
        nOutlier=as.numeric(length(Z.K_outliers_list)),
        mean_IAC=round(as.numeric(mean_IAC),2),
        min_Z.K=round(as.numeric(min_Z.K),2),
        rho=round(as.numeric(rho),2),
        rho_pvalue=as.numeric(signif(rho_pvalue,2)),
        Z.K_outliers
        )
      
      colnames(res) <- c("Group","round","nSamp","nOutlier","mean_IAC","min_Z.K","rho","rho_pvalue","Z.K_outliers")
      
      basic_sampleNetworkIterate_summary[iteration,] <- res
      
      ## if outlier samples then remove them
      if(outlier_running_count >= 1) {
        
        gx_group <- gx_group[,(colnames(gx_group) %in% iac_outlier_samples)==FALSE ]
        group_samples <- colnames(gx_group)
        phenotype_data_group <- phenotype_data_group[group_samples,]
        
        } else {
          gx_group <- gx_group
          group_samples <- colnames(gx_group)
          phenotype_data_group <- phenotype_data_group[group_samples,]
          cat(" You have no outliers!","\r","\n") ;
          write.table(basic_sampleNetworkIterate_summary,file=paste(outfile,".basic_sampleNetworkIterate_summary.csv",sep=""),sep=",",row.names=FALSE,quote=FALSE)
        }
      
      ## number of samples remaining
      n_samp_left <- as.numeric(length(colnames(gx_group)))
      
      cat(" Number of outliers after round [",iteration,"] = [",outlier_running_count,"].  Percentage [",round(outlier_running_count/n_samp,3),"]. Mean IAC [",mean_IAC,"]. Min Z.K [" ,min_Z.K,"]. KvC [",rho,"] [",rho_pvalue ,"] N SAMPLE LEFT [",n_samp_left,"]","\r","\n")
      
      iteration <- iteration + 1; ## increment iteration number
      
      ## keep going until Z.K > -2 ##
      
      is_min_Z.K_bad <- min_Z.K <= -sd_thrs
      is_rho_bad <- abs(rho) < 0.90
      keep_going <- sum(is_min_Z.K_bad,is_rho_bad)
    
      while( iteration >= 1  & keep_going==2  ) {
        
        n_samp_start <- as.numeric(length(colnames(gx_group)))
        
        out <- basic_sampleNetwork(gx_group ,phenotype_data_group, groupsLabel="GROUPS",outfile=paste(outfile,".group.",mgroup,".round.",iteration,sep="" ))
        
        mean_IAC <- round(as.numeric(out$mean_IAC),2);
        min_Z.K <- round(as.numeric(out$min_Z.K),2);
        Z.K_outliers <- out$Z.K_outliers
        rho_pvalue <- as.numeric(signif(out$rho_pvalue,2))
        rho <- round(as.numeric(out$rho),2)
        
        Z.K_outliers_list <- paste(out$Z.K_outliers)
        
        Z.K_outliers_list <- strsplit(Z.K_outliers_list,";")[[1]]
        
        outlier_running_count <- as.numeric(outlier_running_count + length(Z.K_outliers_list));
        
        iac_outlier_samples <- c(Z.K_outliers_list, iac_outlier_samples)
        
        # remove outliers for next round
        gx_group <- gx_group[,(colnames(gx_group) %in% iac_outlier_samples)==FALSE ]
        group_samples <- colnames(gx_group)
        phenotype_data_group <- phenotype_data_group[group_samples,]
        n_samp_left <- as.numeric(length(colnames(gx_group)))
        
        res <- cbind(
          Group=mgroup,
          iteration=as.numeric(iteration),
          n_samp=as.numeric(n_samp_left),
          nOutlier=as.numeric(length(Z.K_outliers_list)),
          mean_IAC=round(as.numeric(mean_IAC),2),
          min_Z.K=round(as.numeric(min_Z.K),2),
          rho=round(as.numeric(rho),2),
          rho_pvalue=as.numeric(signif(rho_pvalue,2)),
          Z.K_outliers
        )
        
        colnames(res) <- c("Group","round","nSamp","nOutlier","mean_IAC","min_Z.K","rho","rho_pvalue","Z.K_outliers")
        
        basic_sampleNetworkIterate_summary[iteration,] <- res
        
        cat(" Number of outliers after round [",iteration,"] = [",outlier_running_count,"].  Percentage [",round(outlier_running_count/n_samp,3),"]. Mean IAC [",mean_IAC,"]. Min Z.K [" ,min_Z.K,"]. KvC [",rho,"] [",rho_pvalue ,"] N SAMPLE LEFT [",n_samp_left,"]","\r","\n")
        
        iteration <- iteration + 1;
        is_min_Z.K_bad <- min_Z.K <= -sd_thrs
        is_rho_bad <- abs(rho) < 0.90
        keep_going <- sum(is_min_Z.K_bad,is_rho_bad)
        
      };
      ##
    }
    ## DONE LOOPING THROUGH GROUPS ##
    
    final_basic_sampleNetworkIterate_summary <- as.data.frame(basic_sampleNetworkIterate_summary)
    
    final_basic_sampleNetworkIterate_summary <- subset(final_basic_sampleNetworkIterate_summary, final_basic_sampleNetworkIterate_summary$mean_IAC!="NA")
    
    write.table(final_basic_sampleNetworkIterate_summary,file=paste(outfile,".group.basic_sampleNetworkIterate_summary.csv",sep=""),sep=",",row.names=FALSE,quote=FALSE)
    
    outlier_samples <- iac_outlier_samples
    
    iac_outlier_samples_df <- as.data.frame(outlier_samples)
    
    colnames(iac_outlier_samples_df) <- c("Sample.ID")
    
    write.table(iac_outlier_samples_df,file=paste(outfile,".group.iac_SampleNetwork_outlier_samples.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
    
    outlier_samples <- list(iac_outlier_samples=outlier_samples)
    
    return(outlier_samples)
    
  }
  


