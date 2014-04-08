# Extras not contained in the R packages.

# Version 02: cleanup of old functions that have not been used in ages, add some new ones.

if (!exists("UseCpmax")) UseCpmax = FALSE;

#--------------------------------------------------------------------------------------------

# Load the requisite libraries

WorkingDirectory = getwd();

#if (!library(MASS, logical.return=TRUE)) { # standard, no need to install
#For some reason, MASS does not seem to be installed on Titan, so we'll try to load it 
# from my own library. If that fails as well, stop.
#  if (!library(MASS, logical.return=TRUE, lib.loc="M:/Work/RLibrary")) stop()
#}   

libraryLocs = c("")

loadPackages = c("MASS", "class", "cluster", "impute", "Hmisc", "survival",
"WGCNA") #, "mclust");

for (pack in loadPackages)
{
  if (!require(pack, character.only = TRUE, quietly = TRUE))
  {
    ok = FALSE;
    for (loc in libraryLocs) ok = ok | require(pack, lib.loc = loc, 
                                               character.only = TRUE, quietly = TRUE);
    if (!ok) stop(paste("Could not find package", pack, 
                        "either in the standard library nor in the personal libraries\n",
                        paste("     ", libraryLocs, collapse = "\n")));
  }
}


#library(MASS);
#library(class)	# standard, no need to install
#library(cluster)	
#library(impute)# install it for imputing missing value
#library(Hmisc)	# install it for the C-index calculations
#library(survival)
#library(fields);
#     
#oldwd = getwd();

if (exists("memory.limit"))
{
  # increase the available memory 
  memory.limit(size=4000)   
}

# kWithinModule = function(TOM, Colors, gene) -> intramodularConnectivity

options(stringsAsFactors = FALSE);

#-----------------------------------------------------------------------------------------------
# Trait selection based on independence and significance
# Assumes that, just like with the PCs, the Traits have the same columns in each dataset (though the
# sample sets need not be the same).

SelectTraits = function(Traits, BranchCut = 0.25, SelectOnSignificance = FALSE, PCs = NULL, 
                        SignifThres = 0.03, Impute = FALSE, verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent);
  if (verbose>0) printFlush(paste(spaces, "SelectTraits: Selecting from ", dim(Traits[[1]]$data)[2],
                                   "traits."));
  No.Sets = length(Traits);
  TDiss = 1-cor(Traits[[1]]$data, use = "pairwise.complete.obs");
  if (No.Sets>1) for (set in 2:No.Sets)
  {
     TDiss = pmax(TDiss, 1-cor(Traits[[set]]$data, use = "pairwise.complete.obs"));
  }
  h = hclust(as.dist(TDiss), method = "average");
  TMods = ModuleNumber(h, CutHeight = BranchCut, MinSize = 1);
  No.TMods = nlevels(as.factor(TMods));
  SelTraits = vector(mode="list", length = No.Sets);
  for (set in 1:No.Sets)
  {
    TData = Traits[[set]]$data;
    TData[is.na(TData)] = 0;
    TPCs = ModulePrinComps1(TData, as.factor(TMods), Impute = Impute, verbose = 0, 
                            GetConformity = TRUE);
    SelTraits[[set]] = list(data = TPCs$PrinComps);
    for (tmod in 1:No.TMods)
    {
      if (sum(TMods==tmod)>1)
      {
        rnk = order(-TPCs$ModuleConformity[TMods==tmod]);
        SelTraits[[set]]$data[, tmod] = (TData[, TMods==tmod])[, rnk[1]];
        names(SelTraits[[set]]$data)[tmod] = (names(TData)[TMods==tmod])[rnk[1]];
      } else {
        SelTraits[[set]]$data[, tmod] = TData[, TMods==tmod];
        names(SelTraits[[set]]$data)[tmod] = names(TData)[TMods==tmod];
      }
    }
  }

  if (verbose>0) printFlush(paste(spaces, "SelectTraits: Clustering led to ", dim(SelTraits[[1]]$data)[2],
                                   "traits."));
 
  if (SelectOnSignificance)
  {
    # Reduce further: calculate cor.tests for each ME with each trait in each set; 
    # keep only traits that have at least one cor.test$p.value below a threshold

    if (is.null(PCs)) stop("PCs must be given when SelectOnSignificance is requested.");

    No.Mods = dim(PCs[[1]]$data)[2];
    if (is.null(No.Mods)) 
       stop("Given PCs do not appear to have the correct structure (vector of list",
            "with \'data\' component being a matrix whose columns are PC vectors");

    No.Traits = dim(SelTraits[[1]]$data)[2];
    
    SelectTrait = rep(FALSE, times = No.Traits);
    
    for (trait in 1:No.Traits)
      for (mod in 1:No.Mods)
      {
        Significant = TRUE;
        for (set in (1:No.Sets))
        {
          ct = cor.test(PCs[[set]]$data[, mod], SelTraits[[set]]$data[, trait]); 
          if (ct$p.value>SignifThres) Significant = FALSE;
        }
        if (Significant) SelectTrait[trait] = TRUE;
      }

    for (set in 1:No.Sets)
    {
      SelTraits[[set]]$data = SelTraits[[set]]$data[, SelectTrait];
    }
    
    # print(paste("No. of selected traits expected by chance:", No.Mods * No.Traits * TraitThres));
  }
  
  # Re-cluster the significant traits for diagnostic purposes
  if (sum(SelectTrait)>1)
  {
    TDiss = 1-cor(SelTraits[[1]]$data, use = "pairwise.complete.obs");
    if (No.Sets>1) for (set in 2:No.Sets)
    {
       TDiss = pmax(TDiss, 1-cor(SelTraits[[set]]$data, use = "pairwise.complete.obs"));
    }
    newh = hclust(as.dist(TDiss), method = "average");
  } else {
    newh = NULL;
  }

  if (verbose>0) 
  {
    printFlush(paste(spaces, "SelectTraits: Selected", sum(SelectTrait), "traits: "));
    printFlush(paste(spaces, paste(names(SelTraits[[1]]$data), collapse = ", ")));
  }
                  

  list(No.SelectedTraits = sum(SelectTrait), Traits = SelTraits, ClusterTree = h, NewClusterTree = newh);
}

#----------------------------------------------------------------------------------------------
# 
# EvaluateClusters
#
#----------------------------------------------------------------------------------------------
# Input: distance matrix and cluster labels. 0 in Labels means the particular point is a
# singleton (not assigned to a cluster).
# Output: a list of cluster quality indicators.

EvaluateClusters = function(DistM, Labels, RefLabels = NULL, Sample = FALSE, SampleProp = 1,
                            verbose = 2, indent = 0)
{

  spaces = indentSpaces(indent);

  d = dim(DistM);
  No.Points = length(Labels);

  if (d[1]!=d[2]) stop("Distance matrix DistM must be square.");

  if (d[1]!=No.Points) stop("Dimension of DistM incompatible with number of cluster labels");

  if (verbose>0) print(paste(spaces, "Evaluating clusters...")); 

  # Calculate indices relating to the internal and external links

  if (Sample)
  {
    if ( (SampleProp<=0) | (SampleProp>=1) )
      stop(paste("Incorrect SampleProp parameter given:", SampleProp));
    No.Sampled = as.integer(No.Points * SampleProp);
    Sampled = sample(x = No.Points, size = No.Sampled);
    DistM = DistM[Sampled, Sampled];
    Labels = Labels[Sampled];
    RefLabels = RefLabels[Sampled];
    No.Points = No.Sampled;
  }

  InternalLinks = matrix(0, nrow = No.Points, ncol = No.Points);
  for (point in 1:No.Points) if (Labels[point]!=0)
    InternalLinks[, point] = ifelse(Labels==Labels[point], Labels, 0);

  IntLinks = as.vector(InternalLinks[upper.tri(InternalLinks)]);
  Dist = as.vector(DistM[upper.tri(DistM)]);
  No.Internal = sum(IntLinks!=0);
  No.External = length(IntLinks) - No.Internal;

  if ( (No.Internal > 0) & (No.External>0) )
  {
     ord = order(Dist, c(1:length(Dist)));

     # ...number of "misplaced" internal and external links

     s1 = Dist[ord[No.Internal]];
     s2 = Dist[ord[No.Internal+1]];
     if (s1==s2)
     {
       MaxInt = s1 + 1e-10;
       MinExt = s1 - 1e-10;
     } else {
       MaxInt = s1;
       MinExt = s2;
     }
     No.Misplaced = sum( ( (IntLinks>0) & (Dist>MaxInt) ) | ( (IntLinks==0) & (Dist<MinExt) ) );

     # ...weight index of the internal and external links

     DistInt = sum(Dist[IntLinks>0]);
     DistExt = sum(Dist[IntLinks==0]);
     DistSm = sum(Dist[ord <= No.Internal]);
     DistLg = sum(Dist[ord > No.Internal]);
     DistIndex = DistInt/DistSm * DistLg/DistExt;
  } else {
     No.Misplaced = 0;
     DistIndex = 0;
  }

  if (!is.null(RefLabels))
  {
    if (length(RefLabels)!=No.Points)
       stop("Length of given RefLabels incompatible with number of points.");
    
    RefInternalLinks = matrix(0, nrow = No.Points, ncol = No.Points);
    for (point in 1:No.Points) if (RefLabels[point]!=0)
      RefInternalLinks[, point] = ifelse(RefLabels==RefLabels[point], RefLabels, 0);

    RefIntLinks = as.vector(RefInternalLinks[upper.tri(RefInternalLinks)]);
 
    No.Agreed = sum( xor(IntLinks>0, RefIntLinks==0));
    No.Disagreed = sum( xor(IntLinks==0, RefIntLinks==0));

    RandIndex = No.Agreed/(No.Agreed+No.Disagreed);

    # Make each unassigned label unique
    XLabels = Labels;
    NUnassigned = sum(XLabels==0);
    start = max(XLabels);
    XLabels[XLabels==0] = start + c(1:NUnassigned);

    # Same for the reference labels 
    XRefLabels = RefLabels;
    NRefUnassigned = sum(XRefLabels==0);
    start = max(XRefLabels);
    XRefLabels[XRefLabels==0] = start + c(1:NRefUnassigned);
 
    AdjRandIndex = adjustedRandIndex(XLabels, XRefLabels);  # This requires package mclust
  } else {
    RandIndex = 0;
    AdjRandIndex = 0;
  }

  list(No.Misplaced = No.Misplaced, DistIndex = DistIndex, RandIndex = RandIndex, 
       AdjRandIndex = AdjRandIndex);
}


# Biweight correlation. The original functions are from
# http://www.unt.edu/benchmarks/archives/2001/december01/rss.htm and that in
# turn seems to be based on Wilcox (1997, page 197).
# The actually useful versions are PL's rewrites that use block code so the
# functions can be used with matrices as inputs. PL's function could be
# improved since it does some calculations unnecessarily twice when y==NULL.

bicov.original<-function(x, y, na.rm = FALSE){
mx <- median(x, na.rm = na.rm)
my <-median(y, na.rm = na.rm)
ux <- abs((x - mx)/(9 * qnorm(0.75) * mad(x, na.rm = na.rm)))
uy <- abs((y - my)/(9 * qnorm(0.75) * mad(y, na.rm = na.rm)))
aval <- ifelse(ux <= 1, 1, 0)
bval <- ifelse(uy <= 1, 1, 0)
top <- sum(aval * (x - mx) * (1 - ux^2)^2 * bval * (y - my) * (1 - uy^2)^2, na.rm = na.rm)
top <- sum(!is.na(x) & !is.na(y)) * top
botx <- sum(aval * (1 - ux^2) * (1 - 5 * ux^2), na.rm = na.rm)
boty <- sum(bval * (1 - uy^2) * (1 - 5 * uy^2), na.rm = na.rm)
bi <- top/(botx * boty)
bi
}

bicor.original<-function(x, y, na.rm = FALSE){
x <-as.numeric(x)
y<-as.numeric(y)
bicov(x,y, na.rm = na.rm)/(sqrt(bicov(x,x, na.rm = na.rm)*bicov(y,y, na.rm = na.rm)))
}

# This function should now work for both vector and matrix x and y (all combinations);
# Big WARNING: The use = pairwise.complete.obs simply works by removing NAs independently in x and
# y!!!!!

# Weighted version of bicov and bicor

bicovw = function(x, xWeight, y = NULL, yWeight = NULL, robustX = TRUE, robustY = TRUE, 
                  use = "all.obs", diagOnly = FALSE)
{
  x = as.matrix(x);
  xWeight = as.matrix(xWeight);
  if (!all.equal(dim(x), dim(xWeight)))
     stop("x and xWeight must have the same dimensions.");

  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)) 
    stop(paste("Unrecognized parameter 'use'. Recognized values are \n",
         "'all.obs', 'pairwise.complete.obs'"))
  na.rm = (na.method==2);
  if (is.null(y)) {
    if (robustX)
    {
      #mx <- apply(x, 2, median, na.rm = na.rm)
      mx = medianw(x, xWeight);
      mxMat = matrix(mx, nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      madxMat = matrix(apply(x, 2, mad, na.rm = na.rm), nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      ux <- as.matrix(abs((x - mxMat)/(9 * qnorm(0.75) * madxMat)))
    } else {
      mx = meanw(x, xWeight)
      mxMat = matrix(mx, nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      ux = matrix(0, nrow = nrow(mxMat), ncol = ncol(mxMat));
    }
    aval  = ifelse(ux <= 1, 1, 0) * xWeight;
    botx <- apply(as.matrix(aval * (1 - ux^2) * (1 - 5 * ux^2)), 2, sum, na.rm = na.rm);
    ux[is.na(ux)] = 1;
    aval[is.na(aval)] = 0;
    if (diagOnly)	
    {
      # diagOnly makes sense particularly for y=NULL because it is used in the
      # normalization in bicor.
      topFact = apply(!is.na(x), 2, sum)
      x[is.na(x)] = mxMat[is.na(x)];
      top <- apply(as.matrix(aval * (x - mxMat) * (1 - ux^2)^2)^2, 2, sum); 
      top <- top * topFact;
      bi <- top/((botx^2) * (1-1/nrow(x)));
    } else { 
      topFact = t(!is.na(x)) %*% !is.na(x) 
      x[is.na(x)] = mxMat[is.na(x)];
      top <- t(as.matrix(aval * (x - mxMat) * (1 - ux^2)^2)) %*% 
               as.matrix((aval * (x - mxMat) * (1 - ux^2)^2))
      top <- top * topFact;
      bi <- top/(botx %o% (botx *(1-1/nrow(x)) ))
    }
    bi
  } else {
    y = as.matrix(y);
    yWeight = as.matrix(yWeight);
    if (robustX)
    {
      #mx <- apply(x, 2, median, na.rm = na.rm)
      mx = medianw(x, xWeight)
      mxMat = matrix(mx, nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      madxMat = matrix(apply(x, 2, mad, na.rm = na.rm), nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      ux <- as.matrix(abs((x - mxMat)/(9 * qnorm(0.75) * madxMat)))
    } else {
      mx = meanw(x, xWeight)
      mxMat = matrix(mx, nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
      ux = matrix(0, nrow = nrow(mxMat), ncol = ncol(mxMat));
    }
    if (robustY)
    {
      #my <- apply(y, 2, median, na.rm = na.rm)
      my = medianw(y, yWeight)
      myMat = matrix(my, nrow = nrow(y), ncol = ncol(y), byrow = TRUE);
      madyMat = matrix(apply(y, 2, mad, na.rm = na.rm), nrow = nrow(y), ncol = ncol(y), byrow = TRUE);
      uy <- as.matrix(abs((y - myMat)/(9 * qnorm(0.75) * madyMat)))
    } else {
      my = meanw(y, yWeight);
      myMat = matrix(my, nrow = nrow(y), ncol = ncol(y), byrow = TRUE);
      uy = matrix(0, nrow = nrow(myMat), ncol = ncol(myMat));
    }
    aval <- ifelse(ux <= 1, 1, 0) * xWeight;
    bval <- ifelse(uy <= 1, 1, 0) * yWeight;
    botx <- apply(as.matrix(aval * (1 - ux^2) * (1 - 5 * ux^2)), 2, sum, na.rm = na.rm)
    boty <- apply(as.matrix(bval * (1 - uy^2) * (1 - 5 * uy^2)), 2, sum, na.rm = na.rm)
    ux[is.na(ux)] = 1;
    uy[is.na(uy)] = 1;
    aval[is.na(aval)] = 0;
    bval[is.na(bval)] = 0;
    topFact = t(!is.na(x)) %*% !is.na(y) 
    x[is.na(x)] = mxMat[is.na(x)];
    y[is.na(y)] = myMat[is.na(y)];
    top <- t(as.matrix(aval * (x - mxMat) * (1 - ux^2)^2)) %*% as.matrix((bval * (y - myMat) * (1 - uy^2)^2))
    top <- top * topFact;
    bi <- top/(botx %o% (boty * (1-1/nrow(x))))
    bi
  }
}

# weighted version of biweight mid-corellation
#
bicorw<-function(x, xWeight, y = NULL, yWeight = NULL, robustX = TRUE, robustY = TRUE, use = "all.obs")
{
  na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method)) 
    stop(paste("Unrecognized parameter 'use'. Recognized values are \n",
         "'all.obs', 'pairwise.complete.obs'"))
  if (!all.equal(dim(x), dim(xWeight)))
    stop("x and xWeight must have the same dimensions");
  if (is.null(y)) 
  {
    bcx = sqrt(bicov(x, robustX = robustX, use = use, diagOnly = TRUE));
    bicov(x, robustX = robustX, use = use)/(bcx %o% bcx);
  } else
  {
    if (!all.equal(dim(y), dim(yWeight)))
       stop("y and yWeight must have the same dimensions");
    bicov(x, y, robustX = robustX, robustY = robustY, use = use)/
                (sqrt(bicov(x, robustX = robustX, use = use, diagOnly = TRUE)) %o%
                 sqrt(bicov(y, robustX = robustY, use = use, diagOnly = TRUE)))
  }
}

meanw = function(x, w, na.rm = FALSE)
{
  x = as.matrix(x);
  w = as.matrix(w);
  if (!all.equal(dim(x), dim(w)))
    stop("Dimensions of x and w must be the same.");

  if (na.rm)
  {
    w[is.na(w)] = 0;
    w[is.na(x)] = 0;
  }
  s = apply(x*w, 2, sum, na.rm = TRUE);
  sw = apply(w, 2, sum, na.rm = TRUE);

  if (sum(sw==0)>0)
    stop("Sum of weights is zero in some columns.")

  s/sw;
}


medianw = function(x, w, na.rm = FALSE)
{
  x = as.matrix(x);
  w = as.matrix(w);
  if (!all.equal(dim(x), dim(w)))
    stop("Dimensions of x and w must be the same.");

  if (na.rm)
  {
    w[is.na(w)] = 0;
    w[is.na(x)] = 0;
  }

  if (sum(w<0)>0) stop("Weights w must be non-negative.");

  order = apply(x, 2, order, na.last = TRUE);
  for (i in 1:ncol(x)) w[, i] = w[order[, i], i];

  nSamples = nrow(w);
  csw = apply(w, 2, cumsum)
  sw = apply(w, 2, sum);
  if (sum(sw==0)>0)
    stop("Some columns of w have no non-zero entries.")
  # csw are cumulative sums of weights
  csw = csw * matrix(nSamples/csw[nSamples, ], nrow = nSamples, ncol = ncol(w), byrow = TRUE);

  midw = nSamples/2;
  # midh will be the index in each column where csw rises above midw. Can be anywehere between 1 and
  # nSamples.
  midh = apply(csw-midw>0, 2, match, x = TRUE)

  # Couldn't figure out a way to do this in block notation.
  medians = rep(0, ncol(x));
  for (i in 1:ncol(x))
    if (midh[i]==1) {
      medians[i] = x[order[1, i], i];
    } else {
      w0 = abs(csw[midh[i]-1]); x0 = x[order[midh[i]-1, i], i];
      w1 = csw[midh[i]]; x1 = x[order[midh[i], i], i];
      medians[i] = (w0*x1 + w1* x0)/(w1+w0);
    }
  medians;
}

propVarExplained = function(datExpr, colors, MEs, corFnc = "cor", corOptions = "use = 'p'")
{
  fc = as.factor(colors);
  mods = levels(fc);
  nMods = nlevels(fc);
  if (nMods!=ncol(MEs))
    stop(paste("Input error: number of distinct 'colors' differs from\n", 
               "the number of module eigengenes given in ME."));

  if (ncol(datExpr)!=length(colors))
    stop("Input error: number of probes (columns) in 'datExpr' differs from the length of goven 'colors'.");

  if (nrow(datExpr)!=nrow(MEs))
    stop("Input error: number of observations (rows) in 'datExpr' and 'MEs' differ.");

  PVE = rep(0, nMods);

  col2MEs = match(mods, substring(names(MEs), 3));

  if (sum(is.na(col2MEs))>0)
    stop("Input error: not all given colors could be matched to names of module eigengenes.");

  for (mod in 1:nMods)
  {
    modGenes = c(1:nGenes)[as.character(colors)==mods[mod]];
    corExpr = parse(text = paste(corFnc, "(datExpr[, modGenes], MEs[, col2MEs[mod]],",
                                 corOptions, ")"));
    PVE[mod] = mean(as.vector(eval(corExpr)^2));
  }

  names(PVE) = paste("PVE", mods, sep = "");
  PVE
}
 

cat("\n");
if (exists("cutreeDynamic", inherits = FALSE, mode = "function"))
{
  cat("\n");
  cat(paste("Warning: a function of the name 'cutreeDynamic' appears to be defined", 
            "in the main environment.\n")); 
  cat(paste("  It is likely that an old version of NetworkFunctions has been loaded.\n"));
  cat(paste("  Please do not use old NetworkFunctions, they are not compatible with this code.\n"));
  cat("\n");
  cat(paste("  The old definition has been removed; if the old functionality of cutreeDynamic is needed, \n",
            "  replace calls to the old cutreeDynamic(...) by cutreeDynamic.1(...).\n",
            "  Arguments to the function can be left unchanged.\n"));
  cat("\n");
  cat(paste("  After loading the dynamicTreeCut package, type help(cutreeDynamic) to see information\n",
            "  about the new dynamic tree cutting function and its use.\n"));
  cat("\n");

  rm("cutreeDynamic", inherits = FALSE);
}

cutreeDynamic.1 = function(hierclust, maxTreeHeight=1, deepSplit=TRUE, minModuleSize=50, 
                           minAttachModuleSize=100, nameallmodules=FALSE, useblackwhite=FALSE)
{
  if ( (minAttachModuleSize!=100) | (nameallmodules!=FALSE) | (useblackwhite!=FALSE) )
    waring("cutreeDynamic.1: Parameters minAttachModuleSize, nameallmodules, useblackwhite are ignored.");
  labels2colors(cutreeDynamicTree(hierclust, maxTreeHeight, deepSplit, minModuleSize));
}

#===================================================================================================
#
# moduleEigenMatrices.R
#
#===================================================================================================


#==========================================================================================================
#
# Module Eigenmatrix functions
#
#==========================================================================================================
#
# input: expr: Expression data of a single set
#	colors: module colors
#	impute: use imputation to replace NA's in data?
#	nPC: optional explicit specification of the number of PCs to keep
#	minVarExpl: the number of kept PCs will be such that at least minVarExpl variance will be
#	explained.
#	align: specification of alignment of the first PC.
#	saveAll: keep all eigenvectors and eigenvalues? Should only be useful for diagnostic purposes.
# output:
# Eigen-information will be stored in a list with the following components: 
# modName: name of module
# d: singular values
# v: singular vectors
# mat: singular matrix ( = sum{ labda_i sv_i sv_i^T })
# NSq: normalization constant = sum lamba_i^2
# varExpl: variance explained by the singular vectors
# averageExpr: average expression of the module as a vector

moduleEigenMatrices=function(expr, colors, excludeGrey = TRUE, 
                             grey = if (is.numeric(colors)) 0 else "grey", 
                             impute = TRUE, nPC = NULL, minVarExpl = 0.9, 
                             align = "along average", saveAll = FALSE, scaleExpr = TRUE,
                             maxNPC = 10, 
                             verbose = 0, indent = 0)
{
  nSamples = nrow(expr);
  spaces = indentSpaces(indent);

  if (verbose>0) 
     printFlush(paste(spaces, "moduleEigenMatrices: Calculating", nlevels(as.factor(colors)), 
                              "module eigenMatrices in given set."));

  if (is.null(expr)) stop("moduleEigenMatrices: Error: expr is NULL. ");
  if (is.null(colors)) stop("moduleEigenMatrices: Error: colors is NULL. ");

  alignRecognizedValues =  c("", "along average");
  if (!is.element(align, alignRecognizedValues))
    stop(paste("moduleEigenMatrices: Error:",
                "parameter align has an unrecognised value:", 
                align, "; Recognized values are ", alignRecognizedValues));

  if ( (minVarExpl < 0) | (minVarExpl > 1))
    stop("Incorrect minVarExpl given. Must be between 0 and 1.");

  modlevels=levels(factor(colors))
  if (excludeGrey) modlevels = modlevels[modlevels != as.character(grey)];

  nMods = length(modlevels);

  if (!is.null(nPC)) maxNPC = nPC;

  eigenMatrices = list(modNames = modlevels, 
                       nPC = rep(NA, nMods), 
                       d = matrix(NA, maxNPC, nMods),
                       PCs = array(NA, dim = c(nSamples, maxNPC, nMods)), 
                       mat = array(NA, dim = c(nSamples, nSamples, nMods)),
                       NSq = rep(NA, nMods),
                       varExpl = matrix(NA, maxNPC, nMods),
                       averageExpr = matrix(NA, nSamples, nMods),
                       nSamples = nSamples);
 # names(PrinComps)=paste("ME",modlevels,sep="")
 # names(aveExpr)=paste("AE",modlevels,sep="")

  #print(paste("nPC:", nPC))
  for(i in c(1:length(modlevels)) ) 
  {
    if (verbose>1) 
      printFlush(paste(spaces, "moduleEigengenes : Working on ME for module", modlevels[i]));
    modulename = modlevels[i]
    restrict1 = as.character(colors)== modulename
    datModule = t(expr[, restrict1])
    if (impute)
    {
        seedSaved = FALSE;
        if (exists(".Random.seed")) {
           saved.seed = .Random.seed;
           seedSaved = TRUE;
        }
        datModule = impute.knn(as.matrix(datModule), k = min(10, nrow(datModule)-1))
        if (seedSaved) .Random.seed = saved.seed;
    }
    if (scaleExpr) datModule=t(scale(t(datModule)));
    svd1=svd(datModule);
    #print(paste("sum(d^2):", sum(svd1$d^2)))
    varExpl= svd1$d^2/sum(svd1$d^2);
    varExplCum = cumsum(varExpl);
    if (is.null(nPC))
    {
      modNPC = min(maxNPC, match(TRUE, minVarExpl <= varExplCum));
      if (is.na(modNPC)) stop ("Internal error: modNPC is NA. Sorry!");
    } else {
      modNPC = min(nPC, length(varExpl));
    }
    #print(paste("modNPC:", modNPC))
    aveExpr = apply(datModule, 2, mean);  #datModule is transposed...
    if (align == "along average")
    {
      if (verbose>4) printFlush(paste(spaces,
                       " .. aligning module eigengene with average expression."))
      if (cor(aveExpr, svd1$v[,1])<0) svd1$v[,1] = -svd1$v[,1];
    }
    modEMat = matrix(0, nrow = nSamples, ncol = nSamples);
    for (pc in 1:modNPC)
    {
      pcMat = svd1$v[, pc] %o%svd1$v[, pc];
      modEMat = modEMat + abs(svd1$d[pc]) * pcMat
    }
    #eigenMatrices[[i]] = list(
    #    modName = paste("ME", modulename, sep = ""),
    #    n = modNPC,
    #    d = svd1$d[1:modNPC],
    #    PCs = scale(svd1$v[, c(1:modNPC)]),
    #    mat = modEMat,
    #    NSq = sum((svd1$d[1:modNPC])^2),
    #    varExpl = varExpl[1:modNPC],
    #    averageExpr = aveExpr,
    #    nSamples = ncol(datModule));
    eigenMatrices$modNames[i] = paste("ME", modulename, sep = "");
    eigenMatrices$nPC[i] = modNPC;
    eigenMatrices$d[c(1:modNPC), i] = svd1$d[1:modNPC];
    eigenMatrices$PCs[, c(1:modNPC), i] = scale(svd1$v[, c(1:modNPC)]);
    eigenMatrices$mat[,,i] = modEMat;
    eigenMatrices$NSq[i] =  sum((svd1$d[1:modNPC])^2);
    eigenMatrices$averageExpr[, i] = aveExpr;
    eigenMatrices$varExpl[1:modNPC, i] = varExpl[1:modNPC];
    #if (saveAll)
    #{
    #  eigenMatrices[[i]]$dAll = svd1$d;
    #  eigenMatrices[[i]]$vAll = svd1$v;
    #}
  }
  class(eigenMatrices) = 'eigenMatrices'
  eigenMatrices;
 
}

subsetMEMs = function(MEMs, subset, drop = TRUE)
{
  nMods = length(MEMs$modNames);
  if (sum(abs(subset)>nMods) > 0)
   stop("Non-existing eigenmatrices selected")

  nNewMods = length(eigenMatrices$modNames[subset]);
  if (nNewMods==0)
    stop("Selection is empty.");

  newMEMs = list(modNames = MEMs$modNames[subset, drop = drop], 
       nPC = MEMs$nPC[subset, drop = drop], 
       d = MEMs$d[, subset, drop = drop],
       PCs = MEMs$PCs[, , subset, drop = drop], 
       mat = MEMs$mat[, , subset, drop = drop],
       NSq = MEMs$NSq[subset, drop = drop],
       varExpl = MEMs$varExpl[, subset, drop = drop],
       averageExpr = MEMs$averageExpr[, subset, drop = drop],
       nSamples = MEMs$nSamples);
  class(newMEMs) = 'eigenMatrices';
  newMEMs;
}
  
corVMat = function(vec, MEMs, vecSelect = NULL, memSelect = NULL, na.rm = FALSE)
{
  vec = scale(as.matrix(vec))
  vecNames = dimnames(vec)[[2]];
  nMods = length(MEMs$nPC);
  nSamples = nrow(vec);
  if (is.null(vecSelect)) vecSelect = c(1:ncol(vec));
  if (is.null(memSelect)) memSelect = c(1:nMods);

  vecSelect = c(1:ncol(vec))[vecSelect];
  memSelect = c(1:nMods)[memSelect];

  cor = matrix(NA, length(vecSelect), length(memSelect));

  for (v in 1:length(vecSelect)) for (m in 1:length(memSelect))
  {
    vvec = vec[, vecSelect[v]];
    if (na.rm)
    {
      keep = c(1:nSamples)[!is.na(vec)];
    } else {
      keep = c(1:nSamples); 
    }
    proj = MEMs$mat[keep, keep, memSelect[m]] %*% vec[keep];
    cor[v, m] = sign(cor(MEMs$PCs[,1,memSelect[m]], vvec, use = 'p')) * 
                 sqrt( sum(proj^2)/sum(vec[keep]^2)/MEMs$NSq[memSelect[m]] );
  }
  dimnames(cor) = list(vecNames[vecSelect], MEMs$modNames[memSelect]);
  cor;
}

# This function may not work!

corMatMat = function(eigenMat1, eigenMat2)
{
  n1 = nrow(eigenMat1$mat);
  n2 = nrow(eigenMat2$mat);
  if (n1!=n2) stop("Matrices in eigenMat1 and eigenMat2 must have the same dimension.");

  Trace = sum(diag(eigenMat1$mat %*% eigenMat2$mat));
  sqrt(Trace/sqrt(eigenMat1$NSq * eigenMat2$NSq));
}

corMats = function(MEMs1, MEMs2 = NULL, select1 = NULL, select2 = NULL)
{
  nMods1 = length(MEMs1$nPC);
  if (is.null(select1)) {
    select1 = c(1:nMods1)
  } else select1 = c(1:nMods1)[select1];
  n1 = length(select1);
  if (is.null(MEMs2))
  {
    cm = matrix(1, n1, n1);
    dimnames(cm) = list(MEMs1$modNames[select1], MEMs1$modNames[select1]);
    for (i in 1:(n1-1)) for (j in (i+1):n1)
    {
      ii = select1[i];
      jj = select1[j];
      Trace = sum(diag(MEMs1$mat[,,ii] %*% MEMs1$mat[,,jj]));
      sign = sign(cov(MEMs1$PC[, 1, ii], MEMs1$PC[, 1, jj]))
      cm[i,j] = sign * sqrt(Trace/sqrt(MEMs1$NSq[ii] * MEMs1$NSq[jj]));
      cm[j,i] = cm[i,j]
    }
  } else
  {
    d1 = dim(MEMs1$mat)[1]
    d2 = dim(MEMs2$mat)[1]
    if (d1!=d2) stop("Eigenmatrices in MEMs1 and MEMs2 have incompatible dimensions.");
    n2 = length(select2);
    cm = matrix(0, n1, n2);
    dimnames(cm) = list(MEMs1$modNames[select1], MEMs2$modNames[select2]);
    for (i in 1:n1) for (j in 1:n2)
    {
      ii = select1[i];
      jj = select2[j];
      Trace = sum(diag(MEMs1$mat[,,ii] %*% MEMs2$mat[,,jj]));
      sign = sign(cov(MEMs1$PC[, 1, ii], MEMs2$PC[, 1, jj]))
      cm[i,j] = sign * sqrt(Trace/sqrt(MEMs1$NSq[ii] * MEMs2$NSq[jj]));
    }
  }
  cm;
}

#--------------------------------------------------------------------------------------
#
# orderMEMs
#
#--------------------------------------------------------------------------------------

orderMEMs = function(MEMs)
{
  nMats = length(MEMs);
  if (nMats==0) stop("Given 'MEMs' is empty.");
  dist = as.dist(1-corMats(MEMs))
  h = hclust(dist, method = "average");
  order = h$order
  subsetMEMs(MEMs, order);
}

checkMultiMEMs = function(multiMEMs)
{
  if (!is.list(multiMEMs))
   stop("Given 'multiMEMs' must be a list.");

  nSets = length(multiMEMs)
  if (nSets==0) stop("Given 'multiMEMs' is empty.");

  for (set in 1:nSets)
  {
    if (class(multiMEMs[[set]])!='eigenMatrices')
      stop("Not all entries in 'multiMEMs' have class 'eigenMatrices'");
  }
  nMods = length(multiMEMs[[1]]$nPC);
  nSamples = rep(NA, nSets)
  for (set in 1:nSets)
  {
    if (length(multiMEMs[[set]]$nPC) != nMods)
      stop('Numbers of modules not consistent across sets');
    nSamples[set] = MEMs[[set]]$nSamples;
  } 
  modNames = MEMs[[1]]$modNames;

  list(nSets = nSets, nMods = nMods, nSamples = nSamples, modNames = modNames);
}

consensusOrderMEMs = function(multiMEMs, order = NULL)
{
  size = checkMultiMEMs(multiMEMs);

  if (is.null(order))
  {
    dist = as.dist(1-corMats(multiMEMs[[1]]))
    if (size$nSets>1)
      for (set in 2:size$nSets)
        dist = pmax(dist, as.dist(1-corMats(multiMEMs[[set]])));
    h = hclust(dist, method = "average");
    order = h$order
  } else {
    if (length(order)!=size$nMods) 
      stop("Number of modules in 'multiMEMs' must match the length of 'order'.");
  }
  ordMEMs =  vector(mode = "list", length = size$nSets);
  for (set in 1:size$nSets)
  {
    ordMEMs[[set]] = subsetMEMs(multiMEMs[[set]], order);
  }
  ordMEMs;
}

#--------------------------------------------------------------------------------------
#
# PlotCorMEMsAndDendros
#
#--------------------------------------------------------------------------------------
# Plots a matrix plot of the correlation of module eigenmatrices. 
# On the diagonal the heatmaps show correlation of MEMs in the
# particular subset; off-diagonal are differences in the correlation matrix. 
# Titles is a vector of titles for the diagonal diagrams; the off-diagonal will have no title
# for now.

PlotCorMEMsAndDendros = function(MEMs, Titles, ColorLabels = TRUE, colors = NULL, IncludeSign = TRUE, 
                      ColoredBarPlot = TRUE, LetterSubPlots = FALSE, Letters = NULL, IncludeGrey = FALSE, 
                      setMargins = TRUE, plotCPMeasure = FALSE, plotMeans = TRUE, CPzlim = c(0,1),
                      printAdj = FALSE, printCPVals = FALSE, CPcex = 0.9, plotErrors = FALSE, 
                      marDendro = NULL,
                      marHeatmap = NULL, PlotDiagAdj = TRUE, invertColors = FALSE,
                      ...)
{
  #Letters = "abcdefghijklmnopqrstuvwxyz";
  if (is.null(Letters)) Letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

  if (is.null(colors)) 
    if (IncludeSign)
    {
      colors = GreenWhiteRed(50);
    } else {
      colors = heat.colors(30);
    }
  size = checkMultiMEMs(MEMs);
  nSets = size$nSets;
  cex = par("cex");
  mar = par("mar");
  par(mfrow = c(nSets+1, nSets));
  par(cex = cex);
  cors = list();
  for (set in 1:nSets) cors[[set]] = corMats(MEMs[[set]]);
  #if (!IncludeGrey)
  #{
  #  for (set in 1:nSets)
  #    PCs[[set]]$data = PCs[[set]]$data[ , substring(names(PCs[[set]]$data),3)!="grey"]
  #}
  letter.ind = 1;
  for (set in 1:nSets)
  {
    #par(cex = StandardCex/1.4);
    par(mar = marDendro);
    labels = size$modNames;
    uselabels = labels # [substring(labels,3)!="grey"];
    disPC = as.dist(1-cors[[set]]);
    clust = hclust(disPC, method = "average");
    if (LetterSubPlots) {
      main = paste(substring(Letters, letter.ind, letter.ind), ". ", Titles[set], sep="");
    } else {
      main = Titles[set];
    }
    #validColors = is.na(match(uselabels, colors()));
    #plotLabels = ifelse(validColors, substring(uselabels[validColors], 3), uselabels[!validColors]);
    plotLabels = uselabels;
    plot(clust, main = main, sub="", xlab="", 
         labels = plotLabels, ylab="", ylim=c(0,1));
    letter.ind = letter.ind + 1;
  }

  for (i.row in (1:nSets))
  {
    for (i.col in (1:nSets))
    {
      letter.ind = i.row * nSets + i.col;
      if (LetterSubPlots) 
      {
         #letter = paste("(", substring(Letters, first = letter.ind, last = letter.ind), ")", sep = "");
         letter = paste( substring(Letters, first = letter.ind, last = letter.ind), ".  ", sep = "");
      } else {
         letter = NULL;
      }
      par(cex = cex);
      if (setMargins) {
        if (is.null(marHeatmap))
        {
          if (ColorLabels) {
            par(mar = c(1,2,3,4)+0.2);
          } else {
            par(mar = c(6,7,3,5)+0.2);
          }
        } else {
          par(mar = marHeatmap);
        }
      }
      No.Modules = size$nMods;
      if (i.row==i.col)
      {
        if (IncludeSign)
        {
          if (PlotDiagAdj) {
           HeatmapWithTextLabels((1+cors[[i.row]])/2, dimnames(cors[[i.row]])[[2]],
                                 dimnames(cors[[i.row]])[[2]],
                                 main=paste(letter, Titles[[i.col]]), InvertColors=invertColors, 
                                 zlim=c(0,1.0),
                                 NumMatrix = if (printAdj) signif((1+cors[[i.row]])/2, 2) else NULL,
                                 cex.Num = 0.6,
                                 ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
          } else {
           HeatmapWithTextLabels(cors[[i.row]], dimnames(cors[[i.row]])[[2]], 
                                 dimnames(cors[[i.row]])[[2]],
                                 main=paste(letter, Titles[[i.col]]), InvertColors=invertColors, 
                                 zlim=c(-1,1.0),
                                 NumMatrix = if (printAdj) signif(cors[[i.row]], 2) else NULL,
                                 cex.Num = 0.6,
                                 ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
          }
        } else {
           HeatmapWithTextLabels(abs(cors[[i.row]]), dimnames(cors[[i.row]])[[2]],
                                 dimnames(cors[[i.row]])[[2]],
                                 main=paste(letter, Titles[[i.col]]), InvertColors=invertColors, 
                                 zlim=c(0,1.0),
                                 NumMatrix = if (printAdj) signif(abs(cors[[i.row]]), 2) else NULL,
                                 cex.Num = 0.6,
                                 ColorLabels = ColorLabels, colors = colors, SetMargins = FALSE, ...);
        }
      } else
      {
        cor1 = cors[[i.col]];
        cor2 = cors[[i.row]];
        cor.dif = (cor1 - cor2)/2;
        d = tanh((cor1 - cor2) / (abs(cor1) + abs(cor2))^2);
        # d = abs(corPC1 - corPC2) / (abs(corPC1) + abs(corPC2));
        dispd = cor.dif;
        if (plotCPMeasure) dispd[upper.tri(d)] = d[upper.tri(d)];
        if (i.row>i.col)
        {
          if (IncludeSign)
          {
            half = as.integer(length(colors)/2);
            if (invertColors) { range = c(1:(half+1));} else {range = c(half:length(colors)); }
            halfColors = colors[range];
          } else {
            halfColors = colors;
          }
          if (printCPVals) {
            printMtx = matrix(paste(".", as.integer((1-abs(dispd))*100), sep = ""), 
                               nrow = nrow(dispd), ncol = ncol(dispd));
            printMtx[printMtx==".100"] = "1";
          } else { 
            printMtx = NULL; 
          }
          if (sum( (1-abs(dispd)<CPzlim[1]) | (1-abs(dispd)>CPzlim[2]) )>0)
            warning("PlotCorPCs: Correlation preservation data out of zlim range!");
          if (plotCPMeasure) {
             HeatmapWithTextLabels(1-abs(dispd), dimnames(cor1)[[2]], dimnames(cor2)[[2]], 
                                   main=paste(letter, "UT: Cor.Pres\nLT: 1-Cor.Diff"),
                                   InvertColors=invertColors, 
                                   ColorLabels = ColorLabels, zlim = CPzlim, colors = halfColors,
                                   SetMargins = FALSE, 
                                   NumMatrix = printMtx, cex.Num = CPcex, ...);
          } else {
             HeatmapWithTextLabels(1-abs(dispd), dimnames(cor1)[[2]], dimnames(cor2)[[2]], 
                                   main=paste(letter, "Preservation"), InvertColors=invertColors, 
                                   ColorLabels = ColorLabels, zlim = CPzlim, colors = halfColors,
                                   SetMargins = FALSE,  NumMatrix= printMtx, cex.Num = CPcex, ...);
          }
        } else {
          if (plotCPMeasure) {
             dp = 1-abs(d);
             method = "Cor.Pres.:";
          } else {
             dp = 1-abs(cor.dif); 
             method = "Preservation:";
          }
          diag(dp) = 0;
          if (plotMeans) {
            sum_dp = mean(dp[upper.tri(dp)]);
            means = apply(dp, 2, sum)/(ncol(dp)-1);
            if (plotErrors) {
               Errors = sqrt( (apply(dp^2, 2, sum)/(ncol(dp)-1) - means^2)/(ncol(dp)-2));
            } else {
               Errors = NULL; 
            }
            BarplotWithTextLabels(means, dimnames(cor1)[[2]], 
                                 main=paste(letter, "D=", signif(sum_dp,2)), 
                                 ylim=c(0,1),
                                 ColorLabels = ColorLabels, Colored = ColoredBarPlot,
                                 SetMargins = FALSE, Errors = Errors, ... )
          } else {
            sum_dp = sum(dp[upper.tri(dp)]);
            BarplotWithTextLabels(dp, names(cor1),
                                 main=paste(letter, method, "sum = ", signif(sum_dp,3)), 
                                 ylim=c(0,dim(dp)[[1]]),
                                 ColorLabels = ColorLabels, Colored = ColoredBarPlot, 
                                 SetMargins = FALSE, ... )
          }
        }
      }
    }
  }
}
  


#=================================================================================================
#
# Dendrogram orderings
#
#=================================================================================================


dendrogramOrderings = function(dendro, maxMerge = 18, verbose = 1, indent = 0)
{

  spaces = indentSpaces(indent)

  if (verbose>0) printFlush(paste(spaces, "Calculating dendrogramOrderings.."));
  nMerge = length(dendro$height);
  nSingl = nMerge + 1;
  mleft = dendro$merge[,1];
  mright = dendro$merge[,2];

  if ((verbose>0) & (nMerge>maxMerge)) 
     printFlush(paste("..FYI: limiting number of flipped merges to", maxMerge)); 

  # First step: get the list of L and R singletons below each merge, i.e. for each merge and each
  # singleton the array sides[, merge, side] gives singletons to the side side of merge merge.
  # The number of those sigletons is recorded in nSides[merge, side]

  sides = array(0, dim = c(nSingl, ncol = nMerge, 2));

  nSides = matrix(0, nrow = nMerge, ncol = 2);

  traceBack = rep(0, times = nMerge);
  curSide = rep(1, times = nMerge);	# 1=left, 2=right, 3=all done.

  nTrace = 1;		# Number of traces includes current merge.
  traceBack[1] = nMerge;
  curMerge = nMerge;

  if (verbose>1) printFlush(paste(spaces, " ..finding `sides' of each merge.."));
  while ((curSide[nMerge]<3) | (nTrace > 1))
  {
    # print(paste("nTrace:", nTrace, "traceBack:", paste(traceBack, collapse = ", ")));
    # print(paste("curMerge:", curMerge, "curSide", paste(curSide, collapse = ", ")));
    if (curSide[curMerge] < 3)
    { 
      # print(paste("Processing side", curSide[curMerge]));
      # Both sides not yet processed.
      if (dendro$merge[curMerge, curSide[curMerge]] < 0)
      {
        # Have a sigleton on the current side. Record the singleton as being on the current side of the
        # current merge and all merges in the traceback.
        # print("Singleton.");
        singleton = -dendro$merge[curMerge, curSide[curMerge]];
        curSide[curMerge] = curSide[curMerge] + 1;
        iTrace = nTrace;
        while (iTrace > 0)
        {
          merge = traceBack[iTrace];
          mergeSide = curSide[merge] - 1;
          i = nSides[merge, mergeSide] + 1;
          sides[i, merge, mergeSide] = singleton;
          nSides[merge, mergeSide] = i;
          iTrace = iTrace - 1;
        }
        # print("nSides after:"); print(nSides);
        # print("sides[,,1]:"); print(sides[,,1]);
        # print("sides[,,2]:"); print(sides[,,2]);
      } else {
        # Record traceback and do the next merge.
        # print("Branch.")
        x = curSide[curMerge];
        curSide[curMerge] = x+1;
        curMerge = dendro$merge[curMerge, x];
        nTrace = nTrace + 1;
        traceBack[nTrace] = curMerge;
      }
    } else {
      # print("Taking one step back.");
      # Both sides were already processed. Return one step back.
      if (nTrace > 1)
      {
        nTrace = nTrace - 1;
        curMerge = traceBack[nTrace];
      } else {
        printFlush("All done, apparently!");
      }
    }
    # print("Press enter"); scan();
  }

  # Generate all orderings by recursively flipping left and right sides of each merge.

  # It is much easier to flip rankings than orderings.

  nRankings = 2**(min(nMerge, maxMerge)-1);
  ranks = matrix(0, nrow = nSingl, ncol = nRankings);

  cuRanking = 1;  # Ranking being currently worked on

  # Convert the dendorgram order into ranks
  dendrRanking = rep(0, nSingl);
  dendrRanking[dendro$order] = c(1:nSingl);

  # Rankings used as a base before flipping for every merge

  baseRankings = matrix(rep(dendrRanking, nMerge), nrow = nSingl, ncol = nMerge);

  flipState = rep(0, nMerge);

  # The very last merge can be left out.
  merge = nMerge-1;

  if (verbose>1) printFlush(paste(spaces, " ..generating rankings.."));
  counter = 0; countedMerge = nMerge - 5;

  while (merge<nMerge) 
  {
    # print(paste("merge:", merge, "flipState:", paste(flipState, collapse = ", ")));
    if (flipState[merge]==0)
    {
      # First pass through this merge.
      flipState[merge] = 1;
      # Copy base order from above
      if (merge<nMerge) baseRankings[,merge] = baseRankings[,merge+1];
      if ((merge>1) & (merge>nMerge-maxMerge+1))
      {
        # Reset everything below
        flipState[1:(merge-1)] = 0;
        # ...and go one step below
        merge = merge-1;
      } else {
        # No merges below: save current ranking
        ranks[,cuRanking] = baseRankings[, merge];
        cuRanking = cuRanking + 1;
        # Perform flip at the merge
        nLeft = nSides[merge, 1];
        nRight = nSides[merge, 2];
        left = sides[1:nLeft, merge, 1];
        right = sides[1:nRight, merge, 2];
        flippedRanking = baseRankings[, merge];
        #print(paste("Flip: left=", paste(left, collapse = ","), "; right=", 
                     #paste(right, collapse = ", ")));
        #print(paste("Ranks before:", paste( flippedRanking, collapse = ", ")));
        flippedRanking[left] = baseRankings[left, merge] + nRight;
        flippedRanking[right] = baseRankings[right, merge] - nLeft;
        #print(paste("Ranks after", paste( flippedRanking, collapse = ", ")));
        # Save the flipped order
        ranks[, cuRanking] = flippedRanking;
        cuRanking = cuRanking + 1;
        flipState[merge] = 2;
      }
    } else if (flipState[merge]==1)
    {
      if ((verbose > 2) & (merge == countedMerge))
      { 
        printFlush(paste(spaces, "  ..counter on counted merge:", counter, "of",
                         2**(nMerge-countedMerge-1)));
        counter = counter + 1;
      }
      flipState[merge] = 2;
      # perform the flip of this merge
      nLeft = nSides[merge, 1];
      nRight = nSides[merge, 2];
      left = sides[1:nLeft, merge, 1];
      right = sides[1:nRight, merge, 2];
      flippedRanking = baseRankings[, merge];
      #print(paste("Flip: left=", paste(left, collapse = ","), "; right=", 
                   #paste(right, collapse = ", ")));
      #print(paste("Ranks before:", paste( flippedRanking, collapse = ", ")));
      flippedRanking[left] = baseRankings[left, merge] + nRight;
      flippedRanking[right] = baseRankings[right, merge] - nLeft;
      #print(paste("Ranks after", paste( flippedRanking, collapse = ", ")));
      # Reset all flips below
      flipState[1:(merge-1)] = 0;
      # Set new base order
      baseRankings[, merge] = flippedRanking;
      # and go to next lower merge
      merge = merge -1;
    } else {
      # Both flips on this merge were performed and saved... go one step back.
      merge = merge + 1;
    }
  }

  collect_garbage();

  if (verbose>1) printFlush(paste(spaces, " ..converting rankings into oderings.."));

  
# The following lines do not work for large dendrograms, presumably because subsetting creates a copy of
# the whole object.

  orders = ranks;
  for (i in 1:nRankings) orders[ranks[, i], i] = c(1:nSingl);

# Replace it with slightly more complicated but block code.
# Want to convert everything into a flat vector.

  #BatchSize = as.integer(1000000/nSingl);

  #nBatch = as.integer( (nRankings-1)/BatchSize ) + 1;
  #for (batch in 1:nBatch)

  #nData = nRankings * nSingl;
  #flatRanks = as.integer( c(0:(nData-1)) / nSingl ) * nSingl + as.vector(ranks);
  #orders = flatRanks;
  #orders[flatRanks] = rep(c(1:nSingl), times = nRankings);
  #dim(orders) = c(nSingl, nRankings);

  list(orders = orders, ranks = ranks);
}

#-----------------------------------------------------------------------------------------------
#
# outlierMeasure
#
#-----------------------------------------------------------------------------------------------

# This function calculates the measure of being an outlier for each sample in each column and averages
# them. It returns a vector of the length = nrow(x) giving the mean measure of outlierishness for each
# sample.   
        
outlierMeasure = function(x)
{
    x = as.matrix(x);
    mx = apply(x, 2, median, na.rm = TRUE)
    mxMat = matrix(mx, nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
    madxMat = matrix(apply(x, 2, mad, na.rm = TRUE), nrow = nrow(x), ncol = ncol(x), byrow = TRUE);
    ux = as.matrix(abs((x - mxMat)/(9 * qnorm(0.75) * madxMat)))
    aval = ifelse(ux <= 1, 1, 0)
    outMeas = ux^2;   # Could multiply by 1-aval as well 
    apply(outMeas, 1, mean, na.rm = TRUE);
}


# Related: this function plots a legend but cleans up the legend rectangle
# first. Do not give it a plot argument, it will cause an error.

legendClean = function(...)
{
  box = legend(..., plot = FALSE)$rect
  rect(box$left, box$top-box$h, box$left +box$w, box$top, col = "white", border = "white");
  legend(...);
}

initLib = function(libname, force = FALSE)
{
  os = R.Version()$os;
  recOS = c("linux-gnu", "mingw32");
  ios = match(os, recOS);
  if (is.na(ios)) stop("Unrecognized OS.");
  exts = c(".so", ".dll");
  filename = paste(libname, exts[ios], sep="");
  if (force && is.loaded(filename)) dyn.unload(filename);
  if (!is.loaded(filename)) dyn.load(filename);
}


# spaste: paste without spaces.

spaste = function(...)
{
  paste(..., sep="");
}


# Plot enrichment barplot of modules in a yes/no indicator variable

plotEnrichmentBarplot = function(indicator, moduleLabels, moduleLevels, moduleColors,
                xlab = "Module number",
                ylab = "-log(p-value)",
                cex.lab = 1,
                ...)
{
   nGenes = length(moduleLabels);
   geneColors = moduleLabels[indicator];
   presentMods = sort(unique(geneColors))
   nPresMods = length(presentMods)
   pValues = rep(NA, nPresMods)
   nIntGenes = rep(0, nPresMods)
   xGeneColors = rep(-1, nGenes)
   xGeneColors[indicator] = geneColors;
   for (imod in 1:nPresMods)
   {
     mod = presentMods[imod]
     testMatrix = matrix( c(
          sum(xGeneColors>-1 & moduleLabels==mod),
          sum(xGeneColors==-1 & moduleLabels==mod),
          sum(xGeneColors>-1 & moduleLabels!=mod),
          sum(xGeneColors==-1 & moduleLabels!=mod)), 2, 2)
     pValues[imod] = fisher.test(testMatrix, alternative = "greater")$p.value;
     nIntGenes[imod] = testMatrix[1,1];
   }

   pValues[pValues==0] = 1e-300;

   max = max(-log10(pValues), -log10(0.01 / nPresMods) + 1);
   min = 0;

   mp = barplot(-log10(pValues), col = moduleColors[match(presentMods, moduleLevels)], names.arg = FALSE,
                xlab = xlab, ylab = ylab, ylim = c(min, max), ...)
   text(x = mp, y = pmax(-log10(pValues) + 1, rep(-log10(0.01 / nPresMods) + 1, length(pValues))),
        labels =  nIntGenes, adj = c(0.5, 0), xpd = TRUE, cex = 0.8 )
   text(x = mp, y = min - (max-min)*0.05, labels =presentMods, adj = c(0.5, 1), xpd = TRUE, cex = cex.lab )
   abline(h = -log10(0.01 / nPresMods), col = "red")
}
 

# Row-wise consensus: must have equal signs
# i.e. consensus across columns, returning one value for each row (for higher-dimensional arrays, returns
# an array of one dimension lower of consensi across the last dimension)
# attempts to ignore NA's 

rowConsensus = function(data)
{
  nDim = length(dim(data))
  if (nDim < 2) stop("'data' must be a matrix or array with at least 2 dimensions.");
  if (nDim==2) data = as.matrix(data);
  dimd = dim(data);
  nc = dimd[length(dimd)];
  sign = apply(sign(data), c(1:(nDim-1)), sum, na.rm = TRUE)
  sign[is.na(sign)] = 0;
  sign[abs(sign) < nc] = 0;
  consensus = apply(abs(data), c(1:(nDim-1)), min, na.rm = TRUE) * sign/nc;
  consensus;
}

plotHistogramAndCF = function(data, nTicks = NULL, ...)
{
   h = hist(data, ...);
   ymax = par("yaxp")[2];
   n = par("yaxp")[3];
   #if (!(n %in% c(2,4,5,10)))
   if (is.null(nTicks)) nTicks = n;
   sum = sum(h$counts);
   lines(h$mids, cumsum(h$counts) / sum * ymax, col="red");
   axis(side = 4, at = c(0:nTicks)/nTicks * ymax, labels = signif(c(0:nTicks)/nTicks, 2))
   addGrid();
   invisible(h);
}

# The [bi]corAndPvalue functions have been moved to WGCNA package

#===============================================================================================
#
# Plot of a set of variables vs. genomic location.
#
#===============================================================================================

.plotQTLtypes = c("lines", "points");

plotQTL = function(data, chr, bp, colors = 1, type = "lines", bpBasedSpacing = FALSE,
                   lty = 1, lwd = 1, ylim = NULL, xLabOffset = 0.04, 
                   gridColor = "grey", grid.lty = 2, ... )
{

  itype = charmatch(type, .plotQTLtypes);
  if (is.na(itype))
    stop(paste("Unrecognized 'type'. Recognized values are", paste(.plotQTLtypes, collapse = ", ")));

  nchr = as.numeric(chr);
  order = order(nchr, bp);
  data = as.matrix(data);
  nCols = ncol(data);
  nSNPs = length(order)
  x = data[order, , drop = FALSE];
  if (bpBasedSpacing)
  {
    chromosomeNumbers = sort(unique(nchr));
    nChromo = length(chromosomeNumbers);
    chromoLength = tapply(bp, nchr, max);
    dividers = cumsum(chromoLength);
    totalLength = sum(chromoLength);
    snpGlobalPositions = bp;
    addLength = 0; index = 1;
    for (ch in chromosomeNumbers)
    {
      snpGlobalPositions[ch==nchr] = bp[ch==nchr] + addLength;
      addLength = addLength + chromoLength[index];
      index = index + 1;
    }
    plotX = snpGlobalPositions;
    breaksX = c(0.5, dividers + 0.5);
    breaks = dividers[-length(dividers)] + 0.5;
  } else {
    plotX = c(1:nSNPs);
    cx = chr[order];
    cx2 = cx[-1];
    breaks = which(cx[-nSNPs]!=cx2) + 0.5;
    breaksX = c(0.5, breaks, nSNPs + 0.5);
  }
  mids = (breaksX[-1] + breaksX[-length(breaksX)])/2;
   
  if (is.null(colors)) colors = matrix(1, nSNPs, nCols);
  if ((length(colors)==nCols) || (length(colors)==1)) 
     colors = matrix(colors, nSNPs, nCols, byrow = TRUE);
  if (length(lty==1)) lty = rep(lty, nCols);
  if (length(lwd==1)) lwd = rep(lwd, nCols);


  if (is.null(ylim))
  {
     min = min(data, na.rm = TRUE);
     max = max(data, na.rm = TRUE);
     ylim = c(min, max);
  }

  if (itype==1)
  {
    plot(plotX, x[, 1], col = colors[, 1], type = "l", lty = lty[1], lwd = lwd[1], 
         ylim = ylim, xaxt = "none", xlab = "Chromosome", xaxs = "i", ...)

    if (nCols > 1) for (col in 2:nCols)
      lines(plotX, x[, col], col = colors[, col], lty = lty[col], lwd = lwd[col]);
  } else {
    plot(plotX, x[, 1], col = colors[, 1], pch = ".",
         ylim = ylim, xaxt = "none", xlab = "Chromosome", xaxs = "i", ...)

    if (nCols > 1) for (col in 2:nCols)
      points(plotX, x[, col], col = colors[, col], pch = ".");
  }


  box = par("usr");
  for (b in breaks)
    lines(x = rep(b, 2), y = box[3:4], col = gridColor, lty = grid.lty );

  text(mids, rep(box[3] - xLabOffset * (box[4]-box[3]), length(mids)), sort(unique(chr)), xpd = TRUE);

  invisible(plotX);
}

distributionOnSNPs2 = function(dataList, snpChr, snpBp, window, filterType = "gaussian", fun = sum, 
                     normalize = FALSE)
{
  nData = length(dataList);
  nSNPs= length(snpChr);
  distribution = matrix(0, nSNPs, nData);
  snpChrLevs = sort(unique(snpchr));

  ft = charmatch(filterType, .filterTypes);
  if (is.na(ft))
    stop(paste("Unrecognized filter type: ", filterType,
               "\nRecognized values are ", paste(.filterTypes, collapse = ", ")))

  for (m in 1:nData)
  {
    sourceChr = dataList[[m]]$chr;
    sourceBp = dataList[[m]]$bp;
    for (c in snpChrLevs)
    {
      chrSnps = snpChr==c;

      targetBp = snpBp[chrSnps];
      sourceBpC = sourceBp[sourceChr==c];
      finiteBp = is.finite(sourceBpC);
      sourceBpC = sourceBpC[finiteBp];

      # outer(x,y,`-`)[i,j] = x[i]-y[j]
      distMat = outer(sourceBpC, targetBp, `-`);
      dataMat = matrix(dataList[[m]]$data[sourceChr==c][finiteBp], nrow(distMat), ncol(distMat));

      if (ft==1)
      {
         weight = exp(-distMat^2/window^2);
      } else {
         weight = matrix(as.numeric(abs(distMat) <= window), nrow(distMat), ncol(distMat));
      }
      d = apply(weight * dataMat, 2, fun);
      if (normalize)
      {
         sw = apply(weight, 2, fun);
         d = d/sw;
         d[sw==0] = NA;
      }
      distribution[chrSnps, m] = d;
    }
  }
  distribution;
}

#================================================================================================
#
# qtlPeaks
#
#================================================================================================

# Find peaks in data (qtl) that have no near higher QTL (within tol). 
# Returns a logical vector that indicates whether
# each entry is a peak or not.

qtlPeaks = function(qtl, chr, bp, minQTL, window, tol = 0.001)
{
  order = order(chr, bp)
  nSNPs = length(qtl);
  qtlL = c(qtl[-1], 0);
  qtlR = c(0, qtl[-nSNPs]);

  chrL = c(chr[-1], chr[nSNPs]);
  chrR = c(chr[1], chr[-nSNPs]);

  peaks = (qtl >= qtlL | chrL != chr) & (qtl >= qtlR | chr != chrR)

  peaks[is.na(peaks)] = FALSE;
  peaks[peaks][qtl[peaks] < minQTL] = FALSE;

  # remove neighboring peaks that have the exact same LOD
  nPeaks2 = sum(peaks);
  delPeaks2 = rep(FALSE, nPeaks2);
  qtl2 = qtl[peaks];
  delPeaks2[qtl2[-1]==qtl2[-nPeaks2]] = TRUE;
  # This looks a bit confusing :)
  peaks[peaks][delPeaks2] = FALSE;

  peakP = qtl[peaks];
  nPeaks = sum(peaks);
  peakchr = chr[peaks];
  peakBp = bp[peaks];

  # Remove peaks that are not highest in their neighborhood
  delPeaks = rep(FALSE, nPeaks);
  for (p in 1:nPeaks)
  {
    near = chr==peakchr[p] & abs(peakBp[p] - bp) <= window
    maxNearQTL = max(qtl[near], na.rm = TRUE);
   if (maxNearQTL > peakP[p] + tol) delPeaks[p] = TRUE;
  }

  peaks[peaks][delPeaks] = FALSE;

  peaks;
}

#===============================================================================================
#
# labeled boxplot
#
#===============================================================================================

labeledBoxplot = function(x, names, srt = 45, adj = c(1, 0.5), namesColor = 1, cex.names = 1, 
                          verbose = FALSE, ...)
{
  x = as.list(x);
  if (verbose)
  {
    lengths = lapply(x, length);
    n = length(x);
    g = rep(1:n, lengths);
    xx = as.vector(unlist(x));
    out = verboseBoxplot(xx, g, names = FALSE, ...);
  } else 
    out = boxplot(x, names = FALSE, ...);
  n = length(x);
  box = par("usr");
  xText = c(1:n);
  yText = rep(box[3] - 0.02 * (box[4] - box[3]), n);

  text(xText, yText, names, srt = srt, adj = adj, col = namesColor, cex = cex.names, xpd = TRUE);
  out;
}


#===============================================================================================
#
# labeled plot
#
#===============================================================================================

labeledPlot = function(x, y, names, srt = 45, adj = c(1, 0.5), namesColor = 1, cex.names = 1, ...)
{
  plot(x, y, xaxt = "n", ...);
  n = length(x);
  box = par("usr");
  xText = c(1:n);
  yText = rep(box[3] - 0.02 * (box[4] - box[3]), n);

  text(xText, yText, names, srt = srt, adj = adj, col = namesColor, cex = cex.names, xpd = TRUE);
}


#===============================================================================================
#
# labeled barplot, version 2 (somewhat different aim than the WGCNA function)
#
#===============================================================================================

labeledBarplot2 = function(x, names, g = NULL, horiz = FALSE, 
                           srt = if (horiz) 0 else 45, 
                           adj = if (horiz) c(1, 0.5) else c(1, 0.5), 
                           namesColor = 1, cex.names = 1, 
                           addGuide = TRUE, guideColor = "grey30", guideLty = 2, ...)
{

  if (!is.null(g)) {
    mp = verboseBarplot(x, g, names.arg = NA, horiz = horiz, ...);
    heights = attr(mp, "height");
  } else {
    mp = barplot(x, names.arg = NA, horiz= horiz, ...);
    heights = x;
    attr(mp, "height") = x;
    attr(mp, "stdErr") = rep(0, length(x));
  }

  box = par("usr");

  if (horiz)
  {
    yText = mp;
    xText = rep(box[1] - 0.01 * (box[2] - box[1]), length(mp));
  } else {
    xText = mp;
    yText = rep(box[3] - 0.02 * (box[4] - box[3]), length(mp));
  }

  text(xText, yText, names, srt = srt, adj = adj, col = namesColor, cex = cex.names, xpd = TRUE);
  if (addGuide) for (i in 1:length(mp))
    if (horiz)
    {
      lines(c(min(0, heights[i]), box[1] - 0.02 * (box[2] - box[1])),  rep(mp[i], 2), 
            col = guideColor, lty = guideLty);
    } else {
      lines(rep(mp[i], 2), c(min(0, heights[i]), box[3] - 0.02 * (box[4] - box[3])), 
            col = guideColor, lty = guideLty);
    }

  invisible(mp);
  
}


#===============================================================================================
#
# Gaussian smoothing filter
#
#===============================================================================================

smoothGauss = function(x,y, xtest, width, normalize = TRUE)
{
  dist = outer(x, xtest, "-");
  weight = exp(-dist^2/(2*width^2));
  yMat = matrix(y, length(x), length(xtest));

  if (normalize)
  {
     val = apply(yMat * weight, 2, sum)/apply(weight, 2, sum)
  } else
     val = apply(yMat * weight, 2, sum)
  val;
}

#===============================================================================================
#
# Convert date in the form month/day/year into number of days
# Accepts the separator (usually /, -, or .) and the position of day, month, and year in the date
#
#===============================================================================================

linDate = function(dateStr, separator= "/", day = 2, month = 1, year = 3)
{
  # Here we assume date specified as month/day/year

  daysInMonth = c(31,28,31,30,31,30,31,31,30,31,30,31);
  sumDIM = c(0, cumsum(daysInMonth));
  daysInYear = sum(daysInMonth);

  split = strsplit(dateStr, split = separator, fixed = TRUE);
  nSubj = length(split);
  date = rep(NA, nSubj);
  for (s in 1:nSubj)
  {
    date[s] = as.numeric(split[[s]][day]) + sumDIM[as.numeric(split[[s]][month])] +
              as.numeric(split[[s]][year]) * daysInYear;
  }
  date;
}

#===============================================================================================
#
# Convert age in the form xxyxxm to number of years. Assumes all months have the same length.
#
#===============================================================================================

linAge = function(ageStr)
{
   # Assume ageStr contains age in the form xxyxxm where xx are numbers and y and m are year and month
   # delimiters

   split = strsplit(ageStr, split = c("y"), fixed = TRUE)

   nSubj = length(split);
   age = rep(NA, nSubj);
   for (s in 1:nSubj)
   {
     age[s] = as.numeric(split[[s]][1]);
     if (length(split[[s]])==2)
     {
        months = strsplit(split[[s]][[2]], split = "m", fixed = TRUE)[[1]];
        age[s] = age[s] + as.numeric(months)/12;
     }
   }
   age;
}

#=========================================================================================================
#
# Simulation of a causal model. From Old/GeneExpressionSimulation/.
# Note: the causal model simulation has been modified to produce a consistent interpretation of the path
# coefficients.
#
#=========================================================================================================

#----------------------------------------------------------------------------
#
# CausalChildren
#
#----------------------------------------------------------------------------
# Note: The returned vector may contain multiple occurences of the same child.

CausalChildren = function(Parents, Cause)
{
  nNodes = dim(Cause)[[1]];

  # print(paste("Length of Parents: ",length(Parents)));
  if (length(Parents)==0) return(NULL);

  Child_ind = apply(as.matrix(abs(Cause[, Parents])), 1, sum)>0;
  if (sum(Child_ind)>0)
  {
     children = c(1:nNodes)[Child_ind] 
  } else {
     children = NULL;
  }
  children;
}

#----------------------------------------------------------------------------
#
# simulateCausal
#
#----------------------------------------------------------------------------
#
# Given a set of causal anchors, this function creates a network of vectors that should satisfy the
# causal relations encoded in the causal matrix Cause, i.e. Cause[j,i] is the causal effect of vector i on
# vector j. (This is Jason's convention.)

# The function starts by initializing all vectors to noise given in the noise specification. (The noise
# can be specified for each vector separately.) Then it runs the standard causal network signal
# propagation and returns the resulting vectors.

simulateCausal = function(Cause, AnchorInd, AnchorVecs, Noise, verbose = 2, indent = 0)
{
  spaces = indentSpaces(indent);

  if (verbose>0) printFlush(paste(spaces, "Creating seed vectors..."));
  nNodes = dim(Cause)[[1]];
  nSamples = dim(AnchorVecs)[[1]];
  
  if (length(AnchorInd)!=dim(AnchorVecs)[[2]])
  {
    printFlush(paste("CreateSeedVectors: Error: Length of AnchorInd must equal",
            "the number of vectors in AnchorVecs."));
    stop();
  }
  if (length(Noise)!=nNodes)
  {
    printFlush(paste("CreateSeedVectors: Error: Length of Noise must equal",
            "the number of nodes as given by the dimension of the Cause matrix."));
    stop();
  }

  # Initialize all node vectors to noise with given standard deviation

  NodeVectors = matrix(0, nrow = nSamples, ncol = nNodes);
  for (i in 1:nNodes)
  {
    NodeVectors[,i] = rnorm(n=nSamples, mean=0, sd=Noise[i]);
  }

  Levels = rep(0, times = nNodes);

  # Calculate levels for all nodes: start from anchors and go through each successive level of children

  level = 0;
  Parents = AnchorInd;
  Children = CausalChildren(Parents = Parents, Cause = Cause);
  if (verbose>1) printFlush(paste(spaces, "..Determining level structure..."));
  while (!is.null(Children))
  {
    # print(paste("level:", level));
    # print(paste("   Parents:", Parents));
    # print(paste("   Children:", Children));
    level = level + 1;
    if ((verbose>1) & (level/10 == as.integer(level/10))) 
          printFlush(paste(spaces, "  ..Detected level", level));
    #printFlush(paste("Detected level", level));
    Levels[Children] = level;
    Parents = Children;
    Children = CausalChildren(Parents = Parents, Cause = Cause);
  }

  HighestLevel = level;

  # Generate the whole network

  if (verbose>1) printFlush(paste(spaces, "..Calculating network..."));
  NodeVectors[,AnchorInd] = NodeVectors[,AnchorInd] + AnchorVecs;
  for (level in (1:HighestLevel))
  {
    if ( (verbose>1) & (level/10 == as.integer(level/10)) ) 
      printFlush(paste(spaces, " .Working on level", level));
    #printFlush(paste("Working on level", level));
    LevelChildren = c(1:nNodes)[Levels==level]
    for (child in LevelChildren) 
    {
      LevelParents = c(1:nNodes)[Cause[child, ]!=0]
      parentEffect = rep(0, nSamples);
      for (parent in LevelParents)
        parentEffect = parentEffect + Cause[child, parent]*scale(NodeVectors[,parent]);
      NodeVectors[, child] = NodeVectors[, child] + parentEffect;
    }
  }

  Nodes = list(Vectors = scale(NodeVectors), Cause = Cause, Levels = Levels, AnchorInd = AnchorInd);
  Nodes;
} 



#===================================================================================================== 
#
# partial correlation of columns of x conditioned on A
#
#=====================================================================================================

matrixPCor = function(x, A, corFnc = "cor", corOptions = "use = 'p'")
{
  cx = eval(parse(text = spaste(corFnc, "(x, ", corOptions, ")")));
  ca = as.numeric(eval(parse(text = spaste(corFnc, "(x, A, ", corOptions, ")"))));

  pc = ( cx - ca %*% t(ca)) /  sqrt( (1-ca^2) %*% t(1-ca^2) )
  pc;
};

zeo = function(x, A, corFnc = "cor", corOptions = "use = 'p'")
{
  n = ncol(x);
  cx = eval(parse(text = spaste(corFnc, "(x, ", corOptions, ")")));
  ca = as.numeric(eval(parse(text = spaste(corFnc, "(x, A, ", corOptions, ")"))));
  caMatCol = matrix(ca, n, n);
  caMatRow = matrix(ca, n, n, byrow = TRUE);
  zeo = (caMatCol - cx * caMatRow) / sqrt( (1-caMatRow^2) * (1-cx^2) );
  diag(zeo) = NA;
  zeo;
}

#===================================================================================================
#
# verbosePairs
#
#===================================================================================================

verbosePairs = function(x, mar1 = c(2,2,1,1), names = colnames(x), cex.names = 3, sample = NULL,
                        cex.lab = 1, cex.main = 1, cex.axis = 1, breaks = 100, corFnc = "cor", ...)
{
  x = as.matrix(x);
  nSets = ncol(x);
  par(mfrow = c(nSets, nSets));
  par(mar = mar1)
  if (is.null(names)) names = c(1:nSets);

  for (s1 in 1:nSets) for (s2 in 1:nSets)
  {
    if (s1==s2)
    {
      #plot(c(0,1), c(0,1), type = "n", xaxt = "n", yaxt = "n");
      #text(0.5, 0.5, adj = c(0.5, 0.5), cex = cex.names, labels = names[s1]);
      if (corFnc=="bicor")
      {
        mean = median(x[, s1], na.rm = TRUE);
        sd = mad(x[, s1], na.rm = TRUE);
        mainLine = spaste("\nmedian: ", signif(mean, 2), ", mad: ", signif(sd, 2));
      } else {
        mean = mean(x[, s1], na.rm = TRUE);
        sd = sd(x[, s1], na.rm = TRUE);
        mainLine = spaste("\nmean: ", signif(mean, 2), ", sd: ", signif(sd, 2));
      }
      par(mar = mar1 + c(0,0,2,0))
      hist(x[, s1], breaks = breaks, main = spaste(names[s1], mainLine), xlab = "", cex.lab = cex.lab,
           cex.axis = cex.axis, cex.main = cex.main);
    } else {
      par(mar = mar1)
      verboseScatterplot(x[, s2], x[, s1], sample = sample, cex.lab = cex.lab,
                         cex.axis = cex.axis, cex.main = cex.main, xlab = "", ylab = "",
                         corFnc = corFnc, ...);
    }
  }
}


#=======================================================================================================
#
# invertEmpiricalDistribution
#
#=======================================================================================================

# Inverts an empirical distribution

invertEmpiricalDistribution = function(data, nBreaks = 100, smoothWindow = 3, plot = FALSE)
{
  if (!is.null(dim(data))) stop("'data' must be a vector.");
  n = length(data);
  h = hist(data, breaks = nBreaks, plot = plot);


  step = h$mids[2] - h$mids[1];
  smooth = smoothGauss(h$mids, h$counts, h$mids, smoothWindow * step);
  smooth[smooth < 0] = 0;
  sums = c(0, cumsum(smooth))/max(cumsum(smooth));

  if (plot) lines(h$mids, smooth, col = "red")
  if (plot) lines(h$breaks, sums * max(h$counts), col = "blue")
  
  y = h$breaks
  x = sums;
  # fit = lm(y~ns(x, df = nBreaks));

  if (FALSE)
  {
    nInvBreaks = nBreaks;  
    invX = seq(from = 0, to=1, length.out = nInvBreaks + 1)
    inverse = predict(fit, newdata = data.frame(x=invX))

    inverse[1] = h$breaks[1];
    inverse[nInvBreaks + 1] = h$breaks[length(h$breaks)];

    if (plot) lines(inverse, invX*max(h$counts), col = "green");

    plot(sums, h$breaks, type = "l", col = 1);
    lines(invX, inverse+1, col = 2);
  }

  #fit;

  list(x = x, y = y);
}

# Make a slightly modified version of the function approx to supress a warning

approx = function (x, y = NULL, xout, method = "linear", n = 50, yleft, 
    yright, rule = 1, f = 0, ties = mean) 
{
    x <- xy.coords(x, y)
    y <- x$y
    x <- x$x
    nx <- length(x)
    method <- pmatch(method, c("linear", "constant"))
    if (is.na(method)) 
        stop("invalid interpolation method")
    stopifnot(is.numeric(rule), (lenR <- length(rule)) >= 1, 
        lenR <= 2)
    if (lenR == 1) 
        rule <- rule[c(1, 1)]
    if (any(na <- is.na(x) | is.na(y))) {
        ok <- !na
        x <- x[ok]
        y <- y[ok]
        nx <- length(x)
    }
    if (!identical(ties, "ordered")) {
        if (length(ux <- unique(x)) < nx) {
#            if (missing(ties)) 
#                warning("collapsing to unique 'x' values")
            y <- as.vector(tapply(y, x, ties))
            x <- sort(ux)
            nx <- length(x)
        }
        else {
            o <- order(x)
            x <- x[o]
            y <- y[o]
        }
    }
    if (nx <= 1) {
        if (method == 1) 
            stop("need at least two non-NA values to interpolate")
        if (nx == 0) 
            stop("zero non-NA points")
    }
    if (missing(yleft)) 
        yleft <- if (rule[1] == 1) 
            NA
        else y[1L]
    if (missing(yright)) 
        yright <- if (rule[2] == 1) 
            NA
        else y[length(y)]
    stopifnot(length(yleft) == 1, length(yright) == 1, length(f) == 
        1)
    if (missing(xout)) {
        if (n <= 0) 
            stop("'approx' requires n >= 1")
        xout <- seq.int(x[1L], x[nx], length.out = n)
    }
    y <- .C("R_approx", as.double(x), as.double(y), as.integer(nx), 
        xout = as.double(xout), as.integer(length(xout)), as.integer(method), 
        as.double(yleft), as.double(yright), as.double(f), NAOK = TRUE, 
        PACKAGE = "stats")$xout
    list(x = xout, y = y)
}

# convenience function

evalInvDist = function(invDist, xnew)
{
  approx(invDist$x, invDist$y, xout = xnew)$y;
}


#===============================================================================================
#
# grey2red
#
#===============================================================================================

grey2red = function(n, base, gamma)
{
  red = seq(from=base^gamma, to=1, length.out = n)^(1/gamma)
  green = blue = seq(from = base^gamma, to=0, length.out = n)^(1/gamma);
  col = rgb(red, green, blue, maxColorValue = 1); 
}


# Example of grey2red:

if (FALSE)
{
  par(mfrow = c(5,1))
  par(mar = c(1,3,1,1))
  n= 100
  barplot(rep(1, n), col = grey2red(n, 0, 1))
  barplot(rep(1, n), col = grey2red(n, 1, 1))
  barplot(rep(1, n), col = grey2red(n, 0.5, 1))
  barplot(rep(1, n), col = grey2red(n, 0.5, 0.2))
  barplot(rep(1, n), col = grey2red(n, 0.5, 5.0))
}



#===============================================================================================
#
# Network plot for generating VisANT-like plots
#
#===============================================================================================

#transform to radians

toRadians = function(angle, degrees = TRUE)
{
  if (degrees) angle = angle * pi/180;
  angle;
}

toDegrees = function(angle, radians = TRUE)
{
  if (radians) angle = angle / pi * 180;
  angle;
}

normalizeAngle = function(angle, degrees = TRUE)
{
  if (degrees)
  {
     step = 360; min = 0;
  } else {
     step = 2*pi; min = 0;
  }

  for (i in 1:length(angle))
  {
    while (angle[i] < min) angle[i] = angle[i] + step;
    while (angle[i] > step + min) angle[i] = angle[i] - step;
  }
  angle;
}

scaleToRange = function(x, min, max, scale = TRUE)
{
  if (scale)
  {
    y = (x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (max - min) + min;
  } else {
    y = x;
    y[!is.na(x)] = min;
  }
  y;
}

# Angles here are assumed in degrees
# This should return the distance from center of a symbol to the center of a label
pointToLabelDist = function(labelWidth, labelHeight, pointRadius, alignAngle, labelAngle, degrees = TRUE)
{
  alignAngle = toRadians(alignAngle, degrees);
  labelAngle = toRadians(labelAngle, degrees);

  # Align and label angles cannot go against each other since that would put the label into the point
  # symbol.
  if (cos(alignAngle - labelAngle) < 0) labelAngle = labelAngle + pi;

  al.m.lab = alignAngle - labelAngle
  ratio = labelHeight/labelWidth;
  epsilon = atan(ratio) * sign( sin(al.m.lab) );
  gamma = al.m.lab - epsilon; 

  if (abs(tan(al.m.lab)) > ratio + 2*pointRadius/labelWidth)
  {
     # "Scenario A": long (width) side touches the circle
     return ( (pointRadius + 0.5 * labelHeight)/abs(sin(al.m.lab)) );
  } else if (abs(tan(al.m.lab)) > ratio) {
     # "Scenario B": corner touches the circle
     a = 0.5 * sqrt(labelWidth^2 + labelHeight^2);
     delta.m.al = asin( a*sin(gamma) / pointRadius );
     kappa = pi - gamma - delta.m.al;
     return ( pointRadius * sin(kappa) / sin(gamma) );
  } else {
     # Scenario C: short (height) side touches
     return ( (pointRadius + 0.5 * labelWidth) / abs(cos(gamma)) );
  }
}

labelPosition = function(labelWidth, labelHeight, pointRadius, alignAngle, labelAngle, degrees = TRUE)
{
  dst = pointToLabelDist(labelWidth, labelHeight, pointRadius, alignAngle, labelAngle, degrees);
  alignAngle = toRadians(alignAngle, degrees);
  c( dst * cos(alignAngle), dst * sin(alignAngle) );
}

drawRectangle = function(center.x, center.y, width, height, angle, degrees = TRUE, ...)
{
  corners = matrix( c(-width/2, -height/2, 
                      width/2, -height/2,
                      width/2, height/2,
                      -width/2, height/2,
                      -width/2, -height/2), 2, 5)
  angle = toRadians(angle, degrees);
  rotMat = matrix(c( cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2);

  points = rotMat %*% corners + matrix(c(center.x, center.y), 2, 5);
  lines(points[1, ], points[2, ], ...);
}

drawEllipse = function(center.x, center.y, width, height = width, angle = 0, startAngle = 0, 
                       stopAngle = 360, degrees = TRUE, nSteps = 100,...)
{
  pixelAngles = toRadians(seq(from=startAngle, to = stopAngle, length.out = nSteps + 1), degrees)
  pixels = matrix( c(width * cos(pixelAngles), 
                     height * sin(pixelAngles)), 2, nSteps + 1, byrow = TRUE);

  angle = toRadians(angle, degrees);
  rotMat = matrix(c( cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2)
  points = rotMat %*% pixels +
            matrix(c(center.x, center.y), 2, nSteps + 1);
  lines(points[1, ], points[2, ], ...);
}



  

# This function needs x and y for the points.

networkPlot = function(
  adjacency,
  labels,
  pos.x, pos.y, 
  alignAngles,
  labelAngles,
  cex.labels,
  cex.points,
  lineAngles = NULL, # Not used for now
  maxAdj = NULL,
  colors = NULL,
  startNewPlot = TRUE,
  plotBox = c(-1, 1, -1, 1),
  variable.line.width = TRUE,
  min.line.width = 1,
  max.line.width = 5,
  nPlotLines = NULL,
  pch = 21,
  labelColors = "black",
  pointColors = "black",
  pointBg = "black",
  xLabelOffset = 0.01,
  yLabelOffset = 0.01,
  degrees = TRUE,
  ...)
{

  if (startNewPlot)
    plot(plotBox[1:2], plotBox[3:4], axes = FALSE, type = "n", xlab = "", ylab = "", ...);

  # plot(c(-1-xMargin,1+xMargin), c(-1-yMargin,1+yMargin), axes = FALSE, type = "n", xlab = "", ylab = "", ...) 
  checkAdjMat(adjacency, min = -1)
  n = length(labels);

  adjx = adjacency
  adjx[is.na(adjx)] = 0;

  if (length(pch)==1) pch = rep(pch, n);
  if (length(labelColors)==1) labelColors = rep(labelColors, n);
  if (length(pointColors)==1) pointColors = rep(pointColors, n);
  if (length(pointBg)==1) pointBg = rep(pointBg, n);
  if (length(xLabelOffset)==1) xLabelOffset = rep(xLabelOffset, n);
  if (length(yLabelOffset)==1) yLabelOffset = rep(yLabelOffset, n);
  if (length(cex.labels)==1) cex.labels = rep(cex.labels, n);
  if (length(cex.points)==1) cex.points = rep(cex.points, n);
  if (length(alignAngles)==1) alignAngles = rep(alignAngles, n);
  if (length(labelAngles)==1) labelAngles = rep(labelAngles, n);

  diag(adjx) = 0;
  diag(adjacency) = 0;
  maxA = max(abs(adjx));
  if (!is.null(maxAdj)) if (maxA<maxAdj) maxA = maxAdj;
  if (sum(adjx < 0) > 0)
  {
     if (is.null(colors)) colors = greenWhiteRed(100);
     adjCol = numbers2colors(adjacency, signed = TRUE, colors = colors, lim = c(-maxA, maxA));
  } else {
     if (is.null(colors)) colors = greenWhiteRed(100)[50:100];
     adjCol = numbers2colors(adjacency, signed = FALSE, colors = colors, lim = c(0, maxA));
  }


  ltA = adjacency;
  diag(ltA) = NA;
  ltA[upper.tri(ltA)] = NA;

  adjOrder = order(c(abs(ltA)))
  rows = row(adjacency)[adjOrder];
  cols = col(adjacency)[adjOrder];

  nLines = n*(n-1)/2;
  if (is.null(nPlotLines)) 
  {
    startLine = 1
  } else { 
    startLine = nLines - nPlotLines + 1;
  }
  for (line in startLine:nLines)
  {
    n1 = rows[line];
    n2 = cols[line];
    a = adjacency[n1, n2];
    normA = abs(a)/maxA;

    w = min.line.width;
    if (variable.line.width)
      w = min.line.width + (max.line.width - min.line.width) * normA;

    #pRadius1 = par("cxy") * cex.points[n1]/35;  # Emprical fudge factor..
    #pRadius2 = par("cxy") * cex.points[n2]/35;
    lineLen = sqrt( (pos.x[n1] - pos.x[n2])^2 + (pos.y[n1] - pos.y[n2])^2);
    x1 = pos.x[n1] #+ pRadius1[1] * (x[n2] - x[n1]) / lineLen
    y1 = pos.y[n1] #+ pRadius1[1] * (y[n2] - y[n1]) / lineLen
    x2 = pos.x[n2] #+ pRadius2[1] * (x[n1] - x[n2]) / lineLen
    y2 = pos.y[n2] #+ pRadius2[1] * (y[n1] - y[n2]) / lineLen

    lines(c(x1,x2),c(y1, y2), lwd = w, col = adjCol[n1, n2]);
  }

  x = pos.x;
  y = pos.y;

  for (node in 1:n)
    points(x[node], y[node], pch = pch[node], cex = cex.points[node], 
           bg = pointBg[node], col = pointColors[node]);

  for (node in 1:n)
  {
    cex = cex.labels[node];
    textWidth = strwidth(labels[node], cex = cex);
    textHeight = strheight(labels[node], cex = cex);
    pRadius = par("cxy") * cex.points[node]/5  ;  # Emprical fudge factor..
    effPointRadius = sqrt(mean(pRadius^2));
    labelShift = labelPosition(textWidth+2*xLabelOffset[node], 
                               textHeight + 2*yLabelOffset[node], effPointRadius, 
                               alignAngles[node], labelAngles[node], degrees);
    #printFlush(paste("Node:", node, " labelShift: ", paste(signif(labelShift, 2), collapse = ", "),
    #               ", width, height: ", paste(signif(c(textWidth, textHeight), 2), collapse = ",")));
    #  Make sure the label will read from left to right
    labAng = toDegrees(normalizeAngle(labelAngles[node], degrees), !degrees);

    if (labAng > 90 & labAng < 270) labAng = labAng - 180;

    text(x[node] + labelShift[1], y[node] + labelShift[2],
         labels = labels[node], adj = c(0.5, 0.5), 
         cex = cex, col = labelColors[node], srt = labAng, xpd = TRUE);
  }

}


# Debug
if (FALSE)
{
labelWidth = textWidth+2*xLabelOffset[node];
labelHeight = textHeight + 2*yLabelOffset[node];
pointRadius = effPointRadius;
alignAngle= alignAngles[node];
labelAngle = labelAngles[node];
}

#===============================================================================================
#
# Circle plot for generating VisANT-like plots
#
#===============================================================================================

circlePlot = function(
  adjacency,
  labels,
  order = NULL,
  maxAdj = NULL,
  colors = NULL,
  startNewPlot = TRUE,
  plotBox = c(-1, 1, -1, 1),
  center = c(0,0), 
  radii = c(0.8, 0.8),
  startAngle = 0,
  variable.cex.labels = TRUE,
  min.cex.labels = 1,
  max.cex.labels = 1.5,
  variable.cex.points = TRUE,
  min.cex.points = 1,
  max.cex.points = 3,
  variable.line.width = TRUE,
  min.line.width = 1,
  max.line.width = 5,
  nPlotLines = NULL,
  pch = 21,
  labelColors = "black",
  pointColors = "black",
  pointBg = "black",
# [x,y]Margin arguments have been taken out. Use radii instead.
  xLabelOffset = 0.01,
  yLabelOffset = 0.01,
  variableLabelAngle = TRUE,
  labelAngleMultiplier = 0.5,
  
  ...)
{

  checkAdjMat(adjacency, min = -1)
  n = length(labels);
  angles = seq(from = startAngle, to = startAngle + 2*pi * (1-1/n), length.out = n);
  x = center[1] + radii[1] * sin(angles);  # This is intentional; top should correspond to angle=0
  y = center[2] + radii[2] * cos(angles);

  adjx = adjacency
  adjx[is.na(adjx)] = 0;
  connectivity = colSums(abs(adjx))-diag(adjx)
  minConn = min(connectivity, na.rm = TRUE);
  maxConn = max(connectivity, na.rm = TRUE);

  if (is.null(order)) order = order(-connectivity, na.last = TRUE)

  if (length(pch)==1) pch = rep(pch, n);
  if (length(labelColors)==1) labelColors = rep(labelColors, n);
  if (length(pointColors)==1) pointColors = rep(pointColors, n);
  if (length(pointBg)==1) pointBg = rep(pointBg, n);
  if (length(xLabelOffset)==1) xLabelOffset = rep(xLabelOffset, n);
  if (length(yLabelOffset)==1) yLabelOffset = rep(yLabelOffset, n);

  oLabs = labels[order]
  oLColors = labelColors[order];
  oPColors = pointColors[order];
  oPBg = pointBg[order];
  oConn = connectivity[order];
  oAdj = adjx[order, order];
  oPch = pch[order];

  alignAngles = normalizeAngle(toDegrees(-angles) + 90); 
  alignAngles[alignAngles > 270] = alignAngles[alignAngles > 270] - 360;
  if (variableLabelAngle)
  {
     labelAngles = ifelse(alignAngles<=90, alignAngles*labelAngleMultiplier, 
                                           180 + (alignAngles - 180)*labelAngleMultiplier)
     alignAngles = labelAngles;
  } else
     labelAngles = 0;

  actualCexPts = scaleToRange(oConn, min.cex.points, max.cex.points, variable.cex.points)
  #actualCexPts = rep(min.cex.points, n);
  #if (variable.cex.points)
  #     actualCexPts = min.cex.points + (max.cex.points - min.cex.points) * 
  #                                (oConn - minConn)/(maxConn - minConn)

  cex.labels = scaleToRange(oConn, min.cex.labels, max.cex.labels, variable.cex.labels);
  #cex.labels = rep(min.cex.labels, n);
  #if (variable.cex.labels)
  #     cex.labels[node] = min.cex.labels + (max.cex.labels - min.cex.labels) *
  #                               (oConn - minConn)/(maxConn - minConn)

  networkPlot(adjacency[order, order],
              oLabs, 
              x, y,
              alignAngles,
              labelAngles,
              cex.labels = cex.labels,
              cex.points = actualCexPts,
              maxAdj = maxAdj,
              colors = colors,
              startNewPlot = startNewPlot,
              variable.line.width = variable.line.width,
              min.line.width = min.line.width,
              max.line.width = max.line.width,
              nPlotLines = nPlotLines,
              pch = oPch,
              labelColors = oLColors,
              pointColors = oPColors,
              pointBg = oPBg,
              xLabelOffset = xLabelOffset[order],
              yLabelOffset = yLabelOffset[order],
              degrees = TRUE, ...);
            
}

# Example of circle plot:

if (FALSE)
{

   sizeGrWindow(8,8)
   par(mfrow = c(1,1));
   nS = 100;
   nn = 30;
   mod = simulateModule(rnorm(nS), nn);
   adjacency = cor(mod)^3;
   
   order = NULL
   
   labels = paste("Gene", c(1:nn));
   circlePlot(adjacency, labels, order, variable.cex.labels = FALSE, radii = c(0.6, 0.6));

   # Plot two circles in one plot

   circlePlot(adjacency, labels, order, variable.cex.labels = FALSE, center = c(-0.5, -0.5), 
              radii = c(0.35, 0.35));

   circlePlot(adjacency, labels, order, startNewPlot = FALSE, 
              variable.cex.labels = FALSE, center = c(0.5, 0.5), 
              radii = c(0.35, 0.35));

   
}


#=======================================================================================================
#
# Hierarchical circle network plot
#
#=======================================================================================================

alternatingAngles = function(n, startAngle, degrees = TRUE)
{
  nm = n-1;
  steps = c(0, rep(c(1:ceil(nm/2)), rep(2, ceil(nm/2))) * rep(c(1, -1), ceil(nm/2)))
  steps = steps[1:n];

  startAngle = toDegrees(startAngle, !degrees);
  angles = normalizeAngle(startAngle + 360/n * steps);
  if (!degrees) angles = toRadians(angles);
  angles 
}


networkPlot.hierarchical = function(
  adjacency,
  labels,
  clusteringAdjacency = adjacency,
  maxAdj = NULL,
  colors = NULL,
  startNewPlot = TRUE,
  plotBox = c(-1, 1, -1, 1),
  center = c(0,0),
  outerRadii = c(0.5, 0.5),

  variableInnerRadius = TRUE,
  minInnerRadii = c(0.1, 0.1),
  maxInnerRadii = c(0.25, 0.25),

  minClusterSize = 3,
  deepSplit = 1,

  angleOfLargestCluster = -90,
  angularSizeOffset = 3,
  variable.cex.labels = TRUE,
  min.cex.labels = 1,
  max.cex.labels = 1.5,
  variable.cex.points = TRUE,
  min.cex.points = 1,
  max.cex.points = 3,
  variable.line.width = TRUE,
  min.line.width = 1,
  max.line.width = 5,
  nPlotLines = NULL,
  pch = 21,
  labelColors = "black",
  pointColors = "black",
  pointBg = "black",
# [x,y]Margin arguments have been taken out. Use radii instead.
  xLabelOffset = 0.01,
  yLabelOffset = 0.01,
  variableLabelAngle = TRUE,
  variableLabelAngleMinSize = 4,
  showSkeleton = FALSE,
  ...)
{
   
  checkAdjMat(adjacency, min = -1)
  n = length(labels);

  adjx = clusteringAdjacency
  adjx[is.na(adjx)] = 0;
  tree = hclust(as.dist(1-adjx), method = "a");

  # Get cluster order around the big circle. 
  # Room for improvement: instead of this could get eigennodes, cluster them and order their dendrogram
  # using dendrogramOrderings 
  clusters.noPam = cutreeDynamic(tree, distM = 1-adjx,
                                 minClusterSize = minClusterSize, 
                                 deepSplit = deepSplit,
                                 pamStage = FALSE);
  clusterOrder = unique(clusters.noPam[tree$order]);
  clusterOrder = clusterOrder[clusterOrder!=0];


  # Here I want all points in a cluster
  clusters = cutreeDynamic(tree, distM = 1-adjx,
                           minClusterSize = minClusterSize,
                           deepSplit = deepSplit,
                           pamStage = TRUE, cutHeight = 2);

  #plotDendroAndColors(tree, cbind(clusters.noPam, clusters))
  clusterSizes = table(clusters)
  nClusters = length(clusterSizes);

  if (nClusters==1)
  {
    circlePlot(
        adjacency,
        labels,
        order = NULL,
        maxAdj = maxAdj,
        colors = colors,
        startNewPlot = startNewPlot,
        plotBox = plotBox,
        center = center,
        radii = radii,
        startAngle = startAngle,
        variable.cex.labels = variable.cex.labels,
        min.cex.labels = min.cex.labels,
        max.cex.labels = max.cex.labels,
        variable.cex.points = variable.cex.points,
        min.cex.points = min.cex.points,
        max.cex.points = max.cex.points,
        variable.line.width = variable.line.width,
        min.line.width = min.line.width,
        max.line.width = max.line.width,
        nPlotLines = nPlotLines,
        pch = pch,
        labelColors = labelColors,
        pointColors = pointColors,
        pointBg = pointBg,
        xLabelOffset = xLabelOffset,
        yLabelOffset = yLabelOffset,
        variableLabelAngle = variableLabelAngle,
        showSkeleton = FALSE,
        ...);
    return ;
  }

  effectiveSizes = clusterSizes + angularSizeOffset 
  if (variableInnerRadius & (max(clusterSizes) > min(clusterSizes))) 
  {
    ratios = effectiveSizes / max(effectiveSizes);
    if (min(ratios) < max(minInnerRadii/maxInnerRadii))
       ratios = scaleToRange(ratios, max(minInnerRadii/maxInnerRadii), 1);
    innerRadii = matrix( ratios, nClusters, 2)  * 
                 matrix( maxInnerRadii, nClusters, 2, byrow = TRUE);
  } else
    innerRadii = matrix( minInnerRadii, nClusters, 2, byrow = TRUE);

  # This may not work well in certain cases but for now will do...
  meanApproxInnerRadius = sqrt(rowMeans(innerRadii^2));

  clusterSlot = match(names(clusterSizes), clusterOrder);

  angleRanges = meanApproxInnerRadius/sum(meanApproxInnerRadius) * 360;
  angleRanges.ordered = angleRanges[ clusterOrder ]
  circleAngles1 = c(0, cumsum(angleRanges.ordered)) + c(angleRanges.ordered/2, 0);
  circleAngles = circleAngles1[1:nClusters] - circleAngles1[2]/2;
  
  largestCluster = which.max(clusterSizes);
  circleAngles = circleAngles - circleAngles[ clusterSlot[largestCluster] ] + angleOfLargestCluster;

  if (showSkeleton) 
  {
     addGrid(v = TRUE)
     abline(h=center[2]);
     abline(v=center[1])
  }

  n = ncol(adjacency);
  x = rep(NA, n);
  y = rep(NA, n);
  kIn = rep(0, n);
  alignAngles = rep(0, n);

  if (showSkeleton) drawEllipse(center[1], center[2], outerRadii[1], outerRadii[2]);
  clusterLevels = names(clusterSizes);
  for (c in 1:nClusters)
  {
    inCluster = clusters==clusterLevels[c];
    kIn[inCluster] = colSums(abs(adjacency[inCluster, inCluster]));
    clusterAngle = circleAngles[ clusterSlot[c] ];
    incOrder = rank(-kIn[inCluster]);
    incAngles = alternatingAngles( sum(inCluster), 180 +clusterAngle)[incOrder];
    clustCenter = center + outerRadii* c(cos(toRadians(clusterAngle)), sin(toRadians(clusterAngle)));
    x[inCluster] = clustCenter[1] + innerRadii[c, 1] * cos(toRadians(incAngles));  
    y[inCluster] = clustCenter[2] + innerRadii[c, 2] * sin(toRadians(incAngles));
    alignAngles[inCluster] = normalizeAngle(incAngles);
    if (showSkeleton) 
      {
        drawEllipse(clustCenter[1], clustCenter[2], innerRadii[c, 1], innerRadii[c, 2]);
        points(clustCenter[1], clustCenter[2], pch = 5)
      }
  }

  cex.points = scaleToRange(kIn, min.cex.points, max.cex.points, variable.cex.points)
  cex.labels = scaleToRange(kIn, min.cex.labels, max.cex.labels, variable.cex.labels);

  alignAngles[alignAngles > 270] = alignAngles[alignAngles > 270] - 360;
  labelAngles = rep(0, n);
  for (c in 1:nClusters)
  {
    inCluster = clusters==clusterLevels[c];
    if (variableLabelAngle & sum(inCluster) >= variableLabelAngleMinSize)
    {
       alignAngles[inCluster] = ifelse(alignAngles[inCluster]<=90, alignAngles[inCluster]/2, 
                                          180 + (alignAngles[inCluster] - 180)/2)
       labelAngles[inCluster] = alignAngles[inCluster];
    } 
  }

  networkPlot(
        adjacency,
        labels,
        pos.x = x, pos.y = y,
        alignAngles = alignAngles,
        labelAngles = labelAngles,
        cex.labels = cex.labels,
        cex.points = cex.points,
        lineAngles = NULL, # Not used for now
        maxAdj = maxAdj,
        colors = colors,
        startNewPlot = startNewPlot & !showSkeleton,
        plotBox = plotBox,
        variable.line.width = variable.line.width,
        min.line.width = min.line.width,
        max.line.width = max.line.width,
        nPlotLines = nPlotLines,
        pch = pch,
        labelColors = labelColors,
        pointColors = pointColors,
        pointBg = pointBg,
        xLabelOffset = xLabelOffset,
        yLabelOffset = yLabelOffset,
        degrees = TRUE,
        ...);
  invisible(list(x = x, y = y, clusters = clusters, tree =tree));
}

#=======================================================================================================
#
# Interval overlap
#
#=======================================================================================================

# Interval overlap: 
#    return -1 if [xl, xh] does not intersect with [yl, yh]
#    return 1 if there's partial overlap between the intervals
#    return 2 if there's partial overlap and the ymid is in the [xl, xh] interval
#    return 3 if [yl, yh] is a subset of [xl, xh]
#    return 4 if [xl, xh] is a subset of [yl, yh]
#    return 5 if [xl, xh] is a subset of [yl, yh] and ymid is in [xl, xh]

intervalOverlap = function(xl, xh, yl, yh, ymid = NULL)
{
  n = length(xl);
  m = length(yl);
  result = matrix(-1, m, n);
  xl = as.double(xl);
  xh = as.double(xh);
  yl = as.double(yl);
  yh = as.double(yh);

  a = outer(yl, xl, `-`) * outer(yl, xh, `-`) <= 0; 
  if (!is.null(ymid))
  {
    ymid = as.double(ymid);
    b = outer(ymid, xl, `-`) * outer(ymid, xh, `-`) <= 0;
  } else 
    b = array(FALSE, dim = dim(a));

  c = outer(yh, xl, `-`) * outer(yh, xh, `-`) <= 0;
  d = outer(yh, xl, `-`) * outer(yl, xl, `-`) <= 0;
  e = outer(yh, xh, `-`) * outer(yl, xh, `-`) <= 0;

  result[xor(a,c)] = 1 + as.numeric(b[xor(a,c)]);
  result[a&c] = 3;
  result[d&e] = 4 + as.numeric(b[d&e]);
  t(result);
}

#========================================================================================================
#
# allocateJobs
#
#========================================================================================================

# Facilitates multi-threading by producing an even allocation of jobs 

allocateJobs = function(nJobs, nThreads)
{
  if (nJobs<nThreads)
    stop("Fewer jobs than threads... please decrease nThreads."); 
  n1 = floor(nJobs/nThreads);
  n2 = nJobs - nThreads*n1;
  allocation = list();
  start = 1;
  for (t in 1:nThreads)
  {
    end = start + n1 - 1 + as.numeric(t<=n2);
    allocation[[t]] = c(start:end);
    start = end+1;
  }

  allocation;
}

#=======================================================================================================
#
# eps - produce eps plots
#
#=======================================================================================================

eps = function(...) { postscript(..., horizontal = FALSE, onefile = FALSE, paper = "special"); }

#========================================================================================================
#
# Array of vertical barplots
#
#=======================================================================================================
    
barplotArray = function(mat, 
                        lim = NULL,
                        commonScale = FALSE, 
                        colors = "grey40",
                        border = "black",
                        barGap = 0,
                        ablines = NULL,
                        abColors = NULL,
                        ablines.lty = NULL,
                        textMat = NULL,
                        cex.text = 0.8,
                        leftThreshold = 0.8,
                        alignRight = TRUE,
                        col.text = 1,
                        textMargin = 0.04/ncol(mat),
                        margin = 0.1, 
                        textMat.minXPosition = 0,
                        main = "",
                        xlab = "", ylab = "",
                        cex.main = 1.4
                        
                        )
{
  nCol = ncol(mat);
  nRow = nrow(mat);

  if (length(colors)==1 | length(colors)==nRow)
  {
    colors = matrix(colors, nRow, nCol)
  } else if (length(colors)==nCol) {
    colors = matrix(colors, nRow, nCol, byrow = TRUE)
  } else if (!isTRUE(all.equal(dim(colors), c(nRow, nCol))))
    stop(spaste("'colors' must be either a single color, a vector of length nrow or ncol(mat),", 
                "or a matrix of the same dimensions as mat."));

  if (length(border)==1 | length(border)==nRow)
  {
    border = matrix(border, nRow, nCol)
  } else if (length(border)==nCol) {
    border = matrix(border, nRow, nCol, byrow = TRUE)
  } else if (!isTRUE(all.equal(dim(border), c(nRow, nCol))))
    stop(spaste("'border' must be either a single color, a vector of length nrow or ncol(mat),", 
                "or a matrix of the same dimensions as mat."));

  plot(c(0,1), c(0,1), xlim = c(0,1), ylim = c(0,1), axes = FALSE, main = main, cex.main = cex.main,
       type = "n", xlab = xlab, ylab = ylab);

  box = par("usr");
  xMin = box[1];
  xMax = box[2];
  yMin = box[3];
  yMax = box[4];

  if (is.null(lim))
  {
    if (commonScale)
    {
      matMax = rep(max(mat, na.rm = TRUE), nCol);
      matMin = rep(min(mat, na.rm = TRUE), nCol);
    } else {
      matMax = apply(mat, 2, max, na.rm = TRUE);
      matMin = apply(mat, 2, min, na.rm = TRUE);
    }
  } else { 
    if (length(lim)==2) lim = matrix(lim, 2, nCol);
    if (!all.equal(dim(lim), c(2, nCol)))
      stop(paste("if 'lim' is given, it must be either a vector of length 2 or a matrix of 2 rows and",
                 "same number of columns as 'mat'"));
    matMax = lim[2, ];
    matMin = lim[1, ];
  }

  matMin[matMin > 0] = 0; 
  mat[is.na(mat)] = 0;

  x0 = matrix(c(0:(nCol-1))/nCol * (xMax - xMin) + xMin, nRow,nCol, byrow = TRUE);
  x1 = matrix(c(1:nCol)/nCol * (xMax - xMin) + xMin, nRow,nCol, byrow = TRUE);

  y0 = matrix(c((nRow-1):0)/nRow * (yMax - yMin)+ yMin, nRow,nCol);
  y1 = matrix(c(nRow:1)/nRow * (yMax - yMin) + yMin, nRow,nCol);

  y1 = y1 - (y1-y0) * barGap;

  if (!is.null(ablines))
  {
    if (is.null(abColors)) abColors = rep(1, length(ablines))
    if (is.null(ablines.lty)) ablines.lty = rep(1, length(ablines))
  }
  
  out = list(xleft = x0, xright = x1, ytop = y0, ybottom = y1,
             x0 = rep(0, nCol),
             box = par("usr"),
             xMid = (x0[1, ] + x1[1, ])/2,
             yMid = (y0[, 1] + y1[, 1])/2);

  for (c in 1:nCol)
  {
    scale = (x1[1, c] - x0[1, c])/(matMax[c] - matMin[c])/(1+margin);
    cx0 = x0[, c] - matMin[c] * scale; # Position of the zero line. Note that matMin by def. cannot positive.
    cx1 = x0[, c] + mat[, c] * scale; # Position of the end of the bar.
    rect(cx0, y1[, c], cx1, y0[, c], col = colors[, c], border = border[, c]);

    if (!is.null(ablines))
      for (al in 1:length(ablines))
        if (ablines[al] >= matMin[c] && ablines[al] <= matMax[c])
          lines(rep(x0[1, c] + ablines[al] * scale, 2), c(y0[nRow, c], y1[1, c]), col = abColors[al],
                lty = ablines.lty[al])
    if (c > 1) lines(x = c(x0[1, c], x0[1, c]), y = c(y0[1, 1], y1[nRow, 1]), col = 1);
    out$x0[c] = cx0[1];
    if (!is.null(textMat))
    {
      onLeft = (cx1-x0[, c])/(x1[, c]-x0[, c]) > leftThreshold
      ytext = (y1[, c]  + y0[, c])/2
      if (alignRight)
      {
        xtext = x1[, c] - textMargin;
        xtext[onLeft] = pmax(cx0, cx1)[onLeft]  - textMargin
        adj.x = rep(1, nRow)
      } else {
        adj.x = ifelse(onLeft, 1, 0);
        xtext = pmax(cx0, cx1) + ifelse(onLeft, -1, 1) * textMargin;
        scaledMinXPosition = x0[1, c] + textMat.minXPosition * scale
        xtext[xtext < scaledMinXPosition + textMargin] = scaledMinXPosition + textMargin
      }
      for (r in 1:nRow)
        text(xtext[r], ytext[r], textMat[r, c], col = col.text, cex = cex.text, adj = c(adj.x[r], 0.5));
    }
  }

  invisible(out);
}


#---------------------------------------------------------------------------------------------------------
# labeledBarplotArray.R
#---------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
#
# .reverseRows = function(Matrix)
#
#--------------------------------------------------------------------------
#


.reverseRows = function(Matrix)
{
  ind = seq(from=dim(Matrix)[1], to=1, by=-1);
  Matrix[ind,];
  #Matrix
}

.reverseVector = function(Vector)
{
  ind = seq(from=length(Vector), to=1, by=-1);
  Vector[ind];
  #Vector
}
  
#--------------------------------------------------------------------------
#
# labeledBarplotArray = function ( Matrix, xLabels, yLabels, ... ) { 
#
#--------------------------------------------------------------------------
# This function plots a heatmap of the specified matrix 
# and labels the x and y axes wit the given labels.
# It is assumed that the number of entries in xLabels and yLabels is consistent 
# with the dimensions in.
# If colorLabels==TRUE, the labels are not printed and instead interpreted as colors --
#  -- a simple symbol with the appropriate color is printed instead of the label.
# The x,yLabels are expected to have the form "..color" as in "MEgrey" or "PCturquoise".
# xSymbol, ySymbols are additional markers that can be placed next to color labels

#barplotArray = function(mat,
#                        lim = NULL,
#                        commonScale = FALSE,
#                        colors = "grey40",
#                        border = "black",
#                        ablines = NULL,
#                        abColors = NULL,
#                        ablines.lty = NULL,
#                        textMat = NULL,
#                        cex.text = 0.8,
#                        leftThreshold = 0.8,
#                        col.text = 1,
#                        textMargin = 0.04/ncol(mat),
#                        margin = 0.1,
#                        main = "",
#                        xlab = "", ylab = ""
#
#                        )


labeledBarplotArray = function (mat, 
                           signed = FALSE,
                           xLabels, yLabels = NULL, 
                           xSymbols = NULL, ySymbols = NULL, 
                           colorLabels = NULL, 
                           xColorLabels = FALSE, yColorLabels = FALSE,
                           checkColorsValid = TRUE,
                           xLabelsPosition = "bottom",
                           xLabelsAngle = 45,
                           xLabelsAdj = 1,
                           xColorWidth = 0.05,
                           yColorWidth = 0.05,
                           colors = if (signed) greenWhiteRed(100) else greenWhiteRed(100)[50:100],
                           invertColors = FALSE, 
                           cex.lab = NULL, 
                           cex.lab.x = cex.lab,
                           cex.lab.y = cex.lab,
                           colors.lab.x = 1,
                           colors.lab.y = 1,
                           ... ) 
{
  if (!is.null(colorLabels)) {xColorLabels = colorLabels; yColorLabels = colorLabels; }
  
  if (is.null(yLabels) & (!is.null(xLabels)) & (dim(mat)[1]==dim(mat)[2])) 
    yLabels = xLabels; 

  if (checkColorsValid)
  {
    xValidColors = !is.na(match(substring(xLabels, 3), colors()));
    yValidColors = !is.na(match(substring(yLabels, 3), colors()));
  } else {
    xValidColors = rep(TRUE, length(xLabels));
    yValidColors = rep(TRUE, length(yLabels));
  }

  if (sum(xValidColors)>0) xColorLabInd = c(1:length(xLabels))[xValidColors]
  if (sum(!xValidColors)>0) xTextLabInd = c(1:length(xLabels))[!xValidColors]

  if (sum(yValidColors)>0) yColorLabInd = c(1:length(yLabels))[yValidColors]
  if (sum(!yValidColors)>0) yTextLabInd = c(1:length(yLabels))[!yValidColors]

  xLabPos = charmatch(xLabelsPosition, c("bottom", "top"));
  if (is.na(xLabPos))
    stop("Argument 'xLabelsPosition' must be (a unique abbreviation of) 'bottom', 'top'");

  if (is.null(colors)) colors = heat.colors(30);
  if (invertColors) colors = .reverseVector(colors);

  # Call the barplotArray function

  labPos = barplotArray(mat, 
                        colors = numbers2colors(mat, colors = colors, signed = signed,
                                                lim = c(0, max(mat))), 
                        ...)

  # Draw the labels

  nxlabels = length(xLabels)
  plotbox = labPos$box;
  xmin = plotbox[1]; xmax = plotbox[2]; ymin = plotbox[3]; yrange = plotbox[4]-ymin;
  ymax = plotbox[4]; xrange = xmax - xmin;

  xspacing = labPos$xMid[2] - labPos$xMid[1];
  yspacing = abs(labPos$yMid[2] - labPos$yMid[1]);

  nylabels = length(yLabels)
  offsetx = par("cxy")[1] / 3;
  offsety = par("cxy")[2] / 3;
  # Transform fractional widths into coordinate widths
  xColW = min(xmax - xmin, ymax - ymin) * xColorWidth;
  yColW = min(xmax - xmin, ymax - ymin) * yColorWidth;
  if (sum(!xValidColors)>0)
  {
    xLabYPos = ifelse(xLabPos==1, ymin - offsety, ymax + offsety)
    if (is.null(cex.lab)) cex.lab = 1;
    text(labPos$xMid[xTextLabInd] , xLabYPos, srt = xLabelsAngle, 
          adj = xLabelsAdj, labels = xLabels[xTextLabInd], xpd = TRUE, cex = cex.lab.x, col = colors.lab.x)
  }
  if (sum(xValidColors)>0)
  {
    baseY = ifelse(xLabPos==1, ymin-offsety-xColW, ymax + offsety + xColW);
    deltaY = ifelse(xLabPos==1, xColW, -xColW);
    rect(xleft = labPos$xMid[xColorLabInd] - xspacing/2, ybottom = baseY,
         xright = labPos$xMid[xColorLabInd] + xspacing/2, ytop = baseY + deltaY,
         density = -1,  col = substring(xLabels[xColorLabInd], 3), 
         border = substring(xLabels[xColorLabInd], 3), xpd = TRUE)
    if (!is.null(xSymbols))
      text ( labPos$xMid[xColorLabInd], baseY - sign(deltaY)* offsety, xSymbols[xColorLabInd], 
             adj = xLabelsAdj, 
             xpd = TRUE, srt = xLabelsAngle, cex = cex.lab.x, col = colors.lab.x);
  }
  if (sum(!yValidColors)>0)
  {
    if (is.null(cex.lab)) cex.lab = 1;
    text(xmin - offsetx, labPos$yMid[yTextLabInd], srt = 0, 
         adj = c(1, 0.5), labels = yLabels[yTextLabInd], xpd = TRUE, cex = cex.lab.y, col = colors.lab.y )
  } 
  if (sum(yValidColors)>0)
  {
    rect(xleft = xmin- yColW - offsetx, ybottom = labPos$yMid[yColorLabInd] - yspacing/2,
         xright = xmin- offsetx, ytop = labPos$yMid[yColorLabInd] + yspacing/2, 
         density = -1,  col = substring(yLabels[yColorLabInd], 3), 
         border = substring(yLabels[yColorLabInd], 3), xpd = TRUE)
    if (!is.null(ySymbols))
      text (xmin- yColW - 2*offsetx, 
            labPos$yMid[yColorLabInd], ySymbols[yColorLabInd], 
            adj = c(1, 0.5), xpd = TRUE, cex = cex.lab.y, col = colors.lab.y);
  }

  axis(1, labels = FALSE, tick = FALSE)
  axis(2, labels = FALSE, tick = FALSE)
  axis(3, labels = FALSE, tick = FALSE)
  axis(4, labels = FALSE, tick = FALSE)
}

#==================================================================================================
#
# fixLabels
#
#=================================================================================================

# Replace (some) spaces in the given labels by newlines so that the length of each line is no more than
# maxCharPerLine


fixLabels = function(labels, maxCharPerLine = 14, split = " ", fixed = TRUE, newsplit = split)
{
  n = length(labels);
  splitX = strsplit(labels, split = split, fixed = fixed);
  newLabels= rep("", n);
  for (l in 1:n)
  {
    nl = "";
    line = "";
    for (s in 1:length(splitX[[l]]))
    {
      newLen = nchar(line) + nchar(splitX [[l]] [s]);
      if (nchar(line) < 5 | newLen < maxCharPerLine)
      {
        nl = paste(nl, splitX[[l]] [s], sep = newsplit)
        line = paste(line, splitX[[l]] [s], sep = newsplit);
      } else {
        nl = paste(nl, splitX[[l]] [s], sep = "\n");
        line = splitX[[l]] [s];
      }
    }
    newLabels[l] = nl;
  }
  substring(newLabels, 2);
}

#==================================================================================================
#
# shortenLabels
#
#=================================================================================================

# Truncate labels at the last space before given maximum length, add ... if the label is shortened.
# maxCharPerLine


shortenLabels = function(labels, maxLength = 25, minLength = 10)
{
  n = length(labels);
  split = strsplit(labels, split = " ", fixed = TRUE);
  newLabels= rep("", n);
  for (l in 1:n)
  {
    nl = "";
    s = 1; len = 0;
    while (s <= length(split[[l]]) && len <= maxLength)
    {
      newLen = len + nchar(split [[l]] [s]);
      if (len < minLength | newLen < maxLength)
      {
        nl = paste(nl, split[[l]] [s]);
        len = nchar(nl);
        s = s+1;
      } else {
        nl = spaste(nl, "...");
        len = maxLength + 1;
      }
    }
    newLabels[l] = nl;
  }
  substring(newLabels, 2);
}


#===============================================================================================
#
# stretchLimits
#
#===============================================================================================

# Calculate minimum and maximum and stretch the limits by a given factor

stretchLimits = function(x, y, rx = 0.03, ry = 0.03)
{
  minx = min(x, na.rm = TRUE);
  maxx = max(x, na.rm = TRUE);
  xlim = c(minx - rx * (maxx-minx), maxx + rx * (maxx-minx));
  miny = min(y, na.rm = TRUE);
  maxy = max(y, na.rm = TRUE);
  ylim = c(miny - ry * (maxy-miny), maxy + ry * (maxy-miny));
  list(xlim = xlim, ylim = ylim);
}


#===============================================================================================
#
# nearestSNP
#
#===============================================================================================

# Find the nearest SNP and the corresponding distance to given gene positions.

nearestSNP = function(snpChr, snpBp, geneChr, geneBp)
{
  snpHasAnno = !is.na(snpChr) & !is.na(snpBp);
  geneHasAnno = !is.na(geneChr) & !is.na(geneBp);

  snpChr[!snpHasAnno] = -131;  # Hopefully no one will ever have chromosome -131 in their data. 
  geneChr[!geneHasAnno] = -132;  # Same for -132
  nGenes = length(geneChr);
  validChr = snpChr[snpChr!=-131]
  chrNames = sort(unique(validChr));
  nChr = length(chrNames);
  nSNPs = length(snpChr);
  nearestSNP = rep(NA, nGenes);
  for (chr in 1:nChr)
  {
    print(paste("chromosome", chrNames[chr]));
    snps = snpChr==chrNames[chr]
    pos = snpBp[snps];
    pos1 = pos[-1];
    breaks = (pos[-length(pos)] + pos1)/2;

    genes = c(1:nGenes)[geneChr==chrNames[chr] ];
    genePos = geneBp[genes];
    nearest = findInterval(genePos, breaks);
    nearestSNP[genes] = c(1:nSNPs)[snps][nearest+1];
  }
  nearestSNPdist = geneBp - snpBp[nearestSNP];
  list(nearestSNP = nearestSNP, nearestSNPdist = nearestSNPdist);
}

#===============================================================================================
#
# peakRegions
#
#===============================================================================================

# This function finds peaks and peak regions (defined as areas around the peak where the function drops less
# than a certain amount). Meant for QTL analysis.


peakRegions = function(qtl, chromo, basePair, minLOD, lodDrop, minPeakDistance)
{
  nSNPs = length(qtl);
  peaks = qtlPeaks(qtl, chromo, basePair, minQTL = minLOD, window = minPeakDistance);
  peakLocs = c(1:nSNPs)[peaks];
  nPeaks = sum(peaks)

  # For each peak, find the nearest base pairs up- and down-stream where the LOD score goes down by 1.

  peakChromo = chromo[peakLocs];
  peakBp = basePair[peakLocs];

  fromBp = rep(0, nPeaks);
  toBp = rep(0, nPeaks);
  fromLoc = rep(0, nPeaks);
  toLoc = rep(0, nPeaks);
  zzzzz = 10123;
  for (p in 1:nPeaks)
  {
    loc = peakLocs[p];
    baseLod = qtl[loc];
    while (loc > 0 && chromo[loc]==peakChromo[p] && qtl[loc] > baseLod - lodDrop)
      loc = loc -1;

    if (loc < 1) loc = 1;
  
    if (chromo[loc]!=peakChromo[p]) fromBp[p] = 0 else fromBp[p] = basePair[loc];
    fromLoc[p] = loc;
  
    loc = peakLocs[p];
    while (loc <= nSNPs && chromo[loc]==peakChromo[p] && qtl[loc] > baseLod - lodDrop)
      loc = loc +1;
  
    if (loc > nSNPs) loc = nSNPs;
    toLoc[p] = loc;
    if (chromo[loc]!=peakChromo[p]) toBp[p] = 1e10 else toBp[p] = basePair[loc];
  }

  list(peakIndicator = peaks, peakIndex = peakLocs, 
       peakChromo = peakChromo, peakBp = peakBp,
       peakStartInd = fromLoc,
       peakEndInd = toLoc,
       peakStartBp = fromBp,
       peakEndBp = toBp);
}


#===============================================================================================
#
# makeCross
#
#===============================================================================================

# prepare cross data
# First column in both genotypes and phenotypes must be sample id
# Second column in phenotypes must be sex

makeCross = function(genotypes, traits, markerNames = NULL, 
                       chromo = NULL, basePair = NULL, fileNameBase = NULL,
                       outlierZ = 3, genoCodes = c(1,2,3), alleles = c("B", "C"))
{

  commonSamples = intersect(genotypes[, 1], traits[, 1]);
  if (length(commonSamples) < 10)
    stop("Something's wrong: less than 10 samples common between genotypes and traits.");
  genotypes = genotypes[match(commonSamples, genotypes[, 1]), ];
  traits = traits[match(commonSamples, traits[, 1]), ];
  
  if (is.null(fileNameBase))
  {
    delete = TRUE
    genoFile = tempfile("genotypes", getwd());
    phenFile = tempfile("phenotypes", getwd());
  } else {
    delete = FALSE
    genoFile = spaste(fileNameBase, "-genotypes.csv");
    phenFile = spaste(fileNameBase, "-phenotypes.csv");
  }

  if (is.null(markerNames))
  {
    # Determine marker names from colnames of genotypes
    markerNames = colnames(genotypes);
    if (is.null(chromo)) # Assume marker names have the form marker.chrN.BpNNNNNN
    {
      split = sapply(strsplit(markerNames[-1], split = ".", fixed = TRUE), I);
      markerNames = split[1, ];
      chromo = as.numeric(substring( split[2, ], 4));
      basePair = as.numeric(substring( split[3, ], 3));
    }
  }

  SNPinfo = rbind( c("", chromo), c("", basePair));
  colnames(SNPinfo) =c(colnames(genotypes)[1], markerNames);
  colnames(genotypes) = colnames(SNPinfo);

  genoForQTL = rbind(SNPinfo, genotypes);

  write.csv(genoForQTL, genoFile, row.names = FALSE, quote = FALSE);

  if (outlierZ > 0)
  {
     stdTraits = scale(apply(traits[, -c(1:2)], 2, as.numeric)); # Do not scale mouse ID and sex:)
     outliers = abs(stdTraits) > outlierZ
     traits2 = as.data.frame(traits);
     traits2[, -c(1:2)][outliers] = NA;
  } else 
     traits2 = traits;
  
  write.csv(traits2, phenFile, row.names = FALSE, quote = FALSE);

  cross = read.cross(format = "csvs",
                     genfile = genoFile,
                     phefile = phenFile,
                     genotypes = genoCodes, alleles = alleles);
   
  if (delete)
  {
    unlink(genoFile);
    unlink(phenFile);
  }
  list(cross = cross, markerNames = markerNames, chromo = chromo, basePair = basePair);
}


#================================================================================================
#
# add x- and y- error bars to a plot
#
#================================================================================================

addErrorBars.xy = function(x, y, sdx, sdy, barWidth = 0.01)
{
  box = par("usr");
  barWidth.x = (box[2] - box[1])*barWidth;
  barWidth.y = (box[4] - box[3])*barWidth;

  n = length(x);

  if (length(y)!=n) stop("Lengths of x and y must be the same.");
  if (length(sdx)!=n) stop("Lengths of x and sdx must be the same.");
  if (length(sdy)!=n) stop("Lengths of x and sdy must be the same.");

  for (i in 1:n) if ( is.finite(x[i]) && is.finite(y[i]))
  {
    if (is.finite(sdx[i]))
    {
      segments(x[i]-sdx[i], y[i], x[i] + sdx[i], y[i]);
      segments(x[i]-sdx[i], y[i] - barWidth.y/2, x[i]-sdx[i], y[i] + barWidth.y/2);
      segments(x[i]+sdx[i], y[i] - barWidth.y/2, x[i]+sdx[i], y[i] + barWidth.y/2);
    }
    if (is.finite(sdy[i]))
    {
      segments(x[i], y[i]-sdy[i], x[i], y[i]+sdy[i]);
      segments(x[i]-barWidth.x/2, y[i] - sdy[i], x[i]+barWidth.x/2, y[i] - sdy[i]);
      segments(x[i]-barWidth.x/2, y[i] + sdy[i], x[i]+barWidth.x/2, y[i] + sdy[i]);
    }
  }
}

#================================================================================================
#
# linearize date
#
#================================================================================================

linearizeDate = function(dates, sep = "-", subtractMin = TRUE, yearPosition = 3, monthPosition = 1)
{
  daysInMonth = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
  add = c(0, cumsum(daysInMonth));
  splitDOB = strsplit(dates, split = sep, fixed = TRUE);
  tab = as.matrix(as.data.frame(splitDOB));

  dayPosition = c(1,2,3)[ -c(yearPosition, monthPosition)];

  month = suppressWarnings(as.numeric(as.character(tab[monthPosition, ])));

  if (any(is.na(month)))
  {
    monthcode = substring(as.character(tab[monthPosition, ]), 1, 3);
    month = match(tolower(monthcode), 
             c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"));
  }

  lin = as.numeric(as.character(tab[yearPosition, ])) * 365 +
        add[month] + as.numeric(as.character(tab[dayPosition, ]));
  lin - as.numeric(subtractMin) * min(lin, na.rm = TRUE);
}


#================================================================================================
#
# verboseMultiplot: overlay several verboseScatterplots
#
#================================================================================================

verboseMultiplot = function(x, y, sample = NULL, corFnc = "cor", corOptions = "use = 'p'", 
    colors = c(1:length(x)),
    main = "", xlab = NA, ylab = NA, cex = 1, cex.axis = 1.5, 
    cex.lab = 1.5, cex.main = 1.5, abline = FALSE, abline.color = colors, 
    abline.lty = 1, corLabel = corFnc, xlim = NULL, ylim = NULL, 
    pch = rep(1, length(x)), ...) 
{
    if (is.na(xlab)) 
        xlab = as.character(match.call(expand.dots = FALSE)$x)
    if (is.na(ylab)) 
        ylab = as.character(match.call(expand.dots = FALSE)$y)
    if (mode(x)!="list") { x = list(x); y = list(y); }
    nSets = length(x);
    if (nSets!=length(y))
      stop("Number of components in 'x' and 'y' must be the same.");
    cor = corp = rep(NA, nSets);
    for (set in 1:nSets)
    {
      corExpr = parse(text = paste(corFnc, "(x[[set]], y[[set]] ", prepComma(corOptions), ")"))
      cor[set] = signif(eval(corExpr), 2)
      corp1 = signif(corPvalueStudent(cor[set], sum(is.finite(x[[set]]) & is.finite(y[[set]]))), 2)
      if (corp1 < 10^(-200)) 
          corp[set] = "<1e-200"
      else corp[set] = spaste("=", corp1);
    }
    if (!is.na(corLabel)) {
        mainX = paste(main, " ", corLabel, "=", cor, ", p", corp, 
            sep = "")
      } else 
         mainX = main

    x1 = unlist(x);
    y1 = unlist(y);
    fin = is.finite(x1) & is.finite(y1);
    if (is.null(xlim))
    {
      minx = min(x1[fin]);
      maxx = max(x1[fin]);
      xlim = c(minx, maxx);
    }
    if (is.null(ylim))
    {
      miny = min(y1[fin]);
      maxy = max(y1[fin]);
      ylim = c(miny, maxy);
    }

    for (set in 1:nSets)
    {
      if (!is.null(sample)) 
      {
        if (length(sample) == 1) {
            sample = sample(length(x[[set]]), sample)
        }
      } else
        sample = c(1:length(x[[set]]));
      if (set==1)
      {
        plot(x[[set]] [sample], y[[set]] [sample], main = mainX, xlab = xlab, 
            ylab = ylab, cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, 
            cex.main = cex.main, xlim = xlim, ylim = ylim, col = colors[set], 
            pch = pch[set], ...)
      } else 
        points(x[[set]] [sample], y[[set]] [sample], col = colors[set], pch = pch[set])
 
      if (abline) 
      {
        fit = lm(y[[set]] ~ x[[set]])
        abline(reg = fit, col = abline.color[set], lty = abline.lty)
      }
    }
}


# Test:

if (FALSE)
{
  x = list(rnorm(20), rnorm(20)+2, rnorm(10)-1);
  y = list(rnorm(20), rnorm(20)+2, rnorm(10)-1);
  y[[1]] = y[[1]] +x[[1]];
  y[[2]] = 0.5 * y[[2]] +x[[2]];
  y[[3]] = y[[3]] - x[[3]];

  verboseMultiplot(x, y, abline = TRUE)
  
}

#=========================================================================================================
#
# Plot module heatmap and eigengene barplot
#
#=========================================================================================================

# This needs finishing

moduleHeatmapAndEigengeneBarplot = function(
   expr, colors = greenWhiteRed(50, gamma = 0.3), geneOrder = NULL, sampleOrder = NULL, 
   setLayout = TRUE, layoutHeights = c(0.6, 0.4), marHeatmap = c(0,3.8,3,1), marBarplot = c(2,3.8,0,1),
   maxNShowGenes = 500, scale = TRUE, main = "" )
{
   if (setLayout) layout(matrix(c(1,2), 2, 1), heights= layoutHeights);
   nGenes = ncol(expr);
   if (nGenes > maxNShowGenes )
   {
     sample = sample(c(1:nGenes), size = maxNShowGenes, replace = FALSE);
     modExpr = expr[, sample];
     nGenes = maxNShowGenes;
   } else
     modExpr = expr;
   if (scale) modExpr = scale(modExpr);
   if (is.null(geneOrder)) geneOrder = hclust(as.dist(1-cor(modExpr, use="p")), method = "average")$order
   if (is.null(sampleOrder)) sampleOrder = hclust(dist(modExpr), method = "average")$order;
   maxZ = max(abs(modExpr), na.rm = TRUE);
   exprColors = numbers2colors(modExpr, signed = TRUE, lim = c(-maxZ, maxZ), colors = colors);
   par(mar = marHeatmap);
   eigengene = moduleEigengenes(modExpr, rep(1, nGenes))$eigengenes[, 1];
   mp = barplot(as.vector(eigengene[sampleOrder]), col = "white", border = "white", axisnames = FALSE,
                main = main, 
                axes = FALSE);
   plotbox = par("usr");
   xstep = mp[2]-mp[1]; xLeft = mp - xstep/2; xRight = mp + xstep/2;

   nrows = ncol(modExpr);
   yRange = plotbox[4]-plotbox[3];
   yBot = plotbox[3] + c(0:(nrows-1)) * yRange/nrows;
   yTop = yBot + yRange/nrows;

   for (sample in 1:nrow(modExpr))
   {
     rect(xleft = rep(xLeft[sample], nrows), xright = rep(xRight[sample], nrows),
          ybottom = yBot, ytop = yTop, col = exprColors[ sampleOrder[sample], geneOrder], 
          border = exprColors[ sampleOrder[sample], geneOrder]);
   }
   par(mar = marBarplot);
   eigengene[is.na(eigengene)] = 0;
   colors = numbers2colors(eigengene, signed = TRUE, colors = colors)
   barplot(as.vector(eigengene[sampleOrder]), col = colors[sampleOrder],
           xlab = "", ylab = "Eigengene expression");
}


#========================================================================================================
#
# Convenient prediction accuracy functions
#
#========================================================================================================

naiveAccuracyFnc = function(response)
{
  max(table(response))/sum(table(response))
}

accuraciesFnc = function(confusion)
{
  c (diag(confusion) / rowSums(confusion), sum(diag(confusion))/sum(confusion) );
}


table2.allLevels = function(x, y, levels.x = sort(unique(x)), levels.y = sort(unique(y)), setNames = FALSE)
{
  nx = length(levels.x);
  ny = length(levels.y);
  t = table(x, y);

  out = matrix(0, nx, ny);
  if (setNames)
  {
    rownames(out) = levels.x;
    colnames(out) = levels.y;
  }
  out[ match(rownames(t), levels.x), match(colnames(t), levels.y) ] = t;
  out;
}

#========================================================================================================
#
# Convenience functions for manipulating multiData structures
#
#========================================================================================================

# Note: many of these function would be simpler to use if I used some sort of class/method technique to keep
# track of the class of each object internally. For example, I could then write a generic function "subset" that
# would work consistently on lists and multiData objects. Similarly, multiData2list would simply become a
# method of as.list, and as.list would be safe to use both on lists and on multiData objects. 

multiData.subset = function(multiData, colIndex = NULL, rowIndex = NULL)
{
  size = checkSets(multiData);
  if (is.null(colIndex)) colIndex = c(1:size$nGenes);
  if (is.null(rowIndex)) rowIndex = lapply(size$nSamples, function(n) {c(1:n)})
  if (length(rowIndex)!=size$nSets) 
    stop("If given, 'rowIndex' must be a list of the same length as 'multiData'.");
  out = list();
  for (set in 1:size$nSets)
    out[[set]] = list(data = multiData[[set]]$data[rowIndex[[set]], colIndex, drop = FALSE]);
  names(out) = names(multiData);
  out;
}

multiData2list = function(multiData)
{
  lapply(multiData, `[[`, 'data');
}

list2multiData = function(data)
{
  out = list();
  for (set in 1:length(data))
    out[[set]] = list(data = data[[set]]);
  names(out) = names(data);
  out;
}

varNames = function(multiData)
{
  colnames(multiData[[1]]$data);
}

multiData.apply = function(multiData, ..., mdaSimplify = TRUE)
{
  size = checkSets(multiData);
  res = lapply(multiData, lapply, ...)
  if (mdaSimplify) return (multiData.simplify(res));

  return(res);
}


multiData.simplify = function(multiData)
{
  len = length(multiData[[1]]$data);
  dim = dim(multiData[[1]]$data);
  simplifiable = TRUE;
  nSets = length(multiData);
  for (set in 1:nSets)
  {
    if (len!=length(multiData[[set]]$data)) simplifiable = FALSE;
    if (!isTRUE(all.equal( dim, dim(multiData[[set]]$data)))) simplifiable = FALSE;
  }
  if (simplifiable)
  {
    if (is.null(dim)) {
       innerDim = len;
       innerNames = names(multiData[[1]]$data);
       if (is.null(innerNames)) innerNames = spaste("X", c(1:len));
    } else {
       innerDim = dim;
       innerNames = dimnames(multiData[[1]]$data);
       if (is.null(innerNames)) 
         innerNames = lapply(innerDim, function(x) {spaste("X", 1:x)})
       nullIN = sapply(innerNames, is.null);
       if (any(nullIN))
         innerNames[nullIN] = lapply(innerDim[nullIN], function(x) {spaste("X", 1:x)})
    }
    setNames = names(multiData);
    if (is.null(setNames)) setNames = spaste("Set_", 1:nSets);
    multiData.s = matrix(NA, prod(innerDim), nSets);
    for (set in 1:nSets)
      multiData.s[, set] = as.vector(multiData[[set]]$data);

    dim(multiData.s) = c(innerDim, nSets);
    if (!is.null(innerNames))
      dimnames(multiData.s) = c (if (is.list(innerNames)) innerNames else list(innerNames), list(setNames));
    return(multiData.s);
  }
  return(multiData);
}

isMultiData = function(x)
{
  !inherits(try(checkSets(x), silent = TRUE), 'try-error');
}

multiData.mapply = function(FUN, ..., MoreArgs = NULL, mdmaSimplify = TRUE, mda.doCollectGarbage = FALSE)
{
  dots = list(...);
  if (length(dots)==0) 
    stop("No arguments were specified. Please type ?multiData.mapply to see the help page.");
  dotLengths = sapply(dots, length);
  if (any(dotLengths!=dotLengths[1]))
    stop(spaste("All arguments to vectorize over must have the same length.\n", 
                "Scalar arguments should be put into the 'MoreArgs' argument.\n",
                "Note: lengths of '...' arguments are: ", paste(dotLengths, collapse = ", ")));
  nArgs = length(dots);
  res = list();
  isMultiSet = sapply(dots, isMultiData);

  FUN = match.fun(FUN);
  nSets = dotLengths[1];
  for (set in 1:nSets)
  {
    localArgs = list();
    for (arg in 1:nArgs)
      localArgs[[arg]] = if (isMultiSet[arg]) dots[[arg]] [[set]] $ data else dots[[arg]] [[set]];
    names(localArgs) = names(dots);
    res[[set]] = list(data = do.call(FUN, c(localArgs, MoreArgs)));
    if (mda.doCollectGarbage) collectGarbage();
  }

  names(res) = names(dots[[1]]);

  if (mdmaSimplify)
    return(multiData.simplify(res));

  return(res);
}


multiData.rbindSelf = function(multiData)
{
  size = checkSets(multiData);
  out = NULL;
  colnames = varNames(multiData);
  for (set in 1:nSets)
  {
    if (!is.null(colnames(multiData[[set]]$data)) && 
        !isTRUE(all.equal(colnames, colnames(multiData[[set]]$data))) )
          colnames(multiData[[set]]$data) = colnames;
    out = rbind(out, multiData[[set]]$data);
  }
  out;
}


#========================================================================================================
#
# Log-transformation of data
#
#========================================================================================================

logTrafo = function(data, base = 2, zeroOrNegativeAction = c("NA", "min"), minCoeff = 0.5)
{
  origDim = dim(data);
  origNames = dimnames(data);
  data = as.matrix(data)
  action = match.arg(zeroOrNegativeAction)
  missing = is.na(data);
  data[data<=0] = NA;
  if (action=="min")
  {
    data[data<=0] = NA;
    mins = apply(data,2,min, na.rm = TRUE);
    minMat = matrix(mins*minCoeff, nrow(data), ncol(data), byrow = TRUE);
    data[is.na(data)] = minMat[is.na(data)];
  }
  out = logb(data, base);
  out[missing] = NA;
  dim(out) = origDim;
  dimnames(out) = origNames;
  out;
}

#======================================================================================================
#
# Convert a multi-level factor variable to a set of binary variables
#
#======================================================================================================

multiLev2numeric = function(x, nameBase, nameSep = ".", minFraction = 0.05, exclude = c("", "NA"))
{
  t = table(x);
  minn = ceil(length(x) * minFraction);
  keepLevels = names(t)[t >= minn & ! (names(t) %in% exclude)];
  if (sum(t<minn) > 0) addOthers = TRUE else addOthers = FALSE;
  nLevels = length(keepLevels)
  nVars = nLevels;

  numer = matrix(NA, length(x), nVars + addOthers)

  for (l in 1:nLevels)
    numer[, l] = as.numeric(x==keepLevels[l]);

  if (addOthers) numer[, nVars + 1] = as.numeric(!(x%in%keepLevels));

  numer[x %in% exclude, ] = NA;

  if (addOthers) keepLevelsX = c(keepLevels, "other") else keepLevelsX = keepLevels;

  colnames(numer) = gsub("[ /]", ".", spaste(nameBase, nameSep, keepLevelsX));

  numer;
}


#=====================================================================================================
#
# removeOutliers (using inter-sample connectivity)
#
#=====================================================================================================

outliers = function(expr, adjacencyFnc = adjacency, adjacencyOptions = list(),
                          Z.min = -3, returnZ = FALSE)
{
  adjacencyFnc = match.fun(adjacency);
  adjacencyOptions$datExpr = t(expr);
  adj = do.call(adjacencyFnc, adjacencyOptions)
  Z = scale(colSums(adj)-1);
  if (returnZ)
  {
     list(outliers = Z<Z.min, Z = Z)
  } else 
     Z<Z.min
}

#=====================================================================================================
#
# goGenes
#
#=====================================================================================================

# Returns a list of gene identifiers in the given GO category.

GOgenes = function(
  termNames = NULL,
  organism = "human",
  ontologies = c("BP", "CC", "MF"),
  evidence = c("IMP", "IGI", "IPI", "ISS", "IDA", "IEA", "TAS", "NAS", "ND", "IC"),
  includeOffspring = TRUE, 
  verbose = 2, indent = 0)
{
   organisms = c("human", "mouse", "rat", "malaria", "yeast", "fly", "bovine", "worm", "canine",
                 "zebrafish", "chicken");
   allEvidence =  c("IMP", "IGI", "IPI", "ISS", "IDA", "IEA", "TAS", "NAS", "ND", "IC");
   allOntologies = c("BP", "CC", "MF");

   spaces = indentSpaces(indent);
   orgInd = pmatch(organism, organisms);
   if (is.na(orgInd))
     stop(paste("Unrecognized 'organism' given. Recognized values are ",
                paste(organisms, collapse = ", ")));

   if (length(evidence)==0)
     stop("At least one valid evidence code must be given in 'evidence'.");
   if (length(ontologies)==0)
     stop("At least one valid ontology code must be given in 'ontology'.");

   evidInd = pmatch(evidence, allEvidence);
   if (sum(is.na(evidInd))!=0)
     stop(paste("Unrecognized 'evidence' given. Recognized values are ",
                paste(allEvidence, collapse = ", ")));

   ontoInd = pmatch(ontologies, allOntologies);
   if (sum(is.na(ontoInd))!=0)
     stop(paste("Unrecognized 'ontologies' given. Recognized values are ",
                paste(allEvidence, collapse = ", ")));

   orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
   orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
   reverseMap = c(rep(".egGO2EG", 4), ".sgdGO2ORF", rep(".egGO2EG", 6))

   missingPacks = NULL;
   packageName = paste("org.", orgCodes[orgInd], orgExtensions[orgInd], ".db", sep="");
   if (!require(packageName, character.only = TRUE))
     missingPacks = c(missingPacks, packageName);

   if (!require(GO.db))
     missingPacks = c(missingPacks, "GO.db");

   if (!is.null(missingPacks))
     stop(paste("Could not load the requisite package(s)",
           paste(missingPacks, collapse = ", "), ". Please install the package(s)."))

   if (verbose > 0)
   {
     printFlush(paste(spaces, "goGenesInCategory: loading annotation data..."));
   }

   egGO = eval(parse(text = paste(packageName, ":::org.", orgCodes[orgInd], orgExtensions[orgInd],
                                  "GO", sep = "")));

   Go2eg = eval(parse(text = paste("AnnotationDbi::as.list(", packageName, ":::org.", orgCodes[orgInd],
                                           reverseMap[orgInd],")", sep = "")));
   nTerms = length(Go2eg);

   goInfo = eval(parse(text = "AnnotationDbi::as.list(GO.db:::GOTERM)"));
   if (length(goInfo) > 0)
   {
      orgGoNames = names(Go2eg);
      dbGoIDs = as.character(sapply(goInfo, GOID));
      dbGoOntologies = as.character(sapply(goInfo, Ontology));
      dbGoTerm = as.character(sapply(goInfo, Term));
   } else {
      dbGoIDs = "";
   }

   goOffSpr = list();
   if (includeOffspring)
   {
     goOffSpr[[1]] = eval(parse(text = "AnnotationDbi::as.list(GO.db:::GOBPOFFSPRING)"));
     goOffSpr[[2]] = eval(parse(text = "AnnotationDbi::as.list(GO.db:::GOCCOFFSPRING)"));
     goOffSpr[[3]] = eval(parse(text = "AnnotationDbi::as.list(GO.db:::GOMFOFFSPRING)"));
   }
   term2info = match(names(Go2eg), names(goInfo));
   orgOntologies = dbGoOntologies[term2info];
   orgIDs = dbGoIDs[term2info];
   orgTerms = dbGoTerm[term2info];

   if (is.null(termNames)) termNames = orgTerms;

   term2orgTerm = match(tolower(termNames), tolower(orgTerms));
   if (any(is.na(term2orgTerm)))
   {
      missingTermNames = term2orgTerm[is.na(term2orgTerm)];
      term2orgID = match(missingTermNames, dbGoIDs);
      if (any(is.na(term2orgID)))
         printFlush(spaste("Warning in GOGenesInCategory: the following terms were found ",
                           "in neither GO names nor GO term IDs:\n",
                           paste(missingTermNames[is.na(term2orgID)], collapse = ", ")));
      term2orgTerm[is.na(term2orgTerm)] = term2orgID;
   }

   keepInputTerms = is.finite(term2orgTerm);
   nInputTerms = length(termNames);

   termCodes = vector(mode="list", length = nInputTerms);
   for (t in 1:nInputTerms) termCodes[[t]] = NA;
   names(termCodes) = termNames;

   termIDs = rep(NA, nInputTerms);
   termIDs[keepInputTerms] = orgIDs[term2orgTerm[keepInputTerms]];

   collectGarbage();
   nExpandLength = 0;
   blockSize = 1000; # For a more efficient concatenating of offspring genes
   nAllInTerm = rep(0, nTerms);

   if (verbose > 0)
   {
      cat(paste(spaces, " ..preparing term lists..."));
      if (nInputTerms > 10) pind = initProgInd();
   }

   for (t in 1:nInputTerms) if (keepInputTerms[t] && !is.na(Go2eg[[ term2orgTerm[t] ]] [[1]]))
   {
      c = term2orgTerm[t];
      te = as.character(names(Go2eg[[c]])); # Term evidence codes
      tc = Go2eg[[c]]; # Gene codes
      if (includeOffspring)
      {
        termOffspring = NULL;
        for (ont in 1:length(goOffSpr))
        {
          term2off = match(names(Go2eg)[c], names(goOffSpr[[ont]]))
          if (!is.na(term2off))
            termOffspring = c(termOffspring, goOffSpr[[ont]][[term2off]]);
        }
        if (length(termOffspring)>0)
        {
           maxLen = blockSize;
           tex = rep("", maxLen);
           tcx = rep("", maxLen);
           ind = 1;
           len = length(te);
           tex[ ind:len ] = te;
           tcx[ ind:len ] = tc;
           ind = len + 1;
           o2go = match(termOffspring, as.character(names(Go2eg)));
           o2go = o2go[is.finite(o2go)]
           if (length(o2go)>0) for (o in 1:length(o2go)) if (!is.na(Go2eg[[o2go[o]]][[1]]))
           {
             #printFlush(paste("Have offspring for term", c, ": ", names(Go2eg)[c], 
             #           Term(goInfo[[term2info[c]]])));
             newc = Go2eg[[o2go[o]]];
             newe = names(newc);
             newl = length(newe);
             if ((len + newl) > maxLen)
             {
               nl = len + newl;
               maxLen = blockSize * ceiling( nl/blockSize);
               tex = c(tex, rep("", maxLen - length(tex)));
               tcx = c(tcx, rep("", maxLen - length(tex)));
               nExpandLength = nExpandLength + 1;
             }
             tex[ind:(len + newl)] = newe;
             tcx[ind:(len + newl)] = newc;
             ind = ind + newl;
             len = len + newl;
           }
           te = tex[1:len];
           tc = tcx[1:len];
        }
      }
      use = is.finite(match(te, evidence));
      termCodes[[t]] = unique(as.character(tc[use]));
      if (nInputTerms > 10 && verbose > 0) pind = updateProgInd(t/nInputTerms, pind);
   }

   if (nInputTerms > 10 && verbose > 0)
   {
      pind = updateProgInd(1, pind);
      printFlush("");
   }

   attr(termCodes, "GOtermID") = termIDs

   termCodes;
}


#====================================================================================================
#
# topEnrichment
#
#====================================================================================================

# Enrichment of top genes in a given list.

topEnrichment = function(ranking, referenceList, 
                         nTopGenes = NULL, geneIDs = NULL, direction = c("increasing", "decreasing"),
                         chiSquare = FALSE,
                         verbose = 1)
{
  direction = match.arg(direction);
  ranking = as.matrix(ranking);
  if (is.null(geneIDs)) geneIDs = rownames(ranking);
  if (is.null(geneIDs)) stop("geneIDs must be specified if 'ranking' does not have row names.");
  if (length(intersect(referenceList, geneIDs))==0)
     stop("There are no common genes between referenceList and geneIDs.");

  if (is.null(nTopGenes)) nTopGenes = c(10, seq(from=20, to=1000, by=20))
  sign = ifelse(direction=="increasing", 1, -1);
  ranking = apply(sign * ranking, 2, rank, na.last = TRUE);

  nRankings = ncol(ranking);
  nNumbers = length(nTopGenes);

  inReference = is.finite(match(geneIDs, referenceList));
  if (length(unique(inReference))==0) 
    stop("referenceList contains either no or all genes in 'ranking'.");

  enrichmentStat = matrix(NA, nNumbers, nRankings);

  if (verbose > 0) pind = initProgInd();
  for (r in 1:nRankings)
    for (n in 1:nNumbers)
    {
      inTop = ranking[, r] <= nTopGenes[n];
      if (length(unique(inTop)) > 1)
      {
        if (chiSquare) {
          enrichmentStat[n, r] = chiSquare.table(table(inReference, inTop));
        } else {
          enrichmentStat[n, r] = -log10(fisher.test(table(inReference, inTop), 
                                      alternative = "greater")$p.value);
        }
      }
      if (verbose > 0) pind = updateProgInd(((r-1)*nNumbers + n)/(nRankings * nNumbers), pind);

    }

  if (verbose > 0) printFlush("");

  colnames(enrichmentStat) = colnames(ranking);
  rownames(enrichmentStat) = spaste("nTopGenes.", nTopGenes);
  attr(enrichmentStat, "nTopGenes") = nTopGenes;
  enrichmentStat;
}

#====================================================================================================
#
# chiSquare.table
#
#====================================================================================================

chiSquare.table = function(tab)
{
  if (sum(tab<5) > 0)
    return(qchisq(fisher.test(tab, alt = "greater")$p.value), df = 1, lower.tail = FALSE);

  n = sum(tab);
  expected = outer(rowSums(tab), colSums(tab))/n;
  chsq = sum( (expected - tab) ^2 / expected );
  chsq;
}

#====================================================================================================
#
# topEnrichment
#
#====================================================================================================

# Enrichment of top genes in a given list.




#====================================================================================================
#
# plotEnrichments
#
#====================================================================================================

plotEnrichments = function(enrichments, counts = attr(enrichments, "nTopGenes"),
                           colors = c(1:ncol(enrichments)),
                           ylim = NULL, pch = c(1:ncol(enrichments)), 
                           lty = rep(1, ncol(enrichments)), 
                           plotLegend = TRUE,
                           leg.position = "bottomright",
                           legend = colnames(enrichments),
                           cex.legend = 1,
                           legend.ncol = 1,
                           ...)
{
  max = max(enrichments, na.rm = TRUE);
  if (is.null(ylim)) ylim = c(0, max)
  plot(counts, enrichments[, 1], col = colors[1], type = 'l', lty = lty[1], ylim = ylim, ...);
  points(counts, enrichments[, 1], col = colors[1], pch = pch[1], bg = colors[1]);
  nPlots = ncol(enrichments);
  if (nPlots > 1) for (c in 2:nPlots)
  {
    lines(counts, enrichments[, c], col = colors[c], lty = lty[c]);
    points(counts, enrichments[, c], col = colors[c], pch = pch[c], bg = colors[c]);
  }

  if (plotLegend)
    legendClean(leg.position, legend = legend, lty = lty, pch = pch, col = colors, pt.bg = colors, 
                cex = cex.legend, ncol = legend.ncol);
}

plotEnrichments.barplot = function(enrichments, counts = attr(enrichments, "nTopGenes"),
                           nTop = as.integer(length(counts)/5),
                           colors = c(1:ncol(enrichments)),
                           ylim = NULL, 
                           legend = colnames(enrichments),
                           barplot = TRUE,
                           top = c("highest", "lowest"),
                           reverseYAxis = NULL,
                           horiz = FALSE,
                           # Ornaments: braced and delimeted groups
                           # Note: ornaments won't work for horizontal barplots
                           braceEdges = NULL,
                           braceText = NULL,
                           braceText.col = 1,
                           braceText.cex = 1,
                           braceSepLine = TRUE,
                           braceSepLine.lty = 2,
                           braceSepLine.col = "darkgrey",
                           braceHeight = 0.05,
                           braceSpace = 0.25,
                           bracePower = 13,
                           braceCol = braceText.col,
                           # Simple groups
                           groupEdges = NULL,
                           groupText = NULL, 
                           groupText.col = 1,
                           groupText.cex = 1,
                           groupSpace = 0.1,
                           ...)
{

  top = match.arg(top);
  if (is.null(reverseYAxis)) reverseYAxis = top=="lowest";

  ordEnr = apply(enrichments, 2, sort, decreasing = (top=="highest"));
  topEnr = ordEnr[1:nTop, ];
  means = colMeans(topEnr);
  stdErr = apply(topEnr, 2, stdErr)
  nBars = ncol(enrichments);

  max = max(means + stdErr, na.rm = TRUE);
  if (is.null(ylim))
  {
    if (barplot) 
    {
       ylim = c(0, max)
    } else {
       ylim = range(topEnr);
       if (reverseYAxis) ylim = rev(ylim);
    }
  }

  yLimExt = 0;

  if (!is.null(braceEdges))
  {
    yLimExt = yLimExt + braceSpace;
    braceStart = 1/(1+braceSpace);
  } else 
    braceSpace = 0;

  if (!is.null(groupEdges))
  {
    yLimExt = yLimExt + groupSpace;
    groupMiddle = (1 + 0.5 * groupSpace)/(1+groupSpace+braceSpace);
    braceStart = (1+groupSpace)/(1+groupSpace + braceSpace);
  } else
    groupSpace = 0;

  ylim[2] = ylim[2] + (ylim[2]-ylim[1]) * yLimExt;
  

  x = c(topEnr);
  g = rep(c(1:nBars), rep(nTop, nBars));

  if (barplot) {
    if (horiz)
    {
      mp = labeledBarplot2(x, g, col = colors,  xlim = ylim, names = legend, horiz = horiz, ...);
    } else {
      mp = labeledBarplot2(x, g, col = colors,  ylim = ylim, names = legend, horiz = horiz, ...);
    }
    yTop = means+stdErr;
  } else {
    bxp = labeledBoxplot(as.list(as.data.frame(topEnr)), names = legend, namesColor = colors, verbose = TRUE, 
                   border = colors, ylim = ylim, ...);
    yTop = apply(bxp$stats, 1, max, na.rm = TRUE);
    mp = c(1:ncol(topEnr));
  }

  addOrnaments(mp, yTop,
    braceEdges = braceEdges,
    braceText = braceText,      
    braceText.col = braceText.col,     
    braceText.cex = braceText.cex,     
    braceSepLine = braceSepLine,   
    braceSepLine.lty = braceSepLine.lty,  
    braceSepLine.col = braceSepLine.col,
    braceHeight = braceHeight,    
    braceSpace = braceSpace,     
    bracePower = bracePower,       
    braceCol = braceCol,
    # Simple groups        
    groupEdges = groupEdges,     
    groupText = groupText,      
    groupText.col = groupText.col,     
    groupText.cex = groupText.cex,     
    groupSpace = groupSpace);

  invisible(list(topEnrichment = topEnr, midpoints = mp, ylim = ylim));
}

# Add ornaments to a barplot/enrichment barplot

addOrnaments = function(
    midPoints,
    heights,
    braceEdges = NULL,
    braceText = NULL,
    braceText.col = 1,
    braceText.cex = 1,
    braceSepLine = TRUE,
    braceSepLine.lty = 2,
    braceSepLine.col = "darkgrey",
    braceHeight = 0.05,
    braceSpace = 0.25,
    bracePower = 13,
    braceCol = braceText.col,
    # Simple groups
    groupEdges = NULL,
    groupText = NULL,
    groupText.col = 1,
    groupText.cex = 1,
    groupSpace = 0.1)
{

  if (!is.null(braceEdges))
  {
    braceStart = 1/(1+braceSpace);
  } else
    braceSpace = 0;

  if (!is.null(groupEdges))
  {
    groupMiddle = (1 + 0.5 * groupSpace)/(1+groupSpace+braceSpace);
    braceStart = (1+groupSpace)/(1+groupSpace + braceSpace);
  } else
    groupSpace = 0;

  nBars = length(midPoints);

  box = par("usr");
  ymin = box[3];
  ymax = box[4]
  yrange = ymax - ymin;
  barSpace = midPoints[2]-midPoints[1];
  if (!is.null(braceEdges))
  {
    braceText.cex = .extend(braceText.cex, nBars);
    braceText.col = .extend(braceText.col, nBars);
    nBraces = ncol(braceEdges);
    for (b in 1:nBraces)
    {
      st = braceEdges[1, b]; en = braceEdges[2, b]; 
      if (!is.null(braceCol))
        brace.horizontal(midPoints[ st ]-0.3*barSpace, midPoints[ en ] + 0.3*barSpace, 
                         ymin + yrange * braceStart, ymin + yrange * (braceStart + braceHeight), 
                         power = bracePower, col = braceCol);
      if (braceSepLine && braceEdges[2, b] < nBars)
         abline(v = (midPoints[en] + midPoints[en +1])/2, col = braceSepLine.col, lty = braceSepLine.lty);
      text((midPoints[st] + midPoints[en])/2, ymax, braceText[b], adj = c(0.5, 1), cex = braceText.cex[b],
           col = braceText.col[b]);
    }
  }

  if (!is.null(groupEdges))
  {
    groupText.col = .extend(groupText.col, nBars)
    groupText.cex = .extend(groupText.cex, nBars)
    nGroups = ncol(groupEdges);
    for (gr in 1:nGroups)
    {
       st = groupEdges[1, gr]; en = groupEdges[2, gr];
       text((midPoints[st] + midPoints[en])/2, 
            max(0, heights[st:en]) + 0.5 * yrange * groupSpace/(1+groupSpace + braceSpace), 
            groupText[gr], adj = c(0.5, 0.5),
            cex = groupText.cex[gr], col = groupText.col[gr]);
    }
  }
    
}


.extend = function(x, n)
{
  nRep = ceil(n/length(x));
  rep(x, nRep)[1:n];
}

#======================================================================================================
#
# Validation success
#
#======================================================================================================

# This form of the function returns the mean observed value and mean rank of observed values for a given
# vector of top number of genes in predicted.
# In other words: for each column of predicted: the predicted values are ranked and nTopGenes is selected;
# the mean of observed and rank(observed) for these values are calculated. 

validationSuccess = function(predicted, observed, nTopGenes, rankPredicted = TRUE,
                             direction = c("increasing", "decreasing"))
{
  direction = match.arg(direction);
  rankSign = ifelse(direction=="decreasing", -1, 1);
  nNGenes = length(nTopGenes)
  nMethods = 2
  methodNames = c("AverageObserved", "AverageObservedRank");

  predicted = as.matrix(predicted);
  nRankings = ncol(predicted);
  if (is.null(colnames(predicted)))
  {  
    rankingNames= spaste("Ranking.", c(1:nRankings));
  } else
    rankingNames = colnames(predicted);

  success = array(NA, dim = c(nNGenes, nRankings, nMethods));
  observedRank = rank(observed * rankSign, na.last = "keep");

  if (rankPredicted) {
    predictedRank = apply(predicted * rankSign, 2, rank);
  } else 
    predictedRank = predicted;
 
  for (ing in 1:nNGenes) for (r in 1:nRankings) 
  {
    ng = nTopGenes[ing];
    keep = predictedRank[, r] <=ng;
    success[ing, r, 1] = mean(observed[keep], na.rm = TRUE);
    success[ing, r, 2] = mean(observedRank[keep], na.rm = TRUE);
  }
  
  dimnames(success) = list(spaste("nTopGenes.", nTopGenes), rankingNames, methodNames);
  success = success[ , , , drop = TRUE];
  attr(success, "nTopGenes") = nTopGenes
  success;
}

# The second form of the function returns only the success measured on correlations but also returns the
# actual values that were averaged. 

validationSuccess.ext = function(predicted, observed, nTopGenes, rankPredicted = TRUE,
                                 direction = c("increasing", "decreasing"))
{
  direction = match.arg(direction);
  rankSign = ifelse(direction=="decreasing", -1, 1);
  nNGenes = length(nTopGenes)

  predicted = as.matrix(predicted);
  nRankings = ncol(predicted);
  if (is.null(colnames(predicted)))
  {
    rankingNames= spaste("Ranking.", c(1:nRankings));
  } else
    rankingNames = colnames(predicted);

  success = array(NA, dim = c(nNGenes, nRankings));

  if (rankPredicted) {
    predictedRank = apply(predicted * rankSign, 2, rank);
  } else
    predictedRank = predicted;

  topObservedValues = matrix(0,max(nTopGenes), nRankings);

  for (r in 1:nRankings)
  {
    for (ing in 1:nNGenes) 
    {
      ng = nTopGenes[ing];
      keep = predictedRank[, r] <=ng;
      success[ing, r] = mean(observed[keep], na.rm = TRUE);
    }
    ng = max(nTopGenes);
    order = order(predicted[, r] * rankSign);
    topObservedValues[, r] = observed[order[1:ng]];
  }

  dimnames(success) = list(spaste("nTopGenes.", nTopGenes), rankingNames);
  colnames(topObservedValues) = rankingNames;
  #success = success[ , , drop = TRUE];
  attr(success, "nTopGenes") = nTopGenes
  list(success = success, topObservedValues = topObservedValues, nTopGenes = nTopGenes);
}

#====================================================================================================#
#
# replace NAs by a finite value
#
#====================================================================================================

replaceNA = function(x, replaceBy=-9)
{
  x[is.na(x)] = replaceBy;
  x;
}

#=====================================================================================================
#
# minWhichMin
#
#=====================================================================================================

# This is a wrapper around my C-level function. 
minWhichMin = function(x)
{
  x = as.matrix(x);
  nc= ncol(x);
  nr = nrow(x);
  min = rep(0, nc);
  which = rep(0, nc);
  whichmin = .C("minWhichMin", as.double(x),
                as.integer(nr), as.integer(nc),
                as.double(min), as.double(which), DUP = FALSE, NAOK = TRUE);
  cbind( min = whichmin[[4]], which = whichmin[[5]] + 1);
}

#====================================================================================================
#
# brace
#
#====================================================================================================

# draw a brace in a plot area
# For now only a horizontal brace

brace.horizontal = function(xLeft, xRight, yBottom, yTop, color=1, lwd=3, fine = 100, power = 9)
{
  xMid = (xLeft + xRight)/2
  xMM1 = (3*xLeft + xRight)/4;
  xMM2 = (xLeft + 3*xRight)/4;
  yMid = (yBottom + yTop)/2;

  xx = seq(from = -1, to=1, length.out = fine);
  yy = xx^power;

  x = c( xMM1 + (xMM1-xLeft) * xx, xMM2 + (xMM2 - xMid) * xx);
  y = yMid + (yMid-yBottom) * c(yy, rev(yy));

  lines(x, y, col = color, lwd = lwd);
} 

#=====================================================================================================
#
# Pretty names for meta-analysis output
#
#=====================================================================================================

prettyNamesForMetaAnalysis = function(names)
{
  translation = matrix( c(
        "Z.equalWeights",  "Stouffer (equal wts)",
        "Z.RootDoFWeights",  "Stouffer (sqrt wts)",
        "Z.DoFWeights", "Stouffer (dof wts)",
        "pValueHighScale.equalWeights", "Scale (equal wts)",
        "pValueHighScale.RootDoFWeights", "Scale (sqrt wts)",
        "pValueHighScale.DoFWeights", "Scale (dof wts)",
        "consensus", "Consensus",
        "weightedAverage.equalWeights", "Mean (equal wts)",
        "weightedAverage.RootDoFWeights", "Mean (sqrt wts)",
        "weightedAverage.DoFWeights", "Mean (dof wts)",
        "meta.Z.equalWeights", "Stouffer (equal wts)",
        "meta.Z.RootDoFWeights", "Stouffer (sqrt wts)",
        "meta.Z.DoFWeights", "Stouffer (dof wts)",
        "pValueHighScale.equalWeights", "Scale (equal wts)",
        "pValueHighScale.RootDoFWeights", "Scale (sqrt wts)",
        "pValueHighScale.DoFWeights", "Scale (dof wts)",
        "pValueHigh.equalWeights", "Rank (equal wts)",
        "pValueHigh.RootDoFWeights", "Rank (sqrt wts)",
        "pValueHigh.DoFWeights", "Rank (dof wts)"
        ), ncol = 2, byrow = TRUE)

  names2trans = match(names, translation[, 1]);
  out = names;
  out[is.finite(names2trans)] = translation[ names2trans[is.finite(names2trans)], 2];
  out;
}

#=====================================================================================================
#
# GOLabels
#
#=====================================================================================================

GOLabels = function(enrichmentTable, moduleCol = "Mod", pValueCol = "p.Bonf", 
                    labelCol = "termName", pValueThreshold = 0.05, numericModLabels = TRUE)
{
  if (!is.numeric(moduleCol))
    moduleCol = match(moduleCol, colnames(enrichmentTable));
  if (is.na(moduleCol)) stop("Undefined module column.")
  
  if (!is.numeric(pValueCol))
    pValueCol = match(pValueCol, colnames(enrichmentTable));
  if (is.na(pValueCol)) stop("Undefined pValue column.")
  
  if (!is.numeric(labelCol))
    labelCol = match(labelCol, colnames(enrichmentTable));
  if (is.na(labelCol)) stop("Undefined label column.")
  
  order = order(enrichmentTable[, moduleCol], enrichmentTable[, pValueCol] )
  enr.ord = enrichmentTable[order, ];
  module = enr.ord[, moduleCol];

  n = length(module);
  bestIndex = c(1, c(2:n)[ module[-1]!=module[-n] ]);

  moduleLabels = module[bestIndex];
  bestP = enr.ord[bestIndex, pValueCol];
  bestLabel = enr.ord[bestIndex, labelCol];

  GOlabels0 = spaste(": ", bestLabel);
  GOlabels0[bestP > pValueThreshold] = "";
  GOlabels = spaste(moduleLabels, GOlabels0);

  out = data.frame(GOLabels = GOlabels, modules = moduleLabels, bestP = bestP, bestLabels = bestLabel);
  if (numericModLabels)
    out = out[order(as.numeric(as.character(moduleLabels))), ]

  out;
}


#=====================================================================================================
#
# prependZeros
#
#=====================================================================================================
# prepend as many zeros as necessary to fill number to a certain width. Assumes an integer input.

prependZeros = function(x, len = max(nchar(x)))
{
  lengths = nchar(x);
  if (len < max(lengths)) stop("Some entries of 'x' are too long.");
  out = as.character(x);
  n = length(x);
  for (i in 1:n) if (lengths[i] < len)
    out[i] = spaste( paste(rep("0", len-lengths[i]), collapse = ""),
                     x[i]);

  out;
}

#=====================================================================================================
#
# entrez2geneSymbol
#
#=====================================================================================================
# convert entrez gene IDs to gene symbols

entrez2geneSymbol = function(entrezCodes, organism = "human")
{
   organisms = c("human", "mouse", "rat", "malaria", "fly", "bovine", "worm", "canine",
                  "zebrafish", "chicken");
   orgInd = pmatch(organism, organisms);
   if (is.na(orgInd))
     stop(paste("Unrecognized 'organism' given. Recognized values are ",
                 paste(organisms, collapse = ", ")));

   orgCodes = c("Hs", "Mm", "Rn", "Pf", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
   orgExtensions = rep(".eg", 10)

   missingPacks = NULL;
   packageName = paste("org.", orgCodes[orgInd], orgExtensions[orgInd], ".db", sep="");
   if (!require(packageName, character.only = TRUE))
     missingPacks = c(missingPacks, packageName);

   if (!is.null(missingPacks))
     stop(paste("Could not load the requisite package(s)",
           paste(missingPacks, collapse = ", "), ". Please install the package(s)."))

   egSymbol = eval(parse(text = paste(packageName, ":::org.", orgCodes[orgInd], orgExtensions[orgInd],
                                  "SYMBOL", sep = "")));

   mapped_entrez = as.numeric(as.character(mappedkeys(egSymbol)));
   mapped_symbols = unlist(as.list(egSymbol[mappedkeys(egSymbol)]));

   entrezCodes = as.matrix(entrezCodes);
   nSets = ncol(entrezCodes);

   symbols = array("NA", dim = dim(entrezCodes));

   for (s in 1:nSets)
   {
     codes2mapped = match(as.numeric(entrezCodes), mapped_entrez);
     fin = is.finite(codes2mapped);
     symbols[fin, s] = mapped_symbols[ codes2mapped[fin] ];
   }

   symbols[, , drop = TRUE];
}

   





