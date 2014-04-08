#===================================================================================================
#
# Attempt a more involved outlier identification and removal
#
#===================================================================================================

robustScale = function(x, var = TRUE)
{
  x = as.matrix(x);
  medians = colQuantileC(x, 0.5);
  if (var)
  {
    mad = apply(x, 2, mad)
    if (any(mad==0) > 0)
    {
       printFlush("Warning in robustScale: some mad = 0. Will substitute standard deviation.");
       sd = apply(x[, mad==0, drop = FALSE], 2, sd);
       mad[mad ==0 ] = sd;
    }
    rsx = (x - matrix(medians, nrow(x), ncol(x), byrow = TRUE)) / 
            matrix(mad, nrow(x), ncol(x), byrow = TRUE)
  } else 
    rsx = (x - matrix(medians, nrow(x), ncol(x), byrow = TRUE))

  rsx;
}


outlierSamples = function( expr, Z)
{
  distM = as.matrix(dist(expr));
  meanDist = robustScale(colSums(distM, na.rm = TRUE))

  return (meanDist > Z);
}

dynamicZ = function(n, dynamic = TRUE, baseZ = 2.5, baseN = 100, minZ = 1.5)
{
  if (dynamic)
  {
    Z = max(minZ, baseZ * log10(n)/log10(baseN));
  } else
    Z = baseZ;

  Z;
}

  
removeOutliersAndBatches = function (
   expr, batchLabels, nDependentZ = TRUE,
   baseZ = 2.5, baseN = 100, minZ = 1.5,
   verbose = 1, indent = 0)
{

  if (any(is.na(expr)))
  {
     printFlush("Warning in removeOutliersAndBatches: imputing missing data.");
     expr = t(impute.knn(t(expr))$data);
  }
  nGenes = ncol(expr);
  nSamples = nrow(expr);
  if (nSamples != length(batchLabels)) 
     stop("Number of samples in expr must equal length of batchLabels.");

  spaces = indentSpaces(indent);

  # Step 1: flag outliers in each sample cluster separately
  outliers = rep(FALSE, nSamples)

  batchNames = sort(unique(batchLabels));
  nBatches = length(batchNames);
  for (b in 1:nBatches)
  {
    batchSamples = c(1:nSamples) [ batchLabels==batchNames[b] ]
    batchSize = length(batchSamples);
    Z = dynamicZ(batchSize, nDependentZ, baseZ, baseN, minZ)
    outliers [ batchSamples] = outlierSamples(expr[batchSamples, ], Z);
                                 
    if (verbose > 0 && sum(outliers[batchSamples]) > 0)
     printFlush(paste(spaces, "Step 1: in batch", batchNames[b], 
                      "used Z:", signif(Z, 2), " to remove samples", 
                      paste(rownames(expr)[batchSamples] [outliers[batchSamples] ], collapse = ", ")));
  }

  # Step 2: use robust standardization to standardize all batches to the same mean. However, leave variance
  # unchanged.

  expr1 = expr[!outliers, ];
  labels1 = batchLabels[!outliers];
  nSamples1 = nrow(expr1);

  batchSizes = table(labels1);
  biggestBatch = names(batchSizes)[ which.max(batchSizes) ];

  bbMedians = colQuantileC(expr1[labels1==biggestBatch, ], 0.5);

  stdExpr = expr1;
  for (b in 1:nBatches)
  {
    batchSamples = c(1:nSamples1) [ labels1==batchNames[b] ]
    stdExpr[batchSamples, ] = robustScale( expr1[batchSamples, ], var = FALSE);
  }
  
  # Step 3: remove outliers from the combined data:

  outliers2 = outlierSamples(stdExpr, dynamicZ(batchSize, nDependentZ, baseZ, baseN, minZ));
  if (verbose > 0 && sum(outliers2) > 0)
     printFlush(paste(spaces, "Step 2: removing samples", 
                      paste(rownames(expr1)[outliers], collapse = ", ")));
  expr2 = stdExpr[!outliers2, ];

  # Step 4: restore the medians of the largest batch into all samples
  expr2 = expr2 + matrix(bbMedians, nrow(expr2), nGenes, byrow = TRUE);

  outliers2.inAll = rep(TRUE, nSamples);
  outliers2.inAll[!outliers] = outliers2;
  # That's the final result.
  list(expr = expr2, batchOutliers = outliers, allOutliers = outliers2.inAll);
}

