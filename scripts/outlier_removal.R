# get outliers for a single PC, specified by name (PC1) or column
# input the principal components
get_outliers_pc <- function(pc, pcs){
  iqr <- IQR(pcs[,pc])
  upper <- quantile(pcs[,pc])[['75%']] + 1.5*iqr
  lower <- quantile(pcs[,pc])[['25%']] - 1.5*iqr
  outliers <- rownames(pcs)[(pcs[,pc]>upper) |(pcs[,pc]<lower) ]
  return(outliers)
}

# identify number of PCs to use
get_num_pcs <- function(pca_res, thresh){
  s <- summary(pca_res)
  # cumvar <- cumsum(s$importance[2,])
  n <- sum(s$importance[2,]>thresh) # the second row gives the percentage of variance explained by each pc
  return(n)
}

# do pca
do_pca <- function(Counts_edgeR){
  data <- t(Counts_edgeR$CPM)
  data <- log2(data+1)
  return(prcomp(data, scale. = TRUE))
  
}
