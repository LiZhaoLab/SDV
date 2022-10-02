
rem_lowexp <- function(counts){
  
  Counts_edgeR <- DGEList(counts)
  # norm factors to account for differences in library suze
  Counts_edgeR <- calcNormFactors(Counts_edgeR, method="upperquartile")
  # calculate CPM
  Counts_edgeR$CPM <- cpm.DGEList(Counts_edgeR)
  # remove genes with no expression in some samples. If not done this might screw up inference
  idx1 <- !apply(Counts_edgeR$CPM, 1, function(x) any(x == 0) )
  # set median threshold to 1 cpm
  idx2 <- rowMedians(Counts_edgeR$CPM) > 1
  Counts_edgeR$counts <- Counts_edgeR$counts[(idx1 & idx2),]
  Counts_edgeR$CPM <- Counts_edgeR$CPM[(idx1 & idx2),]
  
  return(Counts_edgeR)
}
# function to process counts
process_counts <- function(counts, person_run, biomart_info){
  # this effectively removes the gene_id column and only selects individuals of interest
  counts <- counts %>% dplyr::select(person_run$SAMPID)
  # check this
  Counts_edgeR <- rem_lowexp(counts)
  
  # remove individuals whos expression falls more than 1.5 IQR
  # away from the 75th/25th percentiles along any PC that explains more than
  # 5% of the total variance
  pca_res <- do_pca(Counts_edgeR)
  pcs <- as.data.frame(pca_res$x)
  num_pcs <- get_num_pcs(pca_res, thresh = 0.05) 
  to_remove <- map(seq(num_pcs), get_outliers_pc, pcs = pcs) %>% unlist() %>% unique()
  person_run <- person_run %>% filter(!(SAMPID %in% to_remove))
  counts <- counts %>% dplyr::select(person_run$SAMPID)
  
  
  # Recalculate CPM/norm factors and all that after removing bad samples.
  # probably doesn't matter too much
  
  Counts_edgeR <- rem_lowexp(counts)
  
  ofs <- log(Counts_edgeR$samples$lib.size*Counts_edgeR$samples$norm.factors)
  Counts_edgeR$samples$offset <- ofs
  # Remove lowly expressed genes as for these genes technical variations
  # might have a substantial impact on gene noise estimation.

  # specify sample info. Should probably recode this
  Counts_edgeR$samples$sex <- person_run$SEX
  Counts_edgeR$samples$age <- person_run$AGE
  Counts_edgeR$samples$RIN <- person_run$SMRIN
  Counts_edgeR$samples$ISCH <- person_run$SMTSISCH
  Counts_edgeR$samples$PC1 <- person_run$PC1
  Counts_edgeR$samples$PC2 <- person_run$PC2
  Counts_edgeR$samples$PC3 <- person_run$PC3
  Counts_edgeR$samples$PC4 <- person_run$PC4
  Counts_edgeR$samples$PC5 <- person_run$PC5
  # this is important if we dont have all covariates
  #some samples for some reason dont have some covariates
  Counts_edgeR$samples <- na.omit(Counts_edgeR$samples)
  # remove individuals from CPM and counts matrices which are not present in samples df
  Counts_edgeR$CPM <- Counts_edgeR$CPM[,colnames(Counts_edgeR$CPM) %in% rownames(Counts_edgeR$samples)]
  Counts_edgeR$counts <- Counts_edgeR$counts[,colnames(Counts_edgeR$counts) %in% rownames(Counts_edgeR$samples)]
  #only consider genes not on Y
  good_genes <- biomart_info$gene_id[!(biomart_info$chromosome_name %in% c('Y'))]
  Counts_edgeR$CPM <- Counts_edgeR$CPM[gsub("\\..*", "", rownames(Counts_edgeR$CPM)) %in% good_genes,]
  Counts_edgeR$counts <- Counts_edgeR$counts[gsub("\\..*", "", rownames(Counts_edgeR$counts)) %in% good_genes,]
  
  return(Counts_edgeR)
}






# fit gamlss function for 1 gene. Given a row number (i) and a EdgeR counts object 
fit_gamlss <- function(i, Counts_edgeR) {
  # data used in model 
  dat <- data.frame(sex = factor(Counts_edgeR$samples$sex),
                    age = Counts_edgeR$samples$age,
                    y = Counts_edgeR$counts[i,],
                    ofs = Counts_edgeR$samples$offset,
                    RIN = Counts_edgeR$samples$RIN,
                    ISCH = Counts_edgeR$samples$ISCH,
                    PC1 = Counts_edgeR$samples$PC1,
                    PC2 = Counts_edgeR$samples$PC2,
                    PC3 = Counts_edgeR$samples$PC3,
                    PC4 = Counts_edgeR$samples$PC4,
                    PC5 = Counts_edgeR$samples$PC5)
  dat$sex <- relevel(dat$sex, ref = c("female"))
  # specify model formulas
  form_m0_mu <- formula("y ~ sex + age + RIN + ISCH + PC1 + PC2 + PC3 + PC4 + PC5 + offset(ofs)")
  form_m0_sig <- formula(" ~ sex + age + RIN + ISCH + PC1 + PC2 + PC3 + PC4 + PC5")
  form_m1_sig <- formula(" ~ age + RIN + ISCH + PC1 + PC2 + PC3 + PC4 + PC5")
  form_m2_mu <- formula("y ~  age + RIN + ISCH + PC1 + PC2 + PC3 + PC4 + PC5 + offset(ofs)")
  
  # to assist in the fitting, identify starting point for sigma by fitting
  # simple model without sigma parameterized first
  fit_mean_m0 <- gamlss(fo = form_m0_mu,
                        data = dat,
                        family = NBI(), n.cyc = 5000)
  
  s.start_m0 = as.numeric(fitted(fit_mean_m0,"sigma")[1])
  # similarly identify starting point for sigma for model 2
  fit_mean_m2 <- gamlss(fo = form_m2_mu,
                        data = dat,
                        family = NBI(), n.cyc = 5000)
  
  s.start_m2 = as.numeric(fitted(fit_mean_m2,"sigma")[1])
  
  # fit model m0: The most complex model including all terms
  m0 <- tryCatch(
    gamlss(fo = form_m0_mu, 
           sigma.fo = form_m0_sig , data=dat,
           family = NBI(), sigma.start = s.start_m0, n.cyc = 5000, trace = F),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # fit model m1: Reduced model excluding sex from the sigma parameterization 
  m1 <- tryCatch(
    gamlss(fo =form_m0_mu, 
           sigma.fo = form_m1_sig, data=dat,
           family = NBI(), sigma.start = s.start_m0, n.cyc = 5000, trace = F),
    warning= function(w) NULL, error= function(e) NULL
  )
  # 
  # fit model 2: Reduced model excluding sex from the mu parameterization. 
  # Using all terms for sigma
  m2 <- tryCatch(
    gamlss(fo = form_m2_mu, 
           sigma.fo = form_m0_sig, data=dat,
           family = NBI(), sigma.start = s.start_m2, n.cyc = 5000, trace = F),
    warning= function(w) NULL, error= function(e) NULL
  )
  # results dataframe
  res <- data.frame(
    coef.male.sig = NA,
    coef.male.mu = NA,
    p.sig = NA,
    p.mu = NA,
    cpm = NA,
    bias = NA
  )
  
  # Only consider genes where all models could be fit properly
  if(!any(sapply(list(m0, m1, m2), is.null)))
  {
    res$coef.male.sig <- m0$sigma.coefficients[[2]]
    res$coef.male.mu <- m0$mu.coefficients[[2]]
    res$bias <- 'male'
    if(m0$sigma.coefficients[[2]]<0){
      res$bias <- 'female'
    }
    
    # likelihood ratio tests using built in gamlss functions
    tryCatch({res$p.sig = LR.test(m1, m0, print = FALSE)$p.val},
             warning= function(w) NULL, error= function(e) NULL
    )
    tryCatch({res$p.mu = LR.test(m2, m0, print = FALSE)$p.val},
             warning= function(w) NULL, error= function(e) NULL
    )
    # changed from mean
    res$cpm <- median(Counts_edgeR$CPM[i,])
  }
  return(res)
}


