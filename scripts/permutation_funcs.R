# function to run gamlss during permutation. 
# In this case only need to run model m0




fit_gamlss_m0 <- function(gene_i, Counts_edgeR) {
  # data frame used for model fitting  
  dat <- data.frame(sex = factor(Counts_edgeR$samples$sex),
                    age = Counts_edgeR$samples$age,
                    y = Counts_edgeR$counts[gene_i,],
                    ofs = Counts_edgeR$samples$offset,
                    RIN = Counts_edgeR$samples$RIN,
                    ISCH = Counts_edgeR$samples$ISCH,
                    PC1 = Counts_edgeR$samples$PC1,
                    PC2 = Counts_edgeR$samples$PC2,
                    PC3 = Counts_edgeR$samples$PC3,
                    PC4 = Counts_edgeR$samples$PC4,
                    PC5 = Counts_edgeR$samples$PC5)
  dat$sex <- relevel(dat$sex, ref = c("female"))
  
  # formulas used for m0 model
  form_m0_mu <- formula("y ~ sex + age + RIN + ISCH + PC1 + PC2 + PC3 + PC4 + PC5 +  offset(ofs)")
  form_m0_sig <- formula(" ~ sex + age + RIN + ISCH + PC1 + PC2 + PC3 + PC4 + PC5 ")
  # find suitable starting sigma
  fit_mean_m0 <- gamlss(fo = form_m0_mu,
                        data = dat,
                        family = NBI(), n.cyc = 5000)
  
  s.start_m0 = as.numeric(fitted(fit_mean_m0,"sigma")[1])
  # fit model
  m0 <- tryCatch(
    gamlss(fo = form_m0_mu, 
           sigma.fo = form_m0_sig , data=dat,
           family = NBI(), sigma.start = s.start_m0, n.cyc = 5000, trace = F),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # save results
  res <- data.frame(
    coef.male.sig = NA,
    coef.male.mu = NA
  )
  
  # Fitting might fail
  if(!is.null(m0))
  {
    res$coef.male.sig <- m0$sigma.coefficients[[2]]
    res$coef.male.mu <- m0$mu.coefficients[[2]]
  }
  return(res)
}

# function to run 1 permutation and fit gamlss
run_1perm <- function(iter, gene_i, Counts_sig){
  # permute sex labels
  Counts_sig$samples$sex <- permute(Counts_sig$samples$sex)
  out <- fit_gamlss_m0(gene_i, Counts_sig)
  return(out)
}
# function to run all necessary permutations on gene 'i'
run_permutation <- function(gene_i, Counts_sig, num_perm){
  # run permutation on gene 'i'
  # to make reproducible 
  set.seed(gene_i)
  nperm_left <- num_perm
  nperm_comp <- 0
  out <- data.frame(
    coef.male.sig = c(),
    coef.male.mu = c())
  # make sure you obtain 1000 permutations that are not NA
  while (nperm_left > 0){
    perm_res <- map_dfr(seq(nperm_left), 
                        run_1perm, 
                        gene_i = gene_i, 
                        Counts_sig = Counts_sig) %>% 
      na.omit()
    
    out <- rbind(out, perm_res)
    # count number of permutations completed and remaining
    nperm_comp <- dim(out)[1]
    nperm_left <- num_perm - nperm_comp
  }
  out$gene_id <- rownames(Counts_sig$counts)[gene_i]
  return(out)
}

# get empirical p value for each gene
get_empp <- function(gene, gamlss_sig, perm_res){
  # extract observed and simulated coefficients
  obs_p <- gamlss_sig %>%
    filter(gene_id == gene) %>% 
    dplyr::select(coef.male.sig)
  perm <- perm_res %>% 
    filter(gene_id == gene) %>%
    dplyr::select(coef.male.sig)
  # get p value. Fraction of permuted coefficients more extreme
  # than the observed coefficient
  empp <- sum(abs(perm$coef.male.sig) >= abs(obs_p$coef.male.sig))/length(perm$coef.male.sig)
  return(empp)
}