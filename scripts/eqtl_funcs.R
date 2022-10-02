# function for loading eQTL data from a given tissue
# file table is file containing list of necessary files
load_eqtls <- function(tissue, file_table){
  eqtl_file <- file_table$eqtlfile[file_table$tissue == tissue]
  eqtls <- read.delim(paste0('eqtls/', eqtl_file)) %>% 
    mutate(gene_id = gsub("\\..*", "", phenotype_id),
           tissue = tissue) 
 
  return(eqtls)
}

# function to test enrichment/depletion of eGenes among SDV genes using fisher's exact test
# an eGene is defined as a gene with an eQTL at FDR 0.1.
enrich_eGene <- function(tis, gamlss_out, file_table){
  # only consider autosomes
  eqtls <- load_eqtls(tis, file_table)
  tested <- gamlss_out %>% 
    filter(tissue == tis)
  eqtls <- merge(eqtls, tested, by = "gene_id")
  
  frac_eqtl_sdv <- sum((eqtls$meta_bias == 'SDV') &
                           (eqtls$qval <= 0.1))/sum(eqtls$meta_bias == 'SDV')
  # get faction of all tested genes which have eQTLs
  frac_eqtl_all <- sum(eqtls$qval <= 0.1)/dim(eqtls)[1]
  # ratio of observed/expected
  ratio <- frac_eqtl_sdv/frac_eqtl_all
  fet_vec <- c(sum((eqtls$meta_bias == 'SDV') & (eqtls$qval <= 0.1)),
                 sum((eqtls$meta_bias != 'SDV') & (eqtls$qval <= 0.1)),
                 sum((eqtls$meta_bias == 'SDV') & (eqtls$qval > 0.1)),
                 sum((eqtls$meta_bias != 'SDV') & (eqtls$qval > 0.1)))
    
  p <- fisher.test(matrix(fet_vec, nr = 2))$p.value
    
    
  return(data.frame(ratio = ratio,
                      p = p,
                      tissue = tis))
 }



# split gamlss results by tissue
split_gamlss <- function(tis){
  out <- gamlss_out %>% 
    filter(tissue == tis)
  
  return(out)
}


# run a single permutation and determine if random sample is depleted for eQTLs
run_1perm <- function(i, eqtls, expfrac_sdv){
  # i is permutation number
  # eqtls is eqtl results
  # expface_sdv is the expected
  # fraction of SDV genes that will have an eQTL under randomness
  nsdv <- sum(eqtls$meta_bias == 'SDV')
  eqtls$meta_bias <- permute(eqtls$meta_bias)
  obsfrac_sdv <- sum((eqtls$meta_bias == 'SDV') &
                     (eqtls$qval <= 0.1))/nsdv
  # ratio of observed to expected 
  ratio <- obsfrac_sdv/expfrac_sdv
  return(data.frame(ratio = ratio))
}

# run permutation test for a given tissue
eqtl_permute <- function(tis, gamlss_list, eqtl_list, nperm){
  eqtls <- eqtl_list[[tis]]
  
  #eqtls_sig <- eqtls %>% filter(qval<0.1)
  # only consider gamlss output that were tested for eqtls in the first place
  tested <- gamlss_list[[tis]]
  eqtls <- merge(eqtls, tested, by = "gene_id")
  
  # observed number of SDV genes in the tissue
  # calculate the expected number of eQTLs by multiplying the total fraction of genes with eQTLs 
  # by the number of SDV genes in each tissue
  expfrac_sdv <- sum(eqtls$qval <= 0.1)/dim(eqtls)[1]
  # for each tissue run a permutation test
  perm_ratios <- map_dfr(seq(nperm), 
                        run_1perm,
                        eqtls=eqtls,
                        expfrac_sdv = expfrac_sdv)
    
  return(perm_ratios)
    
 }


# permute the SDV label given a tissue
permute_SDVlab <- function(tis, eqtls_sig){
  # permute label
  eqtls_sig <- eqtls_sig %>% 
    filter(tissue == tis) %>%
    mutate(meta_bias = permute(meta_bias))
  
  return(eqtls_sig)
}
# 1 permutation for the SFS spectrum test. 
maf_1perm <- function(i, eqtls_sig){
  
  tissues <- unique(eqtls_sig$tissue)
  # permute labels in each tissue
  eqtls_sig_perm <- map_dfr(tissues, permute_SDVlab, eqtls_sig = eqtls_sig)
  # extract permuted "SDV" genes
  eqtl_sample <- eqtls_sig_perm %>% filter(meta_bias == 'SDV')
  breaks <- seq(0,0.5,by=0.05)
  # bin them
  groups<- cut(eqtl_sample$maf, 
               breaks=breaks, 
               include.lowest=FALSE, 
               right=TRUE)
  # convert to frequency
  freq_bin = summary(groups)/dim(eqtl_sample)[1]
  
  out <- data.frame(freq = freq_bin, bin = names(freq_bin))
  return(out)
  
}

# fishers test to test for the enrichment of low minor allele frequency eQTLs among SDV genes
# enrich_lowMAF <- function(tis, gamlss_out, file_table){
#   
#   eqtls <- load_eqtls(tis, file_table)
#   tested <- gamlss_out %>% 
#     filter(tissue == tis)
#   
#   eqtls <- merge(eqtls, tested, by="gene_id")
#   eqtls_sig <- eqtls %>% filter(qval <= 0.1)
#   # minimum 20 SDV genes with eQTLs to test
#   if(sum(eqtls_sig$bias != 'not SDV') >=  20){
#     
#     # get fraction of SDV genes which have eQTLs
#     lowMAF_sdv <- sum((eqtls_sig$maf <= 0.05) &
#                         (eqtls_sig$meta_bias == 'SDV'))/sum(eqtls_sig$meta_bias == 'SDV')
#     # get faction of all tested genes which have eQTLs
#     lowMAF_all <- sum((eqtls_sig$maf <= 0.05))/dim(eqtls_sig)[1]
#     ratio <- lowMAF_sdv/lowMAF_all
#     exp <- lowMAF_all*sum(eqtls_sig$meta_bias == 'SDV')
#     obs <- sum((eqtls_sig$maf <= 0.05) & (eqtls_sig$meta_bias == 'SDV'))
#     
#     
#     fet_vec <- c(sum((eqtls_sig$maf <= 0.05) & (eqtls_sig$meta_bias == 'SDV')),
#                  sum((eqtls_sig$maf <= 0.05) & (eqtls_sig$meta_bias != 'SDV')),
#                  sum((eqtls_sig$maf > 0.05) & (eqtls_sig$meta_bias == 'SDV')),
#                  sum((eqtls_sig$maf > 0.05) & (eqtls_sig$meta_bias != 'SDV')))
#     
#     p <- fisher.test(matrix(fet_vec, nr=2))$p.value
#     
#     
#     
#     return(data.frame(ratio = ratio,
#                       p = p,
#                       tissue = tis,
#                       obs = obs,
#                       exp = exp))
#   }else{
#     return(data.frame(ratio = NA,
#                       p = NA,
#                       tissue = NA,
#                       obs = NA,
#                       exp = NA))
#   }
# }


# run a single permutation and determine if random sample is depleted for eQTLs

# run_1perm_lmaf <- function(i, eqtls, sdv_num, exp_lmaf){
#   # i is permutation number
#   # eqtls is eqtl results
#   #rows <- sample(dim(eqtls)[1], sdv_num)
#   #rand_sample <- eqtls[rows,]
#   perm_sample <- eqtls %>% 
#     mutate(meta_bias = permute(meta_bias)) %>%
#     filter(meta_bias == 'SDV')
#   obs_lmaf <- sum(perm_sample$maf <= 0.05)
#   ratio <- obs_lmaf/exp_lmaf
# 
#   return(data.frame(ratio = ratio))
# }



# # run permutation test for a given tissue low mAF
# eqtl_permute_lmaf <- function(tis, gamlss_list, eqtl_list, nperm){
#   eqtls <- eqtl_list[[tis]]
#   
#   #eqtls_sig <- eqtls %>% filter(qval<0.1)
#   # only consider gamlss output that were tested for eqtls in the first place
#   tested <- gamlss_list[[tis]] 
#   eqtls <- merge(eqtls, tested, by = "gene_id")
#   eqtls_sig <- eqtls %>% filter(qval < 0.1)
#   
#   
#   if((sum(eqtls$meta_bias == 'SDV') >= 20)){
#     # get fraction of SDV genes which have eQTLs
#     sdv_num <- sum(eqtls_sig$meta_bias == 'SDV')
#     lowMAF_sdv <- sum((eqtls_sig$maf <= 0.05) &
#                         (eqtls_sig$meta_bias == 'SDV'))/sum(eqtls_sig$meta_bias == 'SDV')
#     # get faction of all tested genes which have eQTLs
#     lowMAF_all <- sum((eqtls_sig$maf <= 0.05))/dim(eqtls_sig)[1]
#     exp_lmaf <- sdv_num*lowMAF_all
#     perm_ratio <- map_dfr(seq(nperm), 
#                          run_1perm_lmaf,
#                          eqtls = eqtls_sig,
#                          sdv_num = sdv_num, 
#                          exp_lmaf = exp_lmaf)
#     
#     return(perm_ratio)
#     
#   }else{
#     return(NA)
#   }
# }
# 
# 
# 
# # count observed number of depleted low maf eQTLs
# count_dep_lmaf <- function(tis, gamlss_list, eqtl_list){
#   eqtls <- eqtl_list[[tis]] 
#   #eqtls_sig <- eqtls %>% filter(qval<0.1)
#   # only consider gamlss output that were tested for eqtls in the first place
#   tested <- gamlss_list[[tis]] 
#   eqtls <- merge(eqtls, tested, by = "gene_id")
#   eqtls_sig <- eqtls %>% filter(qval < 0.1)
#   #if(sum(eqtls_sig$meta_bias == 'SDV') >= 20){
#     lowMAF_sdv <- sum((eqtls_sig$maf <= 0.05) &
#                         (eqtls_sig$meta_bias == 'SDV'))/sum(eqtls_sig$meta_bias == 'SDV')
#     # get faction of all tested genes which have eQTLs
#     lowMAF_all <- sum((eqtls_sig$maf <= 0.05))/dim(eqtls_sig)[1]
#     
#     
#     
#     dep <- lowMAF_sdv/lowMAF_all
#     return(dep)
#   #}else{
#   #  return(NA)
#   }
# #}