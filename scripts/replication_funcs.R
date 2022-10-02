library(dplyr)
library(readr)
library(edgeR)
library(gamlss)
library(stringdist)
library(gtools)
#library(biomaRt)
library(parallel)
library(matrixStats)
library(purrr)





extract_counts <- function(tissue, runtable){
  # load counts table
  counts <- as.data.frame(read_csv(paste0('count_tables/', tissue, '.csv')))
  colnames(counts)[1] <- "gene_id"
  rownames(counts) <- counts$gene_id
  # load covariates
  cov_file <- file_table$covfile[file_table$tissue == tissue]
  cov <- read.table(paste0('GTEx_Analysis_v8_eQTL_covariates/', cov_file),
                    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  #cov <- cov %>% filter(ID %in% c('PC1','PC2','PC3','PC4','PC5'))
  cov <- data.frame(SUBJID = colnames(cov)[c(2:dim(cov)[2])],
                    PC1 = as.numeric(cov %>% filter(ID == 'PC1') %>% dplyr::select(c(2:dim(cov)[2]))),
                    PC2 = as.numeric(cov %>% filter(ID == 'PC2') %>% dplyr::select(c(2:dim(cov)[2]))),
                    PC3 = as.numeric(cov %>% filter(ID == 'PC3') %>% dplyr::select(c(2:dim(cov)[2]))),
                    PC4 = as.numeric(cov %>% filter(ID == 'PC4') %>% dplyr::select(c(2:dim(cov)[2]))),
                    PC5 = as.numeric(cov %>% filter(ID == 'PC5') %>% dplyr::select(c(2:dim(cov)[2]))))
  
  # create table linking person to run
  person_run <- runtable %>%
    dplyr::select(SAMPID, SMTSD, SMRIN, SMTSISCH, SMAFRZE) %>%
    mutate(SMTSD = gsub(" ","",SMTSD),
           SUBJID = sub("^([^-]+-[^-]+).*", "\\1", SAMPID))  %>%
    filter(SAMPID %in% colnames(counts))
  
  person_run <- merge(person_run, pheno, by = "SUBJID")
  person_run <- merge(person_run, cov, by = "SUBJID")
  # only leave caucasian individuals
  person_run <- person_run %>%
    filter(RACE == 3)
  set.seed(19)
  
  all_male <- person_run %>% filter(SEX == 'male')
  rm_init <- sample(dim(all_male)[1], round(0.7*dim(all_male)[1]))
  all_female <- person_run %>% filter(SEX == 'female')
  rf_init <- sample(dim(all_female)[1], round(0.7*dim(all_female)[1]))
  
  init_person_run <- rbind(all_male[rm_init,], all_female[rf_init,])
  rep_person_run <- person_run %>% filter(!SUBJID %in% init_person_run$SUBJID)
  out <- list(init_person_run, rep_person_run, counts)
  names(out) <- c("init_pr", "rep_pr", "counts")
  return(out)
}


# rem_lowexp <- function(counts){
#   
#   Counts_edgeR <- DGEList(counts)
#   # norm factors to account for differences in library suze
#   Counts_edgeR <- calcNormFactors(Counts_edgeR, method="upperquartile")
#   # calculate CPM
#   Counts_edgeR$CPM <- cpm.DGEList(Counts_edgeR)
#   # remove genes with no expression in some samples. If not done this might screw up inference
#   idx1 <- !apply(Counts_edgeR$CPM, 1, function(x) any(x == 0) )
#   # set median threshold to 1 cpm
#   idx2 <- rowMedians(Counts_edgeR$CPM) > 1
#   Counts_edgeR$counts <- Counts_edgeR$counts[(idx1 & idx2),]
#   Counts_edgeR$CPM <- Counts_edgeR$CPM[(idx1 & idx2),]
#   
#   return(Counts_edgeR)
# }
# # function to process counts
# process_counts <- function(counts, person_run, biomart_info){
#   # this effectively removes the gene_id column and only selects individuals of interest
#   counts <- counts %>% dplyr::select(person_run$SAMPID)
#   # check this
#   Counts_edgeR <- rem_lowexp(counts)
#   
#   # remove individuals whos expression falls more than 1.5 IQR
#   # away from the 75th/25th percentiles along any PC that explains more than
#   # 5% of the total variance
#   pca_res <- do_pca(Counts_edgeR)
#   pcs <- as.data.frame(pca_res$x)
#   num_pcs <- get_num_pcs(pca_res, thresh = 0.05) 
#   to_remove <- map(seq(num_pcs), get_outliers_pc, pcs = pcs) %>% unlist() %>% unique()
#   person_run <- person_run %>% filter(!(SAMPID %in% to_remove))
#   counts <- counts %>% dplyr::select(person_run$SAMPID)
#   
#   
#   # Recalculate CPM/norm factors and all that after removing bad samples.
#   # probably doesn't matter too much
#   
#   Counts_edgeR <- rem_lowexp(counts)
#   
#   ofs <- log(Counts_edgeR$samples$lib.size*Counts_edgeR$samples$norm.factors)
#   Counts_edgeR$samples$offset <- ofs
#   # Remove lowly expressed genes as for these genes technical variations
#   # might have a substantial impact on gene noise estimation.
#   
#   # specify sample info. Should probably recode this
#   Counts_edgeR$samples$sex <- person_run$SEX
#   Counts_edgeR$samples$age <- person_run$AGE
#   Counts_edgeR$samples$RIN <- person_run$SMRIN
#   Counts_edgeR$samples$ISCH <- person_run$SMTSISCH
#   Counts_edgeR$samples$PC1 <- person_run$PC1
#   Counts_edgeR$samples$PC2 <- person_run$PC2
#   Counts_edgeR$samples$PC3 <- person_run$PC3
#   Counts_edgeR$samples$PC4 <- person_run$PC4
#   Counts_edgeR$samples$PC5 <- person_run$PC5
#   # this is important if we dont have all covariates
#   #some samples for some reason dont have some covariates
#   Counts_edgeR$samples <- na.omit(Counts_edgeR$samples)
#   # remove individuals from CPM and counts matrices which are not present in samples df
#   Counts_edgeR$CPM <- Counts_edgeR$CPM[,colnames(Counts_edgeR$CPM) %in% rownames(Counts_edgeR$samples)]
#   Counts_edgeR$counts <- Counts_edgeR$counts[,colnames(Counts_edgeR$counts) %in% rownames(Counts_edgeR$samples)]
#   #only consider genes not on Y
#   good_genes <- biomart_info$gene_id[!(biomart_info$chromosome_name %in% c('Y'))]
#   Counts_edgeR$CPM <- Counts_edgeR$CPM[gsub("\\..*", "", rownames(Counts_edgeR$CPM)) %in% good_genes,]
#   Counts_edgeR$counts <- Counts_edgeR$counts[gsub("\\..*", "", rownames(Counts_edgeR$counts)) %in% good_genes,]
#   
#   return(Counts_edgeR)
# }


