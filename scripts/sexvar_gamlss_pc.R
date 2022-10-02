library(dplyr)
library(readr)
library(edgeR)
library(gamlss)
library(stringdist)
library(gtools)
#library(biomaRt)
library(parallel)
library(matrixStats)
library(tidyverse)

setwd('/ru-auth/local/home/lzhao/Data_scratch/Khodursky/sexbias_var/')

rm(list=ls())
# scripts containing necessary functions

source('scripts/outlier_removal.R')
source('scripts/gamlss_funcs.R')
source('scripts/permutation_funcs.R')

tissue = commandArgs(trailingOnly=TRUE)

# contains sample information
runtable <- read_delim('v8_phenotype_files/phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt',
                       delim = '\t' )
# exclude some samples. In particular sex-specific samples and
# samples that dont have 2 sexes
#tissues <- gsub(" ", "", unique(runtable$SMTSD))
#tissues <- tissues[!(tissues %in% c("Brain-FrontalCortex(BA9)",
#                                    "Ovary", "Uterus", "Vagina",
#                                    "Bladder", "Prostate", "Testis",
#                                    "FallopianTube","Kidney-Medulla","Cervix-Ectocervix","Cervix-Endocervix",
#                                    "Cells-Leukemiacellline(CML)"))]


pheno <- read_delim('v8_phenotype_files/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt',
                    delim = '\t') %>% 
  dplyr::select(SUBJID, AGE, RACE, SEX)
pheno$SEX[pheno$SEX == 1] <- "male"
pheno$SEX[pheno$SEX == 2] <- "female"

biomart_info <- read.csv('biomart_info_042221.csv', stringsAsFactors = F)


# run gamlss on a given tissue. This is the main function that encompasses
# all of the other functions
run_diffsexvar <- function(tissue, runtable, pheno, biomart_info){
  print(tissue)
  # load counts table
  counts <- as.data.frame(read_csv(paste0('count_tables/', tissue, '.csv')))
  colnames(counts)[1] <- "gene_id"
  rownames(counts) <- counts$gene_id
  # load covariates
  cov_files <- list.files('GTEx_Analysis_v8_eQTL_covariates/')
  cov_file <- cov_files[amatch(paste0(tissue,'v8.covariates.txt'),
                               cov_files,method='lcs',maxDist=1000)]
  print(tissue)
  print(cov_file)
  cov <- read.table(paste0('GTEx_Analysis_v8_eQTL_covariates/', cov_file),
                    header = T, check.names = F)
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
    filter(SAMPID %in% colnames(counts)) %>%
    merge(pheno, by = "SUBJID") %>%
    merge(cov, by = "SUBJID") %>% 
    filter(RACE == 3) # only leave caucasian individuals
  
 
  Counts_edgeR <- process_counts(counts, person_run, biomart_info)
  gene_i <- seq_along(Counts_edgeR$counts[,1])
  
  # To try algorithm for just some genes, change gene_i variable. For example, set gene_i <- c(1:100) 
  # to estimate GAMLSS models for the first hundred genes.
  # this block here runs the initial gamlss fit for every gene in a tissue
  gamlss_df <- mclapply(gene_i, fit_gamlss, Counts_edgeR=Counts_edgeR, mc.cores = 12)
  gamlss_df <- do.call(rbind, gamlss_df) 
  gamlss_df$gene_id <- rownames(Counts_edgeR$counts)
  # remove instances where p.sig can't be calculated
  gamlss_df <- gamlss_df %>% 
    filter(!is.na(p.sig)) %>% 
    mutate(padj.sig = p.adjust(p.sig, method = "BH"), 
                               padj.mu = p.adjust(p.mu, method = "BH"))
  write.csv(gamlss_df, paste0('gamlss_out/',tissue,'_simple_initial_PC_test.csv'),quote=F,row.names = F)
  # initially significant genes
  #sig_genes <- gamlss_df$gene_id[gamlss_df$padj.sig <= 0.05]
  sig_genes <- gamlss_df %>% filter(padj.sig <= 0.05) %>% pull(gene_id)
  Counts_sig <- Counts_edgeR
  # Maintain matrix datatype. That's why drop=FALSE.
  Counts_sig$counts <- Counts_sig$counts[rownames(Counts_sig$counts) %in% sig_genes, , drop=FALSE]
  Counts_sig$CPM <- Counts_sig$CPM[rownames(Counts_edgeR$counts) %in% sig_genes,]
  # run permutations. Use mclapply to speed up
  set.seed(1)
  
  perm_results <- mclapply(seq(length(sig_genes)), 
                           run_permutation, 
                           Counts_sig = Counts_sig, 
                           num_perm = 1000, 
                           mc.cores =12)
  
  perm_results <- do.call(rbind, perm_results)
  write.csv(perm_results, paste0('gamlss_out/',tissue,'_simple_permutation_beta_PC_test.csv'),quote=F,row.names = F)
  # get empirical p values
  gamlss_sig <- gamlss_df %>% filter(padj.sig <= 0.05)
  
  empirical_p <- map_dbl(gamlss_sig$gene_id, get_empp, gamlss_sig, perm_results)
  
  write.csv(data.frame(empp = as.numeric(empirical_p), gene_id=gamlss_sig$gene_id),
            paste0('gamlss_out/',tissue,'_empp_beta_PC_test.csv'),quote=F,row.names = F)
}


run_diffsexvar(tissue, runtable=runtable,pheno=pheno,  biomart_info = biomart_info)
