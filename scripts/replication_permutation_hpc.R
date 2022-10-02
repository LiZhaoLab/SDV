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
rm(list=ls())

setwd('/ru-auth/local/home/lzhao/Data_scratch/Khodursky/sexbias_var/')
i = commandArgs(trailingOnly=TRUE)
i <- as.numeric(i)
#i <- 1
# script containing outlier removing functions
source('scripts/outlier_removal.R')
source('scripts/replication_funcs.R')
source('scripts/permutation_funcs.R')
source('scripts/gamlss_funcs.R')

# contains sample information
runtable <- read_delim('v8_phenotype_files/phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt',
                       delim = '\t', col_types = cols(SMNOTES = col_character(), SMGTC = col_character()) )
# contains individual level phenotype information
pheno <- read_delim('v8_phenotype_files/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt',
                    delim = '\t') %>% 
  dplyr::select(SUBJID, AGE, RACE, SEX)
pheno$SEX[pheno$SEX == 1] <- "male"
pheno$SEX[pheno$SEX == 2] <- "female"

biomart_info <- read.csv('biomart_info_042221.csv', stringsAsFactors = FALSE)

# load file table. This table has information about file names for covariate and eqtl files
file_table <- read.csv('file_table.csv', stringsAsFactors = FALSE)



# run gamlss on a given tissue. This is the main function that encompasses
# all of the other functions
#run_diffsexvar <- function(tissue, runtable, pheno, biomart_info, file_table){
#tissue <- "Esophagus-Mucosa"

run_analysis <- function(gamlss_NB,
                         person_run, 
                         counts, 
                         replication, 
                         rep, 
                         tissue, 
                         i) {
  Counts_edgeR <- process_counts(counts, person_run, biomart_info)
  
  #gene_i <- seq_along(Counts_edgeR$counts[,1])
  
  if (! replication) {
   sig_genes <- gamlss_NB %>%
      filter(padj.sig <= 0.05) %>%
      pull(gene_id)

  } else {
    sig_genes <- gamlss_NB %>%
      filter(p.sig <= 0.05) %>%
      pull(gene_id)
  }
    # only calculate permutations for a block of genes depending on i
    block_size <- ceiling(length(sig_genes)/50)
    #block <- seq((i-1)*block_size + 1 , i*block_size)
    #if (i == 50){
    #  block <- seq((i-1)*block_size + 1 , length(sig_genes))
    #}
    blocks <- split(sig_genes, ceiling(seq_along(sig_genes)/block_size))
    if (length(blocks) >= i){
      sig_genes <- blocks[[i]]
      Counts_sig <- Counts_edgeR
      # Maintain matrix datatype. That's why drop=FALSE.
      Counts_sig$counts <- Counts_sig$counts[rownames(Counts_sig$counts) %in% sig_genes, , drop=FALSE]
      #Counts_sig$CPM <- Counts_sig$CPM[rownames(Counts_edgeR$counts) %in% sig_genes, ]
      # run permutations. Use mclapply to speed up
      set.seed(1)
  
      perm_results <- mclapply(seq(length(sig_genes)),
                               run_permutation,
                               Counts_sig = Counts_sig,
                               num_perm = 1000,
                               mc.cores = 4)
  
      perm_results <- do.call(rbind, perm_results)
      #write.csv(perm_results, paste0('gamlss_rep/',tissue, rep,'_', i, '_permbeta_repanal_test.csv'),
      #          quote=F,row.names = F)
      # get empirical p values
      gamlss_sig <- gamlss_NB %>% filter(gene_id %in% sig_genes)
      empirical_p <- map_dbl(gamlss_sig$gene_id, get_empp, gamlss_sig, perm_results)
      write.csv(data.frame(empp = as.numeric(empirical_p), gene_id = gamlss_sig$gene_id),
                paste0('gamlss_rep/',tissue, rep,'_', i, '_empp_repanal_s19.csv'),quote=F,row.names = F)
    }
  
}


br <- extract_counts("Breast-MammaryTissue", runtable)
gamlss_NB <- read.csv('gamlss_rep/Breast-MammaryTissue_replication0.7_repanals19_test.csv')

run_analysis(gamlss_NB,
             br$rep_pr,
             br$counts,
             replication = TRUE,
             rep = "_replication0.7",
             "Breast-MammaryTissue",
             i = i)

gamlss_NB <- read.csv('gamlss_rep/Breast-MammaryTissue_initial0.7_repanals19_test.csv')

run_analysis(gamlss_NB,
             br$init_pr,
             br$counts,
             replication = FALSE,
             rep = "_initial0.7",
             "Breast-MammaryTissue",
             i = i)

