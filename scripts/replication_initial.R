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

setwd('/ru-auth/local/home/lzhao/Data_scratch/Khodursky/sexbias_var/')
# setwd('~/Desktop/LabWork/sexbias_variability/')
rm(list=ls())
# script containing outlier removing functions
source('scripts/outlier_removal.R')
source('scripts/replication_funcs.R')
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



run_analysis <- function(person_run, counts, replication, rep, tissue) {
 
  Counts_edgeR <- process_counts(counts, person_run)
  
  
  gene_i <- seq_along(Counts_edgeR$counts[,1])
  
  # To try algorithm for just some genes, change gene_i variable. For example, set gene_i <- c(1:100)
  # to estimate GAMLSS models for the first hundred genes.
  # this block here runs the initial gamlss fit for every gene in a tissue
  gamlss_NB <- mclapply(gene_i, fit_gamlss, Counts_edgeR=Counts_edgeR, mc.cores = 12)
  gamlss_NB <- do.call(rbind, gamlss_NB)
  gamlss_NB <- gamlss_NB %>%
    mutate(gene_id = rownames(Counts_edgeR$counts)) %>% # add in gene_ids
    na.omit() %>%
    mutate(padj.sig = p.adjust(p.sig, method = "BH"),  # adjust p values
           padj.mu = p.adjust(p.mu, method = "BH"))
  
  write.csv(gamlss_NB, paste0('gamlss_rep/', tissue, rep,  '_repanals19_test.csv'),quote=F,row.names = F)
  

}



br <- extract_counts("Breast-MammaryTissue", runtable)
run_analysis(br$rep_pr, br$counts, replication = TRUE, rep = "_replication0.7","Breast-MammaryTissue")

run_analysis(br$init_pr, br$counts, replication = FALSE, rep = "_initial0.7", "Breast-MammaryTissue")

