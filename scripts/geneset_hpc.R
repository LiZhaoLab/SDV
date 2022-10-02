library(dplyr)
library(matrixStats)
library(msigdbr)
library(reshape2)
library(gtools)
library(parallel)
library(purrr)
rm(list=ls())


setwd('/ru-auth/local/home/lzhao/Data_scratch/Khodursky/sexbias_var/')


source('scripts/geneset_funcs.R')

gamlss_out <- read.csv('gamlss_sum/gamlss_out_PC.csv', stringsAsFactors = F)
all_genes <- read.csv('ensembl_entrez_map.csv', stringsAsFactors = F)


# identify tissues with at least 20 SDV genes in each class

tissues_20_male <- gamlss_out %>% get_topNtis(kind = 'male high var.', thresh = 20)
tissues_20_female <- gamlss_out %>% get_topNtis(kind = 'female high var.', thresh = 20)



tissues <- unique(c(tissues_20_female, tissues_20_male))

gamlss_list <- map(tissues, split_gamlss, all_genes = all_genes, gamlss_out = gamlss_out)
names(gamlss_list) <- tissues


m2g_h <- get_m2g('H')
# filter genesets that are tissue specific
to_remove <- c("HALLMARK_PANCREAS_BETA_CELLS",
               "HALLMARK_ADIPOGENESIS",
               "HALLMARK_MYOGENESIS",
               "HALLMARK_SPERMATOGENESIS",
               "HALLMARK_BILE_ACID_METABOLISM",
               "HALLMARK_HEME_METABOLISM",
               "HALLMARK_XENOBIOTIC_METABOLISM")
m2g_h <- m2g_h %>% filter(! (gs_name %in% to_remove))


# get p-values in each tissue for each 
male_hall_df <- map(tissues_20_male,
                        enrich_tissue,
                        bias_sex = 'male high var.',
                        gamlss_list = gamlss_list,
                        all_sets = m2g_h,
                        meta = TRUE,
                        permute = FALSE) %>%
  purrr::reduce(function(d1, d2) merge(d1, d2, by = "gs_name"))

# male_hall_df <- purrr::reduce(function(d1, d2) merge(d1, d2, by = "gs_name"),
#                        male_agg_hall)


male_sumlog <- data.frame(stat = apply(male_hall_df[,2:dim(male_hall_df)[2]], 
                                       1, 
                                       get_sumlog)) %>%
  mutate(gs_name = male_hall_df$gs_name)


# number of permutations
n <- 10000

# run permutations
perm_sumlog_resm <-  run_perm_sumlog(n,
                                     tissues_20_male, 
                                     bias_sex = "male high var.",
                                     gamlss_list = gamlss_list,
                                     all_sets = m2g_h,
                                     ncores = 24)

# get empirical p-values
metap_resm <- get_empp_sumlog(male_sumlog, perm_sumlog_resm) %>%
  mutate(type = "male high var.",
         padj = p.adjust(p, method = "BH")) 

# get p-values in each tissue for each 
female_hall_df <- map(tissues_20_female, 
                          enrich_tissue,
                          bias_sex = 'female high var.',
                          gamlss_list = gamlss_list, 
                          all_sets = m2g_h, 
                          meta = TRUE, 
                          permute = FALSE) %>%
  purrr::reduce(function(d1, d2) merge(d1, d2, by = "gs_name"))

#female_hall_df <- purrr::reduce(function(d1, d2) merge(d1, d2, by = "gs_name"), 
#                         female_agg_hall) 


female_sumlog <- data.frame(stat = apply(female_hall_df[,2:dim(female_hall_df)[2]], 1, get_sumlog)) %>%
  mutate(gs_name = female_hall_df$gs_name)

# run permutations
perm_sumlog_resf <-  run_perm_sumlog(n,
                                     tissues_20_female, 
                                     bias_sex = "female high var.",
                                     gamlss_list = gamlss_list,
                                     all_sets = m2g_h,
                                     ncores = 24)

# get empirical p-values
metap_resf <- get_empp_sumlog(female_sumlog, perm_sumlog_resf) %>%
  mutate(type = "female high var.",
         padj = p.adjust(p, method = "BH"))



write.csv(metap_resf, "geneset_res/empp_hallmark_female_090622.csv", 
          quote = F, row.names = F)


write.csv(metap_resm, "geneset_res/empp_hallmark_male_090622.csv", 
          quote = F, row.names = F)



m2g_reg <- get_m2g('C3','TFT:GTRD')

set.seed(15)
female_reg_df <- map(tissues_20_female,
                         enrich_tissue,
                         bias_sex='female high var.',
                         gamlss_list=gamlss_list,
                         all_sets=m2g_reg,
                         meta = TRUE) %>%
  purrr::reduce(function(d1, d2) merge(d1, d2, by = "gs_name"))


#female_reg_df <- purrr::reduce(function(d1, d2) merge(d1, d2, by = "gs_name"),female_agg_reg)


female_sumlog_reg <- data.frame(stat = apply(female_reg_df[,2:dim(female_reg_df)[2]], 1, get_sumlog)) %>%
  mutate(gs_name = female_reg_df$gs_name)


perm_sumlog_resfreg <- run_perm_sumlog(n,
                                       tissues_20_female,
                                       bias_sex = "female high var.",
                                       gamlss_list = gamlss_list,
                                       all_sets = m2g_reg,
                                       ncores = 24)


metap_resf_reg <- get_empp_sumlog(female_sumlog_reg, perm_sumlog_resfreg) %>%
  mutate(padj = p.adjust(p, method = "BH"),
         type = "female high var.")



male_reg_df <- map(tissues_20_male, enrich_tissue,
                       bias_sex='male high var.',
                       gamlss_list=gamlss_list,
                       all_sets=m2g_reg,
                       meta = TRUE,
                       permute = FALSE)  %>%
  purrr::reduce(function(d1, d2) merge(d1, d2, by = "gs_name"))


#male_reg_df <- purrr::reduce(function(d1, d2) merge(d1, d2, by = "gs_name"), male_agg_reg)


male_sumlog_reg <- data.frame(stat = apply(male_reg_df[,2:dim(male_reg_df)[2]], 1, get_sumlog)) %>%
  mutate(gs_name = male_reg_df$gs_name)



perm_sumlog_resmreg <-  run_perm_sumlog(n,
                                        tissues_20_male,
                                        bias_sex = "male high var.",
                                        gamlss_list = gamlss_list,
                                        all_sets = m2g_reg,
                                        ncores = 24)

metap_resm_reg <- get_empp_sumlog(male_sumlog_reg, perm_sumlog_resmreg) %>%
  mutate(padj = p.adjust(p, method = "BH"),
         type = "male high var.")



write.csv(metap_resf_reg, "geneset_res/empp_reg_female_090622.csv", 
          quote = F, row.names = F)

write.csv(metap_resm_reg, "geneset_res/empp_reg_male_090622.csv", 
          quote = F, row.names = F)

