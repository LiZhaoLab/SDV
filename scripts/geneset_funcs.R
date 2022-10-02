
get_topNtis <- function(gamlss_out, kind, thresh){
  out <- gamlss_out %>% 
    group_by(tissue) %>% 
    filter(bias == kind) %>%
    summarize(n = n()) %>%
    filter((n >= thresh)) %>% 
    pull(tissue)
  
  return(out)
}



# get m2g : mapping between gene and set
get_m2g <- function(cat, subcat = NULL){
  m_t2g <- msigdbr(species = "Homo sapiens", 
                   category = cat,
                   subcategory = subcat) %>%
    dplyr::select(gs_name, entrez_gene) %>% distinct()
  # remove gene sets with fewer than 10 genes
  good_sets <- m_t2g %>% group_by(gs_name) %>% summarise(n = n()) %>% filter(n>=10)
  return(m_t2g %>% filter(gs_name %in% good_sets$gs_name))
}

# perform hypergeometric test for a single term

enrich_1term <- function(term,
                         genes,
                         universe,
                         all_sets){
  
  
  # extract gene set from all sets
  gene_set <- all_sets %>% filter(gs_name == term)
  # using the variable names from phyper docs
  # might be useful to consult documentation
  # the number of SDV genes in gene set
  q <- sum(genes %in% gene_set$entrez_gene)

  # the number of expressed genes in gene set
  m <- sum(universe %in% gene_set$entrez_gene)
  # the number of expressed genes which are
  # not in gene set
  n <- sum(!(universe %in% gene_set$entrez_gene))
  # the number of SDV genes
  k <- length(genes)
  # the minus 1 is because we want greater than
  # or EQUAL to. Otherwise it's just greater than
  p <- phyper(q - 1, m, n, k, lower.tail = FALSE)
  # fraction of SDV in a tissue that live in the gene set
  fraction <- q/k
  out <- data.frame(p = p, count = q, fraction = fraction)
  return(out) 
  
}

split_gamlss <- function(tis, all_genes, gamlss_out){
  out <- gamlss_out %>% 
    filter((tissue == tis)) %>%
    merge(all_genes, by = "gene_id") 
  
  return(out)
}
enrich_tissue <- function(tis,
                          bias_sex,
                          gamlss_list,
                          all_sets,
                          meta = FALSE,
                          permute = FALSE){
  
  # extract tissue specific gamlss results.
  # some genes don't have entrez IDs
  # merge with dataframe containing entrezIDs
  # only leave genes found in genesets
  
  
  gamlss_tissue <- gamlss_list[[tis]]
  
  if (!permute){
    genes <- gamlss_tissue %>%
      filter(bias == bias_sex)
    genes <- na.omit(genes$entrez_id)
  } else {
    
    genes <- gamlss_tissue %>%
      mutate(bias = permute(gamlss_tissue$bias)) %>%
      filter(bias == bias_sex)
    genes <- na.omit(genes$entrez_id)
    
  }
  universe <- na.omit(gamlss_tissue$entrez_id) %>% unique()
  
  # only consider genes expressed in given tissue
  all_sets <- all_sets %>% 
    filter(entrez_gene %in% universe)
  # unique list of sets
  sets <- unique(all_sets$gs_name)
  # perform enrichment analysis for tissue
  enrichment_res <- map_dfr(sets,
                            enrich_1term,
                            genes = genes,
                            universe = universe,
                            all_sets = all_sets) %>%
    mutate(gs_name = sets) 
  
  # for meta analysis 
  if(meta){
    enrichment_res[, tis] <- enrichment_res$p
    enrichment_res <- enrichment_res[, c("gs_name", tis)]
    
  }
  
  
  return(enrichment_res)
}
# calculate the sum of log p values. This is the standard test
# statistic for fisher's method
get_sumlog <- function(row){
  row <- na.omit(row)
  return(sum(-2*log(row)))
}
#run a singular permutation
run_perm_1set <- function(set, 
                             tissue_list, 
                             bias_sex,
                             gamlss_list,
                             all_sets, 
                             nperm){
  print(set)
  set.seed(81)
  out <- data.frame(gs_name = set)
  for (i in seq(nperm)){
  
    go_df <- map(tissue_list,
                     enrich_tissue,
                     bias_sex=bias_sex,
                     gamlss_list=gamlss_list, 
                     all_sets=all_sets %>% filter(gs_name == set), 
                     meta = TRUE, 
                     permute = TRUE) %>% purrr::reduce(function(d1, d2) merge(d1, d2, by = "gs_name"))
    
    #go_df <- Reduce(function(d1, d2) merge(d1, d2, by = "gs_name"), agg_go) 
    # check to make sure that number of rows doesnt change on merging
    #stopifnot(dim(agg_go[[1]])[1] == dim(go_df)[1])
    
    sumlog <- data.frame(stat = get_sumlog(go_df[,2:dim(go_df)[2]]), 
                         gs_name = go_df$gs_name) 
    colnames(sumlog)[1] <- paste0("stat", i)
    out <- out %>% merge(sumlog, by="gs_name")
  }
  return(out)
}
# run 'n' permutations
run_perm_sumlog <-  function( 
                             nperm,
                             tissue_list, 
                             bias_sex, 
                             gamlss_list,
                             all_sets,
                             ncores){
  sets <- all_sets$gs_name %>% unique()
  perm_sumlog_res <- mclapply(sets, 
                              run_perm_1set, 
                              nperm = nperm,
                              tissue_list = tissue_list, 
                              bias_sex = bias_sex,
                              gamlss_list = gamlss_list,
                              all_sets = all_sets,
                              mc.cores = ncores) %>% purrr::reduce(function(d1, d2) rbind(d1, d2))
  # merge output
  #perm_sumlog_res <- Reduce(function(d1, d2) merge(d1, d2, by = "gs_name"), 
  #                        perm_sumlog_res) 
  #colnames(perm_sumlog_res)[2:dim(perm_sumlog_res)[2]] <- paste0('stat', seq(n))
  return(perm_sumlog_res)
}
# get empirical p-values given permutation output and observed sumlog values
get_empp_sumlog <- function(obs_sumlog, perm_sumlog_res){
  
  gs_list <- perm_sumlog_res$gs_name
  p <- c()
  for(term in gs_list){
    obs_sumlogt <- obs_sumlog %>% 
      filter(gs_name ==  term) %>% 
      pull(stat)
    row <- perm_sumlog_res %>% 
      filter(gs_name == term) %>%
      dplyr::select(- gs_name) %>%
      as.numeric() 
    p <- c(p, sum(row >= obs_sumlogt)/length(row))
    
  }
  return(data.frame(gs_name = gs_list, p = p))
}
# plot the overrep results in a given tissue using a dotplot
plot_enrichment_tis <- function(enrich, fdr) {
  tis_sig <- enrich %>% 
    filter(padj <= fdr) %>%
    arrange(p) %>%
    mutate(l10 = -log10(p))
  # remove all possible unneeded things in gene set names
  tis_sig$gs_name <- gsub("GO_", '', tis_sig$gs_name)
  tis_sig$gs_name <- gsub("_", ' ', tis_sig$gs_name)
  tis_sig$gs_name <- gsub("HALLMARK ", '', tis_sig$gs_name)
  
  # only plot significant sets
  to_plot <- tis_sig
  to_plot$gs_name <- factor(to_plot$gs_name, levels = rev(to_plot$gs_name))
  
  plot <- ggplot(to_plot, aes(x = fraction, y = gs_name)) + 
    
    theme_classic() + 
    geom_point(mapping = aes(col = l10), size = 5) + 
    scale_colour_gradient(low = "grey",
                          high = "steelblue4",
                          space = "Lab",
                          na.value = "grey50",
                          guide = "colourbar",
                          aesthetics = "colour") + 
    ylab('') + 
    labs(color=expression(-log[10]*p)) + 
    theme(strip.text.y = element_text(size = 10, angle = 0 ,face = "bold")) + 
    ggtitle(tis_sig$tissue)
  return(plot)
}



# plot results for all significant tissues after meta-analysis
plot_meta_heatmap <- function(go_df, metap, fdr, transpose=FALSE){
  # get significant gene sets
  sig_gs <- metap %>% 
    mutate(padj = p.adjust(p, method = "BH")) %>%
    filter(padj <= fdr)
  sig_go <- go_df %>%
    filter(gs_name %in% sig_gs$gs_name) %>%
    mutate(gs_name = gsub("GO_", "", gs_name)) %>%
    mutate(gs_name = gsub("_", " ", gs_name)) %>%
    mutate(gs_name = gsub("TARGET GENES", "", gs_name)) %>%
    mutate(gs_name = gsub("HALLMARK ", "", gs_name))
  
  rownames(sig_go) <- sig_go$gs_name
  sig_go$gs_name <- NULL
  sig_go <- -log10(as.matrix(sig_go))
  plot <- ggplot(melt(sig_go), aes(Var2, Var1, fill=value)) + 
    geom_tile() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_gradient(low = "gold",
                        high = "steelblue4",
                        space = "Lab",
                        guide = "colourbar",
                        aesthetics = "fill",
                        trans=scales::pseudo_log_trans(base = 10),
                        breaks=c(0, 1, 5, 10)) +
    xlab('') + ylab('') + labs(fill=expression(-log[10]*p))+ coord_fixed()
  if(transpose == TRUE){
    plot <- ggplot(melt(sig_go), aes(Var1, Var2, fill=value)) + 
      geom_tile() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_fill_gradient(low = "gold",
                          high = "steelblue4",
                          space = "Lab",
                          guide = "colourbar",
                          aesthetics = "fill",
                          trans=scales::pseudo_log_trans(base = 10),
                          breaks=c(0, 1, 5, 10)) +
      xlab('') + ylab('') + labs(fill=expression(-log[10]*p))+ coord_fixed()
  }
  return(plot)
}
plot_bar <- function(metap, fdr = 0.05, nperm, jama=FALSE){
  
  to_plot <- metap %>%
    mutate(gs_name = gsub("GO_", "", gs_name)) %>%
    mutate(gs_name = gsub("_", " ", gs_name)) %>%
    mutate(gs_name = gsub("TARGET GENES", "", gs_name)) %>%
    mutate(gs_name = gsub("HALLMARK ", "", gs_name),
           type = gsub("var.", "variability", type)) %>%
    mutate(type = gsub(" variability", "\nvariability", type)) %>%
    filter(padj < fdr)
  to_plot$type <- factor(to_plot$type, levels = c("male high\nvariability", "female high\nvariability"))
  
  #to_plot$gs_name <- factor(to_plot$gs_name, levels = to_plot$gs_name)
  #to_plot$lab <- ''
  #to_plot$lab[to_plot$p == 0] <- "§"
  to_plot$p[to_plot$p == 0] <- 1/(nperm + 1*nperm)
  to_plot <- to_plot %>% 
    mutate(l10 = -log10(p), 
           gs_name = reorder_within(gs_name, l10, type)) 
  
  xaxbreaks <- c(seq(0, log10(nperm)), log10(c(nperm + 1*nperm)))
  plot <- ggplot(to_plot, aes(x = l10, y = gs_name, fill = type)) +
    geom_col() + 
    facet_grid(type~., scales = "free", space = "free") + 
    theme_classic() + 
    scale_y_reordered() + 
    xlab(expression(-log[10]*p)) + 
    ylab('')+
    scale_x_continuous(breaks = xaxbreaks, labels = c(seq(0, log10(nperm)), paste0('>', log10(nperm)))) + 
    theme(strip.text.y = element_text(size = 10, angle = 0 ,face = "bold"),
          strip.background = element_blank(),
          legend.position = "none")
  #geom_text(aes(label = lab), hjust=-.5) + 
  
  if(jama){
    plot <- plot + scale_fill_jama() 
  }else{
    plot <- plot + scale_fill_d3()
  }
  return(plot)
  
}
