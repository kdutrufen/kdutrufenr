library(tidyverse)
library(CEMiTool)
library(WGCNA)
library(parallelMap)

# Databases
WikiPathways_2019_Human <- "~/Dropbox/GMT/WikiPathways_2019_Human.txt" # WikiPathways_2019_Human gene sets
WikiPathways_2019_Mouse <- "~/Dropbox/GMT/WikiPathways_2019_Mouse.txt" # WikiPathways_2019_Mouse gene sets
KEGG_2019_Human <- "~/Dropbox/GMT/KEGG_2019_Human.txt" # KEGG_2019_Human gene sets
KEGG_2019_Mouse <- "~/Dropbox/GMT/KEGG_2019_Mouse_Human.txt" # KEGG_2019_Mouse gene sets
LINCS_L1000_Chem_Pert_down <- "~/Dropbox/GMT/LINCS_L1000_Chem_Pert_down_Human.txt" # LINCS_L1000_Chem_Pert_down gene sets
LINCS_L1000_Chem_Pert_up <- "~/Dropbox/GMT/LINCS_L1000_Chem_Pert_up_Human.txt" # LINCS_L1000_Chem_Pert_up gene sets
Reactome_2016 <- "~/Dropbox/GMT/Reactome_2016_Human.txt" # Reactome_2016 gene sets
GO_Molecular_Function_2018 <- "~/Dropbox/GMT/GO_Molecular_Function_2018_Human.txt" # GO_Molecular_Function_2018 gene sets
GO_Cellular_Component_2018 <- "~/Dropbox/GMT/GO_Cellular_Component_2018_Human.txt" # GO_Cellular_Component_2018 gene sets
GO_Biological_Process_2018 <- "~/Dropbox/GMT/GO_Biological_Process_2018_Human.txt" # GO_Biological_Process_2018 gene sets

# ora_analysis <- function(genes, database = KEGG_2019_Human, pvalue_cutoff = 1, qvalue_cutoff = 1, p.adjust_cutoff = 0.05, count_at_least = 3){
#   # library(CEMiTool)
#   library(clusterProfiler)
#   gmt_list <- database %>% CEMiTool::read.gmt()
#   custom_pathways <- gmt_list %>% split(f = .$term) %>% map(dplyr::select, gene) %>% map(pull)
#   allgenes <- gmt_list[, "gene"] %>% unique()
#   enrichment_result <- genes %>% map(clusterProfiler::enricher, pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff, universe = allgenes, TERM2GENE = gmt_list)
#   clProf.df <- enrichment_result %>% map(as.data.frame) %>% bind_rows(.id = "Cluster")
#   clProf.df <- clProf.df %>% filter(p.adjust < p.adjust_cutoff & Count >= count_at_least)
#   geneClusters <- genes
#   compareCluster_output <- new("compareClusterResult", compareClusterResult = clProf.df, geneClusters = geneClusters, fun = "enricher",
#                                .call = match.call(expand.dots = TRUE))
#   return(compareCluster_output)
# }

ora_analysis <- function(genes, database = KEGG_2019_Human, pvalue_cutoff = 1, qvalue_cutoff = 1, p.adjust_cutoff = 0.1, count_at_least = 3){
  library(clusterProfiler)
  gmt_list <- database %>% read.gmt()
  custom_pathways <- gmt_list %>% split(f = .$ont) %>% map(dplyr::select, gene) %>% map(pull)
  allgenes <- gmt_list[, "gene"] %>% unique()
  enrichment_result <- genes %>% map(clusterProfiler::enricher, pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff, universe = allgenes, TERM2GENE = gmt_list)
  clProf.df <- enrichment_result %>% map(as.data.frame) %>% bind_rows(.id = "Cluster")
  clProf.df <- clProf.df %>% filter(p.adjust < p.adjust_cutoff & Count >= count_at_least)
  geneClusters <- genes
  compareCluster_output <- new("compareClusterResult", compareClusterResult = clProf.df, geneClusters = geneClusters, fun = "enricher",
                               .call = match.call(expand.dots = TRUE))
  return(compareCluster_output)
}

gsea_analysis <- function(sorted_genes, database = GO_Biological_Process_2018, p_AdjustMethod = "BH", min_Size = 15, max_Size = 500, n_perm = 10000, n_proc = 0, pvalue_cutoff = 1, qvalue_cutoff = 1, p.adjust_cutoff = 0.1, exponent = 1) {
  library(clusterProfiler)
  library(stats)
  gmt_list <- database %>% clusterProfiler::read.gmt()
  gsea_result <- sorted_genes %>% GSEA(TERM2GENE = gmt_list, verbose = F, nPerm = n_perm, minGSSize = min_Size, maxGSSize = max_Size, pvalueCutoff = pvalue_cutoff, pAdjustMethod = p_AdjustMethod) %>% suppressWarnings()
  geneSet_ID <- gsea_result@result$Description
  gseaplot_results <- map(1:nrow(gsea_result@result), function(i) gsea_result %>% gseaplot(geneSetID = geneSet_ID[i]) + cowplot::draw_figure_label(label = geneSet_ID[i], size = 15)) %>% purrr::set_names(geneSet_ID)
  results <- list(gsea_results = gsea_result, 
                  gseaplot_results = gseaplot_results)
  return(results)
}

go_gsea_analysis <- function(sorted_genes, database = GO_Biological_Process_2018, p_AdjustMethod = "BH", min_Size = 15, max_Size = 500, n_perm = 10000, n_proc = 0, pvalue_cutoff = 1, qvalue_cutoff = 1, p.adjust_cutoff = 0.1, exponent = 1) {
  library(clusterProfiler)
  library(stats)
  gmt_list <- database %>% clusterProfiler::read.gmt()
  gsea_result <- sorted_genes %>% GSEA(TERM2GENE = gmt_list, verbose = F, nPerm = n_perm, minGSSize = min_Size, maxGSSize = max_Size, pvalueCutoff = pvalue_cutoff, pAdjustMethod = p_AdjustMethod) %>% suppressWarnings()
  geneSet_ID <- gsea_result@result$Description
  gseaplot_results <- map(1:nrow(gsea_result@result), function(i) gsea_result %>% gseaplot(geneSetID = geneSet_ID[i]) + cowplot::draw_figure_label(label = geneSet_ID[i], size = 15)) %>% purrr::set_names(geneSet_ID)
  gsea_result@result <- gsea_result@result %>% rownames_to_column()
  gsea_result@result <- gsea_result@result %>% dplyr::select(-c(rowname, ID))
  gsea_result@result <- gsea_result@result %>% separate(col = "Description", into = c("Description", "ID"), sep = " \\(")
  gsea_result@result <- gsea_result@result %>% mutate(ID = ID %>% str_remove(pattern = "\\)"))
  gsea_result@result <- gsea_result@result %>% mutate(ID_1 = ID)
  gsea_result@result <- gsea_result@result %>% column_to_rownames(var = "ID_1")
  gsea_result@result <- gsea_result@result %>% dplyr::select(ID, everything())
  results <- list(gsea_results = gsea_result, 
                  gseaplot_results = gseaplot_results)
  return(results)
}

# parallel_cemitool_pipeline funciton
parallel_cemitool_pipeline <- function(n_cores = detectCores() - 1, cem = cem, database = KEGG_2019_Human) {
  library(CEMiTool)
  library(tidyverse)
  library(parallelMap)
  no_cores <- detectCores() - 1
  parallelStartSocket(no_cores)
  cem <- find_modules(cem, cor_method = "pearson", cor_function = "cor", eps = 0.1, min_ngen = 20, merge_similar = TRUE, diss_thresh = 0.75, network_type = "unsigned", tom_type = "signed", verbose = FALSE, force_beta = TRUE, set_beta = NULL)
  cem <- cem %>% plot_mean_k(title = "Mean connectivity")
  cem <- cem %>% plot_beta_r2()
  cem <- cem %>% plot_profile()
  genes <- cem %>% unclass() %>% attr("selected_genes")
  
  sample_annotation(cem, sample_name_column = "SampleName", class_column = "Class") <- sample_annotation
  cem <- cem %>% mod_gsea()
  cem <- cem %>% plot_gsea()
  # Adding ORA analysis
  gmt_in <- database %>% read_gmt()
  cem <- cem %>% mod_ora(gmt_in)
  # cem@ora <- cem@ora %>% filter(Count >= min_count) %>% filter(p.adjust <= p_adjust_cutoff)
  cem <- cem %>% plot_ora(n = 20, pv_cut = 0.05)
  
  # read interactions
  int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
  int_df <- int_fname %>% read.delim()
  interactions_data(cem) <- int_df # add interactions
  cem <- cem %>% plot_interactions() # generate plot
  
  parallelStop()
  return(cem)
}

genes_in_modules <- function(cem = cem, split_by = split_by, selected_column = selected_column) {
  cem %>%
    module_genes() %>%
    split(f = .[[split_by]]) %>%
    purrr::map(dplyr::select, selected_column) %>%
    purrr::map(pull)
}

convert_some_rownames <- function(countData, conversion_table, from_ct, to_ct) {
  for (i in 1:nrow(countData)) {
    if (is.na(conversion_table[[to_ct]][match(rownames(countData)[i], conversion_table[[from_ct]])]) == FALSE) {
      rownames(countData)[i] <- conversion_table[[to_ct]][match(rownames(countData)[i], conversion_table[[from_ct]])]
    }
  }
  return(countData)
}

estimate_elbow <- function(n_genes = n_genes, dissTOM = dissTOM, deepSplit = deepSplit, gene_names = gene_names) {
  geneTree <- dissTOM %>% as.dist() %>% hclust(method = "average")
  dynamicModstest_list <- map(n_genes, function(i) geneTree %>% cutreeDynamic(distM = dissTOM, deepSplit = deepSplit, pamRespectsDendro = FALSE, minClusterSize = i))
  dynamicModstest <- map(seq_along(dynamicModstest_list), function(i) dynamicModstest_list[[i]] %>% max()) %>% unlist()
  dynamicColorstest <- map(seq_along(dynamicModstest_list), function(i) dynamicModstest_list[[i]] %>% labels2colors())
  dynamicColorstest1 <- map(seq_along(dynamicModstest_list), function(i) data.frame(dynamicColors = dynamicColorstest[[i]], genes = rownames(adjacency)))
  premoduleblockstest <- map(seq_along(dynamicColorstest1), function(i) dynamicColorstest1[[i]] %>% split(f = .$dynamicColors) %>% map(dplyr::select, -dynamicColors) %>% map(pull))
  gene_names <- gene_names
  inModuletest <- map(seq_along(premoduleblockstest), function(i) map(seq_along(premoduleblockstest[[i]]), function(j) gene_names %in% premoduleblockstest[[i]][[j]]))
  # Select the corresponding Topological Overlap
  modTOMtest <- map(seq_along(inModuletest), function(i) map(seq_along(inModuletest[[i]]), function(j) TOM[inModuletest[[i]][[j]], inModuletest[[i]][[j]]]))
  for (i in 1:length(modTOMtest)) {
    for (j in 1:length(modTOMtest[[i]])) {
      dimnames(modTOMtest[[i]][[j]]) <- list(premoduleblockstest[[i]][[j]], premoduleblockstest[[i]][[j]])
    }
  }
  ssd <- function(x) sum((x - mean(x))^2)
  SSD <- map(seq_along(modTOMtest), function(i) map(seq_along(modTOMtest[[i]]), function(j) modTOMtest[[i]][[j]] %>% ssd())) %>% map(unlist) %>% map(sum) %>% unlist()
  dynamicModstest_df <- data.frame(minimum_n_genes = n_genes, n_modules = dynamicModstest, SSD = SSD)
  gg_plot <- dynamicModstest_df %>% ggplot(aes(x = n_modules, y = SSD)) + geom_line() + geom_point() + theme_bw()
  dynamicModstest_df
  gg_plot
}

parallel_estimate_elbow <- function(n_genes = n_genes, dissTOM = dissTOM, deepSplit = deepSplit, gene_names = gene_names) {
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  # Initiate cluster
  parallelStartSocket(no_cores) 
  geneTree <- dissTOM %>% as.dist() %>% hclust(method = "average")
  dynamicModstest_list <- map(n_genes, function(i) geneTree %>% cutreeDynamic(distM = dissTOM, deepSplit = deepSplit, pamRespectsDendro = FALSE, minClusterSize = i))
  dynamicModstest <- map(seq_along(dynamicModstest_list), function(i) dynamicModstest_list[[i]] %>% max()) %>% unlist()
  dynamicColorstest <- map(seq_along(dynamicModstest_list), function(i) dynamicModstest_list[[i]] %>% labels2colors())
  dynamicColorstest1 <- map(seq_along(dynamicModstest_list), function(i) data.frame(dynamicColors = dynamicColorstest[[i]], genes = rownames(adjacency)))
  premoduleblockstest <- map(seq_along(dynamicColorstest1), function(i) dynamicColorstest1[[i]] %>% split(f = .$dynamicColors) %>% map(dplyr::select, -dynamicColors) %>% map(pull))
  gene_names <- gene_names
  inModuletest <- map(seq_along(premoduleblockstest), function(i) map(seq_along(premoduleblockstest[[i]]), function(j) gene_names %in% premoduleblockstest[[i]][[j]]))
  # Select the corresponding Topological Overlap
  modTOMtest <- map(seq_along(inModuletest), function(i) map(seq_along(inModuletest[[i]]), function(j) TOM[inModuletest[[i]][[j]], inModuletest[[i]][[j]]]))
  for (i in 1:length(modTOMtest)) {
    for (j in 1:length(modTOMtest[[i]])) {
      dimnames(modTOMtest[[i]][[j]]) <- list(premoduleblockstest[[i]][[j]], premoduleblockstest[[i]][[j]])
    }
  }
  ssd <- function(x) sum((x - mean(x))^2)
  SSD <- map(seq_along(modTOMtest), function(i) map(seq_along(modTOMtest[[i]]), function(j) modTOMtest[[i]][[j]] %>% ssd())) %>% map(unlist) %>% map(sum) %>% unlist()
  dynamicModstest_df <- data.frame(minimum_n_genes = n_genes, n_modules = dynamicModstest, SSD = SSD)
  gg_plot <- dynamicModstest_df %>% ggplot(aes(x = n_modules, y = SSD)) + geom_line() + geom_point() + theme_bw()
  parallelStop()
  gg_plot
}

estimate_threshold_cut <- function(thresh_cut = thresh_cut, datExpr = datExpr, dynamicColors = dynamicColors, gene_names = gene_names, TOM = TOM) {
  mergetest_list <- map(thresh_cut, function(i) datExpr %>% mergeCloseModules(colors = dynamicColors, cutHeight = i, verbose = 3) %>% .[["colors"]])
  mergetest <- mergetest_list %>% map(table) %>% map(length) %>% unlist()
  mergetest_1 <- map(seq_along(mergetest_list), function(i) data.frame(dynamicColors = mergetest_list[[i]], genes = rownames(adjacency)))
  merge_blockstest <- map(seq_along(mergetest_1), function(i) mergetest_1[[i]] %>% split(f = .$dynamicColors) %>% map(dplyr::select, -dynamicColors) %>% map(pull))
  gene_names <- datExpr %>% colnames()
  inModule_merge_test <- map(seq_along(merge_blockstest), function(i) map(seq_along(merge_blockstest[[i]]), function(j) gene_names %in% merge_blockstest[[i]][[j]]))
  # Select the corresponding Topological Overlap
  modTOM_merge_test <- map(seq_along(inModule_merge_test), function(i) map(seq_along(inModule_merge_test[[i]]), function(j) TOM[inModule_merge_test[[i]][[j]], inModule_merge_test[[i]][[j]]]))
  # Function to calculate the sum of squared deviations (from the mean)
  ssd <- function(x) sum((x - mean(x))^2)
  SSD_merge <- map(seq_along(modTOM_merge_test), function(i) map(seq_along(modTOM_merge_test[[i]]), function(j) modTOM_merge_test[[i]][[j]] %>% ssd())) %>% map(unlist) %>% map(sum) %>% unlist()
  mergetest_df <- data.frame(
    cor_thres = thresh_cut,
    n_modules = mergetest,
    SSD_merge = SSD_merge
  )
  gg_plot1 <- mergetest_df %>% ggplot(aes(x = n_modules, y = SSD_merge)) + geom_line() + geom_point() + theme_bw()
  # gg_plot2 <- mergetest_df %>% ggplot(aes(x = cor_thres, y = SSD_merge)) + geom_line() + geom_point() + theme_bw()
  # return(list(gg_plot1, gg_plot2))
  return(gg_plot1)
}

parallel_estimate_threshold_cut <- function(thresh_cut = thresh_cut, datExpr = datExpr, dynamicColors = dynamicColors, gene_names = gene_names, TOM = TOM) {
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  # Initiate cluster
  parallelStartSocket(no_cores) 
  mergetest_list <- map(thresh_cut, function(i) datExpr %>% mergeCloseModules(colors = dynamicColors, cutHeight = i, verbose = 3) %>% .[["colors"]])
  mergetest <- mergetest_list %>% map(table) %>% map(length) %>% unlist()
  mergetest_1 <- map(seq_along(mergetest_list), function(i) data.frame(dynamicColors = mergetest_list[[i]], genes = rownames(adjacency)))
  merge_blockstest <- map(seq_along(mergetest_1), function(i) mergetest_1[[i]] %>% split(f = .$dynamicColors) %>% map(dplyr::select, -dynamicColors) %>% map(pull))
  gene_names <- datExpr %>% colnames()
  inModule_merge_test <- map(seq_along(merge_blockstest), function(i) map(seq_along(merge_blockstest[[i]]), function(j) gene_names %in% merge_blockstest[[i]][[j]]))
  # Select the corresponding Topological Overlap
  modTOM_merge_test <- map(seq_along(inModule_merge_test), function(i) map(seq_along(inModule_merge_test[[i]]), function(j) TOM[inModule_merge_test[[i]][[j]], inModule_merge_test[[i]][[j]]]))
  # Function to calculate the sum of squared deviations (from the mean)
  ssd <- function(x) sum((x - mean(x))^2)
  SSD_merge <- map(seq_along(modTOM_merge_test), function(i) map(seq_along(modTOM_merge_test[[i]]), function(j) modTOM_merge_test[[i]][[j]] %>% ssd())) %>% map(unlist) %>% map(sum) %>% unlist()
  mergetest_df <- data.frame(
    cor_thres = thresh_cut,
    n_modules = mergetest,
    SSD_merge = SSD_merge
  )
  gg_plot1 <- mergetest_df %>% ggplot(aes(x = n_modules, y = SSD_merge)) + geom_line() + geom_point() + theme_bw()
  # gg_plot2 <- mergetest_df %>% ggplot(aes(x = cor_thres, y = SSD_merge)) + geom_line() + geom_point() + theme_bw()
  # return(list(gg_plot1, gg_plot2))
  parallelStop()
  gg_plot1
}

get_genes_from_pathway <- function(enrichment_df = enrichment_df, pathway = pathway, comparison = comparison){enrichment_df %>% filter(Description == pathway) %>% filter(Cluster == comparison) %>% dplyr::select(Gene_Symbol) %>% pull %>% strsplit(split = "/") %>% unlist}
get_genes_from_pathways <- function(enrichment_df = enrichment_df, pathway = pathway, comparison = comparison) {enrichment_df %>% filter(Description == pathway) %>% filter(Cluster == comparison) %>% dplyr::select(geneID) %>%  pull() %>% strsplit(split = "/") %>% unlist()}

plot_genes_boxplot <- function(df, genes, design, levels, converter, convert_from, convert_to) {
  df <- df %>% cpm() %>% as.data.frame() %>% rownames_to_column() %>% filter(rowname %in% c(genes)) %>% distinct(rowname, .keep_all = TRUE) %>% arrange(rowname) %>% column_to_rownames()  
  if (nrow(df) %% 4 == 0) {
    df_list <- map(1:floor(length(genes) / 4), function(i) df[(4 * i - 3):(4 * i), ])
  } else {
    df_list <- map(1:floor(length(genes) / 4), function(i) df[(4 * i - 3):(4 * i), ])
    df_list[[(length(df_list) + 1)]] <- df[((length(df_list) * 4) + 1):length(genes), ]
  }
  df_list <- df_list %>% map(t) %>% map(function(x) as.data.frame(cbind(x, "Condition" = as.character(design, levels = levels))))
  df_list <- df_list %>% map(function(x) {x %>% gather(Gene, CPM, -Condition)}) %>% map(function(x) {merge(x, converter, by.x = "Gene", by.y = convert_from)})
  appender <- function(string, suffix = converter[[convert_to]][match(string, converter[[convert_from]])]) paste0(string, "\n (", suffix, ")")
  map(df_list, function(x) ggplot(x, aes(x = Condition, y = as.numeric(CPM), color = Condition)) + geom_boxplot(alpha = .5) + geom_jitter(alpha = .5, width = .1) + theme_bw() + scale_x_discrete(limits = levels) + scale_fill_discrete(breaks = levels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab(label = "Read counts (cpm)") + facet_wrap(~Gene, scales = "free", ncol = 2, labeller = as_labeller(appender)))
}

plot_genes_violinplot <- function(df, genes, design, levels, converter, convert_from, convert_to) {
  df <- df %>% cpm() %>% as.data.frame() %>% rownames_to_column() %>% filter(rowname %in% c(genes)) %>% distinct(rowname, .keep_all = TRUE) %>% arrange(rowname) %>% column_to_rownames()  
  if (nrow(df) %% 4 == 0) {
    df_list <- map(1:floor(length(genes) / 4), function(i) df[(4 * i - 3):(4 * i), ])
  } else {
    df_list <- map(1:floor(length(genes) / 4), function(i) df[(4 * i - 3):(4 * i), ])
    df_list[[(length(df_list) + 1)]] <- df[((length(df_list) * 4) + 1):length(genes), ]
  }
  df_list <- df_list %>% map(t) %>% map(function(x) as.data.frame(cbind(x, "Condition" = as.character(design, levels = levels))))
  df_list <- df_list %>% map(function(x) {x %>% gather(Gene, CPM, -Condition)}) %>% map(function(x) {merge(x, converter, by.x = "Gene", by.y = convert_from)})
  appender <- function(string, suffix = converter[[convert_to]][match(string, converter[[convert_from]])]) paste0(string, "\n (", suffix, ")")
  map(df_list, function(x) ggplot(x, aes(x = Condition, y = as.numeric(CPM), color = Condition)) + geom_violin(alpha = .5) + geom_jitter(alpha = .5, width = .1) + theme_bw() + scale_x_discrete(limits = levels) + scale_fill_discrete(breaks = levels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab(label = "Read counts (cpm)") + facet_wrap(~Gene, scales = "free", ncol = 2, labeller = as_labeller(appender)))
}

library(ballgown)
# adapted from ballgown
get_Attribute_Field = function (x, field, attrsep = "; ") 
{
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}

library(stringr)
ora_to_df <- function(ora = ora) {unclass(ora) %>% attr("ora") %>% as.data.frame()}
filter_ora <- function(df = df, p.adjust = p.adjust, Count = Count) {df %>% filter(p.adjust < 0.05) %>% filter(Count > 2)}
add_ensembl <- function(geneID = geneID, convert_table = convert_table, from = from, to = to, split = "/") {geneID %>% strsplit(split = split) %>% map(function(x) as.character(convert_table[[to]][ match(x, convert_table[[from]]) ])) %>% map(function(x) paste(x, collapse = split)) %>% unlist()}
add_ensembl_to_df <- function(df = df, geneID = geneID, convert_table = convert_table, from = from, to = to, split = "/") {df[[geneID]] %>% strsplit(split = split) %>% map(function(x) as.character(convert_table[[to]][ match(x, convert_table[[from]]) ])) %>% map(function(x) paste(x, collapse = split)) %>% unlist()}
genes_to_title <- function(geneID = geneID, split = "/") {geneID %>% strsplit(split = split) %>% map(str_to_title) %>%  map(paste, collapse = split) %>%  unlist()}
order_genes <- function(geneID = geneID, split = "/") {geneID %>% strsplit(split = split) %>% map(sort) %>% map(paste, collapse = split) %>% unlist()}
order_genes_in_df <- function(df = df, geneID = geneID, split = "/") {df[[geneID]] %>% strsplit(split = split) %>% map(sort) %>% map(paste, collapse = split) %>% unlist()}

library(clusterProfiler)
organize_clusterprofiler_results <- function(clusterprofiler_result, convert_table = convert_table){
  library(tidyverse)
  library(formattable)
  clusterprofiler_result@compareClusterResult <- clusterprofiler_result@compareClusterResult %>% mutate(pvalue = pvalue %>% formattable::scientific(digits = 2, format = "e"), p.adjust = p.adjust %>% formattable::scientific(digits = 2, format = "e"), qvalue = qvalue %>% formattable::scientific(digits = 2, format = "e")
  )
  clusterprofiler_result@compareClusterResult <- clusterprofiler_result@compareClusterResult %>% filter(p.adjust < 0.05) %>% filter(Count > 2)
  df <- KEGG_DEGs %>% as.data.frame()
  df <- df %>% mutate(Gene_Ensembl = geneID %>% add_ensembl(convert_table = convert_table, from = "ENTREZID", to = "ENSEMBL") )
  df <- df %>% mutate(Gene_Symbol = geneID %>% add_ensembl(convert_table = convert_table, from = "ENTREZID", to = "SYMBOL") )
  df <- df %>% mutate(Gene_Symbol = Gene_Symbol %>% order_genes() )
  df <- df %>% mutate(Gene_Ensembl = Gene_Symbol %>% add_ensembl(convert_table = convert_table, from = "SYMBOL", to = "ENSEMBL") )
  df <- df %>% mutate(geneID = Gene_Symbol %>% add_ensembl(convert_table = convert_table, from = "SYMBOL", to = "ENTREZID") )
  return(df)
}

# adapted from ballgown
# get_Attribute_Field
get_Attribute_Field <- function(x, field, attrsep = "; ") {
  s <- strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a <- strsplit(atts, split = "=", fixed = TRUE)
    m <- match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv <- a[[m]][2]
    }
    else {
      rv <- as.character(NA)
    }
    return(rv)
  })
}


setGeneric('plot_gsea', function(cem, ...) {
  standardGeneric('plot_gsea')
})

#' @rdname plot_gsea
setMethod('plot_gsea', signature('CEMiTool'),
          function(cem, pv_cut=0.05) {
            if(length(unique(cem@module$modules)) == 0){
              stop("No modules in CEMiTool object! Did you run find_modules()?")
            }
            if(length(cem@enrichment) == 0){
              stop("No GSEA data! Did you run mod_gsea()?")
            }
            if(all(unlist(lapply(cem@enrichment, nrow))) == 0){
              warning("No modules were enriched for any classes. Unable to plot enrichment.")
              return(cem)
            }
            #cem <- get_args(cem, vars=mget(ls()))
            
            stats <- names(cem@enrichment)
            enrichment <- lapply(cem@enrichment, function(stat){
              stat[is.na(stat)] <- 0
              rownames(stat) <- stat[,1]
              stat[,1] <- NULL
              return(stat)
            })
            names(enrichment) <- stats
            
            pval <- enrichment[['padj']]
            nes <- enrichment[['nes']]
            
            pval <- pval[rowSums(pval < pv_cut) >= 1, , drop=FALSE]
            nes <- nes[rownames(pval), , drop=FALSE]
            
            # check if there is any signif. module
            if(nrow(nes) < 0){
              stop("No significant modules found!")
            }
            
            custom_pal <- c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                            "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                            "#D6604D", "#B2182B", "#67001F")
            custom_pal <- colorRampPalette(custom_pal)(200)
            
            nes <- as.matrix(nes)
            pval <- as.matrix(pval)
            nes[which(pval > pv_cut, arr.ind=TRUE)] <- 0
            
            if(nrow(nes) > 2){
              row_order <- rownames(nes)[hclust(dist(nes))$order]
            } else {
              row_order <- rownames(nes)
            }
            
            nes_melted <- reshape2::melt(nes)
            colnames(nes_melted) <- c("Module", "Class", "NES")
            nes_melted$Module <- factor(nes_melted$Module, levels=row_order)
            # nes_melted$Class <- as.character(nes_melted$Class)
            max_abs_nes <- max(abs(nes_melted$NES))
            res <- ggplot(nes_melted, aes_(x=~Class, y=~Module, size=~abs(NES), fill=~NES)) +
              geom_point(color = "white", shape=21) +
              scale_fill_gradientn(colours=custom_pal, space = "Lab",
                                   limits=c(-max_abs_nes, max_abs_nes)) +
              scale_size(range=c(0,9), limits=c(0, NA)) +
              guides(size="none") +
              theme_minimal() +
              theme(panel.grid.major = element_blank()) +
              scale_x_discrete(position = "top")
            res_list <- list(enrichment_plot=res)
            cem@enrichment_plot <- res_list
            
            return(cem)
          })
# Expression patterns in modules
# setGeneric('plot_profile', function(cem, ...) {
#   standardGeneric('plot_profile')
# })

#' @rdname plot_profile
# setMethod('plot_profile', signature('CEMiTool'),
#           function(cem, order_by_class=TRUE, center_func='mean') {
#             if(!tolower(center_func) %in% c("mean", "median")){
#               stop("Invalid center_func type. Valid values are 'mean' and 'median'")
#             }
#             modules <- unique(cem@module$modules)
#             if(is.null(modules)){
#               stop("No modules in this CEMiTool object.")
#             }
#             #vars <- mget(ls())
#             #vars$modules <- NULL
#             #cem <- get_args(cem=cem, vars=vars)
#             
#             modules <- modules[order(as.numeric(stringr::str_extract(modules, "\\d+")))]
#             expr <- expr_data(cem)
#             annot <- sample_annotation(cem)
#             sample_name_column <- cem@sample_name_column
#             class_column <- cem@class_column
#             mod_cols <- mod_colors(cem)
#             plots <- lapply(modules, function(mod){
#               # subsets from expr all genes inside module mod
#               genes <- cem@module[cem@module[,'modules']==mod, 'genes']
#               expr[, 'id'] <- rownames(expr)
#               mod_expr <- data.table::melt(expr[genes,], 'id',
#                                            variable.name='sample',
#                                            value.name='expression')
#               
#               # initialize plot base layer
#               g <- ggplot(mod_expr, aes_(x=~sample, y=~expression))
#               
#               # adds different background colours if annot is provided
#               if (nrow(annot)!=0) {
#                 if (order_by_class) {
#                   # sorts data.frame by class name
#                   annot <- annot[order(annot[, class_column]),]
#                 }
#                 annot[, sample_name_column] <- factor(annot[, sample_name_column],
#                                                       levels=annot[, sample_name_column])
#                 mod_expr[, 'sample'] <- factor(mod_expr[, 'sample'],
#                                                levels=annot[, sample_name_column])
#                 
#                 # y positioning of background tiles
#                 y_pos <- mean(mod_expr[, 'expression'])
#                 
#                 # reinitialize base layer adding background tiles
#                 g <- ggplot(mod_expr, aes_(x=~sample, y=~expression)) +
#                   geom_tile(data=annot, alpha=0.3, height=Inf,
#                             aes(x=get(sample_name_column), y=y_pos,
#                                 fill=as.factor(get(class_column))))
#               }
#               
#               # adding lines
#               g <- g + geom_line(aes_(group=~id), alpha=0.2, colour=mod_cols[mod]) +
#                 stat_summary(aes(group=1), size=1, fun.y=get(tolower(center_func)), geom='line')
#               
#               # custom theme
#               g <- g + theme(plot.title=element_text(lineheight=0.8,
#                                                      face='bold',
#                                                      colour='black',
#                                                      size=15),
#                              axis.title=element_text(face='bold',
#                                                      colour='black',
#                                                      size=15),
#                              axis.text.y=element_text(angle=0,
#                                                       vjust=0.5,
#                                                       size=12),
#                              axis.text.x=element_text(angle=90,
#                                                       vjust=0.5,
#                                                       size=12),
#                              panel.grid=element_blank(),
#                              legend.title=element_blank(),
#                              legend.text=element_text(size = 8),
#                              legend.background=element_rect(fill='gray90',
#                                                             size=0.5,
#                                                             linetype='dotted'),
#                              legend.position='bottom'
#               )
#               # title
#               g <- g + ggtitle(mod)
#               
#               return(g)
#             })
#             names(plots) <- modules
#             cem@profile_plot <- plots
#             return(cem)
#           })
# 
# 
# # plot ora results
# setGeneric('plot_ora', function(cem, ...) {
#   standardGeneric('plot_ora')
# })

#' @rdname plot_ora
# setMethod('plot_ora', signature('CEMiTool'),
#           function(cem, n=20, pv_cut=0.05, ...){
#             if(length(unique(cem@module$modules)) == 0){
#               stop("No modules in CEMiTool object! Did you run find_modules()?")
#             }
#             if(nrow(cem@ora) == 0){
#               stop("No ORA data! Did you run mod_ora()?")
#             }
#             
#             #cem <- get_args(cem=cem, vars=mget(ls()))
#             ora_splitted <- split(cem@ora, cem@ora$Module)
#             mod_cols <- mod_colors(cem)
#             res <- lapply(ora_splitted, function(x){
#               plot_ora_single(head(x, n=n),
#                               pv_cut=pv_cut,
#                               graph_color=mod_cols[unique(x$Module)],
#                               title=unique(x$Module),
#                               ...)
#             })
#             modules <- names(res)
#             modules <- modules[order(as.numeric(stringr::str_extract(modules, "\\d+")))]
#             cem@barplot_ora <- res[modules]
#             return(cem)
#           })

# plot_ora_single <- function(es, ordr_by='p.adjust', max_length=50, pv_cut=0.05,
#                             graph_color="#4169E1", title="Over Representation Analysis"){
#   
#   comsub <- function(x){
#     #split the first and last element by character
#     d_x <- strsplit(x[c(1, length(x))], "")
#     #search for the first not common element, and so, get the last matching one
#     der_com <- match(FALSE, do.call("==", d_x))-1
#     return(substr(x, 1, der_com + 1))
#   }
#   
#   es[, "GeneSet"] <- es[, "ID"]
#   
#   # limits name length
#   ovf_rows <- which(nchar(es[, "GeneSet"]) > max_length) # overflow
#   ovf_data <- es[ovf_rows, "GeneSet"]
#   test <- strtrim(ovf_data, max_length)
#   dupes <- duplicated(test) | duplicated(test, fromLast=TRUE)
#   if(sum(dupes) > 0){
#     test[dupes] <- ovf_data[dupes]
#     test[dupes] <- comsub(test[dupes])
#     max_length <- max(nchar(test))
#   }
#   
#   es[ovf_rows, "GeneSet"] <-  paste0(strtrim(test, max_length), "...")
#   es[, "GeneSet"] <- stringr::str_wrap(es[, "GeneSet"], width = 20)
#   
#   # order bars
#   lvls <- es[order(es[, ordr_by], decreasing=TRUE), "GeneSet"]
#   es[, "GeneSet"] <- factor(es[, "GeneSet"], levels=lvls)
#   
#   es[, "alpha"] <- 1
#   es[es[, ordr_by] > pv_cut, "alpha"] <- 0
#   
#   # Avoid 0's
#   es[es[, ordr_by] > 0.8, ordr_by] <- 0.8
#   my_squish <- function(...){
#     return(scales::squish(..., only.finite=FALSE))
#   }
#   
#   # plot
#   y_axis <- paste('-log10(', ordr_by, ')')
#   pl <- ggplot(es, aes_string(x="GeneSet", y=y_axis, alpha="alpha", fill=y_axis)) +
#     geom_bar(stat="identity") +
#     theme(axis.text=element_text(size=8), legend.title=element_blank()) +
#     coord_flip() +
#     scale_alpha(range=c(0.4, 1), guide="none") +
#     labs(y="-log10(adjusted p-value)", title=title, x="") +
#     geom_hline(yintercept=-log10(pv_cut), colour="grey", linetype="longdash") +
#     scale_fill_gradient(low="gray", high=graph_color, limits=c(2, 5), oob=my_squish)
#   res <- list('pl'=pl, numsig=sum(es[, ordr_by] < pv_cut, na.rm=TRUE))
#   return(res)
# }
# 
# # Adding interactions
# setGeneric('plot_interactions', function(cem, ...) {
#   standardGeneric('plot_interactions')
# })

#' @rdname plot_interactions
# setMethod('plot_interactions', signature('CEMiTool'),
#           function(cem, n=10, ...) {
#             if(length(unique(cem@module$modules)) == 0){
#               stop("No modules in CEMiTool object! Did you run find_modules()?")
#             }
#             if(length(interactions_data(cem)) == 0){
#               stop("No interactions information! Did you run interactions_data()?")
#             }
#             #cem <- get_args(cem, vars=mget(ls()))
#             mod_cols <- mod_colors(cem)
#             mod_names <- names(cem@interactions)
#             mod_names <- mod_names[which(mod_names!="Not.Correlated")]
#             hubs <- get_hubs(cem)
#             zero_ints <- character()
#             zero_ints <- lapply(names(cem@interactions), function(mod){
#               degree <- igraph::degree(cem@interactions[[mod]], normalized=FALSE)
#               if(length(degree) == 0) {
#                 zero_ints <- append(zero_ints, mod)
#               }
#             })
#             zero_ints <- unlist(zero_ints)
#             if(!is.null(zero_ints)){
#               mod_names <- mod_names[which(!(mod_names %in% zero_ints))]
#             }
#             if(length(mod_names) == 0){
#               warning("There are no interactions in the given modules. Please check interactions file.")
#               return(cem)
#             }
#             res <- lapply(mod_names, function(name){
#               plot_interaction(ig_obj=cem@interactions[[name]],
#                                n=n, color=mod_cols[name], name=name,
#                                mod_genes=module_genes(cem, name)$genes,
#                                coexp_hubs=hubs[[name]])
#             })
#             names(res) <- mod_names
#             mod_names_ordered <- mod_names[order(as.numeric(stringr::str_extract(mod_names, "\\d+")))]
#             cem@interaction_plot <- res[mod_names_ordered]
#             return(cem)
#           })

# plot_interaction <- function(ig_obj, n, color, name, mod_genes, coexp_hubs){
#   degrees <- igraph::degree(ig_obj, normalized=FALSE)
#   ig_obj <- igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
#   max_n <- min(n, length(degrees))
#   net_obj <- intergraph::asNetwork(ig_obj)
#   m <- network::as.matrix.network.adjacency(net_obj) # get sociomatrix
#   # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
#   plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(m, NULL))
#   # or get it them from Kamada-Kawai's algorithm:
#   # plotcord <- data.frame(sna::gplot.layout.kamadakawai(m, NULL))
#   colnames(plotcord) <- c("X1","X2")
#   edglist <- network::as.matrix.network.edgelist(net_obj)
#   edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
#   plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
#   plotcord$Degree <- network::get.vertex.attribute(net_obj, "degree")
#   plotcord[, "shouldLabel"] <- FALSE
#   plotcord[, "Hub"] <- ""
#   int_hubs <- names(sort(degrees, decreasing=TRUE))[1:max_n]
#   int_bool <- plotcord[, "vertex.names"] %in% int_hubs
#   plotcord[which(int_bool), "Hub"] <- "Interaction"
#   sel_vertex <- int_hubs
#   if(!missing(coexp_hubs)){
#     coexp_bool <- plotcord[, "vertex.names"] %in% coexp_hubs
#     coexp_and_int <- coexp_bool & int_bool
#     plotcord[which(coexp_bool), "Hub"] <- "Co-expression"
#     plotcord[which(coexp_and_int), "Hub"] <- "Co-expression + Interaction"
#     sel_vertex <- c(sel_vertex, coexp_hubs)
#   }
#   
#   colnames(edges) <-  c("X1","Y1","X2","Y2")
#   #edges$midX  <- (edges$X1 + edges$X2) / 2
#   #edges$midY  <- (edges$Y1 + edges$Y2) / 2
#   plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <- TRUE
#   plotcord$Degree_cut <- cut(plotcord$Degree, breaks=3, labels=FALSE)
#   plotcord$in_mod <- TRUE
#   #mod_genes <- cem@module[cem@module$modules==name,]$genes
#   not_in <- setdiff(plotcord[,"vertex.names"], mod_genes)
#   plotcord[which(plotcord[, "vertex.names"] %in% not_in), "in_mod"] <- FALSE
#   
#   pl <- ggplot(plotcord)  +
#     geom_segment(data=edges, aes_(x=~X1, y=~Y1, xend=~X2, yend=~Y2),
#                  size = 0.5, alpha=0.5, colour="#DDDDDD") +
#     geom_point(aes_(x=~X1, y=~X2, size=~Degree, alpha=~Degree), color=color) +
#     geom_label_repel(aes_(x=~X1, y=~X2, label=~vertex.names, color=~Hub),
#                      box.padding=unit(1, "lines"),
#                      data=function(x){x[x$shouldLabel, ]}) +
#     scale_colour_manual(values=c("Co-expression" = "#005E87",
#                                  "Interaction" = "#540814",
#                                  "Co-expression + Interaction" = "#736E0B")) +
#     labs(title=name) +
#     ggplot2::theme_bw(base_size = 12, base_family = "") +
#     ggplot2::theme(axis.text = ggplot2::element_blank(),
#                    axis.ticks = ggplot2::element_blank(),
#                    axis.title = ggplot2::element_blank(),
#                    legend.key = ggplot2::element_blank(),
#                    panel.background = ggplot2::element_rect(fill = "white",
#                                                             colour = NA),
#                    panel.border = ggplot2::element_blank(),
#                    panel.grid = ggplot2::element_blank())
#   
#   return(pl)
# }

# Multiplot
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots <- length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
      ncol = cols, nrow = ceiling(numPlots / cols)
    )
  }

  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(
        layout.pos.row = matchidx$row,
        layout.pos.col = matchidx$col
      ))
    }
  }
}