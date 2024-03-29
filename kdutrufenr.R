library(tidyverse)
library(CEMiTool)
library(WGCNA)
library(parallelMap)

get_attribute_field = function (x, field, attrsep = "; ") 
{
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

ora_analysis <- function(genes_list, database = KEGG_2019_Human, pvalue_cutoff = 1, qvalue_cutoff = 1, p.adjust_cutoff = 0.1, count_at_least = 3) {
  
  # Load required packages
  library(clusterProfiler)
  library(tidyverse)
  
  # Read the GMT file into a data frame
  gmt_df <- database %>% clusterProfiler::read.gmt()
  
  # Extract the gene sets (terms) from the GMT data frame
  terms_list <- gmt_df %>%
    split(f = .$term) %>%
    purrr::map(pull, gene)
  
  # Get the unique set of genes present in the database
  genes_set <- gmt_df[, "gene"] %>% unique()
  
  # Perform enrichment analysis for each gene set in the genes_list
  enrichment_result <- genes_list %>%
    purrr::map(clusterProfiler::enricher, pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff, universe = genes_set, TERM2GENE = gmt_df)
  
  # Combine the enrichment results into a single data frame
  cl_prof_df <- enrichment_result %>%
    purrr::map(as.data.frame) %>%
    bind_rows(.id = "Cluster")
  
  # Filter the results based on p.adjust and count criteria
  cl_prof_df <- cl_prof_df %>% dplyr::filter(p.adjust < p.adjust_cutoff & Count >= count_at_least)
  
  # Store the input gene clusters
  geneClusters <- genes_list
  
  # Create a result object of class "compareClusterResult"
  compareCluster_output <- new("compareClusterResult",
                               compareClusterResult = cl_prof_df,
                               geneClusters = list(geneClusters),
                               fun = "enricher",
                               .call = match.call(expand.dots = TRUE)
  )
  
  # Return the result object
  return(compareCluster_output)
}

shape_data <- function(ora_results_output) {
  ora_results_output@compareClusterResult <- ora_results_output@compareClusterResult %>%
    mutate_at(.vars = c("pvalue", "p.adjust", "qvalue"), .funs = formatC, digits = 2, format = "e")
  return(ora_results_output)
}

gsea_analysis <- function(database = KEGG_2019_Mouse, sorted_genes, p_value_cutoff = 0.05, p_adjust_cutoff = 0.05, min_size = 10, max_size = 500, p_adjust_method = "BH") {
  
  # Load required packages
  library(clusterProfiler)
  library(tidyverse)
  
  # Read the GMT file into a data frame
  gmt_df <- database %>% clusterProfiler::read.gmt()
  
  # Extract the gene sets (terms) from the GMT data frame
  gmt_list <- gmt_df %>%
    split(f = .$term) %>%
    purrr::map(pull, gene)
  
  # Perform Gene Set Enrichment Analysis (GSEA) using the sorted_genes and the database
  gsea_result <- clusterProfiler::GSEA(
    geneList = sorted_genes, 
    TERM2GENE = gmt_df, 
    verbose = FALSE, 
    minGSSize = min_size, 
    maxGSSize = max_size, 
    pvalueCutoff = p_value_cutoff, 
    pAdjustMethod = p_adjust_method
  ) %>% suppressWarnings()
  
  # Generate enrichment plots for each gene set and store them in a list
  gseaplot_results <- purrr::map(gsea_result@result$Description, function(i) {
    enrichplot::gseaplot2(
      x = gsea_result,
      geneSetID = i,
      title = gsea_result@result$Description[i]
    )
  }) %>% purrr::set_names(gsea_result@result$Description)
  
  # Store the GSEA results and the enrichment plot results in a list
  results <- list(
    gsea_results = gsea_result,
    gseaplot_results = gseaplot_results
  )
  
  # Return the results
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

calculate_TOM <- function(df){
  
  # Check input data
  stopifnot(is.data.frame(df))
  
  # Load the required packages
  library(tidyverse)
  library(WGCNA)
  
  # Scale data
  dat_expr <- df %>%
    t() %>%
    as.data.frame() %>%
    scale()
  
  # Build hierarchical clustering tree
  sample_tree <- dat_expr %>%
    dist() %>%
    hclust(method = "average")
  
  # Determine soft power threshold
  soft_threshold <- dat_expr %>% 
    pickSoftThreshold(powerVector = c(c(2:10), seq(from = 12, to = 30, by = 2)), 
                      verbose = 2, 
                      moreNetworkConcepts = TRUE)
  
  soft_power <- soft_threshold$fitIndices %>%
    filter(SFT.R.sq > 0.8) %>%
    pull(Power) %>%
    min()
  if (is.na(soft_power | is.infinite(soft_power))) {soft_power <- 10}
  
  # Calculate TOM matrix
  dat_expr <- dat_expr %>% apply(MARGIN = c(1, 2), FUN = as.numeric)
  adjacency <- dat_expr %>% adjacency(power = soft_power)
  TOM <- adjacency %>% TOMsimilarity()
  rownames(TOM) <- colnames(TOM) <- df %>% rownames()
  
  return(TOM)
}

# TOM <- WGCNA_tpm_df %>% calculate_TOM()
# 
# minimum_cluster_size <- 1:50
# minimum_cluster_depth <- 4

calculate_ssd <- function(x) sum((x - mean(x))^2)

get_cluster_summary <- function(minimum_cluster_size, TOM, minimum_cluster_depth = 4, number_of_cores = 6, ...) {
  
  # Check input data
  stopifnot(is.numeric(minimum_cluster_size))
  stopifnot(is.numeric(minimum_cluster_depth))
  stopifnot(is.numeric(number_of_cores))
  stopifnot(is.numeric(TOM))
  
  # Load the required packages
  library(dynamicTreeCut)
  library(foreach)
  library(ggeasy)
  
  # Get gene names
  gene_names <- TOM %>% rownames()
  
  # Initiate cluster
  doParallel::registerDoParallel(cores = number_of_cores)
  
  # dissimilarity TOM
  diss_TOM = 1 - TOM
  
  # estimate gene_tree
  gene_tree <- diss_TOM %>% as.dist() %>% hclust(method = "average")
  
  # Compute dynamic cluster count
  dynamic_cluster_list <- foreach(i = minimum_cluster_size) %dopar% 
    dynamicTreeCut::cutreeDynamic(dendro = gene_tree, distM = diss_TOM, deepSplit = minimum_cluster_depth, pamRespectsDendro = FALSE, minClusterSize = i)
  
  # Compute dynamic cluster count
  dynamic_cluster_count <- dynamic_cluster_list %>% 
    purrr::map(max) %>% 
    unlist()
  
  # Compute module membership colors
  module_membership_colors <- dynamic_cluster_list %>% purrr::map(labels2colors)
  
  # Compute module membership colors data frame
  module_membership_colors_df <- foreach(i = seq_along(dynamic_cluster_list)) %dopar% 
    data.frame(dynamic_colors = module_membership_colors[[i]], genes = gene_names)
  
  pre_cluster_blocks <- purrr::map(seq_along(module_membership_colors_df), function(i) {
    module_membership_colors_df[[i]] %>% 
      split(f = .$dynamic_colors) %>% 
      purrr::map(dplyr::select, -dynamic_colors) %>% 
      purrr::map(pull)
  }
  )
  
  in_cluster_membership <- purrr::map(seq_along(pre_cluster_blocks), function(i) purrr::map(seq_along(pre_cluster_blocks[[i]]), function(j) gene_names %in% pre_cluster_blocks[[i]][[j]]))
  
  # Select the corresponding Topological Overlap
  cluster_TOM <- purrr::map(seq_along(in_cluster_membership), function(i) purrr::map(seq_along(in_cluster_membership[[i]]), function(j) TOM[in_cluster_membership[[i]][[j]], in_cluster_membership[[i]][[j]]]))
  
  sum_of_squared_differences <- purrr::map(seq_along(cluster_TOM), function(i) purrr::map(seq_along(cluster_TOM[[i]]), function(j) cluster_TOM[[i]][[j]] %>% calculate_ssd())) %>% purrr::map(unlist) %>% purrr::map(sum) %>% unlist()
  
  cluster_summary <- data.frame(minimum_n_genes = minimum_cluster_size, 
                                n_modules = dynamic_cluster_count, 
                                sum_of_squared_differences = sum_of_squared_differences)
  parallelStop()
  
  return(cluster_summary)
}

# get_elbow_points_indices <- function(x, y, threshold) {
#     first_derivative <- diff(y) / diff(x) 
#     second_derivative <- diff(first_derivative) / diff(x[-1]) 
#     indices <- which(abs(second_derivative) > threshold)
#     return(indices)
#   }

# https://stackoverflow.com/questions/41518870/finding-the-elbow-knee-in-a-curve
# https://stackoverflow.com/questions/14351608/color-one-point-and-add-an-annotation-in-ggplot2
determine_elbow_point <- function(cluster_summary_output, huge_jump_threshold = 1e3) {
  
  # Extract n_modules and sum_of_squared_differences from cluster_summary_output$table
  n_modules <- cluster_summary_output$n_modules
  ssd <- cluster_summary_output$sum_of_squared_differences
  
  # Approximate the points using linear interpolation
  approximated_points <- approx(n_modules, ssd, n = 1000, yleft = min(ssd), yright = max(ssd))
  approx_n_modules <- approximated_points$x
  approx_ssd <- approximated_points$y
  
  # Find the index of the elbow point (i.e., where the SSD curve starts to level off)
  # threshold for huge jump = 1e4
  # elbow_indices <- get_elbow_points_indices(approx_n_modules, approx_ssd, huge_jump_threshold) 
  # elbow_module_index <- approx_n_modules[max(elbow_indices)] %>% ceiling()
  first_derivative <- diff(approx_ssd) / diff(approx_n_modules) 
  second_derivative <- diff(first_derivative) / diff(approx_n_modules[-1]) 
  elbow_index <- which(abs(second_derivative) == max(abs(second_derivative)))
  elbow_module_index <- approx_n_modules[elbow_index] %>% ceiling()
  
  # Find the minimum module size corresponding to the elbow point
  min_module_size <- cluster_summary_output %>%
    dplyr::filter(n_modules == elbow_module_index) %>%
    pull(minimum_n_genes)
  
  # Set the "elbow" column in cluster_summary_output
  cluster_summary_output <- cluster_summary_output %>% 
    mutate(elbow = n_modules == elbow_module_index)
  
  elbow_point <- cluster_summary_output %>% 
    dplyr::filter(elbow == TRUE)
  
  # Generate a plot showing the SSD curve and the elbow point
  elbow_plot <- cluster_summary_output %>% 
    ggplot(aes(x = n_modules, y = sum_of_squared_differences)) + 
    geom_line() + 
    geom_point() + 
    geom_point(data = elbow_point, color = "red") + 
    theme_bw() +
    labs(x = "Number of modules", y = "Sum of squared differences") +
    ggeasy::easy_text_size(size = 20)
  
  # Return a list containing the minimum module size and the elbow plot
  return(list(min_module_size = min_module_size,
              elbow_plot = elbow_plot))
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
