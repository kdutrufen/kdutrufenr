library(tidyverse)
library(CEMiTool)
library(WGCNA)
library(parallelMap)

extract_edges <- function(dend) {
  if(is.leaf(dend)) {
    # Create a 'singleton' edge for each leaf
    return(data.frame(
      from = attr(dend, "label"),    # Use the leaf label directly
      to = attr(dend, "label"),
      height = 0                     # Height is 0 for leaf nodes
    ))
  } else {
    left_edges <- extract_edges(dend[[1]])
    right_edges <- extract_edges(dend[[2]])

    # Create a unique ID for the internal node 
    node_id <- paste0("node_", attr(dend, "members") + 1)

    new_edges <- data.frame(
      from = node_id,
      to = c(
        if (is.leaf(dend[[1]])) attr(dend[[1]], "label") else paste0("node_", attr(dend[[1]], "members") + 1),
        if (is.leaf(dend[[2]])) attr(dend[[2]], "label") else paste0("node_", attr(dend[[2]], "members") + 1)
      ),
      height = attr(dend, "height")
    )

    return(rbind(left_edges, right_edges, new_edges))
  }
}


create_sce_from_seurat <- function(seurat_object){
  
  # Extract counts from the Seurat object, using the 'raw' assay
  counts <- GetAssayData(seurat_object, slot = "counts", assay = "RNA")

  # Create the SingleCellExperiment object
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = seurat_object@meta.data 
    )

  # Add rowData if available in the Seurat object
  if (!is.null(seurat_object@assays$RNA@meta.features)) {
    rowData(sce) <- seurat_object@assays$RNA@meta.features
  }
  
  # Extract the UMAP coordinates from the Seurat object
  umap_coords <- seurat_object@reductions$umap@cell.embeddings %>% as.data.frame() %>% as.matrix()
  # umap_coords <- final_meta_dataset@reductions$integrated_UMAP
  
  # Add the UMAP coordinates to the reducedDims slot of the sce object
  reducedDims(sce) <- SimpleList(UMAP = umap_coords)
  
  return(sce)
}

include_de_results <- function(sce, top_threshold) {
  markers <- sce %>% findMarkers(test = "binom", direction = "up", lfc = 1, pval.type = "any", groups = .$cell_type)
  markers_list <- purrr::map(seq_along(markers), function(i) markers[[i]])
  markers_list <- markers_list %>% purrr::set_names(sce$cell_type %>% levels())
  de_list <- purrr::map(seq_along(markers), function(i) {
    markers_list[[i]] %>%
      as.data.frame() %>%
      dplyr::filter(FDR <= 0.05 & summary.logFC >= 1.0 & Top <= top_threshold)
  }) %>%
    purrr::set_names(sce$cell_type %>% levels())

  de_list <- de_list %>%
    purrr::map(mutate_at, .vars = c("p.value", "FDR"), .funs = formatC, digits = 2, format = "e") %>%
    purrr::map(mutate_if, .predicate = is.numeric, .funs = round, digits = 2)

  de_list <- de_list %>%
    purrr::map(rownames_to_column, var = "gene")

  # Add as a new assay
  metadata(sce)$de_results <- de_list
  
  return(sce)
}
                             
create_imputation_spaghetti_plot <- function(long_df, column_name, title) {
  # Check if column_name exists in the data frame
  if (!(column_name %in% colnames(long_df))) {
    stop(paste("Error: Column", column_name, "not found in the data frame."))
  }

  long_df %>%
    arrange(imputed_data) %>%
    ggplot(aes(x = time, y = !!sym(column_name), color = imputed_data, group = wound_id)) +
    geom_point() +
    scale_color_manual(values = c("yes" = "purple", "no" = "gray")) +
    geom_line(alpha = 0.4) +
    theme_bw() +
    facet_wrap(~treatment) +
    scale_x_discrete(limits = c(0, 2, 5, 10)) +
    ggeasy::easy_text_size(size = 20) +
    ggeasy::easy_x_axis_labels_size(size = 15) +
    ggeasy::easy_y_axis_labels_size(size = 15) +
    labs(x = "Time Post-Bite (Day)", y = bquote("Necrotic Area "(cm^2)), color = "Imputed data", title = title) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title.position = "top"))
}

is_mar <- function(model) {
  # Calculate the p-values from the z-values for each predictor variable
  coefficients <- summary(model)$coefficients[, "Estimate"]
  standard_errors <- summary(model)$coefficients[, "Std. Error"]
  z_values <- coefficients / standard_errors
  p_values <- 2 * (1 - pnorm(abs(z_values)))

  # Compare the p-values with the significance level (e.g., 0.05)
  significance_level <- 0.05
  is_mar <- all(p_values > significance_level)

  # Print the p-values and conclusion
  print("P-Values:")
  print(p_values)

  print("Conclusion:")
  if (is_mar) {
    print("The data is Missing at Random (MAR).")
  } else {
    print("The data is not Missing at Random (not MAR).")
  }
}

create_chord_diagram <- function(stringDB_data, title_text, combined_score_threshold, display_gene_names = F) {
  # mar = c(bottom, left, top, right)
  # par(cex = 1.2, mar = c(4.1, 4.1, 4.1, 4.1))
  circos.clear()
  
  par(cex = 1.2, mar = c(1, 4, 2, 4)) # Smaller margins
  
  circos.par(circle.margin = 0.5)

  stringDB_data %>%
    dplyr::filter(combined_score >= combined_score_threshold) %>%
    dplyr::select(c("#node1", "node2", "combined_score", "color")) %>%
    chordDiagram(
      annotationTrack = "grid",
      grid.col = "lightgrey",
      col = .$color
    )

  # Add the title
  title(title_text, cex.main = 2, font.main = 1) # Adjust cex.main and font.main as needed

  if (display_gene_names) {
    circos.track(track.index = 1, track.height = 2, panel.fun = function(x, y) {
      circos.text(
        CELL_META$xcenter,
        CELL_META$ylim[2],
        CELL_META$sector.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1
      )
    }, bg.border = NA) # here set bg.border to NA is important
  }

  chord_diagram_plot <- recordPlot()
  # circos.clear()
}

ggplot_RLE <- function(data_matrix, design_df, title, y_lim = c(-2, 2), pseudo_count = 1e-9, ...) {
  # Check for missing values
  if (anyNA(data_matrix)) {
    warning("Removing rows with missing values before calculating RLE.")
    data_matrix <- na.omit(data_matrix)
  }
  
  # 1. Calculate median (reference) values 
  ref_values <- apply(data_matrix, 1, median)
  # ref_values <- apply(data_matrix, 2, median)
  
  # 2. Compute log ratios with pseudocount
  pseudocount <- pseudo_count # Choose a small value appropriate for your data
  log_ratio_matrix <- log2(sweep(data_matrix + pseudocount, 1, ref_values + pseudocount, "/"))
  
  # 3. Create data frame for plotting
  log_ratio_pivoted_df <- log_ratio_matrix %>%
    as.data.frame() %>%
    pivot_longer(cols = 1:ncol(.), names_to = "sample", values_to = "log_ratio") %>% 
    left_join(design_df)
  
  # 4. Create the plot
  log_ratio_pivoted_df %>%
    ggplot(aes(x = sample, y = log_ratio, fill = condition)) +
    stat_boxplot(
      geom ='errorbar', 
      linetype = "dashed",
      # linetype = "dotted", 
      width = 0.5
    ) + 
    geom_boxplot(
      outlier.shape = NA,
      fatten = 2,               # Increase fatten to have a larger median line in the box plot
      coef = 0
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") + # Add horizontal line at y = 0
    ylim(y_lim) +
    theme(
      panel.grid = element_blank(), # Remove gridlines
      panel.background = element_blank(), # Set panel background to blank
      panel.border = element_rect(colour = "black", fill = NA, size = 1) # Add this line
    ) +
    theme(
      axis.text.x = element_blank(), # Remove x-axis text labels
      axis.ticks.x = element_blank()
    ) +
    labs(x = "Sample", y = "Relative Log Expression (RLE)", fill = NULL, title = title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggeasy::easy_all_text_size(size = 15) +
    scale_fill_manual(values = unname(polychrome())[1:length(unique(design_df[["condition"]]))]) +
    theme(legend.position = "bottom")
}

create_volcano_plot <- function(df, dot_size = 0.5){
  
  # Check for required columns
  required_cols <- c("logFC", "FDR", "DEG_status")
  if (!all(required_cols %in% colnames(df))) {
    stop(paste("Missing required columns:", paste(setdiff(required_cols, colnames(df)), collapse = ", ")))
  }

  # Enhanced input validation    
  valid_levels <- c("Down", "Not_DEG", "Up")
  if (!all(unique(df$DEG_status) %in% valid_levels)) {
    stop(paste("Invalid `DEG_status` levels. Allowed levels are:", paste(valid_levels, collapse = ", ")))
  }
  
  max_LFC <- df$logFC %>% as.numeric() %>% range() %>% abs() %>% max() %>% ceiling()
  max_FDR <- df$FDR %>% as.numeric() %>% log10() %>% '*'(-1) %>% max() %>% ceiling()
  
  volcano_plot <- df %>%
    ggplot(aes(x = as.numeric(logFC), y = -log10(as.numeric(FDR)), colour = DEG_status)) + 
    geom_point(alpha = 0.4, size = dot_size) +
    theme_bw() +
    xlim(c(-max_LFC, max_LFC)) + 
    ylim(c(0, max_FDR)) +    
    scale_color_manual(values = c("Down" = "#094FED", "Not_DEG" = "black", "Up" = "red"),
                       labels = c("Downregulated DEG", "Not significant", "Upregulated DEG")) +
    # expand_limits(x = 0, y = 0) +  # Ensure axes cross at (0,0)
    labs(x = expression(paste(log[2], " Fold Change")), y = expression(paste("\u2212", log[10], " adjusted p-value"))) +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 5))) + # Change size of dots in legend
    theme(legend.position = "bottom") +
    easy_all_text_size(size = 20) # +
    # scale_x_continuous(labels = function(x) {ifelse(x < 0, paste("\u2212", abs(x)), x)})
  
  return(volcano_plot)
}

create_read_count_barplot <- function(read_count_df, design_df) {
  # Load the required libraries for color palettes and plot customization
  library(pals)
  library(Polychrome)
  
  # Load the required libraries for color palettes and plot customization
  sum_data <- data.frame(
    read_count = c(mean(colSums(read_count_df)), apply(read_count_df, 2, sum)),
    sample = c("Average", colnames(read_count_df)),
    sample_type = c("Average", design_df$condition)
  )
  
  # Generate a color palette based on the unique sample types (conditions)
  col_cell <- unname(polychrome())[design_df$condition %>% as.factor()]
  
  # Create the bar plot
  read_count_barplot <- sum_data %>%
    ggplot(aes(x = sample, y = read_count, fill = sample_type)) +                          # Map sample to x, read_count to y, and sample_type to fill
    geom_bar(colour = "black", stat = "identity") +                                        # Add bars with black borders
    geom_hline(yintercept = 8e+06, colour = "black") +                                     # Horizontal reference line
    theme_bw() +                                                                           # Classic black & white theme
    labs(x = "Samples", y = "Read count", fill = NULL) +                                   # Label axes, remove default fill legend title
    # ggeasy::easy_x_axis_labels_size(size = 10) +                                         # Option to adjust x-axis label size
    ggeasy::easy_y_axis_labels_size(size = 10) +                                           # Adjust y-axis label size
    ggeasy::easy_x_axis_title_size(size = 20) +                                            # Adjust x-axis title size
    ggeasy::easy_y_axis_title_size(size = 20) +                                            # Adjust y-axis title size
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +                           # Option to rotate x-axis labels
    theme(axis.text.x = element_blank(),                                                   # Remove x-axis text labels
          axis.ticks.x = element_blank()) +                                                # Remove x-axis tick marks
    scale_y_continuous(labels = function(x) {                                              # Format y-axis labels
      ifelse(x == 0, "0", paste0(format(x / 1e7, scientific = FALSE), " \u00D7 10\u2077")) # Use Unicode for superscript
    }) +
    scale_fill_manual(values = c("yellow", col_cell), breaks = sum_data$sample_type) +     # Set colors manually
    theme(legend.position = "bottom")                                                      # Place the legend at the bottom
  
  return(read_count_barplot)
}

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
    purrr::map(function(x) x %>% .[["gene"]])
  
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

gsea_analysis <- function(database = KEGG_2019_Mouse, sorted_genes, p_value_cutoff = 0.05, p_adjust_cutoff = 0.05, min_size = 10, max_size = 500, p_adjust_method = "BH", seed = 123) {
  
  # Load required packages
  library(clusterProfiler)
  library(tidyverse)
  
  # Read the GMT file into a data frame
  gmt_df <- database %>% clusterProfiler::read.gmt()
  
  # Extract the gene sets (terms) from the GMT data frame
  gmt_list <- gmt_df %>%
    split(f = .$term) %>%
    purrr::map("gene")
  
  # Perform Gene Set Enrichment Analysis (GSEA) using the sorted_genes and the database
  set.seed(seed)
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
