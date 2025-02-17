library(tidyverse)
library(CEMiTool)
library(WGCNA)
library(parallelMap)

create_gene_tpm_trajectory_dataframe <- function(sce, gene_v) {
  # Input validation
  if (!is(sce, "SingleCellExperiment")) {
    stop("Input 'sce' must be a SingleCellExperiment object.")
  }
  if (!is.character(gene_v) || !all(nzchar(gene_v))) {
    stop("Input 'gene_v' must be a character vector of gene names.")
  }
  if (is.null(sce$pseudo_paths)) {
    stop("Input 'sce' does not contain a 'pseudo_paths' slot.")
  }
  if (!all(gene_v %in% rownames(sce))) {
    missing_genes <- gene_v[!gene_v %in% rownames(sce)]
    warning(paste("The following genes are not found in 'sce':", paste(missing_genes, collapse = ", ")))
    gene_v <- gene_v[gene_v %in% rownames(sce)]
  }
  
  # gene_v <- gene_v %>%
  #   intersect(row.names(sce)) %>%
  #   sort() %>%
  #   unique()

  # Calculate TPM values
  # sce@assays@data$tpm <- log2(cpm(sce@assays@data$counts) / 10 + 1) %>% as(Class = "dgCMatrix")
  sce@assays@data$tpm <- (sce@assays@data$counts %>% edgeR::cpm() %>% "/"(10)) %>% "+"(1) %>% log2() %>% as(Class = "dgCMatrix")

  # Extract and pre-process pseudotime data
  pseudo_paths <- sce$pseudo_paths
  n_lineages <- pseudo_paths %>% ncol()
  average_pseudo_time <- pseudo_paths %>% rowMeans(na.rm = TRUE)

  # Create lineage specific pseudotime list
  lineage_pseudo_time_list <- purrr::map(seq_len(n_lineages), function(i) pseudo_paths[, i]) %>% purrr::set_names(paste0("lineage_", 1:n_lineages))

  # Long format TPM data function (helper)
  to_long_tpm <- function(tpm_matrix, pseudo_time) {
    tpm_matrix[rownames(tpm_matrix) %in% gene_v, ] %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var = "cell") %>%
      mutate(pseudo_time = pseudo_time) %>%
      drop_na(pseudo_time) %>%
      arrange(pseudo_time) %>%
      dplyr::select(-pseudo_time) %>%
      column_to_rownames(var = "cell") %>%
      t() %>%
      as.data.frame()
  }

  # Create lineage specific pseudotime dataframes
  lineage_pseudo_list <- lineage_pseudo_time_list %>% purrr::map(to_long_tpm, tpm_matrix = sce@assays@data$tpm)

  lineage_pseudo_list$average_pseudo <- to_long_tpm(sce@assays@data$tpm, average_pseudo_time)

  return(lineage_pseudo_list)
}


plot_dendro_heatmap_for_paper <- function(df) {
  library(ggdendro)
  library(tidyverse)
  
  cell_names <- names(df)
  
  # Calculate the dendrogram
  dend <- df %>%
      as.matrix() %>%
    dist() %>%
    hclust() %>%
    as.dendrogram()
  
   dend_data <- dend %>% dendro_data()

  # Setup the data, so that the layout is inverted (this is more
  # "clear" than simply using coord_flip())
  segment_data <- with(
    segment(dend_data),
    data.frame(x = y, y = x, xend = yend, yend = xend)
  )

  # Use the dendrogram label data to position the gene labels
  gene_pos_table <- with(
    dend_data$labels,
    data.frame(y_center = x, gene = as.character(label), height = 1)
  )

  # Table to position the samples
  cell_pos_table <- data.frame(cell = cell_names) %>%
    mutate(x_center = (1:nrow(.)), width = 1)
  
  # Neglecting the gap parameters
  heatmap_data <- df %>%
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cell") %>% 
    pivot_longer(cols = names(.)[2]:names(.)[ncol(.)], names_to = "gene", values_to = "tpm") %>% 
    left_join(gene_pos_table) %>%
    left_join(cell_pos_table)

  # Limits for the vertical axes
  gene_axis_limits <- with(
    gene_pos_table,
    c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
  ) +
    0.1 * c(-1, 1) # extra spacing: 0.1
  
  max_expr <- ceiling(max(heatmap_data$tpm))
  
  # Heatmap plot
  plt_hmap <- heatmap_data %>% 
    ggplot(aes(x = x_center, y = y_center, fill = tpm, height = height, width = width)) +
    geom_tile() +
    # scale_fill_gradient2("Expr", high = "red", low = "blue") +
    scale_fill_gradientn(name = "log2(TPM/10 + 1)",
                         limits = c(0, max_expr),
                         colours = colorRampPalette(colors = c("gainsboro", "floralwhite", "yellow", "orange", "red", "darkred"))(50)) +
                         # colours = colorRampPalette(colors = c("blue", "gray", "yellow", "orange", "red", "darkred"))(50)) +
                         # colours = colorRampPalette(colors = c("steelblue1", "yellow", "orangered"))(50)) +
                         # colours = rev(colorRampPalette(colors = RColorBrewer::brewer.pal(n = 5, name = "RdYlBu"))(50))) +
    scale_x_continuous(
      breaks = cell_pos_table$x_center,
      labels = cell_pos_table$cell,
      expand = c(0, 0)
    ) +
    # For the y axis, alternatively set the labels as: gene_position_table$gene
    scale_y_continuous(
      breaks = gene_pos_table[, "y_center"],
      # labels = rep("", nrow(gene_pos_table)),
      labels = heatmap_data %>% arrange(y_center) %>% .[["gene"]] %>% unique(),
      limits = gene_axis_limits,
      expand = c(0, 0),
      position = "right"
    ) +
    labs(x = "", y = "") +
    theme_bw() +
    theme(legend.position = "bottom") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.y = element_blank(),
      # margin: top, right, bottom, and left
      plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm")
      # panel.grid.minor = element_blank()
    ) +
  guides(fill = guide_colorbar(title.position = "top")) # Move title to the top of the legend

  # Dendrogram plot
  plt_dendr <- segment_data %>% 
    ggplot() +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_reverse(expand = c(0, 0.5)) +
    scale_y_continuous(
      breaks = gene_pos_table$y_center,
      # labels = gene_pos_table$Gene,
      labels = NULL,
      limits = gene_axis_limits,
      expand = c(0, 0)
    ) +
    labs(x = "Distance", y = "", colour = "", size = "") +
    theme(panel.grid.minor = element_blank()) +
    theme_void() 
  
  # plot_grid(plt_dendr, plt_hmap, align = "h", rel_widths = c(0.4, 1))
  figure_list <- list(plt_dendr = plt_dendr,
                   plt_hmap = plt_hmap,
                  heatmap_cell_order = cell_pos_table) # Add the cell order to the returned list

}

plot_lineage_trajectories <- function(sce, lineage_column = "Lineage1", color_by = "cell_type", cell_order = NULL) {
  sc_info_df <- data.frame(
    cell_type = sce$cell_type,
    cell_type_color = sce$cell_type_colors,
    batch = as.character(sce$batch),
    cell = colnames(sce)
  ) %>%
    cbind(sce$pseudo_paths)

  # Extract relevant data
  pseudo_paths_df <- sc_info_df %>%
    rowwise() %>%
    mutate(average = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    dplyr::select(cell, lineage_column, cell_type, cell_type_color, batch) %>%
    drop_na(!!sym(lineage_column)) %>%
    arrange(!!sym(lineage_column))

  pseudo_paths_df <- pseudo_paths_df %>%
    mutate(x_center = 1:nrow(.))
  
  # *** KEY CHANGE: Name the color vector ***
  names(pseudo_paths_df$cell_type_color) <- pseudo_paths_df$cell_type  # Names are cell types

    # *** KEY CHANGE: Apply cell order if provided ***
  if (!is.null(cell_order)) {
    pseudo_paths_df <- pseudo_paths_df %>%
      mutate(original_x_center = x_center) %>% # Keep track of original order
      arrange(factor(cell, levels = cell_order)) %>%  # Reorder cells
      mutate(x_center = 1:nrow(.)) # Update x_center
  }
  
  # Create plot
  plot <- pseudo_paths_df %>%
    ggplot(aes(x = x_center, y = 1, fill = !!sym(color_by))) +
    geom_tile() +
    theme_void() +
    scale_x_continuous(
      breaks = pseudo_paths_df$x_center,
      expand = c(0, 0)) + # Set limits explicitly
    # labs(x = NULL, y = NULL, fill = color_by) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))

  # Apply color palette based on color_by variable
  if (color_by == "cell_type") {
    plot <- plot +
      scale_fill_manual(values = pseudo_paths_df$cell_type_color) +
      labs(x = NULL, y = NULL, fill = "Cell types") +
      guides(fill = guide_legend(
        override.aes = list(size = 3), # Adjust legend key size
        nrow = 2,
        label.theme = element_text(margin = margin(r = 20)), # Add right margin to labels
        label.position = "right",
        title.position = "top", # Position title at the top
        title.hjust = 0.5,
        keywidth = unit(0.8, "cm")
      ))
  } else if (color_by == "batch") {
    # Add logic for batch coloring if needed
    plot <- plot +
      scale_fill_manual(values = RColorBrewer::brewer.pal(n = length(unique(pseudo_paths_df$batch)), name = "Set1")) +
      labs(x = NULL, y = NULL, fill = "Embryonic stage")
  } else if (color_by == "pseudotime") {
    # Add logic for batch coloring if needed
    plot <- pseudo_paths_df %>%
      ggplot(aes(x = x_center, y = 1, fill = as.numeric(!!sym(lineage_column)))) +
      geom_tile() +
      theme_void() +
      scale_x_continuous(
      breaks = pseudo_paths_df$x_center,
      expand = c(0, 0)) + # Set limits explicitly
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
      scale_fill_viridis(option = "viridis") +
      labs(x = NULL, y = NULL, fill = paste("Pseudotime", lineage_column)) +
      guides(fill = guide_colorbar(title.position = "top", label.position = "bottom"))  # Set legend title and position
  }
}

plot_lineage_trajectory_figure <- function(sce, lineage_column = "Lineage1", genes_of_interest = NULL) {
  
  new_lineage_column <- lineage_column %>% str_to_lower() %>% str_replace(pattern = "ge", replacement = "ge_")
  
  trajectory_df_list <- create_gene_tpm_trajectory_dataframe(sce = sce, gene_v = genes_of_interest)
  
  # 1. Prepare Heatmap and Dendrogram
  if (new_lineage_column %in% names(trajectory_df_list)) {
    dendro_heatmap_data <- plot_dendro_heatmap_for_paper(trajectory_df_list[[new_lineage_column]])
    dendro <- dendro_heatmap_data$plt_dendr
    heatmap <- dendro_heatmap_data$plt_hmap + theme(legend.position = "none") # Remove heatmap legend here

    # Get the cell order from the heat map
    heatmap_cell_order <- dendro_heatmap_data$heatmap_cell_order$cell
    
    # Get x-axis limits from heatmap
    x_limits <- layer_scales(dendro_heatmap_data$plt_hmap)$x$range$range
    
    # Create trajectory plots
    cell_type_plot <- plot_lineage_trajectories(sce = sce, 
                                                lineage_column = lineage_column, 
                                                color_by = "cell_type", 
                                                cell_order = heatmap_cell_order)
    
  batch_plot <- plot_lineage_trajectories(sce = sce, 
                                          lineage_column = lineage_column, 
                                          color_by = "batch", 
                                          cell_order = heatmap_cell_order)
  
  sling_plot <- plot_lineage_trajectories(sce = sce, 
                                          lineage_column = lineage_column, 
                                          color_by = "pseudotime", 
                                          cell_order = heatmap_cell_order)

    # Extract legends
    dendro_legend <- get_legend(dendro_heatmap_data$plt_hmap)
    cell_type_legend <- get_legend(cell_type_plot)
    batch_legend <- get_legend(batch_plot)
    sling_legend <- get_legend(sling_plot)

     # Remove legends from individual plots
    cell_type_plot <- cell_type_plot + theme(legend.position = "none")
    batch_plot <- batch_plot + theme(legend.position = "none")
    sling_plot <- sling_plot + theme(legend.position = "none")

    # Create combined figure
    figure <- heatmap %>%
      insert_top(plot = sling_plot, height = 1/20) %>%
      insert_top(plot = batch_plot, height = 1/20) %>%
      insert_top(plot = cell_type_plot, height = 1/20) %>%
      insert_left(plot = dendro, width = 1/10)

    figure_legend <- ggdraw(cell_type_legend) /
      (ggdraw(dendro_legend) | ggdraw(batch_legend) | ggdraw(sling_legend))

    return(as.ggplot(figure) / figure_legend)
  } else {
    stop(paste0(lineage_column, "' not found in the data."))
  }
} 
                                         
chemo_div_plot <- function(
    comp_dis_mat = NULL, div_data = NULL, div_prof_data = NULL,
    samp_dis_mat = NULL, group_data = NULL) {
  
  library(ggplot2)
  library(patchwork)
  
  all_plots <- list()
  if (is.data.frame(group_data)) {
    group_data <- as.vector(group_data[, 1])
  }
  if (!is.null(comp_dis_mat)) {
    comp_dis_mat_clust <- stats::hclust(stats::as.dist(comp_dis_mat),
      method = "average"
    )
    comp_dis_mat_clust_dend <- stats::as.dendrogram(comp_dis_mat_clust)
    comp_dis_mat_clust_dend_data <- ggdendro::dendro_data(comp_dis_mat_clust_dend)
    
    comp_dis_mat_tree_plot <- ggplot() +
      geom_segment(
        data = comp_dis_mat_clust_dend_data$segments,
        aes(x = x, y = y, xend = xend, yend = yend)
      ) +
      theme_bw() +
      geom_text(
        data = comp_dis_mat_clust_dend_data$labels,
        aes(x = x, y = y, label = label),
        hjust = -0.1, angle = 0,
        size = 3  # Add this line to decrease text size
      ) +
      scale_y_reverse(limits = c(1, -0.5), breaks = c(1, 0.75, 0.5, 0.25, 0)) +
      ylab("Dissimilarity") +
      ggtitle("") +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5)
      ) +
      coord_flip()
    
    all_plots[["comp_dis_mat_tree_plot"]] <- comp_dis_mat_tree_plot
  }
  if (!is.null(div_data)) {
    div_data_df <- as.data.frame(div_data)
    if (is.null(group_data)) {
      message("No grouping data provided.")
      group_data <- rep("NoGroup", nrow(div_data_df))
    }
    
    div_data_df$group <- group_data
    
    for (i in 1:(ncol(div_data_df) - 1)) {
      all_plots[[paste0("div_plot", colnames(div_data_df)[i])]] <- local({
        i <- i
        current_col <- colnames(div_data_df)[i]
        div_plot <- div_data_df %>% 
          ggplot(aes(x = group, y = .data[[current_col]], fill = group)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(height = 0, width = 0.1, shape = 21) +
          theme_bw() +
          labs(x = "", y = "Functional Hill Diversity", fill = "Season") +
          theme(
            text = element_text(size = 15),
            legend.position = "none"
          )
      })
    }
  }
  if (!is.null(div_prof_data)) {
    if (is.null(group_data)) {
      message("No grouping data provided.")
      group_data <- rep("NoGroup", nrow(div_prof_data$divProf))
    }
    div_prof <- div_prof_data$divProf
    div_prof_mean1 <- stats::aggregate(div_prof,
      by = list(group = group_data),
      mean, na.rm = TRUE
    )
    div_prof_mean2 <- as.data.frame(t(div_prof_mean1[, 2:ncol(div_prof_mean1)]))
    colnames(div_prof_mean2) <- div_prof_mean1$group
    q_all <- seq(
      from = div_prof_data$qMin, to = div_prof_data$qMax,
      by = div_prof_data$step
    )
    div_prof_mean2$q <- q_all
    div_hill_long <- tidyr::pivot_longer(div_prof_mean2, 1:(ncol(div_prof_mean2) -
      1), names_to = "group", values_to = "diversity")
    div_prof_ind <- as.data.frame(t(div_prof))
    div_prof_ind$q <- q_all
    div_hill_long_ind <- tidyr::pivot_longer(div_prof_ind, 1:(ncol(div_prof_ind) -
      1), names_to = "individual", values_to = "diversity")
    div_hill_long_ind$group <- rep(group_data, length(unique(div_prof_ind$q)))
    div_prof_plot <- ggplot() +
      geom_line(
        data = div_hill_long_ind,
        aes(x = q, y = diversity, group = individual, color = group), 
        linewidth = 0.5, alpha = 0.15
      ) +
      theme_bw() +
      geom_line(data = div_hill_long, 
                aes(x = q, y = diversity, color = group), 
                linewidth = 2
                ) +
      labs(x = "Diversity order (q)", y = "Functional Hill Diversity", color = "Season") +
      theme(text = element_text(size = 15))
    all_plots[["div_prof_plot"]] <- div_prof_plot
  }
  if (!is.null(samp_dis_mat)) {
    if (is.null(group_data)) {
      message("No grouping data provided.")
      if (is.matrix(samp_dis_mat)) {
        group_data <- rep("NoGroup", nrow(samp_dis_mat))
      } else {
        group_data <- rep("NoGroup", nrow(samp_dis_mat[[1]]))
      }
    }
    if (is.matrix(samp_dis_mat)) {
      utils::capture.output(nmds <- vegan::metaMDS(samp_dis_mat,
        autotransform = FALSE
      ))
      nmds_coords <- as.data.frame(nmds$points)
      nmds_coords$group <- group_data
      nmds_plot <- nmds_coords %>% 
        ggplot(aes(x = MDS1, y = MDS2, color = group)) +
        theme_bw() +
        geom_point(size = 4, alpha = 0.5) +
        labs(color = "Season") +
        theme(text = element_text(size = 15))
      all_plots[["nmds_plot"]] <- nmds_plot
    } else {
      utils::capture.output(bc_nmds <- vegan::metaMDS(samp_dis_mat$bray_curtis,
        autotransform = FALSE
      ))
      bc_nmds_coords <- as.data.frame(bc_nmds$points)
      bc_nmds_coords$group <- group_data
      bc_nmds_plot <- bc_nmds_coords %>% 
        ggplot(aes(x = MDS1, y = MDS2, color = group)) +
        theme_bw() +
        geom_point(size = 4, alpha = 0.5) +
        theme(text = element_text(size = 15)) +
        ggtitle("Bray-Curtis NMDS")
      all_plots[["bc_nmds_plot"]] <- bc_nmds_plot
      utils::capture.output(gu_nmds <- vegan::metaMDS(samp_dis_mat$gen_uni_frac,
        autotransform = FALSE
      ))
      gu_nmds_coords <- as.data.frame(gu_nmds$points)
      gu_nmds_coords$group <- group_data
      gu_nmds_plot <- gu_nmds_coords %>% 
        ggplot(aes(x = MDS1, y = MDS2, color = group)) +
        theme_bw() +
        geom_point(size = 4, alpha = 0.5) +
        theme(text = element_text(size = 15)) +
        ggtitle("Generalized UniFrac NMDS")
      all_plots[["gu_nmds_plot"]] <- gu_nmds_plot
    }
  }
  
  compound_figure <- (all_plots$comp_dis_mat_tree_plot | (all_plots$div_plotFuncHillDiv / all_plots$nmds_plot / all_plots$div_prof_plot) ) + plot_layout(width = c(3, 1))
  
  return(compound_figure)
  
  # return(gridExtra::grid.arrange(grobs = all_plots, ncol = ceiling(sqrt(length(all_plots)))))
}

create_survival_plot <- function(long_df, time_var, event_var) {
  # 1. Create the survival object
  survival_object <- Surv(time = long_df[[time_var]], event = long_df[[event_var]])
  
  # 2. Fit the Kaplan-Meier model
  survival_fit <- survfit(survival_object ~ treatment, data = long_df)
  
  # 3. Extract survival probabilities and confidence intervals
  survival_data <- surv_summary(survival_fit, data = long_df)
  survival_data_ci <- summary(survival_fit, conf.int = TRUE)
  
  # 4. Calculate proportion of healed wounds and confidence intervals
  survival_data <- survival_data %>%
    mutate(
      prob = surv * 100,
      proportion_healed = 1 - surv,
      upper_ci = 1 - lower,
      lower_ci = 1 - upper
    )
  
  # 5. Perform log-rank test
  logrank_test <- survdiff(survival_object ~ treatment, data = long_df)
  
  # 6. Create the plot
  survival_plot <- survival_data %>%
    ggplot(aes(x = time, y = prob, color = strata)) +
    geom_step() +
    geom_point(data = filter(survival_data, time %in% c(5, 10))) + # Filter data for geom_point
    labs(x = "Time (days)", y = "Wound Persistance (%)", color = "Treatment") +
    scale_color_manual(values = c("darkblue", "orange")) +
    # theme_bw() +
    theme_classic() +
    scale_x_continuous(
      breaks = c(0, 2, 5, 10), # Set desired x-axis breaks
      limits = c(0, 10), # Set x-axis limits between 0 and 10
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      breaks = c(0, 50, 100),
      limits = c(0, 100), # Set y-axis limits between 0 and 1
      expand = expansion(mult = c(0, 0))
    ) +
    annotate("text",
             x = 5, y = 20,
             # label = paste("Log rank ", "P", " = ", round(logrank_test$pvalue, 2)),
             label = bquote(paste("Log rank ", italic(P), " = ", .(round(logrank_test$pvalue, 2)))),  # Add .( ... ) around the p-value
             # label = paste(bquote(paste("Log rank ", italic(P), " = ")), round(logrank_test$pvalue, 2)),  # Use expression()
             size = 7.5
    ) +
    ggeasy::easy_text_size(size = 20) +
    theme(legend.position = "top") +
    guides(
      color = guide_legend(title.position = "left"),
      override.aes = list(shape = c(15, 16)) # 15 = square, 16 = circle +
    ) +
    
    # Add this to remove all grid lines except the major ones on the y-axis
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major.x = element_blank(), # Remove vertical grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      panel.grid.major.y = element_line(
        color = "grey", # Set the color of the line
        linetype = "solid" # You can change this to "dashed" or other linetypes if needed
      )
    ) +
    theme(axis.line = element_line(colour = "black"))
  
  return(survival_plot)
}

combine_gsea_results <- function(gsea_results_list, comparison_names, cluster_labels) {
  # Combine GSEA results into a compareClusterResult object
  
  result_list <- list()
  # Add cluster labels to each result
  for (i in seq_along(gsea_results_list)) {
    result_list[[i]] <- gsea_results_list[[i]]$gsea_results@result
    result_list[[i]]$Cluster <- cluster_labels[i]
  }
  
  # Combine results into a single data frame
  combined_results <- result_list %>%
    bind_rows()
  
  # Create compareClusterResult object
  compare_cluster_obj <- new("compareClusterResult",
                             compareClusterResult = combined_results,
                             geneClusters = gsea_results_list %>%
                               purrr::map("gsea_results") %>%
                               purrr::map("@geneList"),
                             fun = "GSEA")
  
  # Reorder the levels of the 'Cluster' column
  compare_cluster_obj@compareClusterResult$Cluster <- factor(
    compare_cluster_obj@compareClusterResult$Cluster,
    levels = cluster_labels  # Use the provided cluster_labels for ordering
  )
  
  return(compare_cluster_obj)
}

calculate_network_measures_parallel <- function(igraph_object, n_cores = detectCores() - 2) {
  pacman::p_load(centiserve)
  pacman::p_load(igraph)
  pacman::p_load(parallel)
  
  # Get the adjacency matrix of your network
  adjacency_matrix <- as_adjacency_matrix(igraph_object, names = FALSE)
  
  # Calculate the maximum eigenvalue
  max_eigenvalue <- max(eigen(adjacency_matrix)$values)
  
  # Calculate the upper limit for alpha
  alpha_limit <- 1 / max_eigenvalue
  
  # Choose a valid alpha (ensure it's less than alpha_limit and within 0-0.2)
  valid_alpha <- min(0.1, alpha_limit / 2) 
  
  # List of functions to be applied in parallel
  centrality_functions <- list(
    degree_centrality = function(g) igraph::degree(g, normalized = FALSE),
    betweenness_centrality = function(g) igraph::betweenness(g, directed = FALSE, weights = NULL),
    closeness_centrality = function(g) igraph::closeness(g, mode = "all", weights = NULL, normalized = TRUE),
    eigenvector_centrality = function(g) igraph::eigen_centrality(g, weights = NULL)$vector,
    # entropy_centrality = function(g) centiserve::entropy(g, weights = NULL),
    katz_centrality = function(g) centiserve::katzcent(g, alpha = valid_alpha),
    # topological_coefficient = function(g) igraph::topocoefficient(g),
    hubness_score = function(g) igraph::hub_score(g, weights = NULL)$vector
  )
  
  # Apply the functions in parallel
  centrality_results <- mclapply(centrality_functions, function(f) f(igraph_object), mc.cores = n_cores)
  
  # Combine results into a data frame
  centrality_df <- as.data.frame(centrality_results)
  
  return(centrality_df)
}
                                 
calculate_network_measures <- function(igraph_object) {
  pacman::p_load(centiserve)
  pacman::p_load(igraph)
  
  # Get the adjacency matrix of your network
  adjacency_matrix <- as_adjacency_matrix(igraph_object, names = FALSE)
  
  # Calculate the maximum eigenvalue
  max_eigenvalue <- max(eigen(adjacency_matrix)$values)
  
  # Calculate the upper limit for alpha
  alpha_limit <- 1 / max_eigenvalue
  
  # Choose a valid alpha (ensure it's less than alpha_limit and within 0-0.2)
  valid_alpha <- min(0.1, alpha_limit / 2) 
  
  centrality_df <- data.frame(
    # identify hubs
    degree_centrality = igraph::degree(igraph_object, normalized = FALSE),
    # Calculate betweenness centrality
    betweenness_centrality = igraph::betweenness(igraph_object, directed = FALSE, weights = NULL),
    # Calculate closeness centrality
    closeness_centrality = igraph::closeness(igraph_object, mode = "all", weights = NULL, normalized = TRUE),
    # Calculate eigenvector centrality
    eigenvector_centrality = igraph::eigen_centrality(igraph_object, weights = NULL)$vector,
    # Entropy centrality measures centrality of nodes depending on their contribution to the entropy of the graph.
    entropy_centrality = centiserve::entropy(igraph_object, weights = NULL),
    # Katz centrality (Katz Status Index)
    katz_centrality = centiserve::katzcent(igraph_object, alpha = valid_alpha),
    # Topological coefficient quantifies the extent to which a node shares neighbors with other nodes. Nodes with high topological coefficients tend to have neighbors that are also connected to each other.
    topological_coefficient = topocoefficient(igraph_object),
    hubness_score = hub_score(igraph_object, weights = NULL)$vector
  )
  
  return(centrality_df)
}


get_coexpression_hubs <- function(cem) {
  all_degrees <- intramodularConnectivity(
    adjMat = cem@adjacency,
    colors = cem@module$modules
  )
  
  all_degrees <- all_degrees %>%
    rownames_to_column(var = "gene") %>%
    mutate_if(is.numeric, ~ round(., digits = 2)) %>%
    column_to_rownames(var = "gene")
  
  hubness_by_module_df <- all_degrees %>%
    mutate(module = cem@module$modules) %>%
    split(f = .$module) %>%
    purrr::map(function(x) {
      x %>%
        dplyr::select(-module) %>%
        rownames_to_column(var = "gene") %>%
        arrange(desc(kTotal))
    }) %>%
    bind_rows(.id = "module")
  
  return(hubness_by_module_df)
}

                                 
new_plot_network_interaction_v2 <- function(ig_obj, n = 10, color = "blue") {
  # Calculates Degree
  degrees <- ig_obj %>% igraph::degree(normalized = FALSE)
  ig_obj <- ig_obj %>% igraph::set_vertex_attr("degree", value = degrees)
  max_n <- min(n, length(degrees))

  # Get adjacency matrix and coordinates directly from igraph
  adj_matrix <- igraph::as_adjacency_matrix(ig_obj) %>% as.matrix()
  
  # Convert igraph adjacency matrix to network object for sna
  # net_obj <- network::network(adj_matrix, directed = FALSE)
  
  # Get layout coordinates using sna
  # plot_coord <- as.data.frame(sna::gplot.layout.fruchtermanreingold(network::as.matrix.network.adjacency(net_obj), NULL))
  # plot_coord <- as.data.frame(igraph::layout_with_fr(ig_obj)) # Use igraph's layout function
  plot_coord <- as.data.frame(igraph::layout_nicely(ig_obj))
  plot_coord <- plot_coord %>% purrr::set_names("X1", "X2")
  
  # Get edge list from igraph (numeric indices)
  edge_list <- igraph::get.edgelist(ig_obj, names = FALSE) # Get numeric indices

  # Create a mapping from vertex names to indices for efficient lookup
  vertex_index_map <- setNames(seq_along(igraph::V(ig_obj)), igraph::V(ig_obj)$name)

  # Get vertex names and add to plot_coord
  vertex_names <- igraph::V(ig_obj)$name
  plot_coord$vertex_names <- vertex_names
  plot_coord$degree <- degrees
  
  # Identify hubs
  threshold <- quantile(degrees, 0.9)
  int_hubs <- names(sort(degrees, decreasing = TRUE)[1:max_n])
  int_hubs <- int_hubs[degrees[int_hubs] >= threshold]
  plot_coord$hub <- ifelse(plot_coord$vertex_names %in% int_hubs, "Interaction", "")
  
  # Prepare data for ggplot
  plot_coord$degree_cut <- cut(plot_coord$degree, breaks = 3, labels = FALSE)
  
  # Construct edges dataframe using numeric indices
  edges <- data.frame(
    X1 = plot_coord[edge_list[, 1], "X1"],
    Y1 = plot_coord[edge_list[, 1], "X2"],
    X2 = plot_coord[edge_list[, 2], "X1"],
    Y2 = plot_coord[edge_list[, 2], "X2"]
  )

  # Define custom ggplot theme
  network_theme <- ggplot2::theme_bw(base_size = 12, base_family = "") +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  
  # Create plot
  pl <- ggplot(plot_coord) +
    geom_segment(
      data = edges,
      aes(x = X1, y = Y1, xend = X2, yend = Y2),
      linewidth = 0.5, alpha = 0.5, colour = "#DDDDDD"
    ) +
    geom_point(aes(x = X1, y = X2, size = degree, alpha = degree), color = color) +
    ggrepel::geom_label_repel(
      aes(x = X1, y = X2, label = vertex_names, color = hub),
      box.padding = unit(1, "lines"),
      data = ~ subset(.x, vertex_names %in% int_hubs)
    ) +
    labs(title = "KEGG network") +
    network_theme
  
  return(pl)
}


plot_deg_interaction <- function(ig_obj, n = 10, color = "blue", name, degs_by_comparison) {
  degrees <- igraph::degree(ig_obj, normalized = FALSE)
  ig_obj <- igraph::set_vertex_attr(ig_obj, "degree", value = degrees)
  max_n <- min(n, length(degrees))
  net_obj <- intergraph::asNetwork(ig_obj)
  network_adjacency_matrix <- network::as.matrix.network.adjacency(net_obj) # get sociomatrix
  # get coordinates from Fruchterman and Reingold's force-directed placement algorithm.
  plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(network_adjacency_matrix, NULL))
  colnames(plotcord) <- c("X1", "X2")
  edglist <- network::as.matrix.network.edgelist(net_obj)
  edges <- data.frame(plotcord[edglist[, 1], ], plotcord[edglist[, 2], ])
  plotcord$vertex.names <- as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
  plotcord$Degree <- network::get.vertex.attribute(net_obj, "degree")
  plotcord[, "shouldLabel"] <- FALSE
  plotcord[, "Hub"] <- ""
  # Example: Select top 10% of vertices as hubs
  threshold <- quantile(degrees, 0.9)  
  int_hubs <- names(sort(degrees, decreasing = TRUE))[1:max_n]
  int_hubs <- int_hubs[degrees[int_hubs] >= threshold]
  int_bool <- plotcord[, "vertex.names"] %in% int_hubs
  plotcord[which(int_bool), "Hub"] <- "Interaction"
  sel_vertex <- int_hubs
  
  colnames(edges) <- c("X1", "Y1", "X2", "Y2")
  # edges$midX  <- (edges$X1 + edges$X2) / 2
  # edges$midY  <- (edges$Y1 + edges$Y2) / 2
  plotcord[which(plotcord[, "vertex.names"] %in% sel_vertex), "shouldLabel"] <- TRUE
  plotcord$Degree_cut <- cut(plotcord$Degree, breaks = 3, labels = FALSE)
  plotcord$in_mod <- TRUE
  # degs_by_comparison <- cem@module[cem@module$modules==name,]$genes
  not_in <- setdiff(plotcord[, "vertex.names"], degs_by_comparison)
  plotcord[which(plotcord[, "vertex.names"] %in% not_in), "in_mod"] <- FALSE
  
  pl <- ggplot(plotcord) +
    geom_segment(
      data = edges, aes(x = X1, y = Y1, xend = X2, yend = Y2),
      linewidth = 0.5, alpha = 0.5, colour = "#DDDDDD"
    ) +
    geom_point(aes(x = X1, y = X2, size = Degree, alpha = Degree), color = color) +
    geom_label_repel(aes(x = X1, y = X2, label = vertex.names, color = Hub),
                     box.padding = unit(1, "lines"),
                     data = function(x) {
                       x[x$shouldLabel, ]
                     }
    ) +
    scale_colour_manual(values = c(
      "Co-expression" = "#005E87",
      "Interaction" = "#540814",
      "Co-expression + Interaction" = "#736E0B"
    )) +
    labs(title = name) +
    ggplot2::theme_bw(base_size = 12, base_family = "") +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  
  return(pl)
}

get_plot_colors <- function(deg_interaction_list) {
  all_colors <- colors()
  
  # Function to check if a color is a shade of grey
  is_grey <- function(color) {
    rgb_values <- col2rgb(color)
    max_diff <- max(abs(diff(rgb_values)))
    return(max_diff < 20) # Adjust the threshold (20) as needed
  }
  
  # Filter out grey colors
  color_list <- all_colors[!sapply(all_colors, is_grey)]
  
  set.seed(123)
  color_list <- sample(color_list, size = length(deg_interaction_list)) %>% purrr::set_names(names(deg_interaction_list))
  return(color_list)
}


plot_heatmap_and_dendro <- function(df) {
  library(cowplot)
  library(ggdendro)
  library(tidyverse)
  
  sample_names <- names(df)
  # Calculate the dendrogram
  dend <- df %>%
    # column_to_rownames(var = "Gene") %>%
    as.matrix() %>%
    dist() %>%
    hclust() %>%
    as.dendrogram()
  
  dend_data <- dend %>% dendro_data()
  
  # Setup the data, so that the layout is inverted (this is more
  # "clear" than simply using coord_flip())
  segment_data <- with(
    segment(dend_data),
    data.frame(x = y, y = x, xend = yend, yend = xend)
  )
  
  # Use the dendrogram label data to position the gene labels
  gene_pos_table <- with(
    dend_data$labels,
    data.frame(y_center = x, Gene = as.character(label), height = 1)
  )
  
  # Table to position the samples
  sample_pos_table <- data.frame(Sample = sample_names) %>%
    mutate(x_center = (1:nrow(.)), width = 1)
  
  # Neglecting the gap parameters
  heatmap_data <- df %>%
    rownames_to_column("Gene") %>%
    gather(value = "expr", key = "Sample", -Gene) %>%
    # reshape2::melt(value.name = "expr", varnames = c("Gene", "Sample")) %>%
    left_join(gene_pos_table) %>%
    left_join(sample_pos_table)
  
  # Limits for the vertical axes
  gene_axis_limits <- with(
    gene_pos_table,
    c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
  ) +
    0.1 * c(-1, 1) # extra spacing: 0.1
  
  max_expr <- max(abs(floor(min(heatmap_data$expr, na.rm = T))), abs(ceiling(max(heatmap_data$expr, na.rm = T))))
  
  # Heatmap plot
  plt_hmap <- ggplot(
    heatmap_data,
    aes(
      x = x_center, y = y_center, fill = expr,
      height = height, width = width
    )
  ) +
    geom_tile() +
    # scale_fill_gradient2("Expr", high = "red", low = "blue") +
    scale_fill_gradientn("log2FC", 
                         limits = c(-1*max_expr, max_expr),
                         colours = rev(colorRampPalette(colors = c("red", "white", "#094FED"))(50))) +
    scale_x_continuous(
      breaks = sample_pos_table$x_center,
      labels = sample_pos_table$Sample,
      expand = c(0, 0)
    ) +
    # For the y axis, alternatively set the labels as: gene_position_table$gene
    scale_y_continuous(
      breaks = gene_pos_table[, "y_center"],
      # labels = rep("", nrow(gene_pos_table)),
      labels = heatmap_data %>% arrange(y_center) %>% pull(Gene) %>% unique(),
      limits = gene_axis_limits,
      expand = c(0, 0),
      position = "right"
    ) +
    labs(x = "Treatment", y = "") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = rel(1), hjust = 0.5, angle = 0),
      # margin: top, right, bottom, and left
      plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"),
      panel.grid.minor = element_blank()
    ) +
    ggeasy::easy_all_text_size(size = 20) +
    ggeasy::easy_y_axis_labels_size(size = 0) +
    ggeasy::easy_rotate_x_labels(angle = 90)
  
  # Dendrogram plot
  plt_dendr <- ggplot(segment_data) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_reverse(expand = c(0, 0.5)) +
    scale_y_continuous(
      breaks = gene_pos_table$y_center,
      labels = gene_pos_table$Gene,
      # labels = NULL,
      limits = gene_axis_limits,
      expand = c(0, 0)
    ) +
    labs(x = "Distance", y = "", colour = "", size = "") +
    theme(panel.grid.minor = element_blank()) +
    theme_bw() +
    ggeasy::easy_all_text_size(size = 20)
  
  # results <- list(
  # heatmap = plot_grid(plt_dendr, plt_hmap, align = "h", rel_widths = c(1, 1)),
  # heatmap_data = 
  # )
  
  cowplot::plot_grid(plt_dendr, plt_hmap, align = "h", rel_widths = c(0.4, 1))
}

run_cemitool_pipeline <- function(cem,
                                  sample_annotation_df,
                                  database = KEGG_2019_Human,
                                  int_df = NULL,
                                  n_cores = parallel::detectCores() - 1) {
  pacman::p_load(doParallel, 
                 CEMiTool, 
                 tidyverse,
                 parallelMap)
  
  cem@selected_genes <- cem@selected_genes %>% str_to_upper()
  
  # # Set the number of cores to use
  # n_cores <- n_cores  # Or however many cores you want to utilize
  # 
  # # Create a parallel cluster
  # cl <- makeCluster(n_cores)
  # 
  # # Register the parallel backend
  # registerDoParallel(cl)
  
  # # Set parallel backend (using available cores, leaving one free)
  parallelStartSocket(n_cores)
  
  # Module identification
  cem <- cem %>%
    find_modules(
      cor_method = "pearson",
      cor_function = "cor",
      eps = 0.1,
      min_ngen = 20,
      merge_similar = TRUE,
      diss_thresh = 0.75,
      network_type = "unsigned",
      tom_type = "signed",
      verbose = FALSE,
      force_beta = TRUE,
      set_beta = NULL
    )
  
  cem@module <- cem@module %>% 
    mutate(genes = genes %>% str_to_upper())
  # Visualization of module properties (mean connectivity, beta, and gene expression profiles)
  cem <- cem %>% plot_mean_k(title = "Mean connectivity")
  cem <- cem %>% plot_beta_r2()
  cem <- cem %>% plot_profile()
  
  # Extract selected genes for downstream analysis
  genes <- cem %>%
    unclass() %>%
    attr("selected_genes")
  
  # Add sample annotation to the CEMiTool object
  sample_annotation(cem, sample_name_column = "sample_name", class_column = "class") <- sample_annotation_df
  
  # Gene Set Enrichment Analysis (GSEA)
  cem <- cem %>% mod_gsea()
  cem <- cem %>% plot_gsea()
  
  # Over-Representation Analysis (ORA)
  gmt_in <- database %>% read_gmt()
  cem <- cem %>% mod_ora(gmt_in)
  cem <- cem %>% plot_ora(n = 20, pv_cut = 0.05)
  
  # Interactions Analysis (if applicable)
  if (!is.null(int_df)) { # Check if int_df exists
    interactions_data(cem) <- int_df
    genes_by_module <- split(cem@module$genes, cem@module$modules)
    cem@interactions <- genes_by_module %>% purrr::map(function(x) {
      rows <- which(int_df[[1]] %in% x & int_df[[2]] %in% x)
      ig <- igraph::simplify(igraph::graph_from_data_frame(int_df[rows,], directed=FALSE))
      return(ig)
    })
    
  } else {
    message("Skipping interactions analysis: 'int_df' not provided.")
  }
  
  # Stop parallel processing
  parallelStop()
  
  return(cem)
}

get_filtered_deg_set <- function(gsea_results, gseaplot_results, edgeR_topDE, description) {

  # Extract core enrichment based on description
  core_enrichment <- gsea_results@result$core_enrichment[gsea_results@result$Description == description]

  # Uncomment if you need to access gseaplot_results
  # gseaplot_result <- gseaplot_results[[description]] 

  # Get deg_set using str_detect
  deg_set <- str_detect(
    string = core_enrichment, 
    pattern = str_to_upper(edgeR_topDE$mgi_symbol)
  )

  # Filter edgeR_topDE based on deg_set
  filtered_deg_set <- edgeR_topDE %>% 
    dplyr::filter(mgi_symbol %in% edgeR_topDE$mgi_symbol[deg_set])

  return(filtered_deg_set)
}

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
