#' add two fake new pseudotime column
#'
#' @param seu
#' @param trajectoryPseudotimes
#' @param celltype_column
#' @param trajectoryTo
#' @param trajectoryFrom
#'
#' @return
#' @export
#'
#' @examples
pseudotime_fake_order <- function(seu,
                                  trajectoryPseudotimes = c('slingPseudotime_1', 'slingPseudotime_2'),
                                  celltype_column = 'simple_celltype',
                                  trajectoryTo = c('EC', 'EE'),
                                  trajectoryFrom = 'ISC/EB'){
  # add two fake new pseudotime column: fake and fake order to seu@misc$slingshot$PCA$SlingPseudotime for downstream plot in GeneTrendHeatmap
  # only support two trajectory
  # fake: end in ends,  from in middle
  # fake order: based on fake, use the rank instead of the pseudotime
  # trajectoryTo:  corresponding cell type for slingPseudotime_1 and slingPseudotime_2
  # trajectoryTo[2] will be assigned to negative value, and shared trajectoryFrom cells for two trajectory will be randomly assigned.
  df <- seu@misc$slingshot$PCA$SlingPseudotime[,trajectoryPseudotimes]
  df$celltype <- seu@meta.data[,celltype_column]
  df[,trajectoryPseudotimes[1]][df$celltype == trajectoryTo[2]] <- NA
  df[,trajectoryPseudotimes[2]][df$celltype == trajectoryTo[1]] <- NA
  shared_cells <- rownames(df[df$celltype == trajectoryFrom & !is.na(df[,trajectoryPseudotimes[1]]) & !is.na(df[,trajectoryPseudotimes[2]]),])
  set.seed(1)
  selected_cells <- sample(shared_cells, round(length(shared_cells)/2))
  df$slingPseudotime_fake <- df[,trajectoryPseudotimes[1]]
  df$slingPseudotime_fake[df$celltype == trajectoryTo[2]] <- -df[,trajectoryPseudotimes[2]][df$celltype == trajectoryTo[2]]
  df[selected_cells, 'slingPseudotime_fake'] <- -df[selected_cells, trajectoryPseudotimes[2]]
  # table(df$celltype,df$slingPseudotime_fake > 0)
  seu@misc$slingshot$PCA$SlingPseudotime$slingPseudotime_fake <- df$slingPseudotime_fake
  seu@misc$slingshot$PCA$SlingPseudotime$slingPseudotime_fake_order <- rank(df$slingPseudotime_fake, na.last = TRUE)
  seu@misc$slingshot$PCA$SlingPseudotime$slingPseudotime_fake_order[is.na(seu@misc$slingshot$PCA$SlingPseudotime$slingPseudotime_fake)] <- NA
  # df$slingPseudotime_fake_order <- rank(df$slingPseudotime_fake, na.last = TRUE)
  return(seu)
}

#' add cell type group to GeneTrendHeatmap plot
#'
#' @param p
#' @param seu
#' @param title
#' @param celltype_column
#' @param siling_id
#' @param celltype_levels
#' @param layout_heights
#' @param feature_label_size
#'
#' @return
#' @export
#'
#' @examples
GeneTrendHeatmap_withCellGroup <- function(p,
                                           seu,
                                           title = 'slingPseudotme_fake_order',
                                           celltype_column = 'simple_celltype',
                                           siling_id = 'slingPseudotime_fake_order',
                                           celltype_levels = c('ISC/EB', 'EC', 'EE'),
                                           layout_heights = c(40, 1),
                                           feature_label_size = 6){
  library(ggplot2)
  library(patchwork)

  p <- p + scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0.5,  # Adjust midpoint as needed
    limits = c(0, 1)  # Force full range
  ) +
    theme(axis.text.y = element_text(size = feature_label_size))  +
    theme(panel.spacing = unit(0, "pt"),
          plot.margin = margin(0,0,0,0)) +
    labs(title = title)

  heatmap_pseudotime_points <- as.numeric(levels(p@data$variable))

  df <- seu@misc$slingshot$PCA$SlingPseudotime
  df$celltype <- seu@meta.data[,celltype_column]
  # remove cells with NA pseudotime value
  df <- df[!is.na(df[,siling_id]),]

  # 获取原始细胞的伪时序值
  original_pseudotime <- data.frame(pseudotime = df[,siling_id],
                                    cells = rownames(df),
                                    cell_type = df$celltype)

  # ggplot(df, aes(x=slingPseudotime_fake_order)) +
  #   geom_histogram() +
  #   facet_grid(celltype~ .)
  #
  # original_pseudotime <- original_pseudotime[!is.na(original_pseudotime$pseudotime),]
  # table(original_pseudotime$cell_type)
  pillar <- data.frame(
    xmin = min(original_pseudotime$pseudotime),
    xmax = max(original_pseudotime$pseudotime),
    ymin = 0,
    ymax = 0.1
  )

  segments_data <- data.frame(
    x     =  original_pseudotime$pseudotime,
    cell_type = factor(original_pseudotime$cell_type, levels = celltype_levels)
  )

  annotation_bar <- ggplot() +
    geom_rect(data = pillar,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax),
              fill = "gray90", color = "gray50") +
    coord_cartesian(expand = FALSE) +
    geom_segment(data = segments_data,
                 aes(x = x, xend = x,
                     y = pillar$ymin, yend = pillar$ymax,
                     color = cell_type),
                 alpha = 1,
                 linewidth = 1.2) +
    scale_x_continuous(limits = range(heatmap_pseudotime_points)) +
    # scale_color_manual(values = c("red" = "red", "blue" = "blue")) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()) +
    theme(panel.spacing = unit(0, "pt"),
          plot.margin = margin(0, 0, 0, 0)) +
    theme(legend.position = "bottom")  +
    guides(color = guide_legend(nrow = 1))
  annotation_bar

  # 组合图形
  combined_plot <- p / annotation_bar +
    plot_layout(heights = layout_heights)
  return(combined_plot)
}
