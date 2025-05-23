##data cleaning
Draw_box_plot<-function(box,x,width,c,lwd,line_col){
  segments(x, box[2], x, box[3], col = line_col,lwd =lwd)
  segments(x-(width/2), box[2], x+(width/2), box[2], col = line_col,lwd =lwd)
  segments(x-(width/2), box[3], x+(width/2), box[3], col = line_col,lwd =lwd)
  rect(x-width, box[4], x+width, box[5], col = c,lwd =lwd, border = line_col)
  segments(x-width, box[1], x+width, box[1], col = line_col,lwd=2*lwd)}
Means_factor = function(factor, x){
  m = NULL
  for(i1 in c(1:length(levels(factor)))){
    x1 = x[which(factor==levels(factor)[i1])]
    x1 = x1[which(x1!=-1)]
    m = c(m, mean(x1))}
  return(m)}
Medians_factor = function(factor, x){
  m = NULL
  for(i1 in c(1:length(levels(factor)))){
    x1 = x[which(factor==levels(factor)[i1])]
    x1 = x1[which(x1!=-1)]
    m = c(m, median(x1))}
  return(m)}
concat = function(v) {
  res = ""
  for (i in 1:length(v)){res = paste0(res,v[i])}
  res
}
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha)) }

###########
out_dir = "~/Desktop/SpatioEV_manuscript_code/"
batch = "spatial"
##############
file = concat(c(out_dir, "Eigenvectors_Spatial_Spatial.txt"))
file = concat(c(out_dir, "Eigenvectors_fibre_Spatial.txt"))
file = concat(c(out_dir, "Eigenvectors_neighbourhood_Spatial.txt"))
file = concat(c(out_dir, "Imputed_DATA_FINAL_SCALED_Spatial_Spatial.txt"))
file = concat(c(out_dir, "Imputed_DATA_FINAL_SCALED_fibre_Spatial.txt"))
file = concat(c(out_dir, "Imputed_DATA_FINAL_SCALED_neighbourhood_Spatial.txt"))
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
mat = p
heatmap(p)

rownames(mat) == rownames(p)
p = as.data.frame(p)
p$sample_type <- ifelse(grepl("^JRP112|^JRP122|^JRP141", rownames(p)), "RA", "OA")
p$cell_type <- sub(".*_", "", rownames(p))
p$sample <- gsub("_(.*)", "", rownames(p))

data = p
        
# Load required libraries
library(ggplot2)  
library(dplyr)   
library(tidyr)   


#color scheme
dark_set2_colors <- c("RA" = "#1B9E77", "OA" = "#D95F02", "Not Significant" = "#F0F0F0")

# using ANOVA and Benjamini-Hochberg correction for multiple testing
all_cell_type_results <- data %>%
  group_by(cell_type) %>%  # Analyze each cell type separately
  group_map(function(cell_data, cell_info) {
    cell_type_name <- cell_info$cell_type  # Extract cell type name

    df_long <- cell_data %>%
      pivot_longer(
        cols = starts_with("Module_"),
        names_to = "Module",
        values_to = "Expression"
      )

    # Run ANOVA on expression by sample group (RA vs OA)
    anova_results <- df_long %>%
      group_by(Module) %>%
      do(model = aov(Expression ~ sample, data = .)) %>%
      mutate(p.value = summary(model)[[1]][["Pr(>F)"]][1]) %>%  # Extract p-value
      select(Module, p.value)

    # Adjust p-values using Benjamini-Hochberg method
    bh_adjusted_results <- anova_results %>%
      mutate(padj = p.adjust(p.value, method = "BH"))

    # Compute average expression per module per sample group
    # and assign significance direction based on adjusted p-values and group means
    significant_results <- df_long %>%
      group_by(Module, sample) %>%
      summarise(mean_expression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = sample, values_from = mean_expression) %>%
      left_join(bh_adjusted_results, by = "Module") %>%
      mutate(Significance = case_when(
        padj < 0.05 & !is.na(RA) & !is.na(OA) & RA > OA ~ "RA",              # Significantly higher in RA
        padj < 0.05 & !is.na(RA) & !is.na(OA) & OA > RA ~ "OA",              # Significantly higher in OA
        TRUE ~ "Not Significant"                                            # Not significant
      )) %>%
      select(Module, Significance, padj) %>%
      mutate(cell_type = cell_type_name)  # Retain cell type info
  }) %>%
  bind_rows()  # Combine results from all cell types

# Create full grid of Module × Cell Type to ensure completeness in the heatmap
all_combinations <- expand_grid(
  Module = paste0("Module_", 1:32),
  cell_type = unique(data$cell_type)
)


heatmap_matrix_data <- all_combinations %>%
  left_join(all_cell_type_results, by = c("Module", "cell_type")) %>%
  replace_na(list(Significance = "Not Significant", padj = 1)) %>%  
  mutate(Module = factor(Module, levels = paste0("Module_", 1:32))) 

matrix_heatmap_plot <- ggplot(heatmap_matrix_data, aes(y = cell_type, x = Module, fill = Significance)) +
  geom_tile(color = "black") +
  scale_fill_manual(
    values = dark_set2_colors,
    name = "Significantly Higher In"
  ) +
  labs(
    y = "Cell Type",
    x = "Module",
    title = "Module Significance Across Cell Types"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    axis.ticks.y = element_blank() 
  )

# Display the heatmap plot
print(matrix_heatmap_plot)
