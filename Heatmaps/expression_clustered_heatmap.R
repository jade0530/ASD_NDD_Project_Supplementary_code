# Library
library(ggplot2)
library(ggdendro)
library(hrbrthemes)
library(readxl)
library(tidyr)
library(grid)
library(dplyr)
library(readr)
library(RColorBrewer)
library(dendsort)
library(pheatmap)

#Used complex Heatmap instead of pheatmap for more options
library(ComplexHeatmap)
library(viridis)

options(scipen = 0)
# Here to load data
unique_genes <- read_excel('/media/bml/USER_DATA/ASD_CNV_Project/HEATMAP/CNVs_significant_on_cnv_New_V3_heatmap_filled.xlsx', sheet = "Significant_genes")
fantom <- read_csv('/media/bml/USER_DATA/ASD_CNV_Project/HEATMAP/FantomCAT_robust_genes_vs_traits_pvalue_nonzero_04072021_selected.csv')

# Prepare the dataframe
fantom_prep <- fantom
colnames(fantom_prep) <- fantom[1,]
colnames(fantom_prep)[3] <- 'chr_num'
fantom_prep <- fantom_prep[2:nrow(fantom_prep),]
fantom_prep <- subset(fantom_prep, select = -c(2,4,5,6,8))
fantom_prep[, c(4:ncol(fantom_prep))] <- mutate_all(fantom_prep[, c(4:ncol(fantom_prep))] , function(x) as.double(x)) 

unique_genes <- unique_genes[unique_genes$CNV_ID=="chr22:18844632-18876427",]
unique_genes$`'Gene.name'` <- substr(unique_genes$`'Gene.name'`, 1, nchar(unique_genes$`'Gene.name'`)-1)
unique_genes_prep <- unique_genes[,c(2,3)]

fantom_prep_mapped = merge(unique_genes_prep, fantom_prep, by.x = "Gene_ID", by.y = "geneID", all = FALSE)
fantom_prep_mapped_cleaned <- subset(fantom_prep_mapped, select = -c(4,5))
fantom_prep_mapped_cleaned[, 4:ncol(fantom_prep_mapped_cleaned)] = mutate_all(fantom_prep_mapped_cleaned[, 4:ncol(fantom_prep_mapped_cleaned)],function(x) as.double(x)) 

write_csv(fantom_prep_mapped, "fantom_prep_mapped.csv")
rownames(fantom_prep_mapped_cleaned) <- fantom_prep_mapped_cleaned$Gene_ID                                             
fantom_prep_mapped_cleaned_log10 <- mutate_all(fantom_prep_mapped_cleaned, function(x) -log10(x))
fantom_prep_mapped_cleaned_log10 <- mutate_all(fantom_prep_mapped_cleaned_log10, function(x) set_value_to_10(x))

set_value_to_10 <- function(x) {
  for (e in 1:length(x)) {
    if (x[e] > 8.5) {
      x[e] = 8.5
    } 
  }
  return(x)
}

asd_unique = fantom_prep_mapped_cleaned[fantom_prep_mapped_cleaned['Phenotype']=="ASD",]
asd_unique_numeric <- asd_unique[,4:ncol(asd_unique)]
rownames(asd_unique_numeric) <- asd_unique$Gene_ID                                             
asd_unique_numeric_log10 <- mutate_all(asd_unique_numeric, function(x) -log10(x))
asd_unique_numeric_log10 <- mutate_all(asd_unique_numeric_log10, function(x) set_value_to_10(x))


asd_id_dd_unique = fantom_prep_mapped_cleaned[fantom_prep_mapped_cleaned['Phenotype']=="ASD_ID_DD",]
asd_id_dd_unique_numeric <- asd_id_dd_unique[,4:ncol(asd_id_dd_unique)]
rownames(asd_id_dd_unique_numeric) <- asd_id_dd_unique$Gene_ID     
asd_id_dd_unique_numeric_log10 <- mutate_all(asd_id_dd_unique_numeric, function(x) -log10(x))
asd_id_dd_unique_numeric_log10 <- mutate_all(asd_id_dd_unique_numeric_log10, function(x) set_value_to_10(x))


id_dd_unique = fantom_prep_mapped_cleaned[fantom_prep_mapped_cleaned['Phenotype']=="ID_DD",]
id_dd_unique_numeric <- id_dd_unique[,4:ncol(id_dd_unique)]
rownames(id_dd_unique_numeric) <- id_dd_unique$Gene_ID     
id_dd_unique_numeric_log10 <- mutate_all(id_dd_unique_numeric, function(x) -log10(x))
id_dd_unique_numeric_log10 <- mutate_all(id_dd_unique_numeric_log10, function(x) set_value_to_10(x))

id_dd_unique_numeric_log2 <- mutate_all(id_dd_unique_numeric, function(x) -log2(x))
id_dd_unique_numeric_log2 <- mutate_all(id_dd_unique_numeric_log2, function(x) set_value_to_10(x))


scz_unique = fantom_prep_mapped_cleaned[fantom_prep_mapped_cleaned['Phenotype']=="SCZ",]
scz_unique_numeric <- scz_unique[,4:ncol(scz_unique)]
rownames(scz_unique_numeric) <- scz_unique$Gene_ID   
scz_unique_numeric_log10 <- mutate_all(scz_unique_numeric, function(x) -log10(x))
scz_unique_numeric_log10 <- mutate_all(scz_unique_numeric_log10, function(x) set_value_to_10(x))

# build robust dist method
robust_dist = function(x, y) {
  qx = quantile(x, c(0.1, 0.9))
  qy = quantile(y, c(0.1, 0.9))
  l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
  x = x[l]
  y = y[l]
  sqrt(sum((x - y)^2))
}

# BUILD COLOUR PALATTE
palette_blues <- colorRampPalette(colors = c("#B2182B", "#D1E5F0"))(6)
palette_log2 <- colorRampPalette(colors = c("#006633", "#D1E5F0"))(6)

generate_heatmap <- function(mat) {
  stopifnot(!missing(mat))
  h <- Heatmap(
    as.matrix(mat), 
    name = "-log2(p-value)",
    col = rev(palette_log2),
    column_title  = "Tissues",
    row_title = "Genes",
    column_title_side = "bottom",
    show_row_names = FALSE, 
    width = unit(450, "mm"),
    height = unit(150, "mm"),
    #heatmap_height = nrow(asd_id_dd_unique_numeric_log10)*unit(0.5,"mm"),
    column_dend_reorder  = TRUE,
    column_names_rot = 45,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "single",
    column_dend_height = unit(4, "cm"),
    row_dend_width = unit(4, "cm")
  )
  return(h)
}

save_heatmap_pdf <- function(h, filename, width = 50, height = 50) {
  stopifnot(!missing(filename))
  stopifnot(!missing(h))
  pdf(filename, width=width, height=height)
  draw(h,heatmap_legend_side = "right")
  dev.off()
}

h_new <- Heatmap(
        as.matrix(asd_id_dd_unique_numeric_log10), 
        name = "-log10(p-value)",
        col = rev(palette_blues),
        column_title  = "Tissues",
        row_title = "Genes",
        column_title_side = "bottom",
        show_row_names = FALSE, 
        width = unit(450, "mm"),
        height = unit(150, "mm"),
        #heatmap_height = nrow(asd_id_dd_unique_numeric_log10)*unit(0.5,"mm"),
        column_dend_reorder  = TRUE,
        column_names_rot = 45,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "single",
        column_dend_height = unit(4, "cm"),
        row_dend_width = unit(4, "cm")
)

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

h_pheatmap <- pheatmap::pheatmap(
  mat = as.matrix(asd_unique_numeric), 
  name = "FantomCAT Tissue Significances versus Genes",
  color = palette_blues,
  scale = "row",
  show_rownames = FALSE, 
  width = unit(450, "mm"),
  height = unit(150, "mm"),
  #heatmap_height = nrow(asd_id_dd_unique_numeric_log10)*unit(0.5,"mm"),
  angle_col = 45,
  clustering_callback = callback,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "single",
  filename = "asd_id_dd_fantom_pheatmap.pdf"
)

save_pheatmap_pdf(mat=asd_unique_numeric, label=asd_unique, filename="asd_fantom_heatmap_eu_single.pdf", width=30, height=20)
save_pheatmap_pdf(mat=asd_id_dd_unique_numeric, label=asd_id_dd_unique, filename="asd_id_dd_fantom_heatmap.pdf")
save_pheatmap_pdf(mat=id_dd_unique_numeric, label=id_dd_unique, filename="id_dd_fantom_heatmap.pdf")
save_pheatmap_pdf(mat=scz_unique_numeric, label=scz_unique, filename="scz_fantom_heatmap.pdf")


asd_h <- generate_heatmap(asd_unique_numeric_log10)
asd_id_dd_h <- generate_heatmap(asd_id_dd_unique_numeric_log10)
id_dd_h <- generate_heatmap(id_dd_unique_numeric_log10)
id_dd_h_log2 <- generate_heatmap(id_dd_unique_numeric_log2)
scz_h <- generate_heatmap(scz_unique_numeric_log10)


