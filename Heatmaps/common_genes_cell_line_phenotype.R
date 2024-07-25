library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)


setwd("/Users/jadeliang/Desktop/Mona_project_CNV/cell-cell_communication/")

features <- read_csv("./feature.csv")

gene_phenotypic_group_unique <- 
  read_excel("/Users/jadeliang/Desktop/Mona_project_CNV/cell-cell_communication/CNVs_significant_on_cnv_New_V3.xlsx", sheet = "Significant_genes_unique")


gene_phenotypic_group_unique_extract <- gene_phenotypic_group_unique[gene_phenotypic_group_unique$Gene_type=='protein_coding',]
gene_phenotypic_group_unique_extract <- gene_phenotypic_group_unique_extract %>% filter(!grepl("^CAT", Gene_Name))

gene_phenotypic_group_unique_extract <- gene_phenotypic_group_unique_extract[,c("Gene_Name","Phenotype")]


# calculate how many genes significant in each p group versus each cell line cluster 
# x = phenotypes y = cluster fill = number of genes 
features_clusters <- features[,c("features", "clusters")]

# Merge clusters
features_clusters$clusters_merged <- gsub("-[0-9]+$", '', features_clusters$clusters)
features_clusters_dedup <- features_clusters[,c("features", "clusters_merged")]
features_clusters_dedup <- features_clusters_dedup[!duplicated(features_clusters_dedup),]
features_clusters_unique_merged <- merge(gene_phenotypic_group_unique_extract, features_clusters_dedup, by.x = "Gene_Name", by.y = "features")
write_csv(features_clusters_unique_merged,"common_genes_scrna_and_group_wihtout_CAT_100.csv")


common_unique_genes <- features_clusters_unique_merged %>%
  group_by(Phenotype,clusters_merged) %>%
  summarise(n_genes = n(),
            genes=paste(Gene_Name, collapse = ';'))

common_unique_genes$fraction <- 0
# Calculate fractions: no of gees that appears in the cluster specific to phe group/no of all genes specific to the group
# ASD
asd_value <- nrow(gene_phenotypic_group_unique_extract[gene_phenotypic_group_unique_extract$Phenotype=="ASD",])
common_unique_genes$fraction[common_unique_genes$Phenotype=="ASD"] <- round(100*common_unique_genes$n_genes[common_unique_genes$Phenotype=="ASD"]/asd_value,digits = 2)

# ASD_ID
asd_id_value <- nrow(gene_phenotypic_group_unique_extract[gene_phenotypic_group_unique_extract$Phenotype=="ASD_ID",])
common_unique_genes$fraction[common_unique_genes$Phenotype=="ASD_ID"] <- round(100*common_unique_genes$n_genes[common_unique_genes$Phenotype=="ASD_ID"]/asd_id_value,digits = 2)

# ID_DD
id_dd_value <- nrow(gene_phenotypic_group_unique_extract[gene_phenotypic_group_unique_extract$Phenotype=="ID_DD",])
common_unique_genes$fraction[common_unique_genes$Phenotype=="ID_DD"] <- round(100*common_unique_genes$n_genes[common_unique_genes$Phenotype=="ID_DD"]/id_dd_value,digits = 2)

# SCZ
scz_value <- nrow(gene_phenotypic_group_unique_extract[gene_phenotypic_group_unique_extract$Phenotype=="SCZ",])
common_unique_genes$fraction[common_unique_genes$Phenotype=="SCZ"] <- round(100*common_unique_genes$n_genes[common_unique_genes$Phenotype=="SCZ"]/scz_value,digits = 2)


write_csv(common_unique_genes, "common_unique_genes.csv")

# Convert to heatmap suitable table 
common_genes_heatmap <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("clusters", "ASD", "ASD_ID", "ID_DD", "SCZ"))))
common_genes_heatmap_fraction <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("clusters", "ASD", "ASD_ID", "ID_DD", "SCZ"))))

for (i in 1:nrow(common_unique_genes)) {
  clusters <- common_unique_genes$clusters_merged[i]
  asd <- common_unique_genes$n_genes[which(common_unique_genes$Phenotype=='ASD' & common_unique_genes$clusters_merged==clusters)]
  asd_fraction <- common_unique_genes$fraction[which(common_unique_genes$Phenotype=='ASD' & common_unique_genes$clusters_merged==clusters)]
  asd_id <-common_unique_genes$n_genes[which(common_unique_genes$Phenotype=='ASD_ID' & common_unique_genes$clusters_merged==clusters)]
  asd_id_fraction <-common_unique_genes$fraction[which(common_unique_genes$Phenotype=='ASD_ID' & common_unique_genes$clusters_merged==clusters)]
  id_dd <- common_unique_genes$n_genes[which(common_unique_genes$Phenotype=='ID_DD' & common_unique_genes$clusters_merged==clusters)]
  id_dd_fraction <- common_unique_genes$fraction[which(common_unique_genes$Phenotype=='ID_DD' & common_unique_genes$clusters_merged==clusters)]
  scz <- common_unique_genes$n_genes[which(common_unique_genes$Phenotype=='SCZ' & common_unique_genes$clusters_merged==clusters)]
  scz_fraction <- common_unique_genes$fraction[which(common_unique_genes$Phenotype=='SCZ' & common_unique_genes$clusters_merged==clusters)]
  
  if (length(asd)==0) {
    asd <- 0
  }
  if (length(asd_id)==0) {
    asd_id <- 0
  }
  if (length(id_dd)==0) {
    id_dd <- 0
  }
  if (length(scz)==0) {
    scz <- 0
  }
  
  if (length(asd_fraction)<=0) {
    asd_fraction <- 0
  }
  if (length(asd_id_fraction)<=0) {
    asd_id_fraction <- 0
  }
  if (length(id_dd_fraction)<=0) {
    id_dd_fraction <- 0
  }
  if (length(scz_fraction)<=0) {
    scz_fraction <- 0
  }
  res <- data.frame(c(clusters, asd, asd_id, id_dd, scz))
  res <- t(res)
  colnames(res) <- c("clusters", "ASD", "ASD_ID", "ID_DD", "SCZ")
  rownames(res) = NULL
  print(res)
  common_genes_heatmap <- rbind(common_genes_heatmap, res)
  
  res_frac <- data.frame(c(clusters, asd_fraction, asd_id_fraction, id_dd_fraction, scz_fraction))
  print(res_frac)
  res_frac <- t(res_frac)
  colnames(res_frac) <- c("clusters", "ASD", "ASD_ID", "ID_DD", "SCZ")
  rownames(res_frac) = NULL

  common_genes_heatmap_fraction <- rbind(common_genes_heatmap_fraction, res_frac)
}
common_genes_heatmap <- common_genes_heatmap[!duplicated(common_genes_heatmap),]
common_genes_heatmap$ASD_log <- log(common_genes_heatmap$ASD)
common_genes_heatmap_fraction <- common_genes_heatmap_fraction[!duplicated(common_genes_heatmap_fraction),]

common_unique_genes$log <- log(common_unique_genes$fraction/100)

write_csv(common_genes_heatmap_fraction, "common_genes_heatmap_fraction_wihtout_CAT_100.csv")
