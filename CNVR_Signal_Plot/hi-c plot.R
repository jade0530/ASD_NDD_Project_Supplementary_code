library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(karyoploteR)
library(GenomicRanges)
library(ggplotify)
library(ggpubr)
library(cowplot)


# Read the 8 files 
asd_dup <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/HI_C/ASD_Dup.csv")
asd_del <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/HI_C/ASD_Del.csv")
asd_id_dd_dup <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/HI_C/ASD_ID_DD_Del.csv")
asd_id_dd_del <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/HI_C/ASD_ID_DD_Dup.csv")
id_dd_dup <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/HI_C/ID_DD_Dup.csv")
id_dd_del <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/HI_C/ID_DD_Del.csv")
scz_dup <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/HI_C/SCZ_Del.csv")
scz_del <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/HI_C/SCZ_Dup.csv")

# Normalise p-value by min-max method
min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Apply min-max scaling to both duplication and deletion datasets and set up cutoff for p-valueess
# ASD p-value TRD: 6.1E-06
# ASD_ID_DD TRD: 3.2E-24
# ID_DD TRD: 1.3E-55
# SCZ TRD: 8.7E-05

######################################################ASD#####################################################
asd_trd <- 6.1E-06
normalized_asd_trd <- -log2(asd_trd)
asd_dup$median_pos <- (asd_dup$start + asd_dup$end)/2
asd_dup$log2_p_value = -log2(asd_dup$`P-value`)
asd_dup$scaled_p_value <- min_max_scale(asd_dup$log2_p_value)
scaled_asd_dup_trd <- (normalized_asd_trd - min(asd_dup$log2_p_value)) / (max(asd_dup$log2_p_value) - min(asd_dup$log2_p_value))

asd_del$median_pos <- (asd_del$start + asd_del$end)/2
asd_del$log2_p_value = -log2(asd_del$`P-value`)
asd_del$scaled_p_value <- min_max_scale(asd_del$log2_p_value)
scaled_asd_del_trd <- (normalized_asd_trd - min(asd_del$log2_p_value)) / (max(asd_del$log2_p_value) - min(asd_del$log2_p_value))

###################################################ASD_ID_DD###################################################################
asd_id_dd_trd <- 3.2E-24
normalized_asd_id_dd_trd <- -log2(asd_id_dd_trd)
asd_id_dd_dup$median_pos <- (asd_id_dd_dup$start + asd_id_dd_dup$end)/2
asd_id_dd_dup$log2_p_value = -log2(asd_id_dd_dup$`P-value`)
asd_id_dd_dup$scaled_p_value <- min_max_scale(asd_id_dd_dup$log2_p_value)
scaled_asd_id_dd_dup_trd <- (normalized_asd_id_dd_trd - min(asd_id_dd_dup$log2_p_value)) / (max(asd_id_dd_dup$log2_p_value) - min(asd_id_dd_dup$log2_p_value))

asd_id_dd_del$median_pos <- (asd_id_dd_del$start + asd_id_dd_del$end)/2
asd_id_dd_del$log2_p_value = -log2(asd_id_dd_del$`P-value`)
asd_id_dd_del$scaled_p_value <- min_max_scale(asd_id_dd_del$log2_p_value)
scaled_asd_id_dd_del_trd <- (normalized_asd_id_dd_trd - min(asd_id_dd_del$log2_p_value)) / (max(asd_id_dd_del$log2_p_value) - min(asd_id_dd_del$log2_p_value))

####################################################ID_DD#################################################################
# ID_DD TRD: 1.3E-55
id_dd_trd <- 1.3E-55
normalized_id_dd_trd <- -log2(id_dd_trd)

id_dd_dup$median_pos <- (id_dd_dup$start + id_dd_dup$end)/2
id_dd_dup$log2_p_value = -log2(id_dd_dup$`P-value`)
id_dd_dup$scaled_p_value <- min_max_scale(id_dd_dup$log2_p_value)
scaled_id_dd_dup_trd <- (normalized_id_dd_trd - min(id_dd_dup$log2_p_value)) / (max(id_dd_dup$log2_p_value) - min(id_dd_dup$log2_p_value))

id_dd_del$median_pos <- (id_dd_del$start + id_dd_del$end)/2
id_dd_del$log2_p_value = -log2(id_dd_del$`P-value`)
id_dd_del$scaled_p_value <- min_max_scale(id_dd_del$log2_p_value)
scaled_id_dd_del_trd <- (normalized_id_dd_trd - min(id_dd_del$log2_p_value)) / (max(id_dd_del$log2_p_value) - min(id_dd_del$log2_p_value))

###################################################SCZ####################################################################################
# SCZ TRD: 8.7E-05
scz_trd <- 8.7E-05
normalized_scz_trd <- -log2(scz_trd)
scz_dup$median_pos <- (scz_dup$start + scz_dup$end)/2
scz_dup$log2_p_value = -log2(scz_dup$`P-value`)
scz_dup$scaled_p_value <- min_max_scale(scz_dup$log2_p_value)
scaled_scz_dup_trd <- (normalized_scz_trd - min(scz_dup$log2_p_value)) / (max(scz_dup$log2_p_value) - min(scz_dup$log2_p_value))

scz_del$median_pos <- (scz_del$start + scz_del$end)/2
scz_del$log2_p_value = -log2(scz_del$`P-value`)
scz_del$scaled_p_value <- min_max_scale(scz_del$log2_p_value)
scaled_scz_del_trd <- (normalized_scz_trd - min(scz_del$log2_p_value)) / (max(scz_del$log2_p_value) - min(scz_del$log2_p_value))

###################################################PLOTING#############################################
# plot function for each hromosome 
plot_del <- function(chr_num, ymax1) {
  asd_del_prep <- asd_del[asd_del$chr==chr_num,]
  asd_id_dd_del_prep <- asd_id_dd_del[asd_id_dd_del$chr==chr_num,]
  id_dd_del_prep <- id_dd_del[id_dd_del$chr==chr_num,]
  scz_del_prep <- scz_del[scz_del$chr==chr_num,]
  
  # Standardise chromosome name 
  if (chr_num==23) {
    chr_num <- 'X'
  }
  chr_num_alter <- paste('chr', as.character(chr_num), sep = '')

  # Convert datasets into GRanges
  asd_del_data <- toGRanges(data.frame(chr=chr_num_alter, start=asd_del_prep$start, 
                                       end=asd_del_prep$end, y=asd_del_prep$scaled_p_value), genome = "hg19")
  asd_id_dd_del_data <- toGRanges(data.frame(chr=chr_num_alter, start=asd_id_dd_del_prep$start, 
                                             end=asd_id_dd_del_prep$end, y=asd_id_dd_del_prep$scaled_p_value), genome = "hg19")
  id_dd_del_data <- toGRanges(data.frame(chr=chr_num_alter, start=id_dd_del_prep$start, 
                                         end=id_dd_del_prep$end, y=id_dd_del_prep$scaled_p_value), genome = "hg19")
  scz_del_data <- toGRanges(data.frame(chr=chr_num_alter, start=scz_del_prep$start, 
                                       end=scz_del_prep$end, y=scz_del_prep$scaled_p_value), genome = "hg19")

  # Plotting
  chr <- plotKaryotype(genome = "hg19", plot.type = 1, chromosomes = chr_num_alter)
  kpAddBaseNumbers(chr, tick.dist=10000000, tick.len = 10, minor.tick.dist=1000000, cex = 1)
  kpDataBackground(chr, r0=0, r1=1.5, data.panel = 1, col="white")

  kpLines(chr, data = asd_del_data, r0=0, r1=1.5, data.panel = 1, ymin = 0, ymax = ymax1, col="blue3")
  kpLines(chr, data = asd_id_dd_del_data, r0=0, r1=1.5, data.panel = 1, ymin = 0, ymax = ymax1, col="#FF9900")
  kpLines(chr, data = id_dd_del_data, r0=0, r1=1.5, data.panel = 1, ymin = 0, ymax = ymax1, col="#333333")
  kpLines(chr, data = scz_del_data, r0=0, r1=1.5,  data.panel = 1, ymin = 0, ymax = ymax1, col="red2")
  
  # Add horizontal line for 95% percentile
  kpAbline(chr, chr = chr_num_alter, h = scaled_asd_del_trd, col = "blue3", data.panel = 1, lwd = 0.9, lty = 2)
  kpAbline(chr, chr = chr_num_alter, h = scaled_asd_id_dd_del_trd, col = "#FF9900", data.panel = 1, lwd = 0.9, lty = 2)
  kpAbline(chr, chr = chr_num_alter, h = scaled_id_dd_del_trd, col = "#333333", data.panel = 1, lwd = 0.9, lty = 2)
  kpAbline(chr, chr = chr_num_alter, h = scaled_scz_del_trd, col = "red2", data.panel = 1, lwd = 0.9, lty = 2)
  
  kpAxis(chr, data.panel=2, r0=0, r1=1.5, ymin = 0, ymax = ymax1)  
  
  #kpAddLabels(kp, labels="Deletions", r1=1.5, pos = 4, data.panel = 1, font=3)
  
  title(main = chr_num_alter, cex.main = 1.5, line = -1, adj = 0.05)  
}

# Sampe as Deletion
plot_dup <- function(chr_num, ymax1) {
 
  asd_del_prep <- asd_dup[asd_dup$chr==chr_num,]
  asd_id_dd_del_prep <- asd_id_dd_dup[asd_id_dd_dup$chr==chr_num,]
  id_dd_del_prep <- id_dd_dup[id_dd_dup$chr==chr_num,]
  scz_del_prep <- scz_dup[scz_dup$chr==chr_num,]
  
  # Standardise chromosome names  
  if (chr_num==23) {
    chr_num <- 'X'
  }
  chr_num_alter <- paste('chr', as.character(chr_num), sep = '')
  
  # Convert datasets into GRanges
  asd_del_data <- toGRanges(data.frame(chr=chr_num_alter, start=asd_del_prep$start, 
                                       end=asd_del_prep$end, y=asd_del_prep$scaled_p_value), genome = "hg19")
  asd_id_dd_del_data <- toGRanges(data.frame(chr=chr_num_alter, start=asd_id_dd_del_prep$start, 
                                             end=asd_id_dd_del_prep$end, y=asd_id_dd_del_prep$scaled_p_value), genome = "hg19")
  id_dd_del_data <- toGRanges(data.frame(chr=chr_num_alter, start=id_dd_del_prep$start, 
                                         end=id_dd_del_prep$end, y=id_dd_del_prep$scaled_p_value), genome = "hg19")
  scz_del_data <- toGRanges(data.frame(chr=chr_num_alter, start=scz_del_prep$start, 
                                       end=scz_del_prep$end, y=scz_del_prep$scaled_p_value), genome = "hg19")
 
  # Plotting
  chr <- plotKaryotype(genome = "hg19", plot.type = 1, chromosomes = chr_num_alter)
  kpAddBaseNumbers(chr, tick.dist=10000000, tick.len = 10, minor.tick.dist=1000000, cex = 1)
  kpDataBackground(chr, r0=0, r1=1.5, data.panel = 1, col="white")
  
  kpLines(chr, data = asd_del_data, r0=0, r1=1.5, data.panel = 1, ymin = 0, ymax = ymax1, col="blue3")
  kpLines(chr, data = asd_id_dd_del_data, r0=0, r1=1.5, data.panel = 1, ymin = 0, ymax = ymax1, col="#FF9900")
  kpLines(chr, data = id_dd_del_data, r0=0, r1=1.5, data.panel = 1, ymin = 0, ymax = ymax1, col="#333333")
  kpLines(chr, data = scz_del_data, r0=0, r1=1.5,  data.panel = 1, ymin = 0, ymax = ymax1, col="red2")
  
  # Add horizontal line for 95% percentile
  kpAbline(chr, chr = chr_num_alter, h = scaled_asd_dup_trd, col = "blue3", data.panel = 1, lwd = 0.9, lty = 2)
  kpAbline(chr, chr = chr_num_alter, h = scaled_asd_id_dd_dup_trd, col = "#FF9900", data.panel = 1, lwd = 0.9, lty = 2)
  kpAbline(chr, chr = chr_num_alter, h = scaled_id_dd_dup_trd, col = "#333333", data.panel = 1, lwd = 0.9, lty = 2)
  kpAbline(chr, chr = chr_num_alter, h = scaled_scz_dup_trd, col = "red2", data.panel = 1, lwd = 0.9, lty = 2)
  
  kpAxis(chr, data.panel=2, r0=0, r1=1.5, ymin = 0, ymax = ymax1)  
  
  #kpAddLabels(kp, labels="Deletions", r1=1.5, pos = 4, data.panel = 1, font=3)
  
  title(main = chr_num_alter, cex.main = 1.5, line = -1, adj = 0.05)  
}

# Plot each chromosome for Del
#* means normalised p-value is too small or large, need to change scales from 1 to 0.5/1.5 to make the plot visually appropriate
chr1_del <- as.ggplot(expression(plot_del(1, 1)))
chr2_del <- as.ggplot(expression(plot_del(2, 1)))
chr3_del <- as.ggplot(expression(plot_del(3, 0.5))) #*
chr4_del <- as.ggplot(expression(plot_del(4, 1)))
chr5_del <- as.ggplot(expression(plot_del(5, 1)))
chr6_del <- as.ggplot(expression(plot_del(6, 1)))
chr7_del <- as.ggplot(expression(plot_del(7,1)))
chr8_del <- as.ggplot(expression(plot_del(8,1)))
chr9_del <- as.ggplot(expression(plot_del(9,0.5))) #*
chr10_del <- as.ggplot(expression(plot_del(10,1)))
chr11_del <- as.ggplot(expression(plot_del(11,1)))
chr12_del <- as.ggplot(expression(plot_del(12,1)))
chr13_del <- as.ggplot(expression(plot_del(13,0.5))) #*
chr14_del <- as.ggplot(expression(plot_del(14,0.5))) #*
chr15_del <- as.ggplot(expression(plot_del(15,1)))
chr16_del <- as.ggplot(expression(plot_del(16,1)))
chr17_del <- as.ggplot(expression(plot_del(17,1)))
chr18_del <- as.ggplot(expression(plot_del(18,1)))
chr19_del <- as.ggplot(expression(plot_del(19,0.5))) #*
chr20_del <- as.ggplot(expression(plot_del(20,0.5))) #*
chr21_del <- as.ggplot(expression(plot_del(21,1)))
chr22_del <- as.ggplot(expression(plot_del(22,1)))
chr23_del <- as.ggplot(expression(plot_del(23,1)))

# Plot each chromosome for Dup
#* means normalised p-value is too small or large, need to change scales from 1 to 0.5/1.5 to make the plot visually appropriate
chr1_dup <- as.ggplot(expression(plot_dup(1, 1)))
chr2_dup <- as.ggplot(expression(plot_dup(2, 1)))
chr3_dup <- as.ggplot(expression(plot_dup(3, 0.5))) #*
chr4_dup <- as.ggplot(expression(plot_dup(4, 1.5))) #*
chr5_dup <- as.ggplot(expression(plot_dup(5, 0.5))) #*
chr6_dup <- as.ggplot(expression(plot_dup(6, 0.5))) #*
chr7_dup <- as.ggplot(expression(plot_dup(7, 0.5))) #*
chr8_dup <- as.ggplot(expression(plot_dup(8,1)))
chr9_dup <- as.ggplot(expression(plot_dup(9,0.5))) #*
chr10_dup <- as.ggplot(expression(plot_dup(10,1)))
chr11_dup <- as.ggplot(expression(plot_dup(11,1)))
chr12_dup <- as.ggplot(expression(plot_dup(12,0.5))) #*
chr13_dup <- as.ggplot(expression(plot_dup(13,0.25))) #**
chr14_dup <- as.ggplot(expression(plot_dup(14,0.5))) #*
chr15_dup <- as.ggplot(expression(plot_dup(15,1.5))) #*
chr16_dup <- as.ggplot(expression(plot_dup(16,1)))
chr17_dup <- as.ggplot(expression(plot_dup(17,0.5))) #*
chr18_dup <- as.ggplot(expression(plot_dup(18,0.5))) #*
chr19_dup <- as.ggplot(expression(plot_dup(19,0.5))) #*
chr20_dup <- as.ggplot(expression(plot_dup(20,0.5))) #*
chr21_dup <- as.ggplot(expression(plot_dup(21,0.5))) #*
chr22_dup <- as.ggplot(expression(plot_dup(22,1.5))) #*
chr23_dup <- as.ggplot(expression(plot_dup(23,0.5))) #*

# Arrage the plots 
combine_plots_del <- function() {
  plotlist <- list(chr1_del, chr2_del, chr3_del,chr4_del,chr5_del,chr6_del,chr7_del,chr8_del,chr9_del,chr10_del,chr11_del,
                   chr12_del,chr13_del,chr14_del,chr15_del,chr16_del,chr17_del,chr18_del,chr19_del,chr20_del,chr21_del,chr22_del,chr23_del)
  chroms <- plot_grid(plotlist = plotlist, ncol=4, nrow = 6)
  #legend(x = "bottomright", fill = c("#E41A1C", "#377EB8", "#FF7F00", "#984EA3"), cex=1.2, legend = c("ASD", "ASD_ID_DD","ID_DD","SCZ"), box.lty=0)
}

combine_plots_dup <- function() {
  plotlist <- list(chr1_dup, chr2_dup, chr3_dup,chr4_dup,chr5_dup,chr6_dup,chr7_dup,chr8_dup,chr9_dup,chr10_dup,chr11_dup,
                   chr12_dup,chr13_dup,chr14_dup,chr15_dup,chr16_dup,chr17_dup,chr18_dup,chr19_dup,chr20_dup,chr21_dup,chr22_dup,chr23_dup)
  plot_grid(plotlist = plotlist, ncol=4, nrow = 6)
  #legend(x = "bottomright", fill = c("#E41A1C", "#377EB8", "#FF7F00", "#984EA3"), cex=1.2, legend = c("ASD", "ASD_ID_DD","ID_DD","SCZ"), box.lty=0)
}

arranged_chroms <- combine_plots_del()
arrange_chroms_dup <- combine_plots_dup()

# Save Plots
save_plot("arranged_chroms.svg", arranged_chroms, base_height = 500, base_width = 809, unit = "mm")
save_plot("arranged_chroms_dup_.svg", arrange_chroms_dup, base_height = 500, base_width = 809, unit = "mm")
