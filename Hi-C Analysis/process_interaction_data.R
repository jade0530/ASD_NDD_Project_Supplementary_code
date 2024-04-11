library(tidyr)
library(dplyr)
library(readr)
library(readxl)
library(regioneR)
library(stringr)
library(karyoploteR)
library(GenomicRanges)
library(Repitools)
library(ggplot2)
library(ggpubr)

pth <- "/media/bml/USER_DATA/ASD_CNV_Project/interaction_results/"
files <- list.files(path = "/media/bml/USER_DATA/ASD_CNV_Project/interaction_results/")
gene_overlap_10_percent <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/Hi-C_INTERACTION/gene_overlap_result_merged_10_overlap.csv")
gene_overlap_10_percent_small_reg <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/Hi-C_INTERACTION/gene_overlap_result_merged_10_overlap_small_reg.csv")

results_merged <- data.frame()
for (i in files) {
  path <- paste0(pth, i)
  print(path)
  results_merged <- rbind(results_merged, read_csv(file =  path))
}


##########################FUNCTIONS#########################
calculate_metrics <- function(phenotype_result_table) {
  # First separate by each paper group 
  phenotype_result_table_syn <- filter(phenotype_result_table, grepl("Synapse",paper))
  phenotype_result_table_7015 <- filter(phenotype_result_table, grepl("31367015",paper))
  phenotype_result_table_5922 <- filter(phenotype_result_table, grepl("30555922",paper))
  phenotype_result_table_0116 <- filter(phenotype_result_table, grepl("27760116",paper))
  
  # We need total number of reads, total number of interactions 
  phenotype_result_table_syn_general <- group_each_paper_results(phenotype_result_table_syn)
  phenotype_result_table_syn_gene_distance <- calculate_gene_distance(phenotype_result_table_syn, "Synapse")
  phenotype_result_table_syn_grouped <- merge(phenotype_result_table_syn_gene_distance, phenotype_result_table_syn_general, 
                                              by =  c('chr', 'cnv_start', 'cnv_end', 'cnv_type', "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype"), 
                                              all = TRUE)
  print(head(phenotype_result_table_syn_grouped))
  colnames(phenotype_result_table_syn_grouped) =c('chr', 'cnv_start', 'cnv_end', 'cnv_type', "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype", 
                                                  'overlapped_gene_symbol_Synapse', 
                                                  'overlapped_promoter_gene_symbol_Synapse',
                                                  'interactions_Synapse', 'reads_Synapse',
                                                  "hic_distance_Synapse")
  
  phenotype_result_table_7015_grouped <- group_each_paper_results(phenotype_result_table_7015)
 
  phenotype_result_table_7015_gene_distance <- calculate_gene_distance(phenotype_result_table_7015, "7015")
  phenotype_result_table_7015_grouped <- merge(phenotype_result_table_7015_gene_distance, phenotype_result_table_7015_grouped,  
                                               by =  c('chr', 'cnv_start', 'cnv_end', 'cnv_type', "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype"),
                                               all = TRUE)
  
  colnames(phenotype_result_table_7015_grouped) =c('chr', 'cnv_start', 'cnv_end', 'cnv_type', "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype",
                                                   'overlapped_gene_symbol_31367015',
                                                   'overlapped_promoter_gene_symbol_31367015',
                                                   'interactions_31367015', 'reads_31367015',
                                                   "hic_distance_31367015")
  
  phenotype_result_table_5922_grouped <- group_each_paper_results(phenotype_result_table_5922)
  phenotype_result_table_5922_gene_distance <- calculate_gene_distance(phenotype_result_table_5922, "5922")
  phenotype_result_table_5922_grouped <- merge(phenotype_result_table_5922_gene_distance, phenotype_result_table_5922_grouped, 
                                               by =  c('chr', 'cnv_start', 'cnv_end', 'cnv_type', "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype"),
                                               all = TRUE)
  colnames(phenotype_result_table_5922_grouped) =c('chr', 'cnv_start', 'cnv_end', 'cnv_type', "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype", 
                                                   'overlapped_gene_symbol_30555922',
                                                   'overlapped_promoter_gene_symbol_30555922',
                                                   'interactions_30555922', 'reads_30555922', 
                                                   "hic_distance_30555922")

  phenotype_result_table_0116_grouped <- group_each_paper_results(phenotype_result_table_0116)
  phenotype_result_table_0116_gene_distance <- calculate_gene_distance(phenotype_result_table_0116, "0116")
  
  phenotype_result_table_0116_grouped <- merge(phenotype_result_table_0116_gene_distance, phenotype_result_table_0116_grouped, 
                                               by = c('chr', 'cnv_start', 'cnv_end', 'cnv_type', "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype"))
  
  colnames(phenotype_result_table_0116_grouped) =c('chr', 'cnv_start', 'cnv_end', 'cnv_type', "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype", 
                                                   'overlapped_gene_symbol_27760116', 
                                                   'overlapped_promoter_gene_symbol_27760116',
                                                   'interactions_27760116', 'reads_27760116',
                                                   "hic_distance_27760116")
  
  # We need number of interactions for each paper & fraction=nInteractions/totalInteractions
  phenotype_result_table_paper_list <- list(phenotype_result_table_syn_grouped,
                                            phenotype_result_table_7015_grouped,
                                            phenotype_result_table_5922_grouped,
                                            phenotype_result_table_0116_grouped)
  #print(phenotype_result_table_syn_grouped)
  #merge all data frames in list
  phenotype_result_table_paper_all <- Reduce(function(x, y) merge(x, y, by = c('chr', 'cnv_start', 'cnv_end', 'cnv_type', "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype"), 
                                                                  all=TRUE), 
                                             phenotype_result_table_paper_list)
  phenotype_result_table_paper_all$overlapped_genes_symbol_unique <- paste(phenotype_result_table_paper_all$overlapped_gene_symbol_Synapse,
                                                                            phenotype_result_table_paper_all$overlapped_gene_symbol_31367015,
                                                                            phenotype_result_table_paper_all$overlapped_gene_symbol_30555922,
                                                                            phenotype_result_table_paper_all$overlapped_gene_symbol_27760116, sep = ";")
 
  phenotype_result_table_paper_all$overlapped_promoter_genes_symbol_unique <- paste(phenotype_result_table_paper_all$overlapped_promoter_gene_symbol_Synapse,
                                                                            phenotype_result_table_paper_all$overlapped_promoter_gene_symbol_31367015,
                                                                            phenotype_result_table_paper_all$overlapped_promoter_gene_symbol_30555922,
                                                                            phenotype_result_table_paper_all$overlapped_promoter_gene_symbol_27760116, sep = ";")
  
  # phenotype_result_table_paper_all$gene_distance_unique <- paste(phenotype_result_table_paper_all$gene_distance_Synapse,
  #                                                                  phenotype_result_table_paper_all$gene_distance_31367015,
  #                                                                  phenotype_result_table_paper_all$gene_distance_30555922,
  #                                                                  phenotype_result_table_paper_all$gene_distance_27760116, sep = ";")
  #  
  # phenotype_result_table_paper_all$promoter_distance_unique <- paste(phenotype_result_table_paper_all$promoter_distance_Synapse,
  #                                                                      phenotype_result_table_paper_all$promoter_distance_31367015,
  #                                                                      phenotype_result_table_paper_all$promoter_distance_30555922,
  #                                                                      phenotype_result_table_paper_all$promoter_distance_27760116, sep =";")
  phenotype_result_table_paper_all$hic_distance_unique <- paste(phenotype_result_table_paper_all$hic_distance_Synapse,
                                                                     phenotype_result_table_paper_all$hic_distance_31367015,
                                                                     phenotype_result_table_paper_all$hic_distance_30555922,
                                                                     phenotype_result_table_paper_all$hic_distance_27760116, sep =";")
  phenotype_result_table_paper_all <- calculate_fractions(phenotype_result_table_paper_all)

  return(phenotype_result_table_paper_all)
}

calculate_gene_distance <- function(phenotype_result_table, paper_name) {
  # Then separate all genes with the distance 
  print(colnames(phenotype_result_table))
  phenotype_result_table_separated_gene_body <- 
    phenotype_result_table %>%
    separate_longer_delim(c(overlapped_gene_symbol
    ), delim = ';')
  phenotype_result_table_separated_gene_body <- phenotype_result_table_separated_gene_body[,c('chr', 'cnv_start', 'cnv_end', 'cnv_type', "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype", 
                                                                                             "overlapped_gene_symbol")]
  
  phenotype_result_table_separated_promoter <- 
    phenotype_result_table %>%
    separate_longer_delim(c(overlapped_promoter_gene_symbol)
                          ,delim = ';')
  phenotype_result_table_separated_promoter <- phenotype_result_table_separated_promoter[,c('chr', 'cnv_start', 'cnv_end', 'cnv_type', "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype", 
                                                                                              "overlapped_promoter_gene_symbol")]
  # 
  # # Calculate gene distance: if gene_end_cnv_start and gene_start_cnv_end both < 0 or > 0, take the smaller abs value;
  # # else, gene distance is 0. 5 cases
  # phenotype_result_table_separated_gene_body$gene_distance = 'empty'
  # phenotype_result_table_separated_promoter$promoter_distance = 'empty'
  # 
  # for (i in 1:nrow(phenotype_result_table_separated_gene_body)) {
  #   #print(phenotype_result_table_separated[i,])
  #   gsce = phenotype_result_table_separated_gene_body[i,]$overlapped_gene_start_cnv_end
  #   gecs = phenotype_result_table_separated_gene_body[i,]$overlapped_gene_end_cnv_start
  # 
  #   if (gsce == '' && gecs == '') {
  #     #print("ignored")
  #     phenotype_result_table_separated_gene_body[i,]$gene_distance = ""
  #   } else if (as.numeric(gsce) > 0 && as.numeric(gecs) > 0) {
  #     # right
  #     #print("right")
  #     phenotype_result_table_separated_gene_body[i,]$gene_distance = as.character(gsce)
  #   } else if (as.numeric(gsce) > 0 && as.numeric(gecs) > 0) {
  #     #print("left")
  #     phenotype_result_table_separated_gene_body[i,]$gene_distance = as.character(abs(gecs))
  #   } else {
  #     #print("overlap")
  #     phenotype_result_table_separated_gene_body[i,]$gene_distance = "0"
  #   }
  # }
  # 
  # for (i in 1:nrow(phenotype_result_table_separated_promoter)) {
  #   gsce_p = phenotype_result_table_separated_promoter[i,]$overlapped_promoter_gene_start_cnv_end
  #   gecs_p = phenotype_result_table_separated_promoter[i,]$overlapped_promoter_gene_end_cnv_start
  #   if (gsce_p == '' && gecs_p == '') {
  #     phenotype_result_table_separated_promoter[i,]$promoter_distance = ""
  #   } else if (as.numeric(gsce_p) > 0 && as.numeric(gecs_p) > 0) {
  #     # right
  #     phenotype_result_table_separated_promoter[i,]$promoter_distance = as.character(gsce_p)
  #   } else if (as.numeric(gsce_p) < 0 && as.numeric(gecs_p) < 0) {
  #     phenotype_result_table_separated_promoter[i,]$promoter_distance = as.character(abs(as.numeric(gecs_p)))
  #   } else {
  #     phenotype_result_table_separated_promoter[i,]$promoter_distance = "0"
  #   }
  # }
  #phenotype_result_table_separated_gene_body <- phenotype_result_table_separated_gene_body[,c("chr", "cnv_start", "cnv_end", "cnv_type", "no_of_coding_genes", "no_of_noncoding_genes", 
  #                                                                                           "non_coding_cnv", "phenotype", "overlapped_gene_symbol")]

  #phenotype_result_table_separated_promoter <- phenotype_result_table_separated_promoter[,c("chr", "cnv_start", "cnv_end", "cnv_type","no_of_coding_genes", "no_of_noncoding_genes", 
  #                                                                                           "non_coding_cnv", "phenotype", "overlapped_promoter_gene_symbol")]
  #phenotype_result_table_separated_gene_body <- distinct(phenotype_result_table_separated_gene_body)
  #phenotype_result_table_separated_promoter <- distinct(phenotype_result_table_separated_promoter)
  #filename = paste0(paper_name, "_distances.csv")
  #output_raw_file <- merge(phenotype_result_table_separated_gene_body, phenotype_result_table_separated_promoter, by = 
  #                                          c("chr", "cnv_start", "cnv_end", "cnv_type", "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "phenotype"), all = TRUE)
  
  #write_csv(output_raw_file, filename)
  
  phenotype_result_table_grouped_gene_body <- group_each_paper_results_gene_symbol(phenotype_result_table_separated_gene_body)
  phenotype_result_table_grouped_promoter <- group_each_paper_results_promoter_symbol(phenotype_result_table_separated_promoter)
  phenotype_result_table_grouped <- merge(phenotype_result_table_grouped_gene_body, phenotype_result_table_grouped_promoter, by = 
                                            c("chr", "cnv_start", "cnv_end", "cnv_type", "no_of_coding_genes", "no_of_noncoding_genes", 
                                              "non_coding_cnv", "phenotype"), all = TRUE)
  return(phenotype_result_table_grouped)
}
calculate_fractions <- function(phenotype_result_table_paper_merged) {
  phenotype_result_table_paper_merged$interactions_Synapse[is.na(phenotype_result_table_paper_merged$interactions_Synapse)] = 0 
  phenotype_result_table_paper_merged$interactions_31367015[is.na(phenotype_result_table_paper_merged$`interactions_31367015`)] = 0 
  phenotype_result_table_paper_merged$interactions_30555922[is.na(phenotype_result_table_paper_merged$`interactions_30555922`)] = 0 
  phenotype_result_table_paper_merged$interactions_27760116[is.na(phenotype_result_table_paper_merged$`interactions_27760116`)] = 0 
  
  phenotype_result_table_paper_merged$reads_Synapse[is.na(phenotype_result_table_paper_merged$reads_Synapse)] = 0 
  phenotype_result_table_paper_merged$reads_31367015[is.na(phenotype_result_table_paper_merged$reads_31367015)] = 0 
  phenotype_result_table_paper_merged$reads_30555922[is.na(phenotype_result_table_paper_merged$reads_30555922)] = 0 
  phenotype_result_table_paper_merged$reads_27760116[is.na(phenotype_result_table_paper_merged$reads_27760116)] = 0
  phenotype_result_table_paper_merged$`no_of_interactions` <- phenotype_result_table_paper_merged$interactions_Synapse +
    phenotype_result_table_paper_merged$interactions_31367015 + phenotype_result_table_paper_merged$interactions_30555922 + 
    phenotype_result_table_paper_merged$interactions_27760116
  
  phenotype_result_table_paper_merged$frac_Synapse = phenotype_result_table_paper_merged$interactions_Synapse/phenotype_result_table_paper_merged$`no_of_interactions`
  phenotype_result_table_paper_merged$frac_31367015 = phenotype_result_table_paper_merged$interactions_31367015/phenotype_result_table_paper_merged$`no_of_interactions`
  phenotype_result_table_paper_merged$frac_30555922 = phenotype_result_table_paper_merged$interactions_30555922/phenotype_result_table_paper_merged$`no_of_interactions`
  phenotype_result_table_paper_merged$frac_27760116 = phenotype_result_table_paper_merged$interactions_27760116/phenotype_result_table_paper_merged$`no_of_interactions`
  
  phenotype_result_table_paper_merged$`no_of_reads` = phenotype_result_table_paper_merged$reads_Synapse +
    phenotype_result_table_paper_merged$`reads_31367015` + phenotype_result_table_paper_merged$`reads_30555922`+ phenotype_result_table_paper_merged$`reads_27760116`
  
  phenotype_result_table_paper_merged$read_frac_Synapse = phenotype_result_table_paper_merged$reads_Synapse/phenotype_result_table_paper_merged$`no_of_reads`
  phenotype_result_table_paper_merged$read_frac_31367015 = phenotype_result_table_paper_merged$reads_31367015/phenotype_result_table_paper_merged$`no_of_reads`
  phenotype_result_table_paper_merged$read_frac_30555922 = phenotype_result_table_paper_merged$reads_30555922/phenotype_result_table_paper_merged$`no_of_reads`
  phenotype_result_table_paper_merged$read_frac_27760116 = phenotype_result_table_paper_merged$reads_27760116/phenotype_result_table_paper_merged$`no_of_reads`

  return(phenotype_result_table_paper_merged)
}

group_each_paper_results <- function(paper_group) {
  paper_group_grouped <- 
    paper_group %>%
    group_by(chr, cnv_start, cnv_end, cnv_type, no_of_coding_genes, no_of_noncoding_genes, non_coding_cnv, phenotype) %>%
    summarise(`interactions` = n(),
              `no_of_reads` = sum(as.numeric(reads)),
              #overlapped_gene_symbol = str_c(na.omit(overlapped_gene_symbol), collapse=";"),
              #overlapped_promoter_gene_symbol = str_c(na.omit(unique(overlapped_promoter_gene_symbol)), collapse=";"),
              interaction_distance = str_c(na.omit(interaction_distance), collapse=";"))
  return(paper_group_grouped)
}

group_each_paper_results_gene_distance <- function(paper_group) {
  paper_group_grouped <- 
    paper_group %>%
    group_by(chr, cnv_start, cnv_end, cnv_type, no_of_coding_genes, no_of_noncoding_genes, non_coding_cnv, phenotype) %>%
    summarise(
              overlapped_gene_symbol = str_c(na.omit(overlapped_gene_symbol), collapse=" "),
              gene_distance =  str_c(na.omit(gene_distance), collapse=" "))
  return(paper_group_grouped)
}

group_each_paper_results_promoter_distance <- function(paper_group) {
  paper_group_grouped <- 
    paper_group %>%
    group_by(chr, cnv_start, cnv_end, cnv_type, no_of_coding_genes, no_of_noncoding_genes, non_coding_cnv, phenotype) %>%
    summarise(
      overlapped_promoter_gene_symbol = str_c(na.omit(overlapped_promoter_gene_symbol), collapse=" "),        
      promoter_distance =  str_c(na.omit(promoter_distance), collapse=" "))
  return(paper_group_grouped)
}

group_each_paper_results_gene_symbol <- function(paper_group) {
  paper_group_grouped <- 
    paper_group %>%
    group_by(chr, cnv_start, cnv_end, cnv_type, no_of_coding_genes, no_of_noncoding_genes, non_coding_cnv, phenotype) %>%
    summarise(
      overlapped_gene_symbol = str_c(na.omit(unique(overlapped_gene_symbol)), collapse=";"))
  return(paper_group_grouped)
}

group_each_paper_results_promoter_symbol <- function(paper_group) {
  paper_group_grouped <- 
    paper_group %>%
    group_by(chr, cnv_start, cnv_end, cnv_type, no_of_coding_genes, no_of_noncoding_genes, non_coding_cnv, phenotype) %>%
    summarise(
      overlapped_promoter_gene_symbol = str_c(na.omit(unique(overlapped_promoter_gene_symbol)), collapse=";"))
  return(paper_group_grouped)
}

get_final_table <- function(cnvs_with_results) {
  # final format: cnv_coord, noncoding genes, interacting coding genes (promoters), 
  # Hi-C interactions,	# Reads,	Distance with interacting genes,Phenotypes
  # We need to do this through the original cnv list 
  
  asd_cnvs <- asd_cnv_data[,c("chr", "start", "end", "cnv_type")]
  asd_cnvs$phenotype = "ASD"
  asd_id_dd_cnvs <- asd_id_dd_cnv_data[,c("chr", "start", "end", "cnv_type")]
  asd_id_dd_cnvs$phenotype = "ASD_ID_DD"
  id_dd_cnvs <- id_dd_cnv_data[,c("chr", "start", "end", "cnv_type")]
  id_dd_cnvs$phenotype = "ID_DD"
  scz_cnvs <- scz_cnv_data[,c("chr", "start", "end", "cnv_type")]
  scz_cnvs$phenotype = "SCZ"
  total_cnvs <- rbind(asd_cnvs, asd_id_dd_cnvs, id_dd_cnvs, scz_cnvs)
  total_cnvs$cnv_coord <- paste0(total_cnvs$chr, ':',
                                 total_cnvs$start, '-',
                                 total_cnvs$end)
  
  total_cnvs <- total_cnvs[,c('cnv_coord', 'cnv_type', 'phenotype')]
  
  gene_overlap_result <- gene_overlap_10_percent_small_reg
  gene_overlap_result$cnv_coord <- paste0(gene_overlap_result$cnv_chr, ':',
                                          gene_overlap_result$cnv_start, '-',
                                          gene_overlap_result$cnv_end)
  gene_overlap_result <- gene_overlap_result[,c('cnv_coord', 'phenotype', "#_coding_genes", "#_noncoding_genes")]
  
  total_cnv_with_genes <- merge(total_cnvs, gene_overlap_result, by.x = c("cnv_coord","phenotype"),
                                by.y = c("cnv_coord","phenotype"), all = TRUE)
  
  cnvs_with_results$cnv_coord <- paste0(cnvs_with_results$chr, ':',
                                        cnvs_with_results$cnv_start, '-',
                                        cnvs_with_results$cnv_end)
  cnvs_with_results <- cnvs_with_results[,c("cnv_coord", "non_coding_cnv", "phenotype", "no_of_reads", "no_of_interactions", 
                                            "overlapped_genes_symbol_unique", "overlapped_promoter_genes_symbol_unique", 
                                            "hic_distance_unique")]

  
  # Map all the things from cnvs_with_results to total_cnvs
  cnvs_with_results_merged <- merge(total_cnv_with_genes, cnvs_with_results, by.x = c("cnv_coord","phenotype"),
                                    by.y = c("cnv_coord","phenotype"), all = TRUE)
  
  return(cnvs_with_results_merged)  
}

##############################ANALYSIS#############################
results_merged <- distinct(results_merged)

# Map new set of coding/noncoding genes to the dataset 
gene_overlap_10_percent_small_reg_process <- gene_overlap_10_percent_small_reg[,c("cnv_chr", "cnv_start", "cnv_end", "#_coding_genes", "#_noncoding_genes")]
results_merged_gene_overlap_smaller_reg <- merge(results_merged, gene_overlap_10_percent_small_reg_process, 
                                                 by.x = c("chr", "cnv_start", "cnv_end"), by.y = c("cnv_chr", "cnv_start", "cnv_end"), all.x=TRUE)

# Delete the old ones 
results_merged_gene_overlap_smaller_reg <- results_merged_gene_overlap_smaller_reg[,!(names(results_merged_gene_overlap_smaller_reg) %in% c("no_of_coding_genes", "no_of_noncoding_genes"))]

# Rename new ones 
colnames(results_merged_gene_overlap_smaller_reg)[25] <- "no_of_coding_genes"
colnames(results_merged_gene_overlap_smaller_reg)[26] <- "no_of_noncoding_genes"

# Recalculate noncoding CNVs
results_merged_clean <- lapply(results_merged_gene_overlap_smaller_reg, 
                               function(x) gsub("ignore,", "", x))
results_merged_clean <- lapply(results_merged_clean, 
                               function(x) gsub("ignore", "", x))
results_merged_clean <- lapply(results_merged_clean, 
                               function(x) gsub(",", ";", x))
results_merged_clean_df <- data.frame(results_merged_clean)
results_merged_clean_df$no_of_coding_genes <- replace(results_merged_clean_df$no_of_coding_genes, is.na(results_merged_clean_df$no_of_coding_genes), 0)
results_merged_clean_df$no_of_noncoding_genes <- replace(results_merged_clean_df$no_of_noncoding_genes, is.na(results_merged_clean_df$no_of_noncoding_genes), 0)

results_merged_clean_df$non_coding_cnv <- (results_merged_clean_df$no_of_coding_genes==0)


# calculate distance between interactions 
results_merged_clean_df$interaction_distance = as.numeric(results_merged_clean_df$evol_locus_start) - as.numeric(results_merged_clean_df$interacting_locus_end)
write_csv(results_merged_clean_df, "results_merged_clean_df_small_reg.csv")

# Count only overlaps with only one side of HiC
results_merged_one_side <- results_merged_clean_df[results_merged_clean_df$overlap_with_both_hic_side==FALSE,]

results_merged_count_one_side <- calculate_metrics(results_merged_one_side)
# Clean the output 
results_merged_count_one_side_clean <- lapply(results_merged_count_one_side, 
                                              function(x) gsub("NA", "", x))
results_merged_count_one_side_clean[is.na(results_merged_count_one_side_clean)] = ""
results_merged_count_one_side_clean <- lapply(results_merged_count_one_side_clean, 
                                              function(x) gsub(";{2,}", ";", x))
results_merged_count_one_side_clean <- lapply(results_merged_count_one_side_clean, 
                                              function(x) gsub("^;+|;+$","",x))
results_merged_count_one_side_clean <- lapply(results_merged_count_one_side_clean, 
                                              function(x) trimws(x, which = "both", whitespace = "[ ]"))

results_merged_count_one_side_clean <- data.frame(results_merged_count_one_side_clean)

write_csv(results_merged_count_one_side_clean, "hic_interactions_count_one_side_small_reg.csv")



final_table <- get_final_table(results_merged_count_one_side_clean)

# Data cleaning
final_table_clean <- lapply(final_table, 
                            function(x) gsub("NA", "", x))
final_table_clean <- lapply(final_table_clean, 
                            function(x) gsub("^;+|;+$","",x))
final_table_clean <- data.frame(final_table_clean)
colnames(final_table_clean)[4] <- "#_coding_genes"
final_table_clean$`#_coding_genes`[is.na(final_table_clean$`#_coding_genes`)] = 0
colnames(final_table_clean)[5] <- "#_noncoding_genes"
final_table_clean$`#_noncoding_genes`[is.na(final_table_clean$`#_noncoding_genes`)] = 0
final_table_clean$no_of_reads[is.na(final_table_clean$no_of_reads)] = 0
final_table_clean$no_of_interactions[is.na(final_table_clean$no_of_interactions)] = 0
final_table_clean$non_coding_cnv <- (final_table_clean$`#_coding_genes`==0)
final_table_clean[is.na(final_table_clean)] = ''

# Add annotated columns
final_table_clean$gene_overlap = NA
final_table_clean$with_interactions = NA
final_table_clean$with_interactions <- final_table$no_of_interactions>0
final_table_clean$with_interactions[is.na(final_table_clean$with_interactions)] = FALSE


final_table_clean$gene_overlap[final_table_clean$`#_coding_genes`> 0] = 'CNVs overlap with protein coding genes'
final_table_clean$gene_overlap[which(final_table_clean$`#_coding_genes`==0 & final_table_clean$`#_noncoding_genes`> 0)] = 'CNVs overlap with non-coding genes'
final_table_clean$gene_overlap[which(final_table_clean$`#_coding_genes`==0 & final_table_clean$`#_noncoding_genes`==0)] = 'CNVs overlap with non-annotated regions'

final_table_clean <- distinct(final_table_clean)
write_csv(final_table_clean, "cnv_table_annotated_cleaned_small_reg.csv")

# Start annotate interaction table 
results_merged_one_side_process <- results_merged_one_side
# Select all non-coding cnvs
results_merged_one_side_process <- results_merged_one_side_process[results_merged_one_side_process$non_coding_cnv==TRUE,]
results_merged_one_side_process$hic_interactions = NA
results_merged_one_side_process$interact_with_gene <- !results_merged_one_side_process$overlapped_gene_symbol==''
results_merged_one_side_process$interact_with_gene_promoter <- !results_merged_one_side_process$overlapped_promoter_gene_symbol==''
results_merged_one_side_process$hic_interactions[which(results_merged_one_side_process$interact_with_gene_promoter==TRUE)] = 'Interaction of ncCNVs with promoters of genes (coding & noncoding)'
results_merged_one_side_process$hic_interactions[which(results_merged_one_side_process$interact_with_gene==TRUE & results_merged_one_side_process$interact_with_gene_promoter==FALSE)] = 'Interaction of ncCNVs with genes body'
results_merged_one_side_process$hic_interactions[which(results_merged_one_side_process$interact_with_gene==FALSE & results_merged_one_side_process$interact_with_gene_promoter==FALSE)] = 'Interaction of ncCNVs with non-annotated regions'
write_csv(results_merged_one_side_process, "hic_interactions_ncCNVs_table_small_reg.csv")

# # Generate ncCNVs gene distance table - CNVID | GENE NAME | GENE DISTANCE
# # These tables are for Mona to do better analysis 
# nccnv_gene_distance_table <- results_merged_count_one_side[,c("chr", "cnv_start", "cnv_end", "overlapped_genes_symbol_unique")]
# 
# # separate gene names + gene distance
# nccnv_gene_distance_table_separate <- 
#   nccnv_gene_distance_table %>%
#   separate_longer_delim(c(overlapped_genes_symbol_unique,
#                           gene_distance_unique
#   ), delim = ' ') %>%
#   filter(overlapped_genes_symbol_unique != '' & overlapped_genes_symbol_unique != 'NA') %>%
#   distinct()
# 
# write_csv(nccnv_gene_distance_table_separate, "ncCNVs_gene_distances.csv")
# 
# # Same as ncCNVs promoter of genes distance table - CNVID | GENE NAME | PROMOTER DISTANCE
# nccnv_promoter_distance_table <- results_merged_count_one_side[,c("chr", "cnv_start", "cnv_end", "overlapped_promoter_genes_symbol_unique", "promoter_distance_unique")]
# 
# nccnv_promoter_distance_table_separate <- 
#   nccnv_promoter_distance_table %>%
#   separate_longer_delim(c(overlapped_promoter_genes_symbol_unique,
#                           promoter_distance_unique
#   ), delim = ' ') %>%
#   filter(overlapped_promoter_genes_symbol_unique != '' & overlapped_promoter_genes_symbol_unique != 'NA') %>%
#   distinct()
# write_csv(nccnv_promoter_distance_table_separate, "ncCNVs_promoter_distances.csv")

# Prepare for the pie charts
final_table_asd <- final_table_clean[final_table_clean$phenotype=="ASD",]
final_interaction_asd <- results_merged_one_side_process[results_merged_one_side_process$phenotype=="ASD",]
final_interaction_asd <- final_interaction_asd[final_interaction_asd$non_coding_cnv==TRUE,]

final_table_cnvs_asd <- final_table_asd %>%
  group_by(gene_overlap) %>%
  summarise(
    value = n(),
    percentage = 100*n()/nrow(final_table_asd)
  )

final_table_int_asd <- final_interaction_asd %>%
  group_by(hic_interactions) %>%
  summarise(
    value = n(),
    percentage = 100*n()/nrow(final_interaction_asd)
  )

# Pie Chart
asd_cnvs_pie <- ggplot(final_table_cnvs_asd, aes(x="", y=percentage, fill=gene_overlap)) +
  geom_bar(stat="identity", width=1) + 
  scale_fill_manual(values = c("tomato2", "gold", "lemonchiffon")) + 
  geom_text(aes(label = round(percentage, digits = 2)),
            position = position_stack(vjust = 0.5),
            color = 'black',
            size = 7) +
  coord_polar("y", start=0) +
  theme_void() 

asd_int_pie <- ggplot(final_table_int_asd, aes(x="", percentage, fill=hic_interactions)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = c("cornflowerblue", "lightgreen", "mediumpurple1")) + 
  geom_text(aes(label = round(percentage, digits = 2)),
            position = position_stack(vjust = 0.5),
            color = 'black',
            size = 7) +
  coord_polar("y", start=0) +
  theme_void() 

final_table_asd_id_dd <- final_table_clean[final_table_clean$phenotype=="ASD_ID_DD",]
final_interaction_asd_id_dd <- results_merged_one_side_process[results_merged_one_side_process$phenotype=="ASD_ID_DD",]
final_interaction_asd_id_dd <- final_interaction_asd_id_dd[final_interaction_asd_id_dd$non_coding_cnv==TRUE,]

final_table_cnvs_asd_id_dd <- final_table_asd_id_dd %>%
  group_by(gene_overlap) %>%
  summarise(
    value = n(),
    percentage = 100*n()/nrow(final_table_asd_id_dd)
  )

final_table_int_asd_id_dd <- final_interaction_asd_id_dd %>%
  group_by(hic_interactions) %>%
  summarise(
    value = n(),
    percentage = 100*n()/nrow(final_interaction_asd_id_dd)
  )

# Pie Chart
asd_id_dd_cnvs_pie <- ggplot(final_table_cnvs_asd_id_dd, aes(x="", y=percentage, fill=gene_overlap)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = round(percentage, digits = 2)),
            position = position_stack(vjust = 0.5),
            color = 'black',
            size = 7) + 
  scale_fill_manual(values = c("tomato2", "gold", "lemonchiffon")) + 
  coord_polar("y", start=0) +
  theme_void() 

asd_id_dd_int_pie <- ggplot(final_table_int_asd_id_dd, aes(x="", percentage, fill=hic_interactions)) +
  geom_bar(stat="identity", width=1) + 
  scale_fill_manual(values = c("cornflowerblue", "lightgreen", "mediumpurple1")) + 
  geom_text(aes(label = round(percentage, digits = 2)),
            position = position_stack(vjust = 0.5),
            color = 'black',
            size = 7) +
  coord_polar("y", start=0) +
  theme_void() 

final_table_id_dd <- final_table_clean[final_table_clean$phenotype=="ID_DD",]
final_interaction_id_dd <- results_merged_one_side_process[results_merged_one_side_process$phenotype=="ID_DD",]
final_interaction_id_dd <- final_interaction_id_dd[final_interaction_id_dd$non_coding_cnv==TRUE,]

final_table_cnvs_id_dd <- final_table_id_dd %>%
  group_by(gene_overlap) %>%
  summarise(
    value = n(),
    percentage = 100*n()/nrow(final_table_id_dd)
  )

final_table_int_id_dd <- final_interaction_id_dd %>%
  group_by(hic_interactions) %>%
  summarise(
    value = n(),
    percentage = 100*n()/nrow(final_interaction_id_dd)
  )

# Pie Chart
id_dd_cnvs_pie <- ggplot(final_table_cnvs_id_dd, aes(x="", y=percentage, fill=gene_overlap)) +
  geom_bar(stat="identity", width=1) + 
  scale_fill_manual(values = c("tomato2", "gold", "lemonchiffon")) + 
  geom_text(aes(label = round(percentage, digits = 2)),
            position = position_stack(vjust = 0.5),
            color = 'black',
            size = 7) +
  coord_polar("y", start=0) +
  theme_void() 

id_dd_int_pie <- ggplot(final_table_int_id_dd, aes(x="", percentage, fill=hic_interactions)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = c("cornflowerblue", "lightgreen", "mediumpurple1")) + 
  geom_text(aes(label = round(percentage, digits = 2)),
            position = position_stack(vjust = 0.5),
            color = 'black',
            size = 7) +
  coord_polar("y", start=0) +
  theme_void() 

final_table_scz <- final_table_clean[final_table_clean$phenotype=="SCZ",]
final_interaction_scz <- results_merged_one_side_process[results_merged_one_side_process$phenotype=="SCZ",]
final_interaction_scz <- final_interaction_scz[final_interaction_scz$non_coding_cnv==TRUE,]

final_table_cnvs_scz <- final_table_scz %>%
  group_by(gene_overlap) %>%
  summarise(
    value = n(),
    percentage = 100*n()/nrow(final_table_scz)
  )

final_table_int_scz <- final_interaction_scz %>%
  group_by(hic_interactions) %>%
  summarise(
    value = n(),
    percentage = 100*n()/nrow(final_interaction_scz)
  )

# Pie Chart
scz_cnvs_pie <- ggplot(final_table_cnvs_scz, aes(x="", y=percentage, fill=gene_overlap)) +
  geom_bar(stat="identity", width=1) + 
  scale_fill_manual(values = c("tomato2", "gold", "lemonchiffon")) + 
  geom_text(aes(label = round(percentage, digits = 2)),
            position = position_stack(vjust = 0.5),
            color = 'black',
            size = 7) +
  coord_polar("y", start=0) +
  theme_void() 

scz_int_pie <- ggplot(final_table_int_scz, aes(x="", percentage, fill=hic_interactions)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = c("cornflowerblue", "lightgreen", "mediumpurple1")) + 
  geom_text(aes(label = round(percentage, digits = 2)),
            position = position_stack(vjust = 0.5),
            color = 'black',
            size = 7) +
  coord_polar("y", start=0) +
  theme_void() 

cnv_pie <- ggarrange(asd_cnvs_pie, 
          asd_id_dd_cnvs_pie, 
          id_dd_cnvs_pie,
          scz_cnvs_pie,
          labels = c("a) ASD", 
                     "b) ASD_ID_DD", 
                     "c) ID_DD", 
                     "d) SCZ"),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "right") +
    theme(plot.margin = margin(t = 0.4, r = 0.2, b  = 0.4, l = 0.2, unit = "cm"))

  
 

int_pie <- ggarrange(asd_int_pie, 
                     asd_id_dd_int_pie, 
                     id_dd_int_pie,
                     scz_int_pie,
                     labels = c("a) ASD", 
                                "b) ASD_ID_DD", 
                                "c) ID_DD", 
                                "d) SCZ"),
                     ncol = 2, nrow = 2,common.legend = TRUE, legend = "right") + 
  theme(plot.margin = margin(t = 0.4, r = 0.2, b  = 0.4, l = 0.2, unit = "cm")) 
pie_concated <- ggarrange(cnv_pie, int_pie,
                          labels = c("Section B: CNVs Overlap Categories", 
                                     "Section C: Chromatin Interactions of ncCNVs"),
                          ncol = 1, nrow = 2) + theme(plot.margin = margin(t = 0.4, r = 0.2, b  = 0.4, l = 0.2, unit = "cm"))

ggsave("pie_concated.svg", plot = pie_concated, width = 15, height = 25)
