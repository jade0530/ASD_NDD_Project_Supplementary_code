library(tidyr)
library(dplyr)
library(readr)
library(readxl)
library(regioneR)
library(stringr)
library(karyoploteR)
library(GenomicRanges)
library(Repitools)
library(plyranges)

# Load gene regions 
gene_list <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/Hi-C_INTERACTION/genes_list_three_consortiums_hg19.csv")
gene_regions <- toGRanges(data.frame(chr=gene_list$chr, start=gene_list$start, 
                                     end=gene_list$end, gene_symbol=gene_list$gene_symbol, gene_name=gene_list$gene_name, gene_class=gene_list$gene_class), genome = "hg19")
promoter_regions <- toGRanges(data.frame(chr=gene_list$chr, promoter_start=gene_list$promoter_start, 
                                            promoter_end=gene_list$promoter_end, gene_symbol=gene_list$gene_symbol, gene_name=gene_list$gene_name, gene_class=gene_list$gene_class), genome = "hg19")
# load enhancer regions
enhancer_list<- read.csv("/media/bml/USER_DATA/ASD_CNV_Project/Hi-C_INTERACTION/HMMMark_cortex.csv", sep = "\t")
enhancer_regions <-toGRanges(data.frame(chr=enhancer_list$chr, enhancer_start=enhancer_list$start, 
                                        enhancer_end=enhancer_list$end, type=enhancer_list$type), genome = "hg19") 
# Load Hi-C interaction data 
hic_data <- read_csv("/media/bml/USER_DATA/ASD_CNV_Project/Hi-C_INTERACTION/HiC_data_from_papers.csv")
interacting_locus <- toGRanges(data.frame(chr=hic_data$chr, start=hic_data$interacting_locus_start, 
                                     end=hic_data$interacting_locus_end, reads=hic_data$`# reads`, paper=hic_data$Paper), genome = "hg19")
evol_locus <- toGRanges(data.frame(chr=hic_data$chr, start=hic_data$evol_locus_start, 
                                         end=hic_data$evol_locus_end, reads=hic_data$`# reads`, paper=hic_data$Paper), genome = "hg19")

# Load CNV region data 
asd_cnv_data <- read_excel("/media/bml/USER_DATA/ASD_CNV_Project/CNVs_significant_on_cnv_New_V3.xlsx", sheet = "ASD_merged")
asd_cnv_regions <- toGRanges(data.frame(chr=asd_cnv_data$chr, start=asd_cnv_data$start, 
                                        end=asd_cnv_data$end, cnv_type=asd_cnv_data$cnv_type), genome = "hg19")

asd_id_dd_cnv_data <- read_excel("/media/bml/USER_DATA/ASD_CNV_Project/CNVs_significant_on_cnv_New_V3.xlsx", sheet = "ASD_ID_DD_merged")
asd_id_dd_cnv_regions <- toGRanges(data.frame(chr=asd_id_dd_cnv_data$chr, start=asd_id_dd_cnv_data$start, 
                                        end=asd_id_dd_cnv_data$end, cnv_type=asd_id_dd_cnv_data$cnv_type), genome = "hg19")

id_dd_cnv_data <- read_excel("/media/bml/USER_DATA/ASD_CNV_Project/CNVs_significant_on_cnv_New_V3.xlsx", sheet = "ID_DD_merged")
id_dd_cnv_regions <- toGRanges(data.frame(chr=id_dd_cnv_data$chr, start=id_dd_cnv_data$start, 
                                              end=id_dd_cnv_data$end, cnv_type=id_dd_cnv_data$cnv_type), genome = "hg19")

scz_cnv_data <- read_excel("/media/bml/USER_DATA/ASD_CNV_Project/CNVs_significant_on_cnv_New_V3.xlsx", sheet = "SCZ_merged")
scz_cnv_regions <- toGRanges(data.frame(chr=scz_cnv_data$chr, start=scz_cnv_data$start, 
                                              end=scz_cnv_data$end, cnv_type=scz_cnv_data$cnv_type), genome = "hg19")
################################################# Rare and common CNV regions #####################################
common_rare_cnv_regions <- function(query_region, subject_regions) {
  cnv_id <- paste0(seqnames(query_region), ":", 
                   start(query_region), "-", 
                   end(query_region))
  overlapping <- findOverlaps(query_region, subject_regions)
  query_result <- query_region[queryHits(overlapping)]
  subject_result <- subject_regions[subjectHits(overlapping)]

  print(overlapping)
  r1 <- ranges(query_result)
  r2 <- ranges(subject_result)
  
  print(length(countOverlaps(query_result, subject_result)))
  # new GRanges object only with the overlapping intervals
  if (length(countOverlaps(query_result, subject_result)) > 0) {
    print("There is result overlapping regions")
    mcols(subject_result)$common_id = paste0(seqnames(subject_result), ":", 
                                             start(subject_result), "-", 
                                             end(subject_result))
    # Overlapped region start and end 
    gr<-
      GRanges(seqnames = as.character(seqnames(subject_result)),
              ranges = IRanges(start = matrixStats::rowMaxs(as.matrix(data.frame(r1=start(r1), r2=start(r2)))),
                               end = matrixStats::rowMins(as.matrix(data.frame(r1=end(r1), r2=end(r2))))))
    gr_overlap <- gr
    if (countOverlaps(query_region, scz_cnv_regions) > 1) {
      # If there's multiple ranges 
      gr_overlap <- reduce(gr)
      print(gr_overlap)
    }
    if (length(mcols(subject_result)$common_id)>1) {
      common_id <- paste(mcols(subject_result)$common_id, sep = ";")
    } else {
      common_id <- mcols(subject_result)$common_id
    }
    print("calculate difference between the regions")
    diff <- setdiff(gr, query_result, ignore.strand=TRUE)
    print(diff)
    if (nrow(diff_df) == 0) {
      print("The regions are fully overlapped")
      diff_id = 'ignore'
    } else {
      mcols(diff)$diff_id <- paste0(seqnames(diff), ":", 
                                    start(diff), "-", 
                                    end(diff))
      if (length(mcols(diff)$diff_id)>1) {
        diff_id = paste(mcols(diff)$diff_id, sep = ";")
      } else {
        diff_id <- mcols(diff)$diff_id
      }
    }
    res_row = data.frame(cnv_id, common_id, diff_id)
    print(res_row)
    colnames(res_row) <- c("cnv_id", "common_id",
                           "diff_id")
    return(res_row)
  } else {
    print("there's no overlapping regions")
    common_id = 'ignore' 
    diff_id = 'ignore' 
    res_row = data.frame(cnv_id, common_id, diff_id)
    print(res_row)
    return(res_row)
    
  }
  
}

######################################### Run Rare and common CNV regions#####################
# Find overlap regions
common_scz <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("cnv_id", "common_id", "diff_id"))))
for(i in seq(1,nrow(asd_cnv_data))) {
  query_region <- asd_cnv_regions[i:i]
  res <- common_rare_cnv_regions(query_region, scz_cnv_regions)
  common_scz <- rbind(res, common_scz)
}
colnames(common_scz) <- c("cnv_id", "common_id_scz", "diff_id_scz")
common_scz <- common_scz %>%
  group_by(cnv_id) %>%
  summarise(common_id_scz = paste(common_id_scz, sep = ";"))
common_asd_id_dd <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("cnv_id", "common_id", "diff_id"))))
for(i in seq(1,nrow(asd_cnv_data))) {
  query_region <- asd_cnv_regions[i:i]
  res <- common_rare_cnv_regions(query_region, asd_id_dd_cnv_regions)
  common_asd_id_dd <- rbind(res, common_asd_id_dd)
}
colnames(common_asd_id_dd) <- c("cnv_id", "common_id_asd_id_dd", "diff_id_asd_id_dd")
common_asd_id_dd_grouped <- common_asd_id_dd %>%
  group_by(cnv_id) %>%
  summarise(common_id_asd_id_dd_grouped = paste(common_id_asd_id_dd, sep = ";"))

common_id_dd <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("cnv_id", "common_id", "diff_id"))))
for(i in seq(1,nrow(asd_cnv_data))) {
  query_region <- asd_cnv_regions[i:i]
  res <- common_rare_cnv_regions(query_region, id_dd_cnv_regions)
  common_id_dd <- rbind(res, common_asd_id_dd)
}
colnames(common_id_dd) <- c("cnv_id", "common_id_id_dd", "diff_id_id_dd")
common_id_dd <- common_id_dd %>%
  group_by(cnv_id) %>%
  summarise(common_id_id_dd = paste(common_id_id_dd, sep = ";"))
res_asd <- merge(common_asd_id_dd, common_id_dd, by = "cnv_id")
res_asd <- merge(res_asd, common_scz, by = "cnv_id")


###################################################### Gene Analysis ##########################################################
generate_gene_overlapped_table <- function(phenotype, regions) {
  overlapped_genes <- findOverlaps(phenotype, regions)
  cnv_results <- phenotype[queryHits(overlapped_genes)]
  gene_results <- regions[subjectHits(overlapped_genes)]
  
  cnv_results_df <- data.frame(
    chr = as.vector(seqnames(cnv_results)), 
    start = cnv_results@ranges@start,
    end = cnv_results@ranges@start+cnv_results@ranges@width-1, 
    cnv_type = cnv_results@elementMetadata@listData$cnv_type)
  
  gene_results_df <- data.frame(
    chr = as.vector(seqnames(gene_results)), 
    start = gene_results@ranges@start,
    end = gene_results@ranges@start+gene_results@ranges@width-1, 
    enhancer_id = paste0(seqnames(gene_results), ":", 
                         start(gene_results), "-", 
                         end(gene_results)), 
    type = gene_results@elementMetadata@listData$type)

  #gene_results_df$`#_coding_genes` <- nrow(gene_results_df[gene_results_df$gene_class=='protein_coding'])
  #gene_results_df$`#_noncoding_genes` <- nrow(gene_results_df[gene_results_df$gene_class %like% 'noncoding'])

  #result_table_evol <- data.frame(matrix(ncol=18,nrow=0, dimnames=list(NULL, c("chr", "cnv_start", "cnv_end", "cnv_type", "#_coding_genes", "coding_gene_symbols", 
  #                                                                             "coding_gene_coordination", "#_noncoding_genes", 
  #                                                                             "noncoding_gene_symbols", "noncoding_gene_coordination"))))
  cnv_gene_overlap_df <- merge(cnv_results_df, gene_results_df, by = 'row.names', all = TRUE)
  cnv_gene_overlap_df <- cnv_gene_overlap_df[, 2:ncol(cnv_gene_overlap_df)]
  colnames(cnv_gene_overlap_df) <- c("cnv_chr", "cnv_start", "cnv_end", "cnv_type", "enhancer_chr", "enhancer_start", "enhancer_end", "enhancer_id", "enhancer_type")
  
  return(cnv_gene_overlap_df)
}

generate_gene_overlapped_table_10_perc <- function(phenotype) {
  # calculate overlapped bp
  overlap <- gene_regions@ranges@width
  overlapped_genes <- findOverlaps(phenotype, gene_regions)
  cnv_results <- phenotype[queryHits(overlapped_genes)]
  gene_results <- gene_regions[subjectHits(overlapped_genes)]
  
  #/ subset the input ranges to those ranges that have overlaps with the second ranges object
  r1 <- ranges(phenotype[queryHits(overlapped_genes)])
  r2 <- ranges(gene_regions[subjectHits(overlapped_genes)])
  
  #/ new GRanges object only with the overlapping intervals
  gr<-
    GRanges(seqnames = as.character(seqnames(phenotype[queryHits(overlapped_genes)])),
            ranges = IRanges(start = matrixStats::rowMaxs(as.matrix(data.frame(r1=start(r1), r2=start(r2)))),
                             end = matrixStats::rowMins(as.matrix(data.frame(r1=end(r1), r2=end(r2))))))
  
  # overlapped region/gene length 
  # overlapped region = 
  mcols(gene_results)$percentOverlapGenes <- 100 * width(gr) / width(r2)
  mcols(gene_results)$percentOverlapCNVs <- 100 * width(gr) / width(r1)
  mcols(gene_results)$isCNVBigger <- (width(r1) > width(r2))
  
  
  cnv_results_df <- data.frame(
    chr = as.vector(seqnames(cnv_results)), 
    start = cnv_results@ranges@start,
    end = cnv_results@ranges@start+cnv_results@ranges@width-1, 
    cnv_type = cnv_results@elementMetadata@listData$cnv_type)
  
  gene_results_df <- data.frame(
    chr = as.vector(seqnames(gene_results)), 
    start = gene_results@ranges@start,
    end = gene_results@ranges@start+gene_results@ranges@width-1, 
    gene_symbol = gene_results@elementMetadata@listData$gene_symbol,
    gene_name = gene_results@elementMetadata@listData$gene_name,
    gene_class = gene_results@elementMetadata@listData$gene_class,
    percentOverlapGenes = gene_results@elementMetadata@listData$percentOverlapGenes,
    percentOverlapCNVs = gene_results@elementMetadata@listData$percentOverlapCNVs,
    isCNVBigger = gene_results@elementMetadata@listData$isCNVBigger)

  cnv_gene_overlap_df <- merge(cnv_results_df, gene_results_df, by = 'row.names', all = TRUE)
  cnv_gene_overlap_df <- cnv_gene_overlap_df[, 2:ncol(cnv_gene_overlap_df)]
  colnames(cnv_gene_overlap_df) <- c("cnv_chr", "cnv_start", "cnv_end", "cnv_type", "gene_chr", "gene_start", "gene_end", "gene_symbol", "gene_name", "gene_class", "percentOverlapGenes", "percentOverlapCNVs", "isCNVBigger")
  
  return(cnv_gene_overlap_df)
}

###################################################### Interaction Analysis ########################################################

calculate_interacting_interaction_table <- function(phenotype, gene_overlap_result) {
phenotype_overlap_interacting_loci <- findOverlaps(phenotype, interacting_locus)
cnv_results_interacting_loci <- phenotype[queryHits(phenotype_overlap_interacting_loci)]
results_interacting_loci <- interacting_locus[subjectHits(phenotype_overlap_interacting_loci)]

cnv_results_interacting_loci_df <- data.frame(
  chr = as.vector(seqnames(cnv_results_interacting_loci)), 
  start = cnv_results_interacting_loci@ranges@start,
  end = cnv_results_interacting_loci@ranges@start+cnv_results_interacting_loci@ranges@width-1, 
  cnv_type = cnv_results_interacting_loci@elementMetadata@listData$cnv_type)

results_interacting_loci_df <- data.frame(
  chr = as.vector(seqnames(results_interacting_loci)), 
  start = results_interacting_loci@ranges@start,
  end = results_interacting_loci@ranges@start+results_interacting_loci@ranges@width-1, 
  reads = results_interacting_loci@elementMetadata@listData$reads,
  paper = results_interacting_loci@elementMetadata@listData$paper)
                      
interacting_loci_df <- merge(cnv_results_interacting_loci_df, results_interacting_loci_df, by = 'row.names', all = TRUE)
interacting_loci_df <- interacting_loci_df[, 2:ncol(interacting_loci_df)]
interacting_loci_df$is_left_hand_side = TRUE
colnames(interacting_loci_df) <- c('cnv_chr', 'cnv_start', 'cnv_end', 'cnv_type','hic_chr','hic_start', 'hic_end', 'reads', 'paper', 'is_left_hand_side')

result_table_interacting <- data.frame(matrix(ncol=18,nrow=0, dimnames=list(NULL, c("chr", "cnv_start", "cnv_end", "cnv_type", "interacting_locus_start", "interacting_locus_end", 
                                                                                    "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                                                                                    "evol_locus_start", "evol_locus_end", 
                                                                                    "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class"))))
for (i in 1:nrow(interacting_loci_df)) {
  print(interacting_loci_df[i,])
  res <- hic_analysis_function_interacting(interacting_loci_df[i,], gene_overlap_result)
  result_table_interacting <- rbind(result_table_interacting, res) 
}

return(result_table_interacting)
}

calculate_evol_interaction_table <- function(phenotype, gene_overlap_result) {
  phenotype_overlap_evol_loci <- findOverlaps(phenotype, evol_locus)
  cnv_results_evol_loci <- phenotype[queryHits(phenotype_overlap_evol_loci)]
  results_evol_loci <- evol_locus[subjectHits(phenotype_overlap_evol_loci)]
  
  cnv_results_evol_loci_df <- data.frame(
    chr = as.vector(seqnames(cnv_results_evol_loci)), 
    start = cnv_results_evol_loci@ranges@start,
    end = cnv_results_evol_loci@ranges@start+cnv_results_evol_loci@ranges@width-1, 
    cnv_type = cnv_results_evol_loci@elementMetadata@listData$cnv_type)
  
  results_evol_loci_df <- data.frame(
    chr = as.vector(seqnames(results_evol_loci)), 
    start = results_evol_loci@ranges@start,
    end = results_evol_loci@ranges@start+results_evol_loci@ranges@width-1, 
    reads = results_evol_loci@elementMetadata@listData$reads,
    paper = results_evol_loci@elementMetadata@listData$paper)
  
  evol_loci_df <- merge(cnv_results_evol_loci_df, results_evol_loci_df, by = 'row.names', all = TRUE)
  evol_loci_df <- evol_loci_df[, 2:ncol(evol_loci_df)]
  evol_loci_df$is_left_hand_side = FALSE
  colnames(evol_loci_df) <- c('cnv_chr', 'cnv_start', 'cnv_end', 'cnv_type','hic_chr','hic_start', 'hic_end', 'reads', 'paper', 'is_left_hand_side')
  
  result_table_evol <- data.frame(matrix(ncol=18,nrow=0, dimnames=list(NULL, c("chr", "cnv_start", "cnv_end", "cnv_type", "evol_locus_start", "evol_locus_end", 
                                                                               "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                                                                               "interacting_locus_start", "interacting_locus_end", 
                                                                               "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class"))))
  
  for (i in 1:nrow(evol_loci_df)) {
    print(evol_loci_df[i,])
    res <- hic_analysis_function_evol(evol_loci_df[i,],gene_overlap_result)
    result_table_evol <- rbind(result_table_evol, res) 
  }
  return(result_table_evol)
}



hic_analysis_function_evol <- function(r, gene_overlap_result) {
  res <- data.frame(matrix(ncol=25,nrow=0, dimnames=list(NULL, c("chr", "cnv_start", "cnv_end", "cnv_type", "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "evol_locus_start", "evol_locus_end", 
                                                                 "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                                                                 "interacting_locus_start", "interacting_locus_end", 
                                                                 "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                                                                 "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class", 
                                                                 "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end"))))
  r_df <- as.data.frame(r)

  no_coding_genes <- gene_overlap_result[which(gene_overlap_result$cnv_chr %in% r_df$hic_chr &
                                                        gene_overlap_result$cnv_start %in% r_df$cnv_start &
                                                        gene_overlap_result$cnv_end %in% r_df$cnv_end),c("#_coding_genes")]
  
  no_noncoding_genes <- gene_overlap_result[which(gene_overlap_result$cnv_chr %in% r_df$hic_chr &
                                                           gene_overlap_result$cnv_start %in% r_df$cnv_start &
                                                           gene_overlap_result$cnv_end %in% r_df$cnv_end),c("#_noncoding_genes")]
  if (length(no_coding_genes) == 0) {
    no_coding_genes = 0
  }
  if (length(no_noncoding_genes) == 0) {
    no_noncoding_genes = 0
  }
  if (no_coding_genes > 0) {
    # If there's no coding genes -> non_coding_cnv = TRUE
    r_df$noncoding_cnv = FALSE
  } else {
    # If there's coding genes -> non_coding_cnv = FALSE
    r_df$noncoding_cnv = TRUE
  }
  # parse each interacting loci
  e_loci <- as.data.frame(r_df[c('hic_chr', 'hic_start', 'hic_end')])
  
  # Look for the interacting region
  i_loci <- hic_data[which(hic_data$chr %in% e_loci$hic_chr &
                     hic_data$evol_locus_start %in% e_loci$hic_start &
                     hic_data$evol_locus_end %in% e_loci$hic_end),]
  #print(i_loci)
  # Iterate the interacting region
  for (i in 1:nrow(i_loci)) {
    # First, determine if the cnv covers both side 
    single_i_loci <- as.data.frame(i_loci[i,])
    overlapped_other_start <- single_i_loci$interacting_locus_start
    overlapped_other_end <- single_i_loci$interacting_locus_end
    reads <- single_i_loci$`# reads`
    paper <- single_i_loci$Paper
    if (single_i_loci$chr == r_df$cnv_chr &&  (as.integer(r_df$cnv_start) > overlapped_other_end) && (as.integer(e_loci$hic_start) > as.integer(overlapped_other_end))) {
      print("the cnv has one side overlapped")
      overlap_with_both_hic_side <- FALSE

      # Find gene body overlapped with evol regions
      query_range <- toGRanges(data.frame(chr=single_i_loci$chr, start=single_i_loci$interacting_locus_start, 
                                          end=single_i_loci$interacting_locus_end,  reads=single_i_loci$`# reads`, paper=single_i_loci$Paper), genome = "hg19")
      
      gene_body_res <- findOverlaps(query_range, gene_regions)
      gene_body_overlapped <- gene_regions[subjectHits(gene_body_res)]
      cnv_overlapped <- query_range[queryHits(gene_body_res)]

      gene_body_overlapped@elementMetadata@listData$gene_start_cnv_end <- gene_body_overlapped@ranges@start-(cnv_overlapped@ranges@start + cnv_overlapped@ranges@width - 1)
      gene_body_overlapped@elementMetadata@listData$gene_end_cnv_start <- (gene_body_overlapped@ranges@start+gene_body_overlapped@ranges@width-1)-cnv_overlapped@ranges@start
      
      # We need gene distance here 
      if (countOverlaps(query_range, gene_regions)>1) {
        print("genes were overlapped with this region")
        overlapped_gene_name <- paste(gene_body_overlapped@elementMetadata@listData$gene_name, collapse = ',')
        overlapped_gene_symbol <- paste(gene_body_overlapped@elementMetadata@listData$gene_symbol, collapse = ',')
        overlapped_gene_class <- paste(gene_body_overlapped@elementMetadata@listData$gene_class, collapse = ',')
        
        overlapped_gene_end_cnv_start <- paste(gene_body_overlapped@elementMetadata@listData$gene_end_cnv_start, collapse = ',')
        overlapped_gene_start_cnv_end <- paste(gene_body_overlapped@elementMetadata@listData$gene_start_cnv_end, collapse = ',')
        
      } else {
        print("No genes were overlapped with this region")
        overlapped_gene_symbol <- 'ignore'
        overlapped_gene_name <- 'ignore'
        overlapped_gene_class <- 'ignore'
        overlapped_gene_start_cnv_end <- 'ignore'
        overlapped_gene_end_cnv_start <- 'ignore'
      }

      # Find promotor overlapped with evol regions #
      gene_promoter_res <- findOverlaps(query_range, promoter_regions)
      gene_promoter_overlapped <- gene_regions[subjectHits(gene_promoter_res)]
      promoter_cnv_overlapped <- query_range[queryHits(gene_promoter_res)]
      gene_promoter_overlapped@elementMetadata@listData$gene_start_cnv_end <- gene_promoter_overlapped@ranges@start-(promoter_cnv_overlapped@ranges@start+
                                                                                                                       promoter_cnv_overlapped@ranges@width-1)
      gene_promoter_overlapped@elementMetadata@listData$gene_end_cnv_start <- (gene_promoter_overlapped@ranges@start+gene_promoter_overlapped@ranges@width-1)-promoter_cnv_overlapped@ranges@start
      
      
      if (countOverlaps(query_range, promoter_regions)>1) {
            print("noncoding promoters were overlapped with this region")
            overlapped_promoter_gene_symbol <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_symbol, collapse = ',')
            print(overlapped_promoter_gene_symbol)
            overlapped_promoter_gene_name <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_name, collapse = ',')
            overlapped_promoter_gene_class <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_class, collapse = ',')
            overlapped_promoter_gene_start_cnv_end <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_start_cnv_end, collapse = ',')
            overlapped_promoter_gene_end_cnv_start <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_end_cnv_start, collapse = ',')
            
            } else {
            print("No noncoding promoters were overlapped with this region")
            overlapped_promoter_gene_symbol <- 'ignore'
            overlapped_promoter_gene_name <- 'ignore'
            overlapped_promoter_gene_class <- 'ignore'
            overlapped_promoter_gene_start_cnv_end <- 'ignore'
            overlapped_promoter_gene_end_cnv_start <- 'ignore'
      }
      
    } else {
      # Meaning the cnv covers with both side  
      print("the cnv covers with the both sides")
      overlap_with_both_hic_side <- TRUE

      # wouldn't bother calculating the genes 
      overlapped_gene_symbol <- 'ignore'
      overlapped_gene_name <- 'ignore'
      overlapped_gene_class <- 'ignore'
      overlapped_promoter_gene_symbol <- 'ignore'
      overlapped_promoter_gene_name <- 'ignore'
      overlapped_promoter_gene_class <- 'ignore'
      overlapped_gene_end_cnv_start <- 'ignore'
      overlapped_gene_start_cnv_end <- 'ignore'
      overlapped_promoter_gene_start_cnv_end <- 'ignore'
      overlapped_promoter_gene_end_cnv_start <- 'ignore'
    }

    hic_interacted_gene_results <- data.frame(c(r_df$hic_chr, 
                                            r_df$cnv_start,
                                            r_df$cnv_end,
                                            r_df$cnv_type, 
                                            no_coding_genes,
                                            no_noncoding_genes, 
                                            r_df$noncoding_cnv, 
                                            r_df$hic_start,
                                            r_df$hic_end, 
                                            reads, 
                                            paper,
                                            r_df$is_left_hand_side, 
                                            overlap_with_both_hic_side, 
                                            overlapped_other_start, 
                                            overlapped_other_end, 
                                            overlapped_gene_symbol,
                                            overlapped_gene_name, 
                                            overlapped_gene_class, 
                                            overlapped_gene_end_cnv_start,
                                            overlapped_gene_start_cnv_end, 
                                            overlapped_promoter_gene_symbol, 
                                            overlapped_promoter_gene_name, 
                                            overlapped_promoter_gene_class,
                                            overlapped_promoter_gene_end_cnv_start,
                                            overlapped_promoter_gene_start_cnv_end))
    
    print(hic_interacted_gene_results)
    hic_interacted_gene_results = t(hic_interacted_gene_results)
    colnames(hic_interacted_gene_results) <- c("chr", "cnv_start", "cnv_end", "cnv_type", "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv",
                                               "evol_locus_start", "evol_locus_end", 
                   "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                   "interacting_locus_start", "interacting_locus_end", 
                   "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                   "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class", 
                   "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end")
    rownames(hic_interacted_gene_results) = NULL
    
    res[nrow(res) + 1, ] <- hic_interacted_gene_results
  }
  return(res)
}

hic_analysis_function_interacting <- function(r, gene_overlap_result) {
  res <- data.frame(matrix(ncol=25,nrow=0, dimnames=list(NULL, c("chr", "cnv_start", "cnv_end", "cnv_type", 
                                                                 "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", 
                                                                 "interacting_locus_start", "interacting_locus_end", "reads", "paper", 
                                                                 "is_left_hand_side", "overlap_with_both_hic_side", 
                                                                 "evol_locus_start", "evol_locus_end", 
                                                                 "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", 
                                                                 "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                                                                 "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class", 
                                                                 "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end"))))
  r_df <- as.data.frame(r)
  no_coding_genes <- gene_overlap_result[which(gene_overlap_result$cnv_chr %in% r_df$hic_chr &
                                                        gene_overlap_result$cnv_start %in% r_df$cnv_start &
                                                        gene_overlap_result$cnv_end %in% r_df$cnv_end),c("#_coding_genes")]
  
  no_noncoding_genes <- gene_overlap_result[which(gene_overlap_result$cnv_chr %in% r_df$hic_chr &
                                                           gene_overlap_result$cnv_start %in% r_df$cnv_start &
                                                           gene_overlap_result$cnv_end %in% r_df$cnv_end),c("#_noncoding_genes")]
  if (length(no_coding_genes) == 0) {
    no_coding_genes = 0
  }
  if (length(no_noncoding_genes) == 0) {
    no_noncoding_genes = 0
  }
  if (no_coding_genes > 0) {
    # If there's no coding genes -> non_coding_cnv = TRUE
    r_df$noncoding_cnv = FALSE
  } else {
    # If there's coding genes -> non_coding_cnv = FALSE
    r_df$noncoding_cnv = TRUE
  }  
  
  # parse each interacting loci
  i_loci <- as.data.frame(r_df[c('hic_chr', 'hic_start', 'hic_end')])
  
  # Look for the hic interacting region
  e_loci <- hic_data[which(hic_data$chr %in% i_loci$hic_chr &
                             hic_data$interacting_locus_start %in% i_loci$hic_start &
                             hic_data$interacting_locus_end %in% i_loci$hic_end),]
  print(e_loci)
  # Iterate the interacting region
  for (i in 1:nrow(e_loci)) {
    # First, determine if the cnv covers both side 
    single_e_loci <- as.data.frame(e_loci[i,])
    overlapped_other_start <- single_e_loci$evol_locus_start
    overlapped_other_end <- single_e_loci$evol_locus_end
    
    reads <- single_e_loci$`# reads`
    paper <- single_e_loci$Paper
    
    if (single_e_loci$chr == r_df$cnv_chr &&  (as.integer(overlapped_other_start) > as.integer(r_df$cnv_end)) && (as.integer(overlapped_other_start) > as.integer(i_loci$hic_end))) {
      print("the cnv has one side overlapped")
      overlap_with_both_hic_side <- FALSE
      
      # Find gene body overlapped with evol regions
      query_range <- toGRanges(data.frame(chr=single_e_loci$chr, start=single_e_loci$evol_locus_start, 
                                          end=single_e_loci$evol_locus_end,  reads=single_e_loci$`# reads`, paper=single_e_loci$Paper), genome = "hg19")
      
      gene_body_res <- findOverlaps(query_range, gene_regions)
      gene_body_overlapped <- gene_regions[subjectHits(gene_body_res)]
      cnv_overlapped <- query_range[queryHits(gene_body_res)]
      
      gene_body_overlapped@elementMetadata@listData$gene_start_cnv_end <- gene_body_overlapped@ranges@start-(cnv_overlapped@ranges@start + cnv_overlapped@ranges@width - 1)
      gene_body_overlapped@elementMetadata@listData$gene_end_cnv_start <- (gene_body_overlapped@ranges@start+gene_body_overlapped@ranges@width-1)-cnv_overlapped@ranges@start
      
      if (countOverlaps(query_range, gene_regions)>1) {
        print("genes were overlapped with this region")
        overlapped_gene_name <- paste(gene_body_overlapped@elementMetadata@listData$gene_name, collapse = ',')
        overlapped_gene_symbol <- paste(gene_body_overlapped@elementMetadata@listData$gene_symbol, collapse = ',')
        overlapped_gene_class <- paste(gene_body_overlapped@elementMetadata@listData$gene_class, collapse = ',')
        overlapped_gene_end_cnv_start <- paste(gene_body_overlapped@elementMetadata@listData$gene_end_cnv_start, collapse = ',')
        overlapped_gene_start_cnv_end <- paste(gene_body_overlapped@elementMetadata@listData$gene_start_cnv_end, collapse = ',')
        
      } else {
        print("No genes were overlapped with this region")
        overlapped_gene_symbol <- 'ignore'
        overlapped_gene_name <- 'ignore'
        overlapped_gene_class <- 'ignore'
        overlapped_gene_start_cnv_end <- 'ignore'
        overlapped_gene_end_cnv_start <- 'ignore'
      }
      
      # Find promotor overlapped with evol regions #
      gene_promoter_res <- findOverlaps(query_range, promoter_regions)
      gene_promoter_overlapped <- gene_regions[subjectHits(gene_promoter_res)]
      promoter_cnv_overlapped <- query_range[queryHits(gene_promoter_res)]
      gene_promoter_overlapped@elementMetadata@listData$gene_start_cnv_end <- gene_promoter_overlapped@ranges@start-(promoter_cnv_overlapped@ranges@start+
                                                                                                                       promoter_cnv_overlapped@ranges@width-1)
      gene_promoter_overlapped@elementMetadata@listData$gene_end_cnv_start <- (gene_promoter_overlapped@ranges@start+gene_promoter_overlapped@ranges@width-1)-promoter_cnv_overlapped@ranges@start
      
      
      if (countOverlaps(query_range, promoter_regions)>1) {
        print("noncoding promoters were overlapped with this region")
        overlapped_promoter_gene_symbol <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_symbol, collapse = ',')
        print(overlapped_promoter_gene_symbol)
        overlapped_promoter_gene_name <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_name, collapse = ',')
        overlapped_promoter_gene_class <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_class, collapse = ',')
        overlapped_promoter_gene_start_cnv_end <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_start_cnv_end, collapse = ',')
        overlapped_promoter_gene_end_cnv_start <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_end_cnv_start, collapse = ',')
        
      } else {
        print("No noncoding promoters were overlapped with this region")
        overlapped_promoter_gene_symbol <- 'ignore'
        overlapped_promoter_gene_name <- 'ignore'
        overlapped_promoter_gene_class <- 'ignore'
        overlapped_promoter_gene_start_cnv_end <- 'ignore'
        overlapped_promoter_gene_end_cnv_start <- 'ignore'
      }
      
    } else {
      # Meaning the cnv covers with both side  
      print("the cnv covers with the both sides")
      overlap_with_both_hic_side <- TRUE
      
      # wouldn't bother calculating the genes 
      overlapped_gene_symbol <- 'ignore'
      overlapped_gene_name <- 'ignore'
      overlapped_gene_class <- 'ignore'
      overlapped_promoter_gene_symbol <- 'ignore'
      overlapped_promoter_gene_name <- 'ignore'
      overlapped_promoter_gene_class <- 'ignore'
      overlapped_gene_end_cnv_start <- 'ignore'
      overlapped_gene_start_cnv_end <- 'ignore'
      overlapped_promoter_gene_start_cnv_end <- 'ignore'
      overlapped_promoter_gene_end_cnv_start <- 'ignore'
    }
    
    hic_interacted_gene_results <- data.frame(c(r_df$hic_chr, 
                                                r_df$cnv_start,
                                                r_df$cnv_end,
                                                r_df$cnv_type, 
                                                no_coding_genes,
                                                no_noncoding_genes, 
                                                r_df$noncoding_cnv, 
                                                r_df$hic_start,
                                                r_df$hic_end, 
                                                reads, 
                                                paper,
                                                r_df$is_left_hand_side, 
                                                overlap_with_both_hic_side, 
                                                overlapped_other_start, 
                                                overlapped_other_end, 
                                                overlapped_gene_symbol,
                                                overlapped_gene_name, 
                                                overlapped_gene_class, 
                                                overlapped_gene_end_cnv_start,
                                                overlapped_gene_start_cnv_end, 
                                                overlapped_promoter_gene_symbol, 
                                                overlapped_promoter_gene_name, 
                                                overlapped_promoter_gene_class,
                                                overlapped_promoter_gene_end_cnv_start,
                                                overlapped_promoter_gene_start_cnv_end))
    
    print(hic_interacted_gene_results)
    hic_interacted_gene_results = t(hic_interacted_gene_results)
    colnames(hic_interacted_gene_results) <- c("chr", "cnv_start", "cnv_end", "cnv_type", 
                                               "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv",
                                               "evol_locus_start", "evol_locus_end", 
                                               "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                                               "interacting_locus_start", "interacting_locus_end", 
                                               "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", 
                                               "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                                               "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", 
                                               "overlapped_promoter_gene_class", 
                                               "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end")
    rownames(hic_interacted_gene_results) = NULL
    
    res[nrow(res) + 1, ] <- hic_interacted_gene_results
  }
  return(res)
}

count_genes_in_evol_region <- function() {
  # First, determine if the cnv covers both side 
    single_e_loci <- as.data.frame(e_loci[i,])
    overlapped_other_start <- single_e_loci$evol_locus_start
    overlapped_other_end <- single_e_loci$evol_locus_end
    
    reads <- single_e_loci$`# reads`
    paper <- single_e_loci$Paper

    if (single_e_loci$chr == r_df$cnv_chr &&  (as.integer(overlapped_other_start) > as.integer(r_df$cnv_end)) && (as.integer(overlapped_other_start) > as.integer(i_loci$hic_end))) {
      print("the cnv has one side overlapped")
      overlap_with_both_hic_side <- FALSE

      # Find gene body overlapped with evol regions
      query_range <- toGRanges(data.frame(chr=single_e_loci$chr, start=single_e_loci$evol_locus_start, 
                                          end=single_e_loci$evol_locus_end,  reads=single_e_loci$`# reads`, paper=single_e_loci$Paper), genome = "hg19")
      
      gene_body_res <- findOverlaps(query_range, gene_regions)
      gene_body_overlapped <- gene_regions[subjectHits(gene_body_res)]
      cnv_overlapped <- query_range[queryHits(gene_body_res)]
      
      gene_body_overlapped@elementMetadata@listData$gene_start_cnv_end <- gene_body_overlapped@ranges@start-(cnv_overlapped@ranges@start + cnv_overlapped@ranges@width - 1)
      gene_body_overlapped@elementMetadata@listData$gene_end_cnv_start <- (gene_body_overlapped@ranges@start+gene_body_overlapped@ranges@width-1)-cnv_overlapped@ranges@start
      
      if (countOverlaps(query_range, gene_regions)>1) {
        print("genes were overlapped with this region")
        overlapped_gene_name <- paste(gene_body_overlapped@elementMetadata@listData$gene_name, collapse = ',')
        overlapped_gene_symbol <- paste(gene_body_overlapped@elementMetadata@listData$gene_symbol, collapse = ',')
        overlapped_gene_class <- paste(gene_body_overlapped@elementMetadata@listData$gene_class, collapse = ',')
        overlapped_gene_end_cnv_start <- paste(gene_body_overlapped@elementMetadata@listData$gene_end_cnv_start, collapse = ',')
        overlapped_gene_start_cnv_end <- paste(gene_body_overlapped@elementMetadata@listData$gene_start_cnv_end, collapse = ',')
        
      } else {
        print("No genes were overlapped with this region")
        overlapped_gene_symbol <- 'ignore'
        overlapped_gene_name <- 'ignore'
        overlapped_gene_class <- 'ignore'
        overlapped_gene_start_cnv_end <- 'ignore'
        overlapped_gene_end_cnv_start <- 'ignore'
      }
      
      # Find promotor overlapped with evol regions #
      gene_promoter_res <- findOverlaps(query_range, promoter_regions)
      gene_promoter_overlapped <- gene_regions[subjectHits(gene_promoter_res)]
      promoter_cnv_overlapped <- query_range[queryHits(gene_promoter_res)]
      gene_promoter_overlapped@elementMetadata@listData$gene_start_cnv_end <- gene_promoter_overlapped@ranges@start-(promoter_cnv_overlapped@ranges@start+
                                                                                                                       promoter_cnv_overlapped@ranges@width-1)
      gene_promoter_overlapped@elementMetadata@listData$gene_end_cnv_start <- (gene_promoter_overlapped@ranges@start+gene_promoter_overlapped@ranges@width-1)-promoter_cnv_overlapped@ranges@start
      
      
      if (countOverlaps(query_range, promoter_regions)>1) {
        print("noncoding promoters were overlapped with this region")
        overlapped_promoter_gene_symbol <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_symbol, collapse = ',')
        print(overlapped_promoter_gene_symbol)
        overlapped_promoter_gene_name <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_name, collapse = ',')
        overlapped_promoter_gene_class <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_class, collapse = ',')
        overlapped_promoter_gene_start_cnv_end <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_start_cnv_end, collapse = ',')
        overlapped_promoter_gene_end_cnv_start <- paste(gene_promoter_overlapped@elementMetadata@listData$gene_end_cnv_start, collapse = ',')
        
      } else {
        print("No noncoding promoters were overlapped with this region")
        overlapped_promoter_gene_symbol <- 'ignore'
        overlapped_promoter_gene_name <- 'ignore'
        overlapped_promoter_gene_class <- 'ignore'
        overlapped_promoter_gene_start_cnv_end <- 'ignore'
        overlapped_promoter_gene_end_cnv_start <- 'ignore'
      }
      
    } else {
      # Meaning the cnv covers with both side  
      print("the cnv covers with the both sides")
      overlap_with_both_hic_side <- TRUE
      
      # wouldn't bother calculating the genes 
      overlapped_gene_symbol <- 'ignore'
      overlapped_gene_name <- 'ignore'
      overlapped_gene_class <- 'ignore'
      overlapped_promoter_gene_symbol <- 'ignore'
      overlapped_promoter_gene_name <- 'ignore'
      overlapped_promoter_gene_class <- 'ignore'
      overlapped_gene_end_cnv_start <- 'ignore'
      overlapped_gene_start_cnv_end <- 'ignore'
      overlapped_promoter_gene_start_cnv_end <- 'ignore'
      overlapped_promoter_gene_end_cnv_start <- 'ignore'
    }
    
    hic_interacted_gene_results <- data.frame(c(r_df$hic_chr, 
                                                r_df$cnv_start,
                                                r_df$cnv_end,
                                                r_df$cnv_type, 
                                                no_coding_genes,
                                                no_noncoding_genes, 
                                                r_df$noncoding_cnv, 
                                                r_df$hic_start,
                                                r_df$hic_end, 
                                                reads, 
                                                paper,
                                                r_df$is_left_hand_side, 
                                                overlap_with_both_hic_side, 
                                                overlapped_other_start, 
                                                overlapped_other_end, 
                                                overlapped_gene_symbol,
                                                overlapped_gene_name, 
                                                overlapped_gene_class, 
                                                overlapped_gene_end_cnv_start,
                                                overlapped_gene_start_cnv_end, 
                                                overlapped_promoter_gene_symbol, 
                                                overlapped_promoter_gene_name, 
                                                overlapped_promoter_gene_class,
                                                overlapped_promoter_gene_end_cnv_start,
                                                overlapped_promoter_gene_start_cnv_end))
    
    print(hic_interacted_gene_results)
    hic_interacted_gene_results = t(hic_interacted_gene_results)
    colnames(hic_interacted_gene_results) <- c("chr", "cnv_start", "cnv_end", "cnv_type", 
                                               "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv",
                                               "evol_locus_start", "evol_locus_end", 
                                               "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                                               "interacting_locus_start", "interacting_locus_end", 
                                               "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", 
                                               "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                                               "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", 
                                               "overlapped_promoter_gene_class", 
                                               "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end")
    rownames(hic_interacted_gene_results) = NULL
    
    res[nrow(res) + 1, ] <- hic_interacted_gene_results
}
setdiff_cnvs <- function(phenotype_overlap, phenotype) {
  phenotype_overlap_chrs <- phenotype_overlap[,c("cnv_chr", "cnv_start", "cnv_end")]
  colnames(phenotype_overlap_chrs) <- c("chr", "start", "end")
  phenotype_gene_chr <- phenotype[,c("chr", "start", "end")]
  setdiff(phenotype_gene_chr, phenotype_overlap_chrs) 
}

process_cnv_gene_overlap_result <- function(gene_result) {
  gene_result_process <- gene_result
  gene_result_process$gene_coordinations <- paste0(gene_result_process$gene_chr, ':',
                                                   gene_result_process$gene_start, '-',
                                                   gene_result_process$gene_end)
  coding_gene <- gene_result_process[gene_result_process$gene_class=='protein_coding',]
  noncoding_gene <- gene_result_process
  noncoding_gene <- filter(noncoding_gene, grepl("noncoding",gene_result_process$gene_class))
  
  coding_gene_grouped <- coding_gene %>%
    group_by(cnv_chr, cnv_start, cnv_end) %>%
    summarise(`#_coding_genes` = n(),
              coding_gene_symbol = str_c(gene_symbol, collapse=","),
              coding_gene_name = str_c(gene_name, collapse=","),
              coding_gene_coordinations = str_c(gene_coordinations, collapse=",")
    )
  noncoding_gene_grouped <- noncoding_gene %>%
    group_by(cnv_chr, cnv_start, cnv_end) %>%
    summarise(`#_noncoding_genes` = n(),
              noncoding_gene_symbol = str_c(gene_symbol, collapse=","),
              noncoding_gene_name = str_c(gene_name, collapse=","),
              noncoding_gene_coordinations = str_c(gene_coordinations, collapse=",")
    )
  gene_overlap_count <- merge(coding_gene_grouped, noncoding_gene_grouped, by = c("cnv_chr", "cnv_start", "cnv_end"), 
                              all.x = TRUE, all.y = TRUE)
  gene_overlap_count$`#_coding_genes`[is.na(gene_overlap_count$`#_coding_genes`)] = 0
  gene_overlap_count$`#_noncoding_genes`[is.na(gene_overlap_count$`#_noncoding_genes`)] = 0
  return(gene_overlap_count)
}

process_cnv_gene_overlap_result_10 <- function(gene_result) {
  gene_result_process <- gene_result
  gene_result_process_scnv <- gene_result_process[which(gene_result_process$isCNVBigger==FALSE &
                                                  gene_result_process$percentOverlapCNVs>10.0),]
  gene_result_process_sg <- gene_result_process[which(gene_result_process$isCNVBigger==TRUE &
                                                    gene_result_process$percentOverlapGenes>10.0),]
  gene_result_process <- rbind(gene_result_process_scnv, gene_result_process_sg)
  gene_result_process$gene_coordinations <- paste0(gene_result_process$gene_chr, ':',
                                                   gene_result_process$gene_start, '-',
                                                   gene_result_process$gene_end)
  coding_gene <- gene_result_process[gene_result_process$gene_class=='protein_coding',]
  noncoding_gene <- gene_result_process
  noncoding_gene <- filter(noncoding_gene, grepl("noncoding",noncoding_gene$gene_class))
  
  coding_gene_grouped <- coding_gene %>%
    group_by(cnv_chr, cnv_start, cnv_end) %>%
    summarise(`#_coding_genes` = n(),
              coding_gene_symbol = str_c(gene_symbol, collapse=","),
              coding_gene_name = str_c(gene_name, collapse=","),
              coding_gene_coordinations = str_c(gene_coordinations, collapse=",")
    )
  
  noncoding_gene_grouped <- noncoding_gene %>%
    group_by(cnv_chr, cnv_start, cnv_end) %>%
    summarise(`#_noncoding_genes` = n(),
              noncoding_gene_symbol = str_c(gene_symbol, collapse=","),
              noncoding_gene_name = str_c(gene_name, collapse=","),
              noncoding_gene_coordinations = str_c(gene_coordinations, collapse=",")
    )
  gene_overlap_count <- merge(coding_gene_grouped, noncoding_gene_grouped, by = c("cnv_chr", "cnv_start", "cnv_end"), 
                              all.x = TRUE, all.y = TRUE)
  gene_overlap_count$`#_coding_genes`[is.na(gene_overlap_count$`#_coding_genes`)] = 0
  gene_overlap_count$`#_noncoding_genes`[is.na(gene_overlap_count$`#_noncoding_genes`)] = 0
  return(gene_overlap_count)
}

###################################################### Interaction Analysis Enhancer regions ########################################################

calculate_interacting_interaction_table_enhancer <- function(phenotype, gene_overlap_result) {
  phenotype_overlap_interacting_loci <- findOverlaps(phenotype, interacting_locus)
  cnv_results_interacting_loci <- phenotype[queryHits(phenotype_overlap_interacting_loci)]
  results_interacting_loci <- interacting_locus[subjectHits(phenotype_overlap_interacting_loci)]
  
  cnv_results_interacting_loci_df <- data.frame(
    chr = as.vector(seqnames(cnv_results_interacting_loci)), 
    start = cnv_results_interacting_loci@ranges@start,
    end = cnv_results_interacting_loci@ranges@start+cnv_results_interacting_loci@ranges@width-1, 
    cnv_type = cnv_results_interacting_loci@elementMetadata@listData$cnv_type)
  
  results_interacting_loci_df <- data.frame(
    chr = as.vector(seqnames(results_interacting_loci)), 
    start = results_interacting_loci@ranges@start,
    end = results_interacting_loci@ranges@start+results_interacting_loci@ranges@width-1, 
    reads = results_interacting_loci@elementMetadata@listData$reads,
    paper = results_interacting_loci@elementMetadata@listData$paper)
  
  interacting_loci_df <- merge(cnv_results_interacting_loci_df, results_interacting_loci_df, by = 'row.names', all = TRUE)
  interacting_loci_df <- interacting_loci_df[, 2:ncol(interacting_loci_df)]
  interacting_loci_df$is_left_hand_side = TRUE
  colnames(interacting_loci_df) <- c('cnv_chr', 'cnv_start', 'cnv_end', 'cnv_type','hic_chr','hic_start', 'hic_end', 'reads', 'paper', 'is_left_hand_side')
  
  result_table_interacting <- data.frame(matrix(ncol=14,nrow=0, dimnames=list(NULL, c("chr", "cnv_start", "cnv_end", "cnv_type", "interacting_locus_start", "interacting_locus_end", 
                                                                                      "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                                                                                      "evol_locus_start", "evol_locus_end", 
                                                                                      "overlapped_enhancer_id", "overlapped_enhancer_type"))))
  for (i in 1:nrow(interacting_loci_df)) {
    print(interacting_loci_df[i,])
    res <- hic_analysis_function_interacting_enhancer(interacting_loci_df[i,], gene_overlap_result)
    result_table_interacting <- rbind(result_table_interacting, res) 
  }
  
  return(result_table_interacting)
}

calculate_evol_interaction_table_enhancer <- function(phenotype, gene_overlap_result) {
  phenotype_overlap_evol_loci <- findOverlaps(phenotype, evol_locus)
  cnv_results_evol_loci <- phenotype[queryHits(phenotype_overlap_evol_loci)]
  results_evol_loci <- evol_locus[subjectHits(phenotype_overlap_evol_loci)]
  
  cnv_results_evol_loci_df <- data.frame(
    chr = as.vector(seqnames(cnv_results_evol_loci)), 
    start = cnv_results_evol_loci@ranges@start,
    end = cnv_results_evol_loci@ranges@start+cnv_results_evol_loci@ranges@width-1, 
    cnv_type = cnv_results_evol_loci@elementMetadata@listData$cnv_type)
  
  results_evol_loci_df <- data.frame(
    chr = as.vector(seqnames(results_evol_loci)), 
    start = results_evol_loci@ranges@start,
    end = results_evol_loci@ranges@start+results_evol_loci@ranges@width-1, 
    reads = results_evol_loci@elementMetadata@listData$reads,
    paper = results_evol_loci@elementMetadata@listData$paper)
  
  evol_loci_df <- merge(cnv_results_evol_loci_df, results_evol_loci_df, by = 'row.names', all = TRUE)
  evol_loci_df <- evol_loci_df[, 2:ncol(evol_loci_df)]
  evol_loci_df$is_left_hand_side = FALSE
  colnames(evol_loci_df) <- c('cnv_chr', 'cnv_start', 'cnv_end', 'cnv_type','hic_chr','hic_start', 'hic_end', 'reads', 'paper', 'is_left_hand_side')
  
  result_table_evol <- data.frame(matrix(ncol=14,nrow=0, dimnames=list(NULL, c("chr", "cnv_start", "cnv_end", "cnv_type", "evol_locus_start", "evol_locus_end", 
                                                                               "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                                                                               "interacting_locus_start", "interacting_locus_end", 
                                                                               "overlapped_enhancer_id", "overlapped_enhancer_type"))))
  
  for (i in 1:nrow(evol_loci_df)) {
    print(evol_loci_df[i,])
    res <- hic_analysis_function_evol_enhancer(evol_loci_df[i,],gene_overlap_result)
    result_table_evol <- rbind(result_table_evol, res) 
  }
  return(result_table_evol)
}



hic_analysis_function_evol_enhancer <- function(r, gene_overlap_result) {
  res <- data.frame(matrix(ncol=17,nrow=0, dimnames=list(NULL, c("chr", "cnv_start", "cnv_end", "cnv_type", "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", "evol_locus_start", "evol_locus_end", 
                                                                 "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                                                                 "interacting_locus_start", "interacting_locus_end", 
                                                                 "overlapped_enhancer_id", "overlapped_enhancer_type"))))
  r_df <- as.data.frame(r)
  
  no_coding_genes <- gene_overlap_result[which(gene_overlap_result$cnv_chr %in% r_df$hic_chr &
                                                 gene_overlap_result$cnv_start %in% r_df$cnv_start &
                                                 gene_overlap_result$cnv_end %in% r_df$cnv_end),c("#_coding_genes")]
  
  no_noncoding_genes <- gene_overlap_result[which(gene_overlap_result$cnv_chr %in% r_df$hic_chr &
                                                    gene_overlap_result$cnv_start %in% r_df$cnv_start &
                                                    gene_overlap_result$cnv_end %in% r_df$cnv_end),c("#_noncoding_genes")]
  if (length(no_coding_genes) == 0) {
    no_coding_genes = 0
  }
  if (length(no_noncoding_genes) == 0) {
    no_noncoding_genes = 0
  }
  if (no_coding_genes > 0) {
    # If there's no coding genes -> non_coding_cnv = TRUE
    r_df$noncoding_cnv = FALSE
  } else {
    # If there's coding genes -> non_coding_cnv = FALSE
    r_df$noncoding_cnv = TRUE
  }
  # parse each interacting loci
  e_loci <- as.data.frame(r_df[c('hic_chr', 'hic_start', 'hic_end')])
  
  # Look for the interacting region
  i_loci <- hic_data[which(hic_data$chr %in% e_loci$hic_chr &
                             hic_data$evol_locus_start %in% e_loci$hic_start &
                             hic_data$evol_locus_end %in% e_loci$hic_end),]
  #print(i_loci)
  # Iterate the interacting region
  for (i in 1:nrow(i_loci)) {
    # First, determine if the cnv covers both side 
    single_i_loci <- as.data.frame(i_loci[i,])
    overlapped_other_start <- single_i_loci$interacting_locus_start
    overlapped_other_end <- single_i_loci$interacting_locus_end
    reads <- single_i_loci$`# reads`
    paper <- single_i_loci$Paper
    if (single_i_loci$chr == r_df$cnv_chr &&  (as.integer(r_df$cnv_start) > overlapped_other_end) && (as.integer(e_loci$hic_start) > as.integer(overlapped_other_end))) {
      print("the cnv has one side overlapped")
      overlap_with_both_hic_side <- FALSE
      
      # Find enhancer regions overlapped with evol regions
      query_range <- toGRanges(data.frame(chr=single_i_loci$chr, start=single_i_loci$interacting_locus_start, 
                                          end=single_i_loci$interacting_locus_end,  reads=single_i_loci$`# reads`, paper=single_i_loci$Paper), genome = "hg19")
      
      enhancer_res <- findOverlaps(query_range, enhancer_regions)
      enhancer_overlapped <- enhancer_regions[subjectHits(enhancer_res)]
      cnv_overlapped <- query_range[queryHits(enhancer_res)]
      
      if (countOverlaps(query_range, gene_regions)>1) {
        print("enhancer signals were overlapped with this region")
        mcols(enhancer_overlapped)$enhancer_id <- paste0(seqnames(enhancer_overlapped), ":", 
                                                         start(enhancer_overlapped), "-", 
                                                         end(enhancer_overlapped))
        overlapped_enhancer_id <- paste(mcols(enhancer_overlapped)$enhancer_id, collapse = ';')
        overlapped_enhancer_type <- paste(mcols(enhancer_overlapped)$type, collapse = ';')
        
      } else {
        print("No enhancer signals were overlapped with this region")
        overlapped_enhancer_id <- 'ignore'
        overlapped_enhancer_type <- 'ignore'
      }
      
    } else {
      # Meaning the cnv covers with both side  
      print("the cnv covers with the both sides")
      overlap_with_both_hic_side <- TRUE
      
      # wouldn't bother calculating the genes 
      overlapped_enhancer_id <- 'ignore'
      overlapped_enhancer_type <- 'ignore'
    }
    
    hic_interacted_gene_results <- data.frame(c(r_df$hic_chr, 
                                                r_df$cnv_start,
                                                r_df$cnv_end,
                                                r_df$cnv_type, 
                                                no_coding_genes,
                                                no_noncoding_genes, 
                                                r_df$noncoding_cnv, 
                                                r_df$hic_start,
                                                r_df$hic_end, 
                                                reads, 
                                                paper,
                                                r_df$is_left_hand_side, 
                                                overlap_with_both_hic_side, 
                                                overlapped_other_start, 
                                                overlapped_other_end, 
                                                overlapped_enhancer_id,
                                                overlapped_enhancer_type))
    
    print(hic_interacted_gene_results)
    hic_interacted_gene_results = t(hic_interacted_gene_results)
    colnames(hic_interacted_gene_results) <- c("chr", "cnv_start", "cnv_end", "cnv_type", "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv",
                                               "evol_locus_start", "evol_locus_end", 
                                               "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                                               "interacting_locus_start", "interacting_locus_end", 
                                               "overlapped_enhancer_id", "overlapped_enhancer_type")
    rownames(hic_interacted_gene_results) = NULL
    
    res[nrow(res) + 1, ] <- hic_interacted_gene_results
  }
  return(res)
}

hic_analysis_function_interacting_enhancer <- function(r, gene_overlap_result) {
  res <- data.frame(matrix(ncol=17,nrow=0, dimnames=list(NULL, c("chr", "cnv_start", "cnv_end", "cnv_type", 
                                                                 "no_of_coding_genes", "no_of_noncoding_genes", "non_coding_cnv", 
                                                                 "interacting_locus_start", "interacting_locus_end", "reads", "paper", 
                                                                 "is_left_hand_side", "overlap_with_both_hic_side", 
                                                                 "evol_locus_start", "evol_locus_end", 
                                                                 "overlapped_enhancer_id", "overlapped_enhancer_type"))))
  r_df <- as.data.frame(r)
  no_coding_genes <- gene_overlap_result[which(gene_overlap_result$cnv_chr %in% r_df$hic_chr &
                                                 gene_overlap_result$cnv_start %in% r_df$cnv_start &
                                                 gene_overlap_result$cnv_end %in% r_df$cnv_end),c("#_coding_genes")]
  
  no_noncoding_genes <- gene_overlap_result[which(gene_overlap_result$cnv_chr %in% r_df$hic_chr &
                                                    gene_overlap_result$cnv_start %in% r_df$cnv_start &
                                                    gene_overlap_result$cnv_end %in% r_df$cnv_end),c("#_noncoding_genes")]
  if (length(no_coding_genes) == 0) {
    no_coding_genes = 0
  }
  if (length(no_noncoding_genes) == 0) {
    no_noncoding_genes = 0
  }
  if (no_coding_genes > 0) {
    # If there's no coding genes -> non_coding_cnv = TRUE
    r_df$noncoding_cnv = FALSE
  } else {
    # If there's coding genes -> non_coding_cnv = FALSE
    r_df$noncoding_cnv = TRUE
  }  
  
  # parse each interacting loci
  i_loci <- as.data.frame(r_df[c('hic_chr', 'hic_start', 'hic_end')])
  
  # Look for the hic interacting region
  e_loci <- hic_data[which(hic_data$chr %in% i_loci$hic_chr &
                             hic_data$interacting_locus_start %in% i_loci$hic_start &
                             hic_data$interacting_locus_end %in% i_loci$hic_end),]
  print(e_loci)
  # Iterate the interacting region
  for (i in 1:nrow(e_loci)) {
    # First, determine if the cnv covers both side 
    single_e_loci <- as.data.frame(e_loci[i,])
    overlapped_other_start <- single_e_loci$evol_locus_start
    overlapped_other_end <- single_e_loci$evol_locus_end
    
    reads <- single_e_loci$`# reads`
    paper <- single_e_loci$Paper
    
    if (single_e_loci$chr == r_df$cnv_chr &&  (as.integer(overlapped_other_start) > as.integer(r_df$cnv_end)) && (as.integer(overlapped_other_start) > as.integer(i_loci$hic_end))) {
      print("the cnv has one side overlapped")
      overlap_with_both_hic_side <- FALSE
      
      # Find gene body overlapped with evol regions
      query_range <- toGRanges(data.frame(chr=single_e_loci$chr, start=single_e_loci$evol_locus_start, 
                                          end=single_e_loci$evol_locus_end,  reads=single_e_loci$`# reads`, paper=single_e_loci$Paper), genome = "hg19")
      
      enhancer_res <- findOverlaps(query_range, enhancer_regions)
      enhancer_overlapped <- enhancer_regions[subjectHits(enhancer_res)]
      cnv_overlapped <- query_range[queryHits(enhancer_res)]

      if (countOverlaps(query_range, gene_regions)>1) {
        print("enhancer signals were overlapped with this region")
        mcols(enhancer_overlapped)$enhancer_id <- paste0(seqnames(enhancer_overlapped), ":", 
                                                         start(enhancer_overlapped), "-", 
                                                         end(enhancer_overlapped))
        overlapped_enhancer_id <- paste(mcols(enhancer_overlapped)$enhancer_id, collapse = ';')
        overlapped_enhancer_type <- paste(mcols(enhancer_overlapped)$type, collapse = ';')
        
      } else {
        print("No enhancer signals were overlapped with this region")
        overlapped_enhancer_id <- 'ignore'
        overlapped_enhancer_type <- 'ignore'
      }
      
    } else {
      # Meaning the cnv covers with both side  
      print("the cnv covers with the both sides")
      overlap_with_both_hic_side <- TRUE
      
      # wouldn't bother calculating the genes 
      overlapped_enhancer_id <- 'ignore'
      overlapped_enhancer_type <- 'ignore'
    }
    
    hic_interacted_gene_results <- data.frame(c(r_df$hic_chr, 
                                                r_df$cnv_start,
                                                r_df$cnv_end,
                                                r_df$cnv_type, 
                                                no_coding_genes,
                                                no_noncoding_genes, 
                                                r_df$noncoding_cnv, 
                                                r_df$hic_start,
                                                r_df$hic_end, 
                                                reads, 
                                                paper,
                                                r_df$is_left_hand_side, 
                                                overlap_with_both_hic_side, 
                                                overlapped_other_start, 
                                                overlapped_other_end, 
                                                overlapped_enhancer_id,
                                                overlapped_enhancer_type))
    
    print(hic_interacted_gene_results)
    hic_interacted_gene_results = t(hic_interacted_gene_results)
    colnames(hic_interacted_gene_results) <- c("chr", "cnv_start", "cnv_end", "cnv_type", "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv",
                                               "evol_locus_start", "evol_locus_end", 
                                               "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                                               "interacting_locus_start", "interacting_locus_end", 
                                               "overlapped_enhancer_id", "overlapped_enhancer_type")
    rownames(hic_interacted_gene_results) = NULL
    
    res[nrow(res) + 1, ] <- hic_interacted_gene_results
    
    hic_interacted_gene_results <- data.frame(c(r_df$hic_chr, 
                                                r_df$cnv_start,
                                                r_df$cnv_end,
                                                r_df$cnv_type, 
                                                no_coding_genes,
                                                no_noncoding_genes, 
                                                r_df$noncoding_cnv, 
                                                r_df$hic_start,
                                                r_df$hic_end, 
                                                reads, 
                                                paper,
                                                r_df$is_left_hand_side, 
                                                overlap_with_both_hic_side, 
                                                overlapped_other_start, 
                                                overlapped_other_end, 
                                                overlapped_enhancer_id,
                                                overlapped_enhancer_type))
    
    print(hic_interacted_gene_results)
    hic_interacted_gene_results = t(hic_interacted_gene_results)
    colnames(hic_interacted_gene_results) <- c("chr", "cnv_start", "cnv_end", "cnv_type", "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv",
                                               "evol_locus_start", "evol_locus_end", 
                                               "reads", "paper", "is_left_hand_side", "overlap_with_both_hic_side", 
                                               "interacting_locus_start", "interacting_locus_end", 
                                               "overlapped_enhancer_id", "overlapped_enhancer_type")
    rownames(hic_interacted_gene_results) = NULL
    
    res[nrow(res) + 1, ] <- hic_interacted_gene_results
  }
  return(res)
}

process_cnv_gene_overlap_result_10 <- function(gene_result) {
  gene_result_process <- gene_result
  gene_result_process_scnv <- gene_result_process[which(gene_result_process$isCNVBigger==FALSE &
                                                          gene_result_process$percentOverlapCNVs>10.0),]
  gene_result_process_sg <- gene_result_process[which(gene_result_process$isCNVBigger==TRUE &
                                                        gene_result_process$percentOverlapGenes>10.0),]
  gene_result_process <- rbind(gene_result_process_scnv, gene_result_process_sg)
  gene_result_process$gene_coordinations <- paste0(gene_result_process$gene_chr, ':',
                                                   gene_result_process$gene_start, '-',
                                                   gene_result_process$gene_end)
  coding_gene <- gene_result_process[gene_result_process$gene_class=='protein_coding',]
  noncoding_gene <- gene_result_process
  noncoding_gene <- filter(noncoding_gene, grepl("noncoding",noncoding_gene$gene_class))
  
  coding_gene_grouped <- coding_gene %>%
    group_by(cnv_chr, cnv_start, cnv_end) %>%
    summarise(`#_coding_genes` = n(),
              coding_gene_symbol = str_c(gene_symbol, collapse=","),
              coding_gene_name = str_c(gene_name, collapse=","),
              coding_gene_coordinations = str_c(gene_coordinations, collapse=",")
    )
  
  noncoding_gene_grouped <- noncoding_gene %>%
    group_by(cnv_chr, cnv_start, cnv_end) %>%
    summarise(`#_noncoding_genes` = n(),
              noncoding_gene_symbol = str_c(gene_symbol, collapse=","),
              noncoding_gene_name = str_c(gene_name, collapse=","),
              noncoding_gene_coordinations = str_c(gene_coordinations, collapse=",")
    )
  gene_overlap_count <- merge(coding_gene_grouped, noncoding_gene_grouped, by = c("cnv_chr", "cnv_start", "cnv_end"), 
                              all.x = TRUE, all.y = TRUE)
  gene_overlap_count$`#_coding_genes`[is.na(gene_overlap_count$`#_coding_genes`)] = 0
  gene_overlap_count$`#_noncoding_genes`[is.na(gene_overlap_count$`#_noncoding_genes`)] = 0
  return(gene_overlap_count)
}

############################################RUN GENE ANALYSIS###################################
# Run all groups 
asd_cnv_gene_overlap_result <- generate_gene_overlapped_table(asd_cnv_regions)


asd_gene_overlap_count <- process_cnv_gene_overlap_result(asd_cnv_gene_overlap_result)
asd_gene_overlap_count$phenotype = "ASD"

asd_cnv_gene_overlap_result_10_prec <- generate_gene_overlapped_table_10_perc(asd_cnv_regions)
asd_cnv_gene_overlap_count_10_prec <- process_cnv_gene_overlap_result_10(asd_cnv_gene_overlap_result_10_prec)
asd_cnv_gene_overlap_count_10_prec$phenotype = "ASD"

setdiff_cnvs(asd_cnv_gene_overlap_result, asd_cnv_data)
phenotype_overlap_chrs <- asd_id_dd_gene_overlap_count[,c("cnv_chr", "cnv_start", "cnv_end")]
colnames(phenotype_overlap_chrs) <- c("chr", "start", "end")
phenotype_gene_chr <- asd_id_dd_cnv_data[,c("chr", "start", "end")]
setdiff(phenotype_gene_chr, phenotype_overlap_chrs)

setdiff_cnvs(asd_id_dd_gene_overlap_count, asd_id_dd_cnv_data)
setdiff_cnvs(id_dd_gene_overlap_count, id_dd_cnv_data)
setdiff_cnvs(scz_gene_overlap_count, scz_cnv_data)


asd_id_dd_cnv_gene_overlap_result <- generate_gene_overlapped_table(asd_id_dd_cnv_regions)
asd_id_dd_gene_overlap_count <- process_cnv_gene_overlap_result(asd_id_dd_cnv_gene_overlap_result)
asd_id_dd_gene_overlap_count$phenotype <- "ASD_ID_DD"

asd_id_dd_cnv_gene_overlap_result_10_prec <- generate_gene_overlapped_table_10_perc(asd_id_dd_cnv_regions)
asd_id_dd_gene_overlap_count_10_prec <- process_cnv_gene_overlap_result_10(asd_id_dd_cnv_gene_overlap_result_10_prec)
asd_id_dd_gene_overlap_count_10_prec$phenotype = "ASD_ID_DD"

id_dd_cnv_gene_overlap_result <- generate_gene_overlapped_table(id_dd_cnv_regions)
id_dd_gene_overlap_count <- process_cnv_gene_overlap_result(id_dd_cnv_gene_overlap_result)
id_dd_gene_overlap_count$phenotype <- "ID_DD"

id_dd_cnv_gene_overlap_result_10_prec <- generate_gene_overlapped_table_10_perc(id_dd_cnv_regions)
id_dd_gene_overlap_count_10_prec <- process_cnv_gene_overlap_result_10(id_dd_cnv_gene_overlap_result_10_prec)
id_dd_gene_overlap_count_10_prec$phenotype = "ID_DD"

scz_gene_overlap_result <- generate_gene_overlapped_table(scz_cnv_regions)
scz_gene_overlap_count <- process_cnv_gene_overlap_result(scz_gene_overlap_result)
scz_gene_overlap_count$phenotype <- "SCZ"

scz_gene_overlap_result_10_prec <- generate_gene_overlapped_table_10_perc(scz_cnv_regions)
scz_gene_overlap_count_10_prec <- process_cnv_gene_overlap_result_10(scz_gene_overlap_result_10_prec)
scz_gene_overlap_count_10_prec$phenotype <- "SCZ"

# Merge the 4? 
gene_overlap_result_merged <- rbind(asd_gene_overlap_count, asd_id_dd_gene_overlap_count, id_dd_gene_overlap_count, scz_gene_overlap_count)
write_csv(gene_overlap_result_merged, "/media/bml/USER_DATA/ASD_CNV_Project/Hi-C_INTERACTION/gene_overlap_result_merged.csv")

gene_overlap_result_merged_10 <- rbind(asd_cnv_gene_overlap_count_10_prec, asd_id_dd_gene_overlap_count_10_prec,
                                       id_dd_gene_overlap_count_10_prec, scz_gene_overlap_count_10_prec)
write_csv(gene_overlap_result_merged_10, "/media/bml/USER_DATA/ASD_CNV_Project/Hi-C_INTERACTION/gene_overlap_result_merged_10_overlap_small_reg.csv")

##########################################RUN ENHANCER OVERLAP###############################
asd_cnv_enhancer_overlap_result <- generate_gene_overlapped_table(asd_cnv_regions, enhancer_regions)
asd_cnv_enhancer_overlap_result$phenotype = "ASD"

asd_id_dd_cnv_enhancer_overlap_result <- generate_gene_overlapped_table(asd_id_dd_cnv_regions, enhancer_regions)
asd_id_dd_cnv_enhancer_overlap_result$phenotype = "ASD_ID_DD"

id_dd_cnv_enhancer_overlap_result <- generate_gene_overlapped_table(id_dd_cnv_regions, enhancer_regions)
id_dd_cnv_enhancer_overlap_result$phenotype = "ID_DD"

scz_cnv_enhancer_overlap_result <- generate_gene_overlapped_table(scz_cnv_regions, enhancer_regions)
scz_cnv_enhancer_overlap_result$phenotype = "SCZ"

enhancer_overlap_result_merged <- rbind(asd_cnv_enhancer_overlap_result, asd_id_dd_cnv_enhancer_overlap_result, 
                                        id_dd_cnv_enhancer_overlap_result, scz_cnv_enhancer_overlap_result)
write_csv(enhancer_overlap_result_merged, "enhancer_overlap_result_merged.csv")
############################################RUN INTERACTION ANALYSIS##############################################
asd_evol_locus_result <- calculate_evol_interaction_table(asd_cnv_regions)
asd_evol_locus_result_distinct <- distinct(asd_evol_locus_result)

colnames(asd_evol_locus_result_distinct) <- c("chr", "cnv_start", "cnv_end", "cnv_type", 
                                              "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv", 
                                              "evol_locus_start", "evol_locus_end", "reads", "paper", 
                                              "is_left_hand_side", "overlap_with_both_hic_side", 
                                              "interacting_locus_start", "interacting_locus_end", 
                                              "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", 
                                              "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                                              "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class", 
                                              "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end")
asd_evol_locus_result_distinct$phenotype = "ASD"
asd_evol_locus_interaction_table <- calculate_metrics(asd_evol_locus_result_distinct)
write_csv(asd_evol_locus_interaction_table, 'asd_evol_locus_interaction_table.csv')
write_csv(asd_evol_locus_result_distinct, 'asd_evol_locus_result_distinct.csv')

asd_id_dd_evol_locus_result <- calculate_evol_interaction_table(asd_id_dd_cnv_regions)
asd_id_dd_evol_locus_result_distinct <- distinct(asd_id_dd_evol_locus_result)
colnames(asd_id_dd_evol_locus_result_distinct) <- c("chr", "cnv_start", "cnv_end", "cnv_type", 
                                              "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv", 
                                              "evol_locus_start", "evol_locus_end", "reads", "paper", 
                                              "is_left_hand_side", "overlap_with_both_hic_side", 
                                              "interacting_locus_start", "interacting_locus_end", 
                                              "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", 
                                              "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                                              "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class", 
                                              "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end")
asd_id_dd_evol_locus_result_distinct$phenotype = "ASD_ID_DD"
asd_id_dd_evol_locus_interaction_table <- calculate_metrics(asd_id_dd_evol_locus_result_distinct)
write_csv(asd_id_dd_evol_locus_result_distinct, 'asd_id_dd_evol_locus_result_distinct.csv')
write_csv(asd_id_dd_evol_locus_interaction_table, 'asd_id_dd_evol_locus_interaction_table.csv')

scz_evol_locus_result <- calculate_evol_interaction_table(scz_cnv_regions)
scz_evol_locus_result_distinct <- distinct(scz_evol_locus_result)
colnames(scz_evol_locus_result_distinct) <- c("chr", "cnv_start", "cnv_end", "cnv_type", 
                                                    "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv", 
                                                    "evol_locus_start", "evol_locus_end", "reads", "paper", 
                                                    "is_left_hand_side", "overlap_with_both_hic_side", 
                                                    "interacting_locus_start", "interacting_locus_end", 
                                                    "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", 
                                                    "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                                                    "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class", 
                                                    "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end")
scz_evol_locus_result_distinct$phenotype = "SCZ"
scz_evol_locus_interaction_table <- calculate_metrics(scz_evol_locus_result_distinct)
write_csv(id_dd_evol_locus_interaction_table, 'id_dd_evol_locus_interaction_table.csv')
write_csv(scz_evol_locus_result_distinct, 'scz_evol_locus_result_distinct.csv')

id_dd_evol_locus_result <- calculate_evol_interaction_table(id_dd_cnv_regions)
id_dd_evol_locus_result_distinct <- distinct(id_dd_evol_locus_result)
colnames(id_dd_evol_locus_result_distinct) <- c("chr", "cnv_start", "cnv_end", "cnv_type", 
                                              "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv", 
                                              "evol_locus_start", "evol_locus_end", "reads", "paper", 
                                              "is_left_hand_side", "overlap_with_both_hic_side", 
                                              "interacting_locus_start", "interacting_locus_end", 
                                              "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", 
                                              "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                                              "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class", 
                                              "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end")
id_dd_evol_locus_result_distinct$phenotype = "ID_DD"

id_dd_evol_locus_interaction_table <- calculate_metrics(asd_id_dd_evol_locus_result_distinct)
write_csv(id_dd_evol_locus_interaction_table, 'id_dd_evol_locus_interaction_table.csv')
write_csv(id_dd_evol_locus_result_distinct, 'id_dd_evol_locus_result_distinct.csv')


asd_id_dd_interacting_locus_result <- calculate_interacting_interaction_table(asd_id_dd_cnv_regions)
asd_id_dd_interacting_locus_result_distinct <- distinct(asd_id_dd_interacting_locus_result)
colnames(asd_id_dd_interacting_locus_result_distinct) <- c("chr", "cnv_start", "cnv_end", "cnv_type", 
                                              "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv", 
                                              "evol_locus_start", "evol_locus_end", "reads", "paper", 
                                              "is_left_hand_side", "overlap_with_both_hic_side", 
                                              "interacting_locus_start", "interacting_locus_end", 
                                              "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", 
                                              "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                                              "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class", 
                                              "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end")
asd_id_dd_interacting_locus_result_distinct$phenotype = "ASD_ID_DD"
if (nrow(asd_id_dd_interacting_locus_result_distinct > 0)) {
  asd_id_dd_interacting_locus_interaction_table <- calculate_metrics(asd_id_dd_interacting_locus_result_distinct)
  write_csv(asd_id_dd_interacting_locus_interaction_table, 'asd_id_dd_interacting_locus_interaction_table.csv')
}
write_csv(asd_id_dd_interacting_locus_result_distinct, 'asd_id_dd_interacting_locus_result_distinct.csv')


id_dd_interacting_locus_result <- calculate_interacting_interaction_table(id_dd_cnv_regions)
id_dd_interacting_locus_result_distinct <- distinct(id_dd_interacting_locus_result)
colnames(id_dd_interacting_locus_result_distinct) <- c("chr", "cnv_start", "cnv_end", "cnv_type", 
                                                           "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv", 
                                                           "evol_locus_start", "evol_locus_end", "reads", "paper", 
                                                           "is_left_hand_side", "overlap_with_both_hic_side", 
                                                           "interacting_locus_start", "interacting_locus_end", 
                                                           "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", 
                                                           "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                                                           "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class", 
                                                           "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end")

if (nrow(id_dd_interacting_locus_result_distinct > 0)) {
  id_dd_interacting_locus_interaction_table <- calculate_metrics(id_dd_interacting_locus_result_distinct)
  write_csv(id_dd_interacting_locus_interaction_table, 'id_dd_interacting_locus_interaction_table.csv')
}
write_csv(id_dd_interacting_locus_result_distinct, 'id_dd_interacting_locus_result_distinct.csv')

scz_interacting_locus_result <- calculate_interacting_interaction_table(scz_cnv_regions)
scz_interacting_locus_result_distinct <- distinct(scz_interacting_locus_result)
colnames(scz_interacting_locus_result_distinct) <- c("chr", "cnv_start", "cnv_end", "cnv_type", 
                                                       "no_of_coding_genes", "no_of_noncoding_genes", "noncoding_cnv", 
                                                       "evol_locus_start", "evol_locus_end", "reads", "paper", 
                                                       "is_left_hand_side", "overlap_with_both_hic_side", 
                                                       "interacting_locus_start", "interacting_locus_end", 
                                                       "overlapped_gene_symbol", "overlapped_gene_name", "overlapped_gene_class", 
                                                       "overlapped_gene_end_cnv_start", "overlapped_gene_start_cnv_end", 
                                                       "overlapped_promoter_gene_symbol","overlapped_promoter_gene_name", "overlapped_promoter_gene_class", 
                                                       "overlapped_promoter_gene_end_cnv_start", "overlapped_promoter_gene_start_cnv_end")
if (nrow(scz_interacting_locus_result_distinct > 0)) {
  scz_interacting_locus_interaction_table <- calculate_metrics(scz_interacting_locus_result_distinct)
  write_csv(scz_interacting_locus_interaction_table, 'scz_interacting_locus_interaction_table.csv')
}
write_csv(scz_interacting_locus_result_distinct, 'scz_interacting_locus_result_distinct.csv')

# Merge them ALL
hic_interaction_results <- rbind(asd_interacting_locus_result_distinct, asd_evol_locus_result_distinct)



###################RUN INTERACTION ANALYSIS WITH 10% OVERLAPPED RESULT###########################
asd_evol_locus_result_10 <- calculate_evol_interaction_table(asd_cnv_regions, gene_overlap_result_merged_10)
asd_evol_locus_result_distinct_10 <- distinct(asd_evol_locus_result_10)
asd_evol_locus_result_distinct_10$phenotype = "ASD"
write_csv(asd_evol_locus_result_distinct_10, "asd_evol_locus_result_distinct_10.csv")

asd_interacting_locus_result_10 <- calculate_interacting_interaction_table(asd_cnv_regions, gene_overlap_result_merged_10)
asd_interacting_locus_result_distinct_10 <- distinct(asd_interacting_locus_result_10)
asd_interacting_locus_result_distinct_10$phenotype = "ASD"
write_csv(asd_interacting_locus_result_distinct_10, "asd_interacting_locus_result_distinct_10.csv")


asd_id_dd_evol_locus_result_10 <- calculate_evol_interaction_table(asd_id_dd_cnv_regions, gene_overlap_result_merged_10)
asd_id_dd_evol_locus_result_distinct_10 <- distinct(asd_id_dd_evol_locus_result_10)
asd_id_dd_evol_locus_result_distinct_10$phenotype = "ASD_ID_DD"
write_csv(asd_id_dd_evol_locus_result_distinct_10, "asd_id_dd_evol_locus_result_distinct_10.csv")

asd_id_dd_cnv_regions_first <- asd_id_dd_cnv_regions[1:13]

asd_id_dd_cnv_regions_second <- asd_id_dd_cnv_regions[46:69]


asd_id_dd_interacting_locus_first_result_10 <- calculate_interacting_interaction_table(asd_id_dd_cnv_regions_first, gene_overlap_result_merged_10)
asd_id_dd_interacting_locus_first_result_distinct_10 <- distinct(asd_id_dd_interacting_locus_first_result_10)
asd_id_dd_interacting_locus_first_result_distinct_10$phenotype = "ASD_ID_DD"

asd_id_dd_interacting_locus_4669_result_10 <- calculate_interacting_interaction_table(asd_id_dd_cnv_regions_second, gene_overlap_result_merged_10)
asd_id_dd_interacting_locus_4669_result_distinct_10 <- distinct(asd_id_dd_interacting_locus_4669_result_10)
asd_id_dd_interacting_locus_4669_result_distinct_10$phenotype = "ASD_ID_DD"
write_csv(asd_id_dd_interacting_locus_4669_result_10, "asd_id_dd_interacting_locus_4669_result_distinct_10.csv")

id_dd_evol_locus_result_10 <- calculate_evol_interaction_table(id_dd_cnv_regions, gene_overlap_result_merged_10)
id_dd_evol_locus_result_distinct_10 <- distinct(id_dd_evol_locus_result_10)
id_dd_evol_locus_result_distinct_10$phenotype = "ID_DD"
write_csv(id_dd_evol_locus_result_distinct_10, "id_dd_evol_locus_result_distinct_10.csv")

##HAVE PROBLEMS TO CAUSE RSESSION QUIT
id_dd_interacting_locus_result_10 <- calculate_interacting_interaction_table(id_dd_cnv_regions, gene_overlap_result_merged_10)
id_dd_interacting_locus_result_distinct_10 <- distinct(id_dd_interacting_locus_result_10)
id_dd_interacting_locus_result_distinct_10$phenotype = "ID_DD"
write_csv(id_dd_interacting_locus_result_distinct_10, "id_dd_interacting_locus_result_distinct_10.csv")


scz_evol_locus_result_10 <- calculate_evol_interaction_table(scz_cnv_regions, gene_overlap_result_merged_10)
scz_evol_locus_result_distinct_10 <- distinct(scz_evol_locus_result_10)
scz_evol_locus_result_distinct_10$phenotype = "SCZ"
write_csv(scz_evol_locus_result_distinct_10, "scz_evol_locus_result_distinct_10.csv")


scz_interacting_locus_result_10 <- calculate_interacting_interaction_table(scz_cnv_regions, gene_overlap_result_merged_10)
scz_interacting_locus_result_distinct_10 <- distinct(scz_interacting_locus_result_10)
scz_interacting_locus_result_distinct_10$phenotype = "SCZ"
write_csv(scz_interacting_locus_result_distinct_10, "scz_interacting_locus_result_distinct_10.csv")




###################RUN INTERACTION ANALYSIS ENHANCER WITH 10% OVERLAPPED RESULT###########################
asd_evol_locus_result_10_enhancer <- calculate_evol_interaction_table_enhancer(asd_cnv_regions, gene_overlap_result_merged_10)
asd_evol_locus_result_distinct_10_enhancer <- distinct(asd_evol_locus_result_10_enhancer)
asd_evol_locus_result_distinct_10_enhancer$phenotype = "ASD"
write_csv(asd_evol_locus_result_distinct_10_enhancer, "asd_evol_locus_result_distinct_10_enhancer.csv")

asd_interacting_locus_result_10_enhancer <- calculate_interacting_interaction_table_enhancer(asd_cnv_regions, gene_overlap_result_merged_10)
asd_interacting_locus_result_distinct_10_enhancer <- distinct(asd_interacting_locus_result_10_enhancer)
asd_interacting_locus_result_distinct_10_enhancer$phenotype = "ASD"
write_csv(asd_interacting_locus_result_distinct_10_enhancer, "asd_interacting_locus_result_distinct_10_enhancer.csv")

asd_id_dd_evol_locus_result_10_enhancer <- calculate_evol_interaction_table_enhancer(asd_id_dd_cnv_regions, gene_overlap_result_merged_10)
asd_id_dd_evol_locus_result_distinct_10_enhancer <- distinct(asd_id_dd_evol_locus_result_10_enhancer)
asd_id_dd_evol_locus_result_distinct_10_enhancer$phenotype = "ASD_ID_DD"
write_csv(asd_id_dd_evol_locus_result_distinct_10_enhancer, "asd_id_dd_evol_locus_result_distinct_10_enhancer.csv")

asd_id_dd_cnv_regions_0105 <- asd_id_dd_cnv_regions[1:5]
asd_id_dd_cnv_regions_0610 <- asd_id_dd_cnv_regions[6:10]
asd_id_dd_cnv_regions_1120 <- asd_id_dd_cnv_regions[11:20]
asd_id_dd_cnv_regions_2130 <- asd_id_dd_cnv_regions[21:30]
asd_id_dd_cnv_regions_3140 <- asd_id_dd_cnv_regions[31:40]
asd_id_dd_cnv_regions_4150 <- asd_id_dd_cnv_regions[41:50]
asd_id_dd_cnv_regions_5160 <- asd_id_dd_cnv_regions[51:60]
asd_id_dd_cnv_regions_6169 <- asd_id_dd_cnv_regions[61:69]

asd_id_dd_interacting_locus_result_10_0105 <- calculate_interacting_interaction_table_enhancer(asd_id_dd_cnv_regions_0105, gene_overlap_result_merged_10)
asd_id_dd_interacting_locus_result_10_0105_distinct <- distinct(asd_id_dd_interacting_locus_result_10_0105)
asd_id_dd_interacting_locus_result_10_0105_distinct$phenotype = "ASD_ID_DD"
write_csv(asd_id_dd_interacting_locus_result_10_0105_distinct, "asd_id_dd_interacting_locus_result_10_0105_distinct.csv")

asd_id_dd_interacting_locus_result_10_0610 <- calculate_interacting_interaction_table_enhancer(asd_id_dd_cnv_regions_0610, gene_overlap_result_merged_10)
asd_id_dd_interacting_locus_result_10_0610_distinct <- distinct(asd_id_dd_interacting_locus_result_10_0610)
asd_id_dd_interacting_locus_result_10_0610_distinct$phenotype = "ASD_ID_DD"
write_csv(asd_id_dd_interacting_locus_result_10_0610_distinct, "asd_id_dd_interacting_locus_result_10_0610_distinct.csv")

asd_id_dd_interacting_locus_result_10_1120 <- calculate_interacting_interaction_table_enhancer(asd_id_dd_cnv_regions_1120, gene_overlap_result_merged_10)
asd_id_dd_interacting_locus_result_10_1120_distinct <- distinct(asd_id_dd_interacting_locus_result_10_1120)
asd_id_dd_interacting_locus_result_10_1120_distinct$phenotype = "ASD_ID_DD"
write_csv(asd_id_dd_interacting_locus_result_10_1120_distinct, "asd_id_dd_interacting_locus_result_10_1120_distinct.csv")

asd_id_dd_interacting_locus_result_10_2130 <- calculate_interacting_interaction_table_enhancer(asd_id_dd_cnv_regions_2130, gene_overlap_result_merged_10)
asd_id_dd_interacting_locus_result_10_2130_distinct <- distinct(asd_id_dd_interacting_locus_result_10_2130)
asd_id_dd_interacting_locus_result_10_2130_distinct$phenotype = "ASD_ID_DD"
write_csv(asd_id_dd_interacting_locus_result_10_2130_distinct, "asd_id_dd_interacting_locus_result_10_2130_distinct.csv")

asd_id_dd_interacting_locus_result_10_3140 <- calculate_interacting_interaction_table_enhancer(asd_id_dd_cnv_regions_3140, gene_overlap_result_merged_10)
asd_id_dd_interacting_locus_result_10_3140_distinct <- distinct(asd_id_dd_interacting_locus_result_10_3140)
asd_id_dd_interacting_locus_result_10_3140_distinct$phenotype = "ASD_ID_DD"
write_csv(asd_id_dd_interacting_locus_result_10_3140_distinct, "asd_id_dd_interacting_locus_result_10_3140_distinct.csv")

asd_id_dd_interacting_locus_result_10_4150 <- calculate_interacting_interaction_table_enhancer(asd_id_dd_cnv_regions_4150, gene_overlap_result_merged_10)
asd_id_dd_interacting_locus_result_10_4150_distinct <- distinct(asd_id_dd_interacting_locus_result_10_4150)
asd_id_dd_interacting_locus_result_10_4150_distinct$phenotype = "ASD_ID_DD"
write_csv(asd_id_dd_interacting_locus_result_10_4150_distinct, "asd_id_dd_interacting_locus_result_10_4150_distinct.csv")

asd_id_dd_interacting_locus_result_10_5160 <- calculate_interacting_interaction_table_enhancer(asd_id_dd_cnv_regions_5160, gene_overlap_result_merged_10)
asd_id_dd_interacting_locus_result_10_5160_distinct <- distinct(asd_id_dd_interacting_locus_result_10_5160)
asd_id_dd_interacting_locus_result_10_5160_distinct$phenotype = "ASD_ID_DD"
write_csv(asd_id_dd_interacting_locus_result_10_5160_distinct, "asd_id_dd_interacting_locus_result_10_5160_distinct.csv")

asd_id_dd_interacting_locus_result_10_6169 <- calculate_interacting_interaction_table_enhancer(asd_id_dd_cnv_regions_6169, gene_overlap_result_merged_10)
asd_id_dd_interacting_locus_result_10_6169_distinct <- distinct(asd_id_dd_interacting_locus_result_10_6169)
asd_id_dd_interacting_locus_result_10_6169_distinct$phenotype = "ASD_ID_DD"
write_csv(asd_id_dd_interacting_locus_result_10_6169_distinct, "asd_id_dd_interacting_locus_result_10_6169_distinct.csv")

id_dd_evol_locus_result_10_enhancer <- calculate_evol_interaction_table_enhancer(id_dd_cnv_regions, gene_overlap_result_merged_10)
id_dd_evol_locus_result_distinct_10_enhancer <- distinct(id_dd_evol_locus_result_10_enhancer)
id_dd_evol_locus_result_distinct_10_enhancer$phenotype = "ID_DD"
write_csv(id_dd_evol_locus_result_distinct_10_enhancer, "id_dd_evol_locus_result_distinct_10_enhancer.csv")

##HAVE PROBLEMS TO CAUSE RSESSION QUIT
id_dd_cnv_0150 <- id_dd_cnv_regions[1:50]
id_dd_cnv_51100 <- id_dd_cnv_regions[51:100]
id_dd_cnv_101137 <- id_dd_cnv_regions[101:137]

id_dd_interacting_locus_result_0150_enhancer <- calculate_interacting_interaction_table(id_dd_cnv_0150, gene_overlap_result_merged_10)
id_dd_interacting_locus_result_distinct_0150_enhancer <- distinct(id_dd_interacting_locus_result_0150_enhancer)
id_dd_interacting_locus_result_distinct_0150_enhancer$phenotype = "ID_DD"
write_csv(id_dd_interacting_locus_result_distinct_0150_enhancer, "id_dd_interacting_locus_result_distinct_0150_enhancer.csv")

id_dd_interacting_locus_result_51100_enhancer <- calculate_interacting_interaction_table(id_dd_cnv_51100, gene_overlap_result_merged_10)
id_dd_interacting_locus_result_distinct_51100_enhancer <- distinct(id_dd_interacting_locus_result_51100_enhancer)
id_dd_interacting_locus_result_distinct_51100_enhancer$phenotype = "ID_DD"
write_csv(id_dd_interacting_locus_result_distinct_51100_enhancer, "id_dd_interacting_locus_result_distinct_51100_enhancer.csv")

id_dd_interacting_locus_result_101137_enhancer <- calculate_interacting_interaction_table(id_dd_cnv_101137, gene_overlap_result_merged_10)
id_dd_interacting_locus_result_distinct_101137_enhancer <- distinct(id_dd_interacting_locus_result_101137_enhancer)
id_dd_interacting_locus_result_distinct_101137_enhancer$phenotype = "ID_DD"
write_csv(id_dd_interacting_locus_result_distinct_101137_enhancer, "id_dd_interacting_locus_result_distinct_101137_enhancer.csv")

scz_evol_locus_result_10_enhancer <- calculate_evol_interaction_table_enhancer(scz_cnv_regions, gene_overlap_result_merged_10)
scz_evol_locus_result_distinct_10_enhancer <- distinct(scz_evol_locus_result_10_enhancer)
scz_evol_locus_result_distinct_10_enhancer$phenotype = "SCZ"
write_csv(scz_evol_locus_result_distinct_10_enhancer, "scz_evol_locus_result_distinct_10_enhancer.csv")


scz_interacting_locus_result_10_enhancer <- calculate_interacting_interaction_table_enhancer(scz_cnv_regions, gene_overlap_result_merged_10)
scz_interacting_locus_result_distinct_10_enhancer <- distinct(scz_interacting_locus_result_10_enhancer)
scz_interacting_locus_result_distinct_10_enhancer$phenotype = "SCZ"
write_csv(scz_interacting_locus_result_distinct_10_enhancer, "scz_interacting_locus_result_distinct_10_enhancer.csv")


